"""Select representative sequences from MSA using MinHash clustering.

This module includes LSH code modified from CATCH (https://github.com/broadinstitute/catch),
released under the MIT License. Copyright (c) 2018 Hayden Metsky, Broad Institute.
"""

import argparse
import hashlib
import heapq
import math
import random
import sys
import time
from collections import defaultdict

import numpy as np
from Bio import SeqIO
from scipy.cluster import hierarchy

from qprimer_designer.utils import parse_params, get_probe_params


def register(subparsers):
    """Register the pick-representatives subcommand."""
    parser = subparsers.add_parser(
        "pick-representatives",
        help="Select representative sequences from MSA",
        description="""
Select a representative subset of sequences from a multiple sequence alignment
for downstream primer design using MinHash-based approximate clustering.
""",
    )
    parser.add_argument("--in", dest="input", required=True, help="Input aligned FASTA (MSA)")
    parser.add_argument("--out", dest="output", required=True, help="Output FASTA of representatives")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.add_argument("--name", required=True, help="Target name")
    parser.add_argument("--probe-regions", dest="probe_regions", default=None,
                       help="Output TSV of conserved probe regions (optional)")
    parser.set_defaults(func=run)


# --- LSH Classes (from CATCH) ---

class MinHashFamily:
    """LSH family using MinHash for Jaccard distance."""

    # BUG FIX: Added seed parameter for reproducibility.
    # Original __init__ did not accept a seed, and make_h() used random.randint()
    # without seeding, causing non-reproducible clustering results.
    # Original:
    # def __init__(self, kmer_size, N=1, use_fast_str_hash=False):
    #     self.kmer_size = kmer_size
    #     self.N = N
    #     self.use_fast_str_hash = use_fast_str_hash
    def __init__(self, kmer_size, N=1, use_fast_str_hash=False, seed=None):
        self.kmer_size = kmer_size
        self.N = N
        self.use_fast_str_hash = use_fast_str_hash
        self.seed = seed
        self._rng = random.Random(seed)  # Instance-specific RNG for reproducibility

    def make_h(self):
        p = 2**31 - 1
        # BUG FIX: Use instance-specific RNG instead of global random.
        # Original used random.randint() which is non-deterministic.
        # Original:
        # a = random.randint(1, p)
        # b = random.randint(0, p)
        a = self._rng.randint(1, p)
        b = self._rng.randint(0, p)

        if self.use_fast_str_hash:
            def kmer_hash(x):
                x_hash = abs(hash(x))
                return (a * x_hash + b) % p
        else:
            def kmer_hash(x):
                x_hash = int(hashlib.md5(x.encode('utf-8')).hexdigest(), 16)
                return (a * x_hash + b) % p

        def h(s):
            assert self.kmer_size <= len(s)
            num_kmers = len(s) - self.kmer_size + 1

            def kmer_hashes():
                num_yielded = 0
                while num_yielded < self.N:
                    for i in range(num_kmers):
                        kmer = s[i:(i + self.kmer_size)]
                        num_yielded += 1
                        yield kmer_hash(kmer)

            if self.N == 1:
                return tuple([min(kmer_hashes())])
            else:
                return tuple(sorted(heapq.nsmallest(self.N, kmer_hashes())))

        return h

    def P1(self, dist):
        return 1.0 - dist

    def estimate_jaccard_dist(self, hA, hB):
        hA_i, hB_i = 0, 0
        intersect_count = 0
        union_count = 0

        while hA_i < len(hA) and hB_i < len(hB):
            if union_count == self.N:
                break
            elif hA[hA_i] < hB[hB_i]:
                hA_i += 1
                union_count += 1
            elif hA[hA_i] > hB[hB_i]:
                hB_i += 1
                union_count += 1
            else:
                intersect_count += 1
                union_count += 1
                hA_i += 1
                hB_i += 1

        similarity = float(intersect_count) / union_count if union_count > 0 else 0
        return 1.0 - similarity


class HashConcatenation:
    """Concatenated hash functions (AND constructions)."""

    def __init__(self, family, k, join_as_str=False):
        self.family = family
        self.k = k
        self.join_as_str = join_as_str
        self.hs = [family.make_h() for _ in range(k)]

    def g(self, x):
        if self.join_as_str:
            # BUG FIX: h(x) returns a tuple, not a string, so str.join() fails.
            # Convert each hash result to string before joining.
            # Original:
            # return ''.join(h(x) for h in self.hs)
            return ''.join(str(h(x)) for h in self.hs)
        else:
            return tuple([h(x) for h in self.hs])


class NearNeighborLookup:
    """Support for approximate near neighbor lookups."""

    def __init__(self, family, k, dist_thres, dist_fn, reporting_prob, hash_idx=None, join_concat_as_str=False):
        self.family = family
        self.k = k
        self.dist_thres = dist_thres
        self.dist_fn = dist_fn
        self.hash_idx = hash_idx

        P1 = self.family.P1(dist_thres)
        if P1 == 1.0:
            self.num_tables = 1
        else:
            self.num_tables = math.log(1.0 - reporting_prob, 1.0 - math.pow(P1, k))
            self.num_tables = int(math.ceil(self.num_tables))

        self.hashtables = []
        self.hashtables_g = []
        for j in range(self.num_tables):
            g = HashConcatenation(self.family, self.k, join_as_str=join_concat_as_str)
            self.hashtables.append(defaultdict(set))
            self.hashtables_g.append(g)

    def add(self, pts):
        for j in range(self.num_tables):
            ht = self.hashtables[j]
            g = self.hashtables_g[j].g
            for p in pts:
                p_key = p[self.hash_idx] if self.hash_idx is not None else p
                ht[g(p_key)].add(p)

    def query(self, q, max_n=None):
        neighbors = set()
        q_key = q[self.hash_idx] if self.hash_idx is not None else q

        for j in range(self.num_tables):
            q_hash = self.hashtables_g[j].g(q_key)
            for p in self.hashtables[j].get(q_hash, set()):
                if p == q:
                    continue
                p_key = p[self.hash_idx] if self.hash_idx is not None else p
                dist = self.dist_fn(q_key, p_key)
                if dist <= self.dist_thres:
                    neighbors.add(p)
                if max_n and len(neighbors) >= max_n:
                    return neighbors
        return neighbors


# --- Clustering logic ---

def find_low_gap_window(seqs, window_size=200):
    """Find window position with minimum gap content."""
    if not seqs:
        return 0

    seq_len = len(seqs[0])
    best_pos = 0
    min_gaps = float('inf')

    # BUG FIX: Off-by-one error - original code missed the last valid window position.
    # For a 30bp sequence with window_size=20, last valid position is 10 (0-indexed).
    # Original range(0, 10) checked 0-9, missing position 10.
    # Fixed range(0, 11) checks 0-10 (all valid positions).
    # Original:
    # for pos in range(0, max(1, seq_len - window_size)):
    for pos in range(0, max(1, seq_len - window_size + 1)):
        gaps = sum(s[pos:pos + window_size].count('-') for s in seqs)
        if gaps < min_gaps:
            min_gaps = gaps
            best_pos = pos

    return best_pos


def trim_to_window(seqs, window_start, window_size):
    """Trim sequences to window and remove gaps."""
    trimmed = []
    for s in seqs:
        window = s[window_start:window_start + window_size]
        trimmed.append(window.replace('-', ''))
    return trimmed


def find_conserved_probe_regions(
    aligned_seqs,
    window_start,
    window_size,
    probe_len_max=28,
    max_variable_sites=2,
    conservation_threshold=0.8,
):
    """Find conserved regions within the design window suitable for probe design.

    Scans the MSA design window for sub-windows of probe_len_max where the
    number of variable columns is at most max_variable_sites.

    Args:
        aligned_seqs: List of aligned sequences (with gaps).
        window_start: Start position of the design window in MSA coordinates.
        window_size: Size of the design window.
        probe_len_max: Window size for scanning (use max probe length).
        max_variable_sites: Maximum variable columns allowed per window.
        conservation_threshold: Column conservation fraction below which
            a column is considered "variable".

    Returns:
        List of (start, end, num_variable) tuples in ungapped coordinates
        relative to the design window consensus.
    """
    from collections import Counter

    seq_len = len(aligned_seqs[0])
    win_end = min(window_start + window_size, seq_len)

    # Compute per-column conservation within the design window
    is_variable = []
    non_gap_columns = []  # Track which columns have non-gap bases
    for col_idx in range(window_start, win_end):
        bases = [s[col_idx].upper() for s in aligned_seqs if s[col_idx] != '-']
        if not bases:
            is_variable.append(True)
            non_gap_columns.append(False)
            continue
        non_gap_columns.append(True)
        cnt = Counter(bases)
        most_common_frac = cnt.most_common(1)[0][1] / len(bases)
        is_variable.append(most_common_frac < conservation_threshold)

    # Build list of conservation flags for non-gap columns only
    ungapped_variable = []
    for i, has_base in enumerate(non_gap_columns):
        if has_base:
            ungapped_variable.append(is_variable[i])

    if len(ungapped_variable) < probe_len_max:
        return []

    # Find windows where variable sites <= max_variable_sites
    conserved_starts = set()
    var_count = sum(ungapped_variable[:probe_len_max])

    for start in range(len(ungapped_variable) - probe_len_max + 1):
        if start > 0:
            var_count -= ungapped_variable[start - 1]
            var_count += ungapped_variable[start + probe_len_max - 1]
        if var_count <= max_variable_sites:
            conserved_starts.add(start)

    if not conserved_starts:
        return []

    # Merge overlapping windows into contiguous regions
    sorted_starts = sorted(conserved_starts)
    regions = []
    reg_start = sorted_starts[0]
    reg_end = sorted_starts[0] + probe_len_max

    for s in sorted_starts[1:]:
        if s <= reg_end:
            reg_end = s + probe_len_max
        else:
            # Count variable sites in this merged region
            n_var = sum(ungapped_variable[reg_start:reg_end])
            regions.append((reg_start, reg_end, n_var))
            reg_start = s
            reg_end = s + probe_len_max

    n_var = sum(ungapped_variable[reg_start:reg_end])
    regions.append((reg_start, reg_end, n_var))

    return regions


def cluster_sequences(seqs, kmer_size=10, signature_size=100, dist_threshold=0.3):
    """Cluster sequences using MinHash and hierarchical clustering."""
    if len(seqs) <= 1:
        return list(range(len(seqs)))

    family = MinHashFamily(kmer_size=kmer_size, N=signature_size)
    h = family.make_h()

    signatures = [h(s) for s in seqs]
    n = len(seqs)

    # Build condensed distance matrix
    dists = []
    for i in range(n):
        for j in range(i + 1, n):
            d = family.estimate_jaccard_dist(signatures[i], signatures[j])
            dists.append(d)

    if not dists:
        return list(range(len(seqs)))

    linkage = hierarchy.linkage(dists, method='average')
    clusters = hierarchy.fcluster(linkage, t=dist_threshold, criterion='distance')

    return clusters


def select_medoids(seqs, clusters, signatures, family):
    """Select medoid (most central sequence) from each cluster."""
    cluster_to_seqs = defaultdict(list)
    for i, c in enumerate(clusters):
        cluster_to_seqs[c].append(i)

    medoids = []
    for c, indices in cluster_to_seqs.items():
        if len(indices) == 1:
            medoids.append(indices[0])
        else:
            # Find sequence with minimum average distance to others
            min_avg_dist = float('inf')
            medoid = indices[0]
            for i in indices:
                dists = [family.estimate_jaccard_dist(signatures[i], signatures[j])
                         for j in indices if j != i]
                avg_dist = np.mean(dists) if dists else 0
                if avg_dist < min_avg_dist:
                    min_avg_dist = avg_dist
                    medoid = i
            medoids.append(medoid)

    return sorted(medoids)


def run(args):
    """Run the pick-representatives command."""
    params = parse_params(args.param_file)
    num_reps = int(params.get("NUM_REPRESENTATIVE_SEQS", 5))
    window_size = int(params.get("DESIGN_WINDOW", 200))
    kmer_size = int(params.get("MINHASH_KMER_SIZE", 10))

    print(f"Selecting representative sequences from {args.input}...")
    start_time = time.time()

    # Load sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    seqs = [str(r.seq) for r in records]
    names = [r.id for r in records]

    if len(seqs) <= num_reps:
        # Fewer sequences than requested representatives - return all
        with open(args.output, "w") as out:
            for r in records:
                out.write(f">{r.id}\n{str(r.seq).replace('-', '')}\n")
        # Write empty probe regions (no MSA-based filtering needed)
        if getattr(args, 'probe_regions', None):
            with open(args.probe_regions, "w") as f:
                f.write("start\tend\tnum_variable_sites\n")
        print(f"Wrote {len(seqs)} sequences (all) to {args.output}")
        return

    # Find low-gap window and trim
    window_start = find_low_gap_window(seqs, window_size)
    trimmed = trim_to_window(seqs, window_start, window_size)

    # Filter out empty sequences
    valid_indices = [i for i, s in enumerate(trimmed) if len(s) >= kmer_size]
    if not valid_indices:
        # Fall back to first num_reps sequences
        with open(args.output, "w") as out:
            for i in range(min(num_reps, len(records))):
                r = records[i]
                out.write(f">{r.id}\n{str(r.seq).replace('-', '')}\n")
        if getattr(args, 'probe_regions', None):
            with open(args.probe_regions, "w") as f:
                f.write("start\tend\tnum_variable_sites\n")
        print(f"Wrote {min(num_reps, len(records))} sequences to {args.output}")
        return

    valid_trimmed = [trimmed[i] for i in valid_indices]

    # Cluster and select medoids
    family = MinHashFamily(kmer_size=kmer_size, N=100)
    h = family.make_h()
    signatures = [h(s) for s in valid_trimmed]

    clusters = cluster_sequences(valid_trimmed, kmer_size=kmer_size)
    medoid_local = select_medoids(valid_trimmed, clusters, signatures, family)

    # Map back to original indices
    medoid_indices = [valid_indices[i] for i in medoid_local]

    # Limit to num_reps
    if len(medoid_indices) > num_reps:
        medoid_indices = medoid_indices[:num_reps]

    # Write output
    with open(args.output, "w") as out:
        for i in medoid_indices:
            r = records[i]
            out.write(f">{r.id}\n{str(r.seq).replace('-', '')}\n")

    # Compute conserved probe regions if requested
    if getattr(args, 'probe_regions', None):
        probe_params = get_probe_params(params)
        regions = find_conserved_probe_regions(
            aligned_seqs=seqs,
            window_start=window_start,
            window_size=window_size,
            probe_len_max=probe_params["len_max"],
            max_variable_sites=probe_params["max_mismatches"],
            conservation_threshold=probe_params["conservation_threshold"],
        )
        with open(args.probe_regions, "w") as f:
            f.write("start\tend\tnum_variable_sites\n")
            for start, end, n_var in regions:
                f.write(f"{start}\t{end}\t{n_var}\n")
        if regions:
            total_bp = sum(e - s for s, e, _ in regions)
            print(f"Found {len(regions)} conserved probe region(s) covering {total_bp} bp")
        else:
            print(
                f"[ERROR] No conserved probe regions found with current settings:\n"
                f"  PROBE_CONSERVATION_THRESHOLD = {probe_params['conservation_threshold']}\n"
                f"  PROBE_MAX_MISMATCHES = {probe_params['max_mismatches']}\n"
                f"  Try lowering PROBE_CONSERVATION_THRESHOLD (e.g. 0.6) or "
                f"increasing PROBE_MAX_MISMATCHES (e.g. 3-4) in params.txt.",
                file=sys.stderr,
            )
            sys.exit(1)

    runtime = time.time() - start_time
    print(f"Wrote {len(medoid_indices)} representative sequences to {args.output} ({runtime:.1f} sec)")
