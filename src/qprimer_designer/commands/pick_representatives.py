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
    parser.add_argument("--window-seqs", dest="window_seqs", default=None,
                       help="Output FASTA of ALL sequences trimmed to the design window (optional)")
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

def _column_entropy(seqs, col_idx):
    """Compute Shannon entropy for a single MSA column (ignoring gaps)."""
    bases = [s[col_idx].upper() for s in seqs if s[col_idx] != '-']
    if not bases:
        return 0.0
    n = len(bases)
    counts = {}
    for b in bases:
        counts[b] = counts.get(b, 0) + 1
    ent = 0.0
    for c in counts.values():
        p = c / n
        if p > 0:
            ent -= p * np.log2(p)
    return ent


def find_best_design_window(seqs, window_size=500, primer_len=20, min_gc=0.35):
    """Find window position that maximizes conservation-weighted primer coverage.

    For each column, compute: gap-free fraction × identity fraction.
    This rewards positions where many sequences are gap-free AND share the
    same base. For each primer position (primer_len bp), the score is the
    mean of these per-column scores. The design window score is the sum
    across all primer positions in the window.

    Only considers windows where the consensus GC content >= min_gc.
    Falls back to no GC filter if no window passes.
    """
    if not seqs:
        return 0

    seq_len = len(seqs[0])
    n_seqs = len(seqs)

    # Pre-compute per-column gap/N flags as numpy array (n_seqs x seq_len)
    gap_matrix = np.array([[c in ('-', 'N', 'n') for c in s] for s in seqs], dtype=np.int8)

    # Per-column conservation-weighted score:
    # gap_free_frac × identity_frac (both relative to n_seqs)
    from collections import Counter as _Counter
    col_scores = np.zeros(seq_len, dtype=np.float64)
    consensus = []
    for col in range(seq_len):
        col_bases = [s[col].upper() for s in seqs if s[col] not in ('-', 'N', 'n')]
        if col_bases:
            cnt = _Counter(col_bases)
            most_common_count = cnt.most_common(1)[0][1]
            consensus.append(cnt.most_common(1)[0][0])
            gap_free_frac = len(col_bases) / n_seqs
            identity_frac = most_common_count / len(col_bases)
            col_scores[col] = gap_free_frac * identity_frac
        else:
            consensus.append('N')
            col_scores[col] = 0.0

    # For each primer position, compute mean col_score across primer_len columns
    max_primer_start = seq_len - primer_len + 1
    if max_primer_start <= 0:
        return 0

    cum_col = np.concatenate([[0.0], np.cumsum(col_scores)])
    primer_scores = (cum_col[primer_len:max_primer_start + primer_len] -
                     cum_col[:max_primer_start]) / primer_len  # shape: (max_primer_start,)

    # Score each design window by sum of primer_scores
    n_primer_positions = window_size - primer_len + 1
    if n_primer_positions <= 0:
        return 0

    cum_ps = np.concatenate([[0.0], np.cumsum(primer_scores)])
    max_window_start = min(max_primer_start - n_primer_positions + 1, seq_len - window_size + 1)
    if max_window_start <= 0:
        return 0

    window_scores = cum_ps[n_primer_positions:max_window_start + n_primer_positions] - cum_ps[:max_window_start]

    # GC filtering
    gc_arr = np.array([1.0 if b in 'GC' else 0.0 for b in consensus])
    cum_gc = np.concatenate([[0.0], np.cumsum(gc_arr)])

    # Find best window with GC filter
    best_score = -1.0
    best_pos = -1
    best_gc = 0.0
    for w in range(max_window_start):
        gc = (cum_gc[w + window_size] - cum_gc[w]) / window_size
        if gc < min_gc:
            continue
        score = window_scores[w]
        if score > best_score:
            best_score = score
            best_pos = w
            best_gc = gc

    max_possible = n_primer_positions  # each primer position can score at most 1.0
    if best_pos >= 0:
        print(f"  Design window: pos {best_pos}, conservation score {best_score:.1f}/{max_possible}, GC={best_gc:.3f}")
        return best_pos

    # Fallback: no GC filter
    best_pos = int(np.argmax(window_scores))
    best_score = window_scores[best_pos]
    fb_gc = (cum_gc[best_pos + window_size] - cum_gc[best_pos]) / window_size
    print(f"  WARNING: No window with GC>={min_gc}. Fallback: pos {best_pos}, "
          f"conservation score {best_score:.1f}/{max_possible}, GC={fb_gc:.3f}")
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
    min_cluster_frac = float(params.get("MIN_CLUSTER_FRAC", 0.05))
    window_size = int(params.get("DESIGN_WINDOW", 500))
    kmer_size = int(params.get("MINHASH_KMER_SIZE", 10))

    print(f"Selecting representative sequences from {args.input}...")
    start_time = time.time()

    # Load sequences
    records = list(SeqIO.parse(args.input, "fasta"))
    seqs = [str(r.seq) for r in records]
    names = [r.id for r in records]

    primer_len = int(params.get("PRIMER_LEN_MIN", 20))

    if len(seqs) <= 2:
        # Fewer sequences than requested representatives - return all (window-trimmed)
        window_start = find_best_design_window(seqs, window_size, primer_len)
        trimmed = trim_to_window(seqs, window_start, window_size)
        with open(args.output, "w") as out:
            for r, t in zip(records, trimmed):
                out.write(f">{r.id}\n{t}\n")
        # Write empty probe regions (no MSA-based filtering needed)
        if getattr(args, 'probe_regions', None):
            with open(args.probe_regions, "w") as f:
                f.write("start\tend\tnum_variable_sites\n")
        # Window seqs = same as output when all seqs are representatives
        if getattr(args, 'window_seqs', None):
            with open(args.window_seqs, "w") as out:
                for r, t in zip(records, trimmed):
                    out.write(f">{r.id}\n{t}\n")
        print(f"Wrote {len(seqs)} sequences (all, window-trimmed) to {args.output}")
        return

    # Find best design window and trim
    window_start = find_best_design_window(seqs, window_size, primer_len)
    trimmed = trim_to_window(seqs, window_start, window_size)

    # Filter out empty sequences and N-rich sequences (>10% N)
    def _is_valid_trimmed(s):
        if len(s) < kmer_size:
            return False
        n_count = s.upper().count('N')
        return n_count / len(s) <= 0.1
    valid_indices = [i for i, s in enumerate(trimmed) if _is_valid_trimmed(s)]
    if not valid_indices:
        # Fall back to first num_reps sequences (window-trimmed)
        with open(args.output, "w") as out:
            for i in range(min(num_reps, len(records))):
                r = records[i]
                out.write(f">{r.id}\n{trimmed[i]}\n")
        if getattr(args, 'probe_regions', None):
            with open(args.probe_regions, "w") as f:
                f.write("start\tend\tnum_variable_sites\n")
        print(f"Wrote {min(num_reps, len(records))} sequences to {args.output}")
        return

    valid_trimmed = [trimmed[i] for i in valid_indices]

    # Subsample if too many sequences for O(n²) clustering
    MAX_CLUSTER_SEQS = 500
    if len(valid_trimmed) > MAX_CLUSTER_SEQS:
        rng = random.Random(42)
        sample_local = sorted(rng.sample(range(len(valid_trimmed)), MAX_CLUSTER_SEQS))
        sampled_trimmed = [valid_trimmed[i] for i in sample_local]
        sampled_valid_indices = [valid_indices[i] for i in sample_local]
        print(f"  Subsampled {len(valid_trimmed)} -> {MAX_CLUSTER_SEQS} sequences for clustering")
    else:
        sampled_trimmed = valid_trimmed
        sampled_valid_indices = valid_indices

    # Cluster and select medoids
    family = MinHashFamily(kmer_size=kmer_size, N=100)
    h = family.make_h()
    signatures = [h(s) for s in sampled_trimmed]

    clusters = cluster_sequences(sampled_trimmed, kmer_size=kmer_size)
    medoid_local = select_medoids(sampled_trimmed, clusters, signatures, family)

    # Compute cluster sizes and fractions
    from collections import Counter as _Counter
    cluster_sizes = _Counter(clusters)
    total_seqs_clustered = len(clusters)

    # Map medoid local index → cluster label → fraction
    medoid_cluster_frac = {}
    for ml in medoid_local:
        cluster_label = clusters[ml]
        frac = cluster_sizes[cluster_label] / total_seqs_clustered
        medoid_cluster_frac[ml] = round(frac, 4)

    # Map back to original indices
    medoid_indices = [sampled_valid_indices[i] for i in medoid_local]

    # Sort medoids by cluster size (largest first)
    sorted_pairs = sorted(
        zip(medoid_local, medoid_indices),
        key=lambda pair: medoid_cluster_frac[pair[0]],
        reverse=True,
    )
    medoid_local = [p[0] for p in sorted_pairs]
    medoid_indices = [p[1] for p in sorted_pairs]

    # Select medoids whose cluster represents >= min_cluster_frac of sequences
    selected_ml = []
    selected_idx = []
    excluded = []
    for ml, oi in zip(medoid_local, medoid_indices):
        frac = medoid_cluster_frac[ml]
        if frac >= min_cluster_frac:
            selected_ml.append(ml)
            selected_idx.append(oi)
        else:
            excluded.append((records[oi].id, frac))
    # Always keep at least 1 medoid (the largest cluster)
    if not selected_ml:
        selected_ml = [medoid_local[0]]
        selected_idx = [medoid_indices[0]]
    medoid_local = selected_ml
    medoid_indices = selected_idx

    total_coverage = sum(medoid_cluster_frac[ml] for ml in medoid_local)

    # Print cluster info
    print(f"  Clusters: {len(cluster_sizes)}, selected {len(medoid_indices)} medoids "
          f"(coverage={total_coverage:.4f}, threshold={min_cluster_frac})")
    for ml, oi in zip(medoid_local, medoid_indices):
        print(f"    {records[oi].id}: cluster_frac={medoid_cluster_frac[ml]:.4f}")
    if excluded:
        print(f"  Excluded {len(excluded)} small clusters (<{min_cluster_frac}):")
        for name, frac in excluded:
            print(f"    {name}: cluster_frac={frac:.4f}")

    # Build consensus from design window (majority-rule, non-gap bases)
    consensus_seq = []
    for col in range(window_start, min(window_start + window_size, len(seqs[0]))):
        bases = [s[col].upper() for s in seqs if s[col] not in ('-', 'N', 'n')]
        if bases:
            consensus_seq.append(_Counter(bases).most_common(1)[0][0])
    consensus_seq = ''.join(consensus_seq)

    # Write output (window-trimmed sequences, not full length)
    # Medoids + consensus, with cluster_frac in header
    with open(args.output, "w") as out:
        for ml, i in zip(medoid_local, medoid_indices):
            r = records[i]
            frac = medoid_cluster_frac[ml]
            out.write(f">{r.id} cluster_frac={frac}\n{trimmed[i]}\n")
        out.write(f">{args.name}_consensus cluster_frac=0\n{consensus_seq}\n")

    # Write ALL sequences trimmed to window (for bowtie2 indexing in quick mode)
    # Include all non-empty sequences, not just valid_indices (no N% or length filter)
    if getattr(args, 'window_seqs', None):
        with open(args.window_seqs, "w") as out:
            count = 0
            for i, r in enumerate(records):
                if len(trimmed[i]) > 0:
                    out.write(f">{r.id}\n{trimmed[i]}\n")
                    count += 1
        print(f"Wrote {count} window-trimmed sequences to {args.window_seqs}")

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
                f"[WARNING] No conserved probe regions found with current settings:\n"
                f"  PROBE_CONSERVATION_THRESHOLD = {probe_params['conservation_threshold']}\n"
                f"  PROBE_MAX_MISMATCHES = {probe_params['max_mismatches']}\n"
                f"  Probe design will be skipped for this target.\n"
                f"  To enable probes, try lowering PROBE_CONSERVATION_THRESHOLD (e.g. 0.6) or "
                f"increasing PROBE_MAX_MISMATCHES (e.g. 3-4) in params.txt.",
                file=sys.stderr,
            )

    runtime = time.time() - start_time
    n_written = len(medoid_indices) + 1  # medoids + consensus
    print(f"Wrote {n_written} representative sequences ({len(medoid_indices)} medoids + 1 consensus) to {args.output} ({runtime:.1f} sec)")
