"""Quick mode primer design with iterative batching and early stopping."""

import argparse
import ast
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from sklearn.preprocessing import MultiLabelBinarizer
from torch.utils.data import DataLoader

from Bio.SeqUtils import gc_fraction

from qprimer_designer.commands.generate import generate_primers_multi, generate_primers_single
from qprimer_designer.commands.prepare_input import run as run_prepare_input
from qprimer_designer.external import compute_batch_dimer_dg
from qprimer_designer.external.bowtie import build_index, find_bowtie2
from qprimer_designer.models import load_models, PcrDataset, FEATURE_COLUMNS
from qprimer_designer.utils import (
    encode_batch_parallel,
    get_tm,
    has_homopolymer,
    reverse_complement_dna,
    parse_params,
    get_primer_params,
    get_probe_params,
    sanitize_iupac,
    WOBBLE_W_PRIMER,
    WOBBLE_W_PROBE,
    wobble_mismatch_count_cols,
    wobble_mismatch_count_gapped,
)


def register(subparsers):
    """Register the quick-design subcommand."""
    parser = subparsers.add_parser(
        "quick-design",
        help="Quick primer design with iterative batching",
        description="""
Quick mode for primer design. Generates all primer candidates, then
evaluates batches of ~100 pairs with early stopping when coverage and
activity thresholds are met. Designed for multi-sequence targets where
full pipeline would take too long.
""",
    )
    parser.add_argument("--target", required=True, help="Original target FASTA (all sequences)")
    parser.add_argument("--selected", default=None,
                       help="Selected/representative target FASTA (required for single-window mode)")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.add_argument("--name", required=True, help="Target name")
    parser.add_argument("--msa", default=None,
                       help="Aligned MSA FASTA (enables multi-region amplicon scoring)")
    parser.add_argument("--target-window", dest="target_window", default=None,
                       help="Window-trimmed target FASTA (legacy single-window mode)")
    parser.add_argument("--tmpdir", default=None,
                       help="Directory for intermediate files (default: next to output CSV)")
    parser.add_argument("--init-fa", dest="init_fa", default=None,
                       help="Output FASTA with all evaluated primers (pipeline integration)")
    parser.add_argument("--init-feat", dest="init_feat", default=None,
                       help="Output features file for evaluated primers (pipeline integration)")
    parser.add_argument("--probe-mode", dest="probe_mode", action="store_true",
                       help="Enable probe-first design mode (find conserved probe regions first)")
    parser.add_argument("--probe-fa", dest="probe_fa", default=None,
                       help="Output FASTA for probe candidates (probe mode)")
    parser.add_argument("--probe-feat", dest="probe_feat", default=None,
                       help="Output features CSV for probe candidates (probe mode)")
    parser.add_argument("--probe-csv", dest="probe_csv", default=None,
                       help="Output probe mapping CSV (probe mode, replaces bowtie2 alignment)")
    parser.add_argument("--probe-pair-csv", dest="probe_pair_csv", default=None,
                       help="Output probe pair assignment CSV (probe mode, records which probe goes with which primer pair)")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.set_defaults(func=run)


PRIMERS_PER_BATCH = 100  # primers per batch (fwd+rev combined)
CONSENSUS_QUOTA_FRAC = 0.2  # consensus gets 20% of batch slots
MAX_BATCHES = 5
TOP_PAIRS = 1000  # top primer pairs to evaluate (after Tm/GC/dG filtering)
N_POSITIONS = 100  # amplicon positions to score (multi-region mode)
N_VARIANTS = 10   # primer variants per position per direction

# Primer scoring uses the more lenient wobble weights (ML validates further)
_WOBBLE_W = WOBBLE_W_PRIMER
_COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def _best_primer_base(template_freqs):
    """Pick primer base with highest effective frequency against template."""
    best_base, best_score = 'N', -1.0
    for pb in 'ATGC':
        score = sum(_WOBBLE_W.get((pb, tb), 0) * f for tb, f in template_freqs.items())
        if score > best_score:
            best_score = score
            best_base = pb
    return best_base, best_score


def _column_sense_freqs(aligned_seqs, col):
    """Get base frequency distribution at an MSA column (sense strand).

    Only counts unambiguous bases (A, T, G, C). IUPAC codes and gaps are skipped.
    """
    from collections import Counter
    bases = [s[col].upper() for s in aligned_seqs
             if s[col].upper() in ('A', 'T', 'G', 'C')]
    if not bases:
        return {}
    counts = Counter(bases)
    total = len(bases)
    return {b: c / total for b, c in counts.items()}


def generate_primer_variants(aligned_seqs, start_col, primer_len, strand, n_variants=7,
                              _freq_cache=None):
    """Generate consensus + wobble-optimized primer variants at a position.

    For fwd: primer bases = sense strand bases, template = antisense.
    For rev: primer bases = complement of sense, template = sense.
             Final sequence is reversed (RC construction).

    Args:
        _freq_cache: optional dict mapping col -> sense_freqs to avoid recomputing.

    Returns list of primer sequences (5'->3'), up to n_variants.
    """
    consensus_bases = []
    wobble_bases = []
    gains = []

    for i in range(primer_len):
        col = start_col + i
        if _freq_cache is not None and col in _freq_cache:
            sense_freqs = _freq_cache[col]
        else:
            sense_freqs = _column_sense_freqs(aligned_seqs, col)
        if not sense_freqs:
            consensus_bases.append('N')
            wobble_bases.append('N')
            gains.append((i, 0.0))
            continue

        if strand == 'fwd':
            cons_base = max(sense_freqs, key=sense_freqs.get)
            template_freqs = {_COMPLEMENT[b]: f for b, f in sense_freqs.items()}
        else:
            sens_cons = max(sense_freqs, key=sense_freqs.get)
            cons_base = _COMPLEMENT[sens_cons]
            template_freqs = sense_freqs

        consensus_bases.append(cons_base)

        wob_base, wob_score = _best_primer_base(template_freqs)
        wobble_bases.append(wob_base)

        cons_score = sum(_WOBBLE_W.get((cons_base, tb), 0) * f
                         for tb, f in template_freqs.items())
        gains.append((i, wob_score - cons_score))

    # Variable positions sorted by gain descending
    variable = [(i, g) for i, g in gains if wobble_bases[i] != consensus_bases[i]]
    variable.sort(key=lambda x: -x[1])

    seen = set()
    variants = []

    def _add(bases):
        key = tuple(bases)
        if key not in seen:
            seen.add(key)
            variants.append(list(bases))

    # 1. All consensus
    _add(consensus_bases)
    # 2. All wobble
    _add(wobble_bases)
    # 3. Single substitutions
    for pos, _ in variable:
        single = list(consensus_bases)
        single[pos] = wobble_bases[pos]
        _add(single)
        if len(variants) >= n_variants:
            break
    # 4. Progressive accumulation (consensus → wobble one by one)
    if len(variants) < n_variants:
        cur = list(consensus_bases)
        for pos, _ in variable:
            cur = list(cur)
            cur[pos] = wobble_bases[pos]
            _add(cur)
            if len(variants) >= n_variants:
                break
    # 5. Reverse progressive (wobble → consensus one by one)
    if len(variants) < n_variants:
        cur = list(wobble_bases)
        for pos, _ in reversed(variable):
            cur = list(cur)
            cur[pos] = consensus_bases[pos]
            _add(cur)
            if len(variants) >= n_variants:
                break

    result = []
    for v in variants[:n_variants]:
        if 'N' in v:
            continue
        if strand == 'rev':
            result.append(''.join(reversed(v)))
        else:
            result.append(''.join(v))
    return result


def _best_probe_base(template_freqs):
    """Pick probe base with highest effective frequency against template.

    Uses WOBBLE_W_PROBE (stricter: G-T/A-G = 0.50) since probes have no ML layer.
    """
    best_base, best_score = 'N', -1.0
    for pb in 'ATGC':
        score = sum(WOBBLE_W_PROBE.get((pb, tb), 0) * f
                    for tb, f in template_freqs.items())
        if score > best_score:
            best_score = score
            best_base = pb
    return best_base, best_score


def generate_probe_variants(aligned_seqs, msa_cols, strand, n_variants=5,
                             _freq_cache=None):
    """Generate consensus + wobble-optimized probe variants at a position.

    Same diversification logic as generate_primer_variants but uses
    WOBBLE_W_PROBE (stricter wobble weights).

    Args:
        aligned_seqs: list of MSA row strings
        msa_cols: list of MSA column indices (may be non-contiguous after gap trimming)
        strand: 'fwd' or 'rev'
        n_variants: max number of variants to generate
        _freq_cache: optional dict mapping col -> sense_freqs

    Returns list of probe sequences (5'->3'), up to n_variants.
    """
    probe_len = len(msa_cols)
    consensus_bases = []
    wobble_bases = []
    gains = []

    for i, col in enumerate(msa_cols):
        if _freq_cache is not None and col in _freq_cache:
            sense_freqs = _freq_cache[col]
        else:
            sense_freqs = _column_sense_freqs(aligned_seqs, col)
        if not sense_freqs:
            consensus_bases.append('N')
            wobble_bases.append('N')
            gains.append((i, 0.0))
            continue

        if strand == 'fwd':
            cons_base = max(sense_freqs, key=sense_freqs.get)
            template_freqs = {_COMPLEMENT[b]: f for b, f in sense_freqs.items()}
        else:
            sens_cons = max(sense_freqs, key=sense_freqs.get)
            cons_base = _COMPLEMENT[sens_cons]
            template_freqs = sense_freqs

        consensus_bases.append(cons_base)

        wob_base, wob_score = _best_probe_base(template_freqs)
        wobble_bases.append(wob_base)

        cons_score = sum(WOBBLE_W_PROBE.get((cons_base, tb), 0) * f
                         for tb, f in template_freqs.items())
        gains.append((i, wob_score - cons_score))

    # Variable positions sorted by gain descending
    variable = [(i, g) for i, g in gains if wobble_bases[i] != consensus_bases[i]]
    variable.sort(key=lambda x: -x[1])

    seen = set()
    variants = []

    def _add(bases):
        key = tuple(bases)
        if key not in seen:
            seen.add(key)
            variants.append(list(bases))

    # 1. All consensus
    _add(consensus_bases)
    # 2. All wobble
    _add(wobble_bases)
    # 3. Single substitutions
    for pos, _ in variable:
        single = list(consensus_bases)
        single[pos] = wobble_bases[pos]
        _add(single)
        if len(variants) >= n_variants:
            break
    # 4. Progressive accumulation
    if len(variants) < n_variants:
        cur = list(consensus_bases)
        for pos, _ in variable:
            cur = list(cur)
            cur[pos] = wobble_bases[pos]
            _add(cur)
            if len(variants) >= n_variants:
                break
    # 5. Reverse progressive
    if len(variants) < n_variants:
        cur = list(wobble_bases)
        for pos, _ in reversed(variable):
            cur = list(cur)
            cur[pos] = consensus_bases[pos]
            _add(cur)
            if len(variants) >= n_variants:
                break

    result = []
    for v in variants[:n_variants]:
        if 'N' in v:
            continue
        if strand == 'rev':
            result.append(''.join(reversed(v)))
        else:
            result.append(''.join(v))
    return result


def score_amplicon_positions(aligned_seqs, primer_len, min_amp_len, max_amp_len,
                             min_gc=0.35, max_gc=0.65,
                             gap_threshold=0.50, min_primer_score=0.80,
                             top_n=200, out_csv=None):
    """Score amplicon positions across the full MSA using gap-aware conservation.

    1. Columns with gap fraction > gap_threshold are skipped.
    2. Per-column identity_frac = most_common_count / non_gap_count.
    3. Consensus is built only from kept (non-gap-dominant) columns.
    4. 20-mer positions are scored by mean conservation within contiguous
       blocks of kept columns (no gap_score multiplication).
    5. GC filter on consensus.

    Returns:
        positions: sorted list of (fwd_start_msa, rev_start_msa, score)
        consensus_aligned: gap-aware consensus string (full MSA length,
            '-' at skipped columns)
        kept_cols: boolean array indicating which MSA columns are kept
    """
    seq_len = len(aligned_seqs[0])
    n_seqs = len(aligned_seqs)

    # Vectorize: convert to numpy array for fast column operations
    _msa = np.array([list(s.upper()) for s in aligned_seqs], dtype='U1')

    # Per-column base counts using vectorized operations
    identity = np.zeros(seq_len, dtype=np.float64)
    gap_frac = np.zeros(seq_len, dtype=np.float64)
    consensus_bases = []

    _is_gap = (_msa == '-') | (_msa == 'N')
    _gap_counts = _is_gap.sum(axis=0)
    gap_frac = _gap_counts / n_seqs

    # Count each base per column
    _base_A = (_msa == 'A').sum(axis=0)
    _base_T = (_msa == 'T').sum(axis=0)
    _base_G = (_msa == 'G').sum(axis=0)
    _base_C = (_msa == 'C').sum(axis=0)
    _n_bases = _base_A + _base_T + _base_G + _base_C

    for col in range(seq_len):
        if _n_bases[col] > 0 and gap_frac[col] <= gap_threshold:
            counts = {'A': int(_base_A[col]), 'T': int(_base_T[col]),
                       'G': int(_base_G[col]), 'C': int(_base_C[col])}
            most_common_base = max(counts, key=counts.get)
            consensus_bases.append(most_common_base)
            sense_freqs = {b: c / n_seqs for b, c in counts.items() if c > 0}
            anti_freqs = {_COMPLEMENT.get(b, b): f for b, f in sense_freqs.items()}
            _, fwd_score = _best_primer_base(anti_freqs)
            _, rev_score = _best_primer_base(sense_freqs)
            identity[col] = max(fwd_score, rev_score)
        else:
            consensus_bases.append('-')
            identity[col] = 0.0
    del _msa, _is_gap, _gap_counts, _base_A, _base_T, _base_G, _base_C, _n_bases

    # Boolean mask: kept columns (not gap-dominant)
    kept_cols = np.array([b != '-' for b in consensus_bases], dtype=bool)
    n_kept = kept_cols.sum()
    n_skipped = seq_len - n_kept
    print(f"  Gap-aware filtering: {n_kept} columns kept, {n_skipped} skipped "
          f"(gap fraction > {gap_threshold})")

    # Build consensus aligned string (full MSA length, '-' at skipped positions)
    consensus_aligned = ''.join(consensus_bases)

    # Identify contiguous blocks of kept columns for valid primer placement
    # A valid 20-mer must sit entirely within a contiguous block of kept columns.
    # Use cumsum of kept_cols to efficiently check: for a window [start, start+primer_len),
    # all columns must be kept → sum of kept_cols in window == primer_len
    max_start = seq_len - primer_len + 1
    if max_start <= 0:
        return [], consensus_aligned, kept_cols

    cum_kept = np.concatenate([[0], np.cumsum(kept_cols.astype(np.int32))])
    kept_in_window = cum_kept[primer_len:max_start + primer_len] - cum_kept[:max_start]
    valid_pos = (kept_in_window == primer_len)  # True if all columns in window are kept

    # Per-primer-position: mean conservation (only for valid positions)
    cum_id = np.concatenate([[0.0], np.cumsum(identity)])
    mean_cons = (cum_id[primer_len:max_start + primer_len] - cum_id[:max_start]) / primer_len

    # Primer scores = mean conservation (no gap_score multiplication)
    primer_scores = mean_cons.copy()
    # Zero out invalid positions (those spanning gap-dominant columns)
    primer_scores[~valid_pos] = 0.0

    # Save per-position scores CSV
    if out_csv:
        import csv
        half_w = primer_len // 2
        with open(out_csv, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['msa_position', 'conservation', 'valid', 'primer_score'])
            for pos in range(len(primer_scores)):
                center = pos + half_w
                w.writerow([center, f'{mean_cons[pos]:.6f}',
                           '1' if valid_pos[pos] else '0',
                           f'{primer_scores[pos]:.6f}'])
        print(f"  Saved position scores: {out_csv} ({len(primer_scores)} positions, "
              f"{valid_pos.sum()} valid)")

    # GC content per primer position (from consensus, only for valid positions)
    gc_arr = np.array([1.0 if b in 'GC' else 0.0 for b in consensus_bases])
    cum_gc = np.concatenate([[0.0], np.cumsum(gc_arr)])
    primer_gc = (cum_gc[primer_len:max_start + primer_len] - cum_gc[:max_start]) / primer_len

    # For each fwd_start, find best rev within valid amplicon range
    min_gap = 10
    results = []
    for fwd_start in range(len(primer_scores)):
        if not valid_pos[fwd_start]:
            continue
        fwd_score = primer_scores[fwd_start]
        if fwd_score < min_primer_score:
            continue
        fwd_gc = primer_gc[fwd_start]
        if fwd_gc < min_gc or fwd_gc > max_gc:
            continue

        rev_lo = fwd_start + min_amp_len - primer_len
        rev_hi = fwd_start + max_amp_len - primer_len + 1
        rev_lo = max(0, rev_lo)
        rev_hi = min(len(primer_scores), rev_hi)
        if rev_lo >= rev_hi:
            continue

        # Filter rev positions: must be valid and pass score/GC
        rev_valid_mask = valid_pos[rev_lo:rev_hi]
        rev_score_mask = primer_scores[rev_lo:rev_hi] >= min_primer_score
        rev_slice = primer_scores[rev_lo:rev_hi].copy()
        rev_gc_slice = primer_gc[rev_lo:rev_hi]
        gc_mask = ((rev_gc_slice >= min_gc) & (rev_gc_slice <= max_gc)
                   & rev_valid_mask & rev_score_mask)
        if not gc_mask.any():
            continue
        rev_slice[~gc_mask] = -1.0

        best_rev_local = np.argmax(rev_slice)
        best_rev_idx = rev_lo + best_rev_local
        if rev_slice[best_rev_local] < 0:
            continue

        amp_score = fwd_score + primer_scores[best_rev_idx]
        results.append((fwd_start, best_rev_idx, amp_score))

    results.sort(key=lambda x: -x[2])

    # Deduplicate: skip positions within min_gap of a higher-scoring one
    filtered = []
    used_fwd = []
    for fwd_start, rev_start, score in results:
        if any(abs(fwd_start - u) < min_gap for u in used_fwd):
            continue
        used_fwd.append(fwd_start)
        filtered.append((fwd_start, rev_start, score))

    return filtered[:top_n], consensus_aligned, kept_cols


def extract_primers_at_positions(aligned_seqs, positions, primer_len,
                                 min_tm, max_tm, max_gc, min_dg,
                                 n_variants=N_VARIANTS, max_pri_len=None):
    """Generate wobble-optimized primer variants at scored positions.

    For each position, generates up to n_variants fwd and rev primers
    at each length from primer_len to max_pri_len, mixing consensus
    and wobble-tolerant bases.  Filters Tm/GC/dG.

    Returns list of primer pair dicts and features dict.
    """
    if max_pri_len is None:
        max_pri_len = primer_len
    primer_lengths = list(range(primer_len, max_pri_len + 1))

    primers = []
    features = {}

    fwd_candidates = []  # (pos_idx, var_idx, seq, tm, gc, pos_score)
    rev_candidates = []

    # Pre-compute column frequencies for all columns needed across all positions
    seq_len = len(aligned_seqs[0])
    needed_cols = set()
    for fwd_start, rev_start, _ in positions:
        for plen in primer_lengths:
            for i in range(plen):
                if fwd_start + i < seq_len:
                    needed_cols.add(fwd_start + i)
                if rev_start + i < seq_len:
                    needed_cols.add(rev_start + i)
    freq_cache = {}
    for col in needed_cols:
        freq_cache[col] = _column_sense_freqs(aligned_seqs, col)

    for pos_idx, (fwd_start, rev_start, pos_score) in enumerate(positions):
        for plen in primer_lengths:
            if fwd_start + plen > seq_len or rev_start + plen > seq_len:
                continue
            fwd_vars = generate_primer_variants(aligned_seqs, fwd_start, plen, 'fwd', n_variants,
                                                 _freq_cache=freq_cache)
            for vi, fwd_seq in enumerate(fwd_vars):
                if 'N' in fwd_seq:
                    continue
                tm = get_tm(fwd_seq)
                gc = gc_fraction(fwd_seq)
                if min_tm <= tm <= max_tm and gc <= max_gc / 100.0:
                    fwd_candidates.append((pos_idx, vi, fwd_seq, tm, gc, pos_score))

            rev_vars = generate_primer_variants(aligned_seqs, rev_start, plen, 'rev', n_variants,
                                                 _freq_cache=freq_cache)
            for vi, rev_seq in enumerate(rev_vars):
                if 'N' in rev_seq:
                    continue
                tm = get_tm(rev_seq)
                gc = gc_fraction(rev_seq)
                if min_tm <= tm <= max_tm and gc <= max_gc / 100.0:
                    rev_candidates.append((pos_idx, vi, rev_seq, tm, gc, pos_score))

    # Batch self-dG
    all_seqs = [c[2] for c in fwd_candidates] + [c[2] for c in rev_candidates]
    if not all_seqs:
        return [], {}

    dg_values = compute_batch_dimer_dg([(s, s) for s in all_seqs])
    fwd_dgs = dg_values[:len(fwd_candidates)]
    rev_dgs = dg_values[len(fwd_candidates):]

    fwd_passed = {}
    for (pos_idx, vi, seq, tm, gc, pos_score), dg in zip(fwd_candidates, fwd_dgs):
        if dg < min_dg:
            continue
        fwd_passed[(pos_idx, vi)] = {
            'seq': seq, 'tm': tm, 'gc': gc, 'dg': dg, 'pos_score': pos_score,
        }
        features[seq] = {'len': len(seq), 'Tm': round(tm, 2), 'GC': round(gc, 4), 'dG': round(dg, 2)}

    rev_passed = {}
    for (pos_idx, vi, seq, tm, gc, pos_score), dg in zip(rev_candidates, rev_dgs):
        if dg < min_dg:
            continue
        rev_passed[(pos_idx, vi)] = {
            'seq': seq, 'tm': tm, 'gc': gc, 'dg': dg, 'pos_score': pos_score,
        }
        features[seq] = {'len': len(seq), 'Tm': round(tm, 2), 'GC': round(gc, 4), 'dG': round(dg, 2)}

    # Build all fwd×rev pairs at each position
    candidate_pairs = []
    for pos_idx, (fwd_start, rev_start, pos_score) in enumerate(positions):
        pos_fwd = [(vi, info) for (pi, vi), info in fwd_passed.items() if pi == pos_idx]
        pos_rev = [(vi, info) for (pi, vi), info in rev_passed.items() if pi == pos_idx]
        for fwd_vi, _ in pos_fwd:
            for rev_vi, _ in pos_rev:
                candidate_pairs.append((pos_idx, fwd_vi, rev_vi,
                                        fwd_start, rev_start, pos_score))

    # Batch cross-dimer dG
    if candidate_pairs:
        cross_seqs = [(fwd_passed[(p[0], p[1])]['seq'], rev_passed[(p[0], p[2])]['seq'])
                      for p in candidate_pairs]
        cross_dgs = compute_batch_dimer_dg(cross_seqs)
    else:
        cross_dgs = []

    for (pos_idx, fwd_vi, rev_vi, fwd_start, rev_start, pos_score), cross_dg in zip(candidate_pairs, cross_dgs):
        if cross_dg < min_dg:
            continue
        fwd_info = fwd_passed[(pos_idx, fwd_vi)]
        rev_info = rev_passed[(pos_idx, rev_vi)]
        primers.append({
            'pos_idx': pos_idx,
            'seq_idx': 0,
            'seq_id': f"var_f{fwd_vi}_r{rev_vi}",
            'is_consensus': (fwd_vi == 0 and rev_vi == 0),
            'fwd_start_msa': fwd_start,
            'rev_start_msa': rev_start,
            'fwd_seq': fwd_info['seq'],
            'rev_seq': rev_info['seq'],
            'fwd_tm': fwd_info['tm'],
            'rev_tm': rev_info['tm'],
            'fwd_gc': fwd_info['gc'],
            'rev_gc': rev_info['gc'],
            'fwd_dg_self': fwd_info['dg'],
            'rev_dg_self': rev_info['dg'],
            'pair_dg_cross': cross_dg,
            'score': pos_score,
        })

    return primers, features


def group_into_batches(primers, primer_len, max_ref_span=500):
    """Group primers by MSA position proximity, ensuring each batch's bowtie2
    reference stays compact (within max_ref_span MSA columns).

    Returns list of (batch_primers, min_fwd_start, max_rev_end) tuples.
    """
    if not primers:
        return []

    # Sort by fwd position
    sorted_primers = sorted(primers, key=lambda p: p['fwd_start_msa'])

    batches = []
    current_batch = [sorted_primers[0]]
    batch_min_fwd = sorted_primers[0]['fwd_start_msa']

    for p in sorted_primers[1:]:
        rev_end = p['rev_start_msa'] + primer_len
        # Check if adding this primer would make the reference too wide
        if rev_end - batch_min_fwd > max_ref_span:
            # Finalize current batch
            max_rev = max(pp['rev_start_msa'] + primer_len for pp in current_batch)
            batches.append((current_batch, batch_min_fwd, max_rev))
            current_batch = [p]
            batch_min_fwd = p['fwd_start_msa']
        else:
            current_batch.append(p)

    # Finalize last batch
    if current_batch:
        max_rev = max(pp['rev_start_msa'] + primer_len for pp in current_batch)
        batches.append((current_batch, batch_min_fwd, max_rev))

    return batches


def extract_window_sequences(aligned_seqs, seq_ids, win_start, win_end):
    """Extract window region from aligned sequences, strip gaps.

    Returns list of (seq_id, ungapped_seq) tuples.
    """
    result = []
    for seq_id, seq in zip(seq_ids, aligned_seqs):
        window = seq[win_start:win_end]
        ungapped = window.replace('-', '')
        if len(ungapped) > 0:
            result.append((seq_id, ungapped))
    return result


def _build_origin_map(target_seqs, for_filt, rev_filt,
                      step, min_pri_len, max_pri_len, min_amp_len):
    """Track which representative sequence each filtered primer originated from.

    Returns origin_map: dict mapping primer_seq → set of origin indices.
    The last index (len(target_seqs)-1) is the consensus sequence.
    """
    origin_map = {}  # primer_seq → set of origin indices
    for idx, tseq in enumerate(target_seqs):
        fwd, rev = generate_primers_single(tseq, step, min_pri_len, max_pri_len, min_amp_len)
        for seq in fwd:
            if seq in for_filt:
                origin_map.setdefault(seq, set()).add(idx)
        for seq in rev:
            if seq in rev_filt:
                origin_map.setdefault(seq, set()).add(idx)
    return origin_map


def _compute_position_bins(for_filt, n_bins):
    """Create n_bins overlapping windows across forward primer positions.

    Uses stride = range/n_bins and width = stride*2, so adjacent bins
    overlap by ~50%. This prevents primer pairs near bin boundaries
    from being missed.

    Returns list of (bin_start, bin_end) tuples.
    """
    positions = sorted(for_filt.values())
    if not positions:
        return []
    min_pos = positions[0]
    max_pos = positions[-1]
    total_range = max_pos - min_pos + 1
    stride = max(1, total_range / n_bins)
    width = stride * 2  # 50% overlap between adjacent bins
    return [(min_pos + i * stride, min_pos + i * stride + width)
            for i in range(n_bins)]


def _select_batch_primers(for_filt, rev_filt, features, origin_map,
                          n_origins, origin_weights, bin_start, bin_end,
                          min_amp_len, max_amp_len,
                          max_per_direction=PRIMERS_PER_BATCH):
    """Select primers from a position bin with proportional origin allocation.

    Each origin gets a quota proportional to its cluster_frac (weight).
    For each origin, forward primers are selected first, then for each
    forward primer at least one valid reverse primer is guaranteed.

    Args:
        origin_weights: dict mapping origin_idx → weight (cluster_frac).
            Consensus (last origin) gets equal share with smallest medoid cluster.

    Returns (fwd_list, rev_list) of primer sequences.
    """
    # Group fwd candidates by primary origin
    fwd_by_origin = {i: [] for i in range(n_origins)}
    for seq, pos in for_filt.items():
        if not (bin_start <= pos < bin_end):
            continue
        origins = origin_map.get(seq, set())
        if not origins:
            continue
        primary = min(origins)
        fwd_by_origin[primary].append(seq)

    # Group rev candidates by primary origin
    rev_min = bin_start + min_amp_len
    rev_max = bin_end + max_amp_len
    rev_by_origin = {i: [] for i in range(n_origins)}
    for seq, pos in rev_filt.items():
        if not (rev_min <= pos <= rev_max):
            continue
        origins = origin_map.get(seq, set())
        if not origins:
            continue
        primary = min(origins)
        rev_by_origin[primary].append(seq)

    # Sort each origin's candidates by rep_count descending
    for origin_idx in range(n_origins):
        fwd_by_origin[origin_idx].sort(
            key=lambda s: -features.get(s, {}).get("rep_count", 1))
        rev_by_origin[origin_idx].sort(
            key=lambda s: -features.get(s, {}).get("rep_count", 1))

    # Compute per-origin quota proportional to weights
    active_origins = [i for i in range(n_origins)
                      if fwd_by_origin[i] and rev_by_origin[i]]
    if not active_origins:
        # Fallback: try origins with at least fwd OR rev
        all_fwd = [s for i in range(n_origins) for s in fwd_by_origin[i]]
        all_rev = [s for i in range(n_origins) for s in rev_by_origin[i]]
        all_fwd.sort(key=lambda s: -features.get(s, {}).get("rep_count", 1))
        all_rev.sort(key=lambda s: -features.get(s, {}).get("rep_count", 1))
        return all_fwd[:max_per_direction], all_rev[:max_per_direction]

    total_weight = sum(origin_weights.get(i, 1.0) for i in active_origins)
    quotas = {}
    for i in active_origins:
        w = origin_weights.get(i, 1.0)
        quotas[i] = max(1, int(round(max_per_direction * w / total_weight)))

    # Select fwd+rev pairs per origin
    fwd_selected = []
    rev_selected = []
    rev_selected_set = set()

    for origin_idx in active_origins:
        quota = quotas[origin_idx]
        fwd_cands = fwd_by_origin[origin_idx]
        rev_cands = rev_by_origin[origin_idx]

        n_fwd_take = min(quota, len(fwd_cands))
        fwd_taken = fwd_cands[:n_fwd_take]
        fwd_selected.extend(fwd_taken)

        # For each fwd, ensure at least one valid rev exists in selection
        # Take up to quota rev primers, prioritizing those that form valid
        # amplicons with the selected fwd primers
        n_rev_take = min(quota, len(rev_cands))
        for rev_seq in rev_cands[:n_rev_take]:
            if rev_seq not in rev_selected_set:
                rev_selected.append(rev_seq)
                rev_selected_set.add(rev_seq)

    # Fill remaining slots from any origin
    remaining_fwd = max_per_direction - len(fwd_selected)
    remaining_rev = max_per_direction - len(rev_selected)

    if remaining_fwd > 0 or remaining_rev > 0:
        used_fwd = set(fwd_selected)
        used_rev = rev_selected_set
        leftover_fwd = []
        leftover_rev = []
        for i in range(n_origins):
            leftover_fwd.extend(s for s in fwd_by_origin[i] if s not in used_fwd)
            leftover_rev.extend(s for s in rev_by_origin[i] if s not in used_rev)
        leftover_fwd.sort(key=lambda s: -features.get(s, {}).get("rep_count", 1))
        leftover_rev.sort(key=lambda s: -features.get(s, {}).get("rep_count", 1))
        fwd_selected.extend(leftover_fwd[:remaining_fwd])
        rev_selected.extend(leftover_rev[:remaining_rev])

    return fwd_selected, rev_selected


def _write_batch_fasta_and_features(fwd_seqs, rev_seqs, features, name, batch_idx, tmpdir):
    """Write FASTA and features CSV for a batch of primers.

    Deduplicates sequences so each unique seq gets one FASTA entry and one
    name for bowtie2.  Returns (fasta_path, features_path, fwd_name_to_seq,
    rev_name_to_seq) where the dicts map assigned_name -> seq for every input
    entry (including duplicates).
    """
    fasta_path = Path(tmpdir) / f"batch_{batch_idx}.fa"
    feat_path = Path(tmpdir) / f"batch_{batch_idx}.feat"

    fwd_name_to_seq = {}  # pname -> seq  (all entries)
    rev_name_to_seq = {}
    fwd_seen = {}  # seq -> canonical pname (first seen)
    rev_seen = {}
    feat_rows = []

    with open(fasta_path, "w") as f:
        for i, seq in enumerate(fwd_seqs):
            pname = f"{name}_q{batch_idx}_{i+1}_f"
            fwd_name_to_seq[pname] = seq
            if seq not in fwd_seen:
                fwd_seen[seq] = pname
                f.write(f">{pname}\n{seq}\n")
                feat = features.get(seq, {})
                feat_rows.append({
                    "pname": pname, "pseq": seq, "forrev": "f",
                    "len": feat.get("len", len(seq)),
                    "Tm": feat.get("Tm", 0), "GC": feat.get("GC", 0),
                    "dG": feat.get("dG", 0),
                })

        for i, seq in enumerate(rev_seqs):
            pname = f"{name}_q{batch_idx}_{i+1}_r"
            rev_name_to_seq[pname] = seq
            if seq not in rev_seen:
                rev_seen[seq] = pname
                f.write(f">{pname}\n{seq}\n")
                feat = features.get(seq, {})
                feat_rows.append({
                    "pname": pname, "pseq": seq, "forrev": "r",
                    "len": feat.get("len", len(seq)),
                    "Tm": feat.get("Tm", 0), "GC": feat.get("GC", 0),
                    "dG": feat.get("dG", 0),
                })

    pd.DataFrame(feat_rows).to_csv(feat_path, index=False)
    return fasta_path, feat_path, fwd_name_to_seq, rev_name_to_seq


def _align_batch(fasta_path, index_prefix, n_targets, threads, tmpdir, batch_idx):
    """Align batch primers with bowtie2 and parse SAM into mapped format.

    Returns path to mapped file.
    """
    sam_path = Path(tmpdir) / f"batch_{batch_idx}.sam"
    mapped_path = Path(tmpdir) / f"batch_{batch_idx}.mapped"

    k = 50000
    cmd = [
        "bowtie2",
        "-x", str(index_prefix),
        "-U", str(fasta_path),
        "-f",
        "-p", str(threads),
        "-k", str(k),
        "--mp", "2,2", "--np", "2",
        "--rdg", "4,4", "--rfg", "4,4",
        "-L", "6", "-N", "1",
        "--score-min", "L,-0.6,-0.6",
        "--no-hd", "--no-unal",
        "-S", str(sam_path),
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    if not sam_path.exists() or sam_path.stat().st_size == 0:
        mapped_path.touch()
        return mapped_path

    # Parse SAM with sam2pairwise
    parsed_path = Path(tmpdir) / f"batch_{batch_idx}.parsed"
    with open(sam_path) as fin, open(parsed_path, "w") as fout:
        result = subprocess.run(
            ["sam2pairwise"], input=fin.read(), capture_output=True, text=True
        )
        fout.write(result.stdout)

    if parsed_path.stat().st_size == 0:
        mapped_path.touch()
        return mapped_path

    # Parse the sam2pairwise output and combine with SAM columns
    # sam2pairwise outputs 4-line groups: header, pseq, match, tseq
    sam_lines = []
    with open(sam_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # columns: pname, flag, tname, pos
            sam_lines.append((parts[0], parts[1], parts[2], parts[3]))

    parsed_lines = []
    with open(parsed_path) as f:
        lines = f.read().splitlines()

    # Group into 4-line blocks
    n_blocks = len(lines) // 4
    with open(mapped_path, "w") as out:
        for i in range(n_blocks):
            # pseq = line 1 (0-indexed: 4*i + 1), match = line 2, tseq = line 3
            pseq = lines[4 * i + 1]
            match = lines[4 * i + 2]
            tseq = lines[4 * i + 3]
            if i < len(sam_lines):
                pname, flag, tname, pos = sam_lines[i]
                out.write(f"{pname}\t{flag}\t{tname}\t{pos}\t{pseq}\t{tseq}\t{match}\n")

    return mapped_path


def _apply_coverage_filter(mapped_path, n_targets, tmpdir, batch_idx, cov_frac=0.95):
    """Filter primers covering less than cov_frac of target sequences. Returns filtered mapped path."""
    filtered_path = Path(tmpdir) / f"batch_{batch_idx}.mapped.filt"
    min_cov = int(n_targets * cov_frac)

    if mapped_path.stat().st_size == 0:
        filtered_path.touch()
        return filtered_path

    df = pd.read_csv(mapped_path, sep="\t", header=None,
                     names=["pname", "flag", "tname", "pos", "pseq", "tseq", "match"])

    # Count unique targets per primer
    cov = df.groupby("pname")["tname"].nunique()
    good_primers = set(cov[cov >= min_cov].index)

    if not good_primers:
        filtered_path.touch()
        return filtered_path

    df[df["pname"].isin(good_primers)].to_csv(
        filtered_path, sep="\t", header=False, index=False
    )
    return filtered_path


def _run_evaluate(input_path, ref_path, scaler, classifier, regressor, device, threads,
                   n_targets_global=None):
    """Run ML evaluation on prepared input. Returns DataFrame with scored pairs.

    Args:
        n_targets_global: if provided, use this as denominator for coverage
            instead of the batch-local target count.

    Returns:
        (res, eval_re, eval_full, eval_cl): res is the scored summary, eval_re is per-target
        activity, eval_full is the raw per-row ML output (for .eval.full file),
        eval_cl is the classifier per-target table (for .eval.cl file).
    """
    if not input_path.exists() or input_path.stat().st_size == 0:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    tnames = [s.id for s in SeqIO.parse(ref_path, "fasta")]
    n_for_coverage = n_targets_global if n_targets_global is not None else len(tnames)

    clstbl_list, regtbl_list = [], []
    full_chunks = []

    for chunk in pd.read_csv(input_path, chunksize=20000):
        if chunk.empty:
            continue
        chunk["targets"] = chunk["targets"].apply(ast.literal_eval)

        inps_fe = chunk[FEATURE_COLUMNS]
        inps_fe = scaler.transform(inps_fe)
        inps_se = chunk[["pseq_f", "tseq_f", "pseq_r", "tseq_r"]]
        inps_se = encode_batch_parallel(inps_se, threads)
        dataset = PcrDataset(inps_se, inps_fe, np.array([0] * len(chunk)))
        loader = DataLoader(dataset, batch_size=256, shuffle=False)

        predict_cls, predict_reg = [], []
        with torch.no_grad():
            for seq_in, fea_in, _ in loader:
                seq_in = seq_in.to(device).float()
                fea_in = fea_in.to(device).float()
                out_cls = classifier(fea_in, seq_in)
                out_reg = regressor(fea_in, seq_in)
                if len(seq_in) == 1:
                    predict_cls.append(np.array([out_cls.squeeze().detach().cpu().numpy()]))
                    predict_reg.append(np.array([out_reg.squeeze().detach().cpu().numpy()]))
                else:
                    predict_cls.append(out_cls.squeeze().detach().cpu().numpy())
                    predict_reg.append(out_reg.squeeze().detach().cpu().numpy())

        predict_cls = np.concatenate(predict_cls)
        predict_reg = np.round(np.concatenate(predict_reg), decimals=3)

        chunk.loc[:, "classifier"] = predict_cls
        chunk.loc[:, "regressor"] = predict_reg
        full_chunks.append(chunk.copy())

        mlb = MultiLabelBinarizer()
        onehot = mlb.fit_transform(chunk["targets"])
        target_cols = list(mlb.classes_)

        for label, lst in zip(["classifier", "regressor"], [clstbl_list, regtbl_list]):
            targets_df = pd.DataFrame(onehot, columns=target_cols, index=chunk.index)
            targets_df = targets_df.mul(chunk[label], axis=0)
            evaltbl = pd.concat([chunk[["pname_f", "pname_r"]], targets_df], axis=1)
            agg_dict = {c: "max" for c in evaltbl.columns[2:]}
            evaltbl = evaltbl.groupby(["pname_f", "pname_r"]).agg(agg_dict).reset_index()
            lst.append(evaltbl)

    if not clstbl_list:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    eval_full = pd.concat(full_chunks, ignore_index=True) if full_chunks else pd.DataFrame()

    clstbl = pd.concat(clstbl_list, ignore_index=True)
    regtbl = pd.concat(regtbl_list, ignore_index=True).reindex(columns=clstbl.columns)
    agg_dict = {c: "max" for c in regtbl.columns[2:]}
    regtbl = regtbl.reset_index().groupby(["pname_f", "pname_r"]).agg(agg_dict)
    clstbl = clstbl.reset_index().groupby(["pname_f", "pname_r"]).agg(agg_dict)

    coverage = (clstbl > 0.5).sum(axis=1).reset_index(name="coverage")
    coverage["coverage"] = coverage["coverage"] / n_for_coverage

    active_mask = clstbl > 0.5
    active_reg = regtbl.where(active_mask)
    activity = active_reg.sum(axis=1).reset_index(name="activity")
    activity["activity"] = activity["activity"] / n_for_coverage

    res = coverage.merge(activity, on=["pname_f", "pname_r"])
    res["activity"] = res["activity"].fillna(0)
    res["score"] = res["coverage"] * res["activity"]
    res = res.sort_values("activity", ascending=False)

    # Build per-target classifier table (for .eval.cl file)
    eval_cl = clstbl.fillna(0).round(3).reset_index()

    # Build per-target activity table (active_reg with 0 for inactive)
    eval_re = active_reg.fillna(0).reset_index()
    eval_re = eval_re.round(4)

    return res.round(4), eval_re, eval_full, eval_cl


def run(args):
    """Run quick-design command."""
    start_time = time.time()
    params = parse_params(args.param_file)
    primer_params = get_primer_params(params)

    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    min_dg = primer_params["min_dg"]

    cov_min = float(params.get("QUICK_COV_MIN", 1.0))
    act_min = float(params.get("QUICK_ACT_MIN", 1.0))
    min_pairs = int(params.get("QUICK_MIN_PAIRS", 5))

    # Route based on mode
    if args.msa and getattr(args, 'probe_mode', False):
        return _run_primers_first(args, params, primer_params, cov_min, act_min, min_pairs, start_time)
    elif args.msa:
        return _run_multi_region(args, params, primer_params, cov_min, act_min, min_pairs, start_time)
    else:
        return _run_single_window(args, params, primer_params, cov_min, act_min, min_pairs, start_time)


############################################
# Probe-first quick design
############################################


def _find_probe_regions_msa(aligned_seqs, probe_params, gap_threshold=0.50,
                            fallback_top_n=5, fallback_min_coverage=0.70):
    """Scan gap-filtered MSA for conserved probe regions.

    Uses conventional conservation (most-common-base fraction) to identify
    regions where variable columns <= max_mismatches within a probe_len_max
    window.

    If no strictly conserved regions are found, falls back to a per-sequence
    scoring approach: for each window, counts how many sequences match the
    consensus with <= max_mismatches, and returns the top windows by coverage.

    Returns:
        List of (msa_start, msa_end, n_variable) tuples in MSA coordinates.
        For fallback regions, n_variable is set to max_mismatches.
    """

    seq_len = len(aligned_seqs[0])
    n_seqs = len(aligned_seqs)
    probe_len_max = probe_params["len_max"]
    max_var = probe_params["max_mismatches"]
    cons_thresh = probe_params["conservation_threshold"]

    # Per-column: compute conservation and identify valid (non-gap-dominant) columns
    is_variable = np.zeros(seq_len, dtype=bool)
    is_valid = np.zeros(seq_len, dtype=bool)

    _msa = np.array([list(s.upper()) for s in aligned_seqs], dtype='U1')
    _is_gap = (_msa == '-') | (_msa == 'N')
    _gap_counts = _is_gap.sum(axis=0)

    _base_A = (_msa == 'A').sum(axis=0)
    _base_T = (_msa == 'T').sum(axis=0)
    _base_G = (_msa == 'G').sum(axis=0)
    _base_C = (_msa == 'C').sum(axis=0)
    _n_bases = _base_A + _base_T + _base_G + _base_C

    # Consensus base per column (for fallback scoring)
    consensus_arr = np.full(seq_len, 'N', dtype='U1')
    for col in range(seq_len):
        gap_frac = _gap_counts[col] / n_seqs
        if gap_frac > gap_threshold or _n_bases[col] == 0:
            continue
        is_valid[col] = True
        counts = [int(_base_A[col]), int(_base_T[col]),
                  int(_base_G[col]), int(_base_C[col])]
        most_common = max(counts)
        bases = ['A', 'T', 'G', 'C']
        consensus_arr[col] = bases[counts.index(most_common)]
        conservation = most_common / int(_n_bases[col])
        if conservation < cons_thresh:
            is_variable[col] = True

    if seq_len < probe_len_max:
        del _msa, _is_gap, _gap_counts, _base_A, _base_T, _base_G, _base_C, _n_bases
        return []

    cum_valid = np.concatenate([[0], np.cumsum(is_valid.astype(np.int32))])
    cum_var = np.concatenate([[0], np.cumsum(is_variable.astype(np.int32))])

    # --- Strict pass: all columns valid, variable count <= max_var ---
    conserved_starts = []
    for start in range(seq_len - probe_len_max + 1):
        end = start + probe_len_max
        n_valid = cum_valid[end] - cum_valid[start]
        if n_valid < probe_len_max:
            continue
        n_var = cum_var[end] - cum_var[start]
        if n_var <= max_var:
            conserved_starts.append(start)

    if conserved_starts:
        del _msa, _is_gap, _gap_counts, _base_A, _base_T, _base_G, _base_C, _n_bases
        return _merge_windows(conserved_starts, probe_len_max, cum_var)

    # --- Fallback: score each valid window by per-sequence coverage ---
    print("  No strictly conserved probe regions found. "
          "Searching for best-available regions...")

    # Find all valid windows (no gap-dominant columns)
    valid_starts = []
    for start in range(seq_len - probe_len_max + 1):
        end = start + probe_len_max
        n_valid = cum_valid[end] - cum_valid[start]
        if n_valid < probe_len_max:
            continue
        valid_starts.append(start)

    if not valid_starts:
        del _msa, _is_gap, _gap_counts, _base_A, _base_T, _base_G, _base_C, _n_bases
        return []

    # For each valid window, count sequences matching consensus with <= max_var mismatches
    # Build per-column mismatch matrix: 1 where seq differs from consensus
    mismatch = (_msa != consensus_arr[np.newaxis, :]).astype(np.int32)
    # Gaps don't count as mismatches (they're filtered by valid columns)
    mismatch[_is_gap] = 0
    cum_mm = np.concatenate([np.zeros((n_seqs, 1), dtype=np.int32),
                             np.cumsum(mismatch, axis=1)], axis=1)

    del _msa, _is_gap, _gap_counts, _base_A, _base_T, _base_G, _base_C, _n_bases, mismatch

    scored = []
    for start in valid_starts:
        end = start + probe_len_max
        # Per-sequence mismatch count in this window
        mm_counts = cum_mm[:, end] - cum_mm[:, start]
        n_matching = int(np.sum(mm_counts <= max_var))
        coverage = n_matching / n_seqs
        if coverage >= fallback_min_coverage:
            scored.append((start, coverage))

    if not scored:
        return []

    # Sort by coverage descending, pick top N non-overlapping windows
    scored.sort(key=lambda x: -x[1])
    selected = []
    for start, cov in scored:
        end = start + probe_len_max
        # Check no overlap with already selected
        if any(not (end <= s or start >= e) for s, e, _, _ in selected):
            continue
        selected.append((start, end, max_var, cov))
        if len(selected) >= fallback_top_n:
            break

    if not selected:
        return []

    # Sort by position
    selected.sort()
    regions = [(s, e, nv) for s, e, nv, _ in selected]

    for s, e, nv, cov in selected:
        print(f"    Fallback region: MSA cols {s}-{e} "
              f"({e-s} bp, covers {cov:.1%} of sequences with ≤{max_var} mismatches)")

    return regions


def _merge_windows(conserved_starts, probe_len_max, cum_var):
    """Merge overlapping conserved windows into contiguous regions."""
    regions = []
    reg_start = conserved_starts[0]
    reg_end = conserved_starts[0] + probe_len_max

    for s in conserved_starts[1:]:
        if s <= reg_end:
            reg_end = s + probe_len_max
        else:
            n_var = int(cum_var[reg_end] - cum_var[reg_start])
            regions.append((reg_start, reg_end, n_var))
            reg_start = s
            reg_end = s + probe_len_max

    n_var = int(cum_var[reg_end] - cum_var[reg_start])
    regions.append((reg_start, reg_end, n_var))

    return regions


def _score_probes_by_coverage(probes, aligned_seqs, max_mismatches,
                              max_indels=0, seq_ids=None):
    """Score each probe by fraction of sequences matching within thresholds.

    Uses wobble-aware mismatch counting: G-T and A-G pairs contribute 0.20
    instead of 1.0, and gap positions are counted as indels.

    Adds 'coverage' field to each probe dict.
    When seq_ids is provided, also adds 'covered_seq_ids' (set of matching IDs).
    """
    n_seqs = len(aligned_seqs)

    for probe in probes:
        probe_seq = probe['seq']
        # Use the strand the probe was designed for (thermodynamics validated)
        if probe['strand'] == 'rev':
            sense_seq = reverse_complement_dna(probe_seq)
        else:
            sense_seq = probe_seq

        msa_cols = probe.get('msa_cols')
        start = probe['msa_start']
        end = probe['msa_end']

        n_matching = 0
        covered = set()
        for idx, seq_row in enumerate(aligned_seqs):
            row_upper = seq_row.upper()
            if msa_cols is not None:
                mm, indels = wobble_mismatch_count_cols(
                    sense_seq, row_upper, msa_cols)
            else:
                mm, indels = wobble_mismatch_count_gapped(
                    sense_seq, row_upper, start, end)
            if mm <= max_mismatches and indels <= max_indels:
                n_matching += 1
                if seq_ids is not None:
                    covered.add(seq_ids[idx])

        probe['coverage'] = n_matching / n_seqs
        if seq_ids is not None:
            probe['covered_seq_ids'] = covered


def _select_top_probes(probes, top_n=5, min_distance=30):
    """Select top N probes at distinct positions, ranked by coverage.

    Ensures selected probes are at least min_distance apart in MSA coordinates.
    """
    # Sort by coverage descending, then by Tm closer to midpoint as tiebreaker
    sorted_probes = sorted(probes, key=lambda p: (-p['coverage'], -p['tm']))

    selected = []
    for probe in sorted_probes:
        mid = (probe['msa_start'] + probe['msa_end']) / 2
        if any(abs(mid - (s['msa_start'] + s['msa_end']) / 2) < min_distance
               for s in selected):
            continue
        selected.append(probe)
        if len(selected) >= top_n:
            break

    return selected


def _write_probe_mapping_csv(probes, aligned_seqs, seq_ids, out_path,
                              max_mismatches, max_indels=0):
    """Write probe mapping CSV directly from MSA data.

    Uses wobble-aware mismatch counting (G-T/A-G = 0.20 penalty) and
    correctly records gap positions as indels.
    Columns: probe_name, probe_seq, target_id, start_pos, orientation, mismatches, indels
    """
    mappings = []
    for probe in probes:
        start = probe['msa_start']
        end = probe['msa_end']
        msa_cols = probe.get('msa_cols')
        probe_seq = probe['seq']
        strand = probe.get('strand', 'fwd')
        orientation = '+' if strand == 'fwd' else '-'

        if strand == 'rev':
            sense_seq = reverse_complement_dna(probe_seq)
        else:
            sense_seq = probe_seq

        for i, seq_id in enumerate(seq_ids):
            seq_row = aligned_seqs[i].upper()
            if msa_cols is not None:
                mm, indels = wobble_mismatch_count_cols(
                    sense_seq, seq_row, msa_cols)
            else:
                mm, indels = wobble_mismatch_count_gapped(
                    sense_seq, seq_row, start, end)

            if mm > max_mismatches or indels > max_indels:
                continue

            ungapped_pos = sum(1 for c in seq_row[:start] if c != '-')

            mappings.append({
                'probe_name': probe['name'],
                'probe_seq': probe_seq,
                'target_id': seq_id,
                'start_pos': ungapped_pos,
                'orientation': orientation,
                'mismatches': round(mm, 2),
                'indels': indels,
            })

    df = pd.DataFrame(mappings)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"  Saved probe mapping: {out_path} ({len(df)} mappings)")


def _generate_probes_from_msa(aligned_seqs, regions, probe_params):
    """Generate probe candidates from conserved MSA regions.

    Builds consensus at each region, tiles probes of len_min to len_max on
    both strands, filters by Tm/GC/homopolymer/5'G/self-dimer.

    Returns:
        probes: list of dicts {name, seq, region_idx, msa_start, msa_end, tm, gc, dg, strand}
        probe_features: dict mapping probe_name -> {pseq, len, Tm, GC, dG}
    """
    from collections import Counter

    len_min = probe_params["len_min"]
    len_max = probe_params["len_max"]
    min_tm = probe_params["min_tm"]
    max_tm = probe_params["max_tm"]
    max_gc = probe_params["max_gc"] / 100.0
    homopolymer_max = probe_params["homopolymer_max"]
    avoid_5g = probe_params["avoid_5prime_g"]
    min_dg = probe_params["min_dg"]

    # Trim gap-frequent columns from MSA for probe consensus
    n_seqs = len(aligned_seqs)
    seq_len = len(aligned_seqs[0])
    _msa = np.array([list(s.upper()) for s in aligned_seqs], dtype='U1')
    _gap_freq = (_msa == '-').sum(axis=0) / n_seqs

    gap_threshold = 0.5
    kept_probe_cols = np.where(_gap_freq <= gap_threshold)[0]
    n_trimmed = seq_len - len(kept_probe_cols)
    if n_trimmed > 0:
        print(f"  Trimmed {n_trimmed}/{seq_len} gap-frequent columns "
              f"(>{gap_threshold:.0%} gaps) for probe consensus")
    del _msa

    # Precompute base frequency cache for kept columns (used by generate_probe_variants)
    freq_cache = {}
    for msa_col in kept_probe_cols:
        freq_cache[int(msa_col)] = _column_sense_freqs(aligned_seqs, int(msa_col))

    n_probe_variants = int(probe_params.get("n_variants", 5))

    candidates = []
    for reg_idx, (reg_start, reg_end, _) in enumerate(regions):
        # Convert MSA region boundaries to trimmed coordinates
        trim_start = np.searchsorted(kept_probe_cols, reg_start, side='left')
        trim_end = np.searchsorted(kept_probe_cols, reg_end, side='left')
        trim_region_len = trim_end - trim_start

        for probe_len in range(len_min, len_max + 1):
            for offset in range(trim_region_len - probe_len + 1):
                # MSA column indices for this probe window
                msa_cols = [int(kept_probe_cols[trim_start + offset + j])
                            for j in range(probe_len)]

                # Check for N in consensus at these positions
                if any(not freq_cache.get(c) for c in msa_cols):
                    continue

                msa_start = msa_cols[0]
                msa_end = msa_cols[-1] + 1

                for strand in ('fwd', 'rev'):
                    # Generate wobble-diversified variants
                    variant_seqs = generate_probe_variants(
                        aligned_seqs, msa_cols, strand,
                        n_variants=n_probe_variants, _freq_cache=freq_cache)

                    for seq in variant_seqs:
                        # Fast filters
                        if avoid_5g and seq[0] == 'G':
                            continue
                        if has_homopolymer(seq, homopolymer_max):
                            continue
                        tm = get_tm(seq)
                        if tm < min_tm or tm > max_tm:
                            continue
                        gc = gc_fraction(seq)
                        if gc > max_gc:
                            continue

                        candidates.append({
                            'seq': seq,
                            'region_idx': reg_idx,
                            'msa_start': msa_start,
                            'msa_end': msa_end,
                            'msa_cols': msa_cols,
                            'tm': round(tm, 2),
                            'gc': round(gc * 100, 1),
                            'strand': strand,
                        })

    if not candidates:
        return [], {}

    # Batch self-dimer dG check
    dimer_pairs = [(c['seq'], c['seq']) for c in candidates]
    dg_values = compute_batch_dimer_dg(dimer_pairs)

    probes = []
    probe_features = {}
    seen_seqs = set()
    for i, cand in enumerate(candidates):
        dg = dg_values[i] if dg_values[i] is not None else 0.0
        if dg < min_dg:
            continue
        cand['dg'] = round(dg, 2)
        # Deduplicate by sequence
        if cand['seq'] in seen_seqs:
            continue
        seen_seqs.add(cand['seq'])
        name = f"probe_r{cand['region_idx']}_{len(probes)+1}"
        cand['name'] = name
        probes.append(cand)
        probe_features[name] = {
            'pseq': cand['seq'],
            'len': len(cand['seq']),
            'Tm': cand['tm'],
            'GC': cand['gc'],
            'dG': cand['dg'],
        }

    return probes, probe_features


def _find_flanking_primer_positions(selected_probes, identity, kept_cols,
                                     primer_len, min_amp_len, max_amp_len,
                                     buffer, min_primer_score=0.70,
                                     min_gc=0.30, max_gc=0.65,
                                     top_per_probe=50):
    """Find valid primer positions flanking each selected probe.

    For each individual probe (~25bp), searches for fwd positions upstream
    and rev positions downstream, constrained by amplicon length and buffer.
    Scores by wobble-aware primer conservation.

    Args:
        selected_probes: list of probe dicts with 'msa_start', 'msa_end', 'name'

    Returns:
        List of (fwd_start, rev_start, score, probe_idx) tuples,
        where probe_idx indexes into selected_probes.
    """
    seq_len = len(identity)

    # Pre-compute valid primer positions (contiguous kept columns of primer_len)
    cum_kept = np.concatenate([[0], np.cumsum(kept_cols.astype(np.int32))])
    max_start = seq_len - primer_len + 1
    if max_start <= 0:
        return []
    kept_in_window = cum_kept[primer_len:max_start + primer_len] - cum_kept[:max_start]
    valid_pos = (kept_in_window == primer_len)

    # Mean conservation per primer position
    cum_id = np.concatenate([[0.0], np.cumsum(identity)])
    primer_scores = np.zeros(max_start)
    primer_scores[:] = (cum_id[primer_len:max_start + primer_len] - cum_id[:max_start]) / primer_len
    primer_scores[~valid_pos] = 0.0

    all_positions = []
    for probe_idx, probe in enumerate(selected_probes):
        probe_start = probe['msa_start']
        probe_end = probe['msa_end']

        # Fwd zone: fwd_start + primer_len + buffer <= probe_start
        fwd_max = probe_start - primer_len - buffer
        # Rev zone: rev_start >= probe_end + buffer
        rev_min = probe_end + buffer

        if fwd_max < 0 or rev_min >= max_start:
            continue

        # Collect valid fwd positions
        fwd_candidates = []
        for fwd in range(max(0, fwd_max - max_amp_len), fwd_max + 1):
            if fwd < 0 or fwd >= max_start:
                continue
            if primer_scores[fwd] < min_primer_score:
                continue
            fwd_candidates.append((fwd, primer_scores[fwd]))

        if not fwd_candidates:
            continue

        # Sort by score descending, then by proximity to probe (closest first)
        fwd_candidates.sort(key=lambda x: (-x[1], -x[0]))
        fwd_candidates = fwd_candidates[:top_per_probe]

        for fwd, fwd_score in fwd_candidates:
            # Valid rev range for this fwd
            amp_rev_lo = max(rev_min, fwd + min_amp_len - primer_len)
            amp_rev_hi = min(max_start, fwd + max_amp_len - primer_len + 1)
            if amp_rev_lo >= amp_rev_hi:
                continue

            # Find top rev positions in range (not just the best)
            rev_slice = primer_scores[amp_rev_lo:amp_rev_hi].copy()
            rev_valid = valid_pos[amp_rev_lo:amp_rev_hi]
            rev_score_mask = rev_slice >= min_primer_score
            combined = rev_valid & rev_score_mask
            if not combined.any():
                continue
            rev_slice[~combined] = -1.0
            # Take top 3 rev positions per fwd
            top_rev_idxs = np.argsort(rev_slice)[::-1][:3]
            for ri in top_rev_idxs:
                if rev_slice[ri] < 0:
                    break
                rev = amp_rev_lo + int(ri)
                combined_score = fwd_score + primer_scores[rev]
                all_positions.append((fwd, rev, combined_score, probe_idx))

    # Sort by score descending
    all_positions.sort(key=lambda x: -x[2])
    return all_positions


def _check_probe_primer_dimers(results_df, probes, name_to_seq, min_dg=-6.0):
    """Check probe-primer cross-dimer dG for passing pairs.

    For each primer pair in results, checks dG of each probe against
    both fwd and rev primers. Keeps only probes where both dG > min_dg.

    Returns:
        dict mapping (pname_f, pname_r) -> [compatible probe dicts]
    """
    if results_df.empty or not probes:
        return {}

    # Build all dimer check pairs
    check_list = []  # (pair_key, probe, 'fwd'/'rev')
    for _, row in results_df.iterrows():
        pf, pr = row['pname_f'], row['pname_r']
        fwd_seq = name_to_seq.get(pf)
        rev_seq = name_to_seq.get(pr)
        if not fwd_seq or not rev_seq:
            continue
        for probe in probes:
            check_list.append(((pf, pr), probe, fwd_seq, rev_seq))

    if not check_list:
        return {}

    # Build pairs for batch dG: (probe_seq, fwd_seq), (probe_seq, rev_seq)
    dimer_pairs = []
    for pair_key, probe, fwd_seq, rev_seq in check_list:
        dimer_pairs.append((probe['seq'], fwd_seq))
        dimer_pairs.append((probe['seq'], rev_seq))

    dg_values = compute_batch_dimer_dg(dimer_pairs)

    # Process results
    assignments = {}
    for i, (pair_key, probe, _, _) in enumerate(check_list):
        dg_fwd = dg_values[2 * i] if dg_values[2 * i] is not None else 0.0
        dg_rev = dg_values[2 * i + 1] if dg_values[2 * i + 1] is not None else 0.0
        if dg_fwd > min_dg and dg_rev > min_dg:
            assignments.setdefault(pair_key, []).append(probe)

    return assignments


def _run_probe_first(args, params, primer_params, cov_min, act_min, min_pairs, start_time):
    """Probe-first quick mode: find conserved probe regions, then flanking primers."""
    from Bio import AlignIO

    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    min_dg = primer_params["min_dg"]
    primer_len = primer_params["min_pri_len"]
    min_gc = primer_params.get("min_gc", 30) / 100.0 if primer_params.get("min_gc") else 0.30
    max_gc_frac = primer_params["max_gc"] / 100.0

    top_pairs = int(params.get("QUICK_TOP_PAIRS", TOP_PAIRS))
    n_variants = int(params.get("QUICK_N_VARIANTS", N_VARIANTS))

    probe_params = get_probe_params(params)
    buffer = int(params.get("PROBE_AMPLICON_BUFFER", 20))

    print(f"Quick mode (probe-first): {args.name}")
    print(f"  Thresholds: coverage >= {cov_min}, activity >= {act_min}")
    print(f"  Probe: len {probe_params['len_min']}-{probe_params['len_max']}, "
          f"Tm {probe_params['min_tm']}-{probe_params['max_tm']}, "
          f"conservation >= {probe_params['conservation_threshold']}")

    # Step 1: Load MSA
    print("Step 1: Loading MSA...")
    aln = AlignIO.read(args.msa, "fasta-pearson")
    raw_seqs = [str(r.seq) for r in aln]
    all_seq_ids = [r.id for r in aln]
    n_raw_cols = len(raw_seqs[0])

    # Pre-filter gap-dominant columns
    gap_thresh = 0.50
    n_seqs_raw = len(raw_seqs)
    _msa_arr = np.array([list(s) for s in raw_seqs], dtype='U1')
    _gap_counts = np.sum((_msa_arr == '-') | (_msa_arr == 'N') | (_msa_arr == 'n'), axis=0)
    keep_mask = (_gap_counts / n_seqs_raw) <= gap_thresh
    kept_indices = np.where(keep_mask)[0]
    _kept_arr = _msa_arr[:, kept_indices]
    aligned_seqs = [''.join(row) for row in _kept_arr]
    del _msa_arr, _kept_arr
    n_kept = int(keep_mask.sum())
    print(f"  MSA: {len(aligned_seqs)} seqs, {n_raw_cols} columns → "
          f"{n_kept} after removing {n_raw_cols - n_kept} gap-dominant columns")

    # Step 2: Find conserved probe regions
    # Cap region width so flanking primers fit within amplicon constraints
    print("Step 2: Scanning for conserved probe regions...")
    probe_regions = _find_probe_regions_msa(aligned_seqs, probe_params)
    if not probe_regions:
        print("WARNING: No probe regions found (even with relaxed fallback). "
              "Falling back to primer-only mode.")
        print("  Suggestions: increase PROBE_MAX_MISMATCHES or lower fallback_min_coverage")
        if args.probe_fa:
            Path(args.probe_fa).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_fa).touch()
        if args.probe_feat:
            Path(args.probe_feat).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_feat).touch()
        return _run_multi_region(args, params, primer_params, cov_min, act_min, min_pairs, start_time)

    print(f"  Found {len(probe_regions)} conserved probe regions:")
    for i, (s, e, nv) in enumerate(probe_regions):
        print(f"    Region {i+1}: MSA cols {s}-{e} ({e-s} bp, {nv} variable sites)")

    # Step 3: Generate probe candidates across all regions
    print("Step 3: Generating probe candidates...")
    probes, probe_features = _generate_probes_from_msa(aligned_seqs, probe_regions, probe_params)
    if not probes:
        print("WARNING: No probes passed filters. Falling back to primer-only mode.")
        if args.probe_fa:
            Path(args.probe_fa).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_fa).touch()
        if args.probe_feat:
            Path(args.probe_feat).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_feat).touch()
        return _run_multi_region(args, params, primer_params, cov_min, act_min, min_pairs, start_time)

    print(f"  Generated {len(probes)} probe candidates")

    # Step 3b: Score probes by per-sequence conservation and select top N
    max_mm = probe_params["max_mismatches"]
    max_indel = probe_params["max_indels"]
    top_n_probes = int(params.get("PROBE_TOP_N", 50))
    print(f"  Scoring probes by sequence coverage (≤{max_mm} mismatches, ≤{max_indel} indels)...")
    _score_probes_by_coverage(probes, aligned_seqs, max_mm, max_indel,
                             seq_ids=all_seq_ids)

    selected_probes = _select_top_probes(probes, top_n=top_n_probes, min_distance=30)
    print(f"  Selected {len(selected_probes)} probes at distinct positions:")
    for i, p in enumerate(selected_probes):
        print(f"    Probe {i+1}: {p['name']} ({p['seq'][:20]}...) "
              f"MSA {p['msa_start']}-{p['msa_end']}, "
              f"coverage={p['coverage']:.1%}, Tm={p['tm']}")

    # Update probe_features to only include selected probes
    selected_names = {p['name'] for p in selected_probes}
    probe_features = {k: v for k, v in probe_features.items() if k in selected_names}

    # Build probes_by_probe_idx for dimer checking later
    probes_by_probe_idx = {i: [p] for i, p in enumerate(selected_probes)}

    # Write probe mapping CSV directly from MSA (replaces build_index + align_probes + parse_probe_mapping)
    if getattr(args, 'probe_csv', None):
        _write_probe_mapping_csv(
            selected_probes, aligned_seqs, all_seq_ids,
            args.probe_csv, max_mm, max_indel,
        )

    # Step 4: Compute primer conservation and find flanking positions
    print("Step 4: Scoring flanking primer positions...")
    scored, consensus_aligned, kept_cols = score_amplicon_positions(
        aligned_seqs, primer_len, min_amp_len, max_amp_len,
        min_gc=min_gc, max_gc=max_gc_frac,
        top_n=0,
    )

    # Compute wobble-aware identity per column
    seq_len = len(aligned_seqs[0])
    _msa = np.array([list(s.upper()) for s in aligned_seqs], dtype='U1')
    _base_A = (_msa == 'A').sum(axis=0)
    _base_T = (_msa == 'T').sum(axis=0)
    _base_G = (_msa == 'G').sum(axis=0)
    _base_C = (_msa == 'C').sum(axis=0)
    _n_bases = _base_A + _base_T + _base_G + _base_C
    _gap_counts_col = ((_msa == '-') | (_msa == 'N')).sum(axis=0)
    n_seqs = len(aligned_seqs)

    identity = np.zeros(seq_len, dtype=np.float64)
    for col in range(seq_len):
        if _n_bases[col] > 0 and (_gap_counts_col[col] / n_seqs) <= gap_thresh:
            counts = {'A': int(_base_A[col]), 'T': int(_base_T[col]),
                      'G': int(_base_G[col]), 'C': int(_base_C[col])}
            sense_freqs = {b: c / n_seqs for b, c in counts.items() if c > 0}
            anti_freqs = {_COMPLEMENT.get(b, b): f for b, f in sense_freqs.items()}
            _, fwd_score = _best_primer_base(anti_freqs)
            _, rev_score = _best_primer_base(sense_freqs)
            identity[col] = max(fwd_score, rev_score)
    del _msa, _base_A, _base_T, _base_G, _base_C, _n_bases, _gap_counts_col

    positions_with_probe = _find_flanking_primer_positions(
        selected_probes, identity, kept_cols,
        primer_len, min_amp_len, max_amp_len, buffer,
        min_gc=min_gc, max_gc=max_gc_frac,
    )
    if not positions_with_probe:
        print("ERROR: No valid flanking primer positions found.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    # Map probe_idx for each primer pair (a position can serve multiple probes)
    primer_probe_map = {}  # (fwd, rev) -> set of probe_idxs
    positions = [(f, r, s) for f, r, s, _ in positions_with_probe]
    for f, r, s, probe_idx in positions_with_probe:
        primer_probe_map.setdefault((f, r), set()).add(probe_idx)

    print(f"  Found {len(positions)} flanking primer positions across "
          f"{len(set(pi for _, _, _, pi in positions_with_probe))} probes")

    # Step 5: Generate primer variants
    print("Step 5: Generating wobble-optimized primer variants...")
    primers, features = extract_primers_at_positions(
        aligned_seqs, positions, primer_len,
        primer_params["min_tm"], primer_params["max_tm"],
        primer_params["max_gc"], min_dg,
        n_variants=n_variants,
        max_pri_len=primer_params["max_pri_len"],
    )

    if not primers:
        print("ERROR: No primers passed Tm/GC/dG filters.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    # Propagate probe_idxs to primers (set of all probes this position serves)
    for p in primers:
        key = (p['fwd_start_msa'], p['rev_start_msa'])
        p['probe_idxs'] = primer_probe_map.get(key, set())

    n_unique_pos = len(set(p['pos_idx'] for p in primers))
    print(f"  Passed filters: {len(primers)} pairs ({n_unique_pos} positions)")

    # Take top N by score
    primers.sort(key=lambda p: -p['score'])
    if len(primers) > top_pairs:
        primers = primers[:top_pairs]
    print(f"  Selected top {len(primers)} pairs")

    # Step 6: Group into batches
    print("Step 6: Grouping into batches...")
    batches = group_into_batches(primers, primer_len)
    for i, (batch, min_fwd, max_rev) in enumerate(batches):
        print(f"  Batch {i+1}: {len(batch)} pairs, MSA {min_fwd}-{max_rev}")

    # Step 7: Load ML models
    print("Step 7: Loading ML models...")
    scaler, classifier, regressor, device = load_models()
    torch.set_num_threads(args.threads)

    # Setup tmpdir
    if args.tmpdir:
        tmpdir = str(Path(args.tmpdir) / "_intermediate")
    else:
        tmpdir = str(Path(args.out).parent)
    os.makedirs(tmpdir, exist_ok=True)

    # Step 8: Per-batch evaluation loop
    all_results = []
    all_eval_re = []
    all_eval_full = []
    all_eval_cl = []
    name_to_seq = {}
    all_probe_assignments = {}
    pair_to_probe_idx = {}  # (pname_f, pname_r) -> set of probe_idxs

    try:
        for batch_idx, (batch_primers, min_fwd, max_rev) in enumerate(batches):
            batch_start = time.time()
            print(f"\nBatch {batch_idx+1}/{len(batches)}: {len(batch_primers)} pairs, MSA {min_fwd}-{max_rev}")

            # Extract bowtie2 reference
            margin = 20
            win_start = max(0, min_fwd - margin)
            win_end = max_rev + margin
            all_window = extract_window_sequences(aligned_seqs, all_seq_ids, win_start, win_end)
            if not all_window:
                print("  No sequences in region. Skipping.")
                continue

            win_target_fa = Path(tmpdir) / f"batch_{batch_idx}_targets.fa"
            with open(win_target_fa, "w") as f:
                for sid, seq in all_window:
                    f.write(f">{sid}\n{seq}\n")
            n_targets = len(all_window)

            # Build bowtie2 index
            index_prefix = str(Path(tmpdir) / f"batch_{batch_idx}_bt2idx")
            print("  Building bowtie2 index...", end=" ", flush=True)
            build_index(str(win_target_fa), index_prefix, threads=args.threads)
            print("done.")

            # Write batch FASTA and features
            fwd_seqs = [p['fwd_seq'] for p in batch_primers]
            rev_seqs = [p['rev_seq'] for p in batch_primers]
            fasta_path, feat_path, fwd_n2s, rev_n2s = _write_batch_fasta_and_features(
                fwd_seqs, rev_seqs, features, args.name, batch_idx, tmpdir
            )
            name_to_seq.update(fwd_n2s)
            name_to_seq.update(rev_n2s)

            # Map primer pair names to probe_idx
            for i, p in enumerate(batch_primers):
                pname_f = f"{args.name}_q{batch_idx}_{i+1}_f"
                pname_r = f"{args.name}_q{batch_idx}_{i+1}_r"
                pair_to_probe_idx[(pname_f, pname_r)] = p.get('probe_idxs', set())

            # Align
            print("  Aligning...", end=" ", flush=True)
            mapped_path = _align_batch(
                fasta_path, index_prefix, n_targets, args.threads, tmpdir, batch_idx
            )
            print("done.")

            if mapped_path.stat().st_size == 0:
                print("  No alignments found. Skipping.")
                continue

            # Coverage filter
            filtered_path = _apply_coverage_filter(mapped_path, n_targets, tmpdir, batch_idx, cov_frac=0.90)
            if filtered_path.stat().st_size == 0:
                print("  No primers pass coverage filter. Skipping.")
                continue

            # Prepare ML input
            print("  Preparing ML input...", end=" ", flush=True)
            input_path = Path(tmpdir) / f"batch_{batch_idx}.input"
            prep_args = argparse.Namespace(
                mapped=str(filtered_path),
                ml_input=str(input_path),
                reference=str(win_target_fa),
                param_file=args.param_file,
                reftype="on",
                pri_features=str(feat_path),
                prev="",
            )
            try:
                run_prepare_input(prep_args)
            except SystemExit:
                pass
            print("done.")

            if not input_path.exists() or input_path.stat().st_size == 0:
                print("  No valid pairs formed. Skipping.")
                continue

            # ML evaluate
            print("  Evaluating...", end=" ", flush=True)
            res, eval_re, eval_full, eval_cl = _run_evaluate(
                input_path, str(win_target_fa), scaler, classifier, regressor, device, args.threads,
                n_targets_global=len(aligned_seqs),
            )
            print("done.")

            if res.empty:
                print("  No results. Skipping.")
                continue

            # Probe-primer dimer check
            batch_probe_idxs = set()
            for p in batch_primers:
                batch_probe_idxs.update(p.get('probe_idxs', set()))
            batch_probes = []
            for pi in batch_probe_idxs:
                batch_probes.extend(probes_by_probe_idx.get(pi, []))

            if batch_probes:
                print(f"  Checking probe-primer dimers ({len(batch_probes)} probes)...", end=" ", flush=True)
                assignments = _check_probe_primer_dimers(res, batch_probes, name_to_seq, min_dg=min_dg)
                all_probe_assignments.update(assignments)
                n_with_probes = sum(1 for _, row in res.iterrows()
                                   if (row['pname_f'], row['pname_r']) in assignments)
                print(f"{n_with_probes}/{len(res)} pairs have compatible probes.")

            all_results.append(res)
            all_eval_re.append(eval_re)
            if not eval_full.empty:
                all_eval_full.append(eval_full)
            if not eval_cl.empty:
                all_eval_cl.append(eval_cl)
            batch_time = time.time() - batch_start

            passing = res[(res["coverage"] >= cov_min) & (res["activity"] >= act_min)]
            total_passing = sum(
                len(r[(r["coverage"] >= cov_min) & (r["activity"] >= act_min)])
                for r in all_results
            )
            top = res.iloc[0]
            print(f"  Top pair: coverage={top['coverage']:.3f}, activity={top['activity']:.3f}, score={top['score']:.3f}")
            print(f"  Pairs meeting thresholds: {len(passing)} (total: {total_passing}/{min_pairs})")
            print(f"  Batch time: {batch_time:.1f}s")

            if total_passing >= min_pairs:
                print(f"\n  SUCCESS! Found {total_passing} pairs meeting thresholds (target: {min_pairs}).")
                break

        # Write outputs
        probe_data = {
            'probes': probes,
            'features': probe_features,
            'assignments': all_probe_assignments,
            'pair_to_probe_idx': pair_to_probe_idx,
            'probes_by_probe_idx': probes_by_probe_idx,
        }
        _write_output(args, all_results, all_eval_re, name_to_seq,
                      cov_min, act_min, len(batches), tmpdir, min_dg=min_dg,
                      features=features, probe_data=probe_data,
                      all_eval_full=all_eval_full, all_eval_cl=all_eval_cl,
                      pair_to_probe_idx=pair_to_probe_idx)

    finally:
        print(f"  Temp files kept at: {tmpdir}")

    runtime = time.time() - start_time
    print(f"\nProbe-first quick mode completed in {runtime:.1f}s")


def _run_primers_first(args, params, primer_params, cov_min, act_min, min_pairs, start_time):
    """Primers-first quick mode: design primers globally, then find probes in amplicon regions."""
    from Bio import AlignIO

    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    min_dg = primer_params["min_dg"]
    primer_len = primer_params["min_pri_len"]
    min_gc = primer_params.get("min_gc", 30) / 100.0 if primer_params.get("min_gc") else 0.30
    max_gc_frac = primer_params["max_gc"] / 100.0

    top_pairs = int(params.get("QUICK_TOP_PAIRS", TOP_PAIRS))
    n_positions = int(params.get("QUICK_N_POSITIONS", N_POSITIONS))
    n_variants = int(params.get("QUICK_N_VARIANTS", N_VARIANTS))

    probe_params = get_probe_params(params)
    max_mm = probe_params["max_mismatches"]
    max_indel = probe_params["max_indels"]

    print(f"Quick mode (primers-first + probe): {args.name}")
    print(f"  Thresholds: coverage >= {cov_min}, activity >= {act_min}")
    print(f"  Positions: {n_positions}, variants/direction: {n_variants}")
    print(f"  Probe: len {probe_params['len_min']}-{probe_params['len_max']}, "
          f"Tm {probe_params['min_tm']}-{probe_params['max_tm']}, "
          f"max_mm {max_mm}, max_indels {max_indel}")

    # Step 1: Load MSA and score amplicon positions (identical to _run_multi_region)
    print("Step 1: Scoring amplicon positions across MSA...")
    aln = AlignIO.read(args.msa, "fasta-pearson")
    raw_seqs = [str(r.seq) for r in aln]
    all_seq_ids = [r.id for r in aln]
    n_raw_cols = len(raw_seqs[0])

    gap_thresh = 0.50
    n_seqs_raw = len(raw_seqs)
    _msa_arr = np.array([list(s) for s in raw_seqs], dtype='U1')
    _gap_counts = np.sum((_msa_arr == '-') | (_msa_arr == 'N') | (_msa_arr == 'n'), axis=0)
    keep_mask = (_gap_counts / n_seqs_raw) <= gap_thresh
    kept_indices = np.where(keep_mask)[0]
    _kept_arr = _msa_arr[:, kept_indices]
    aligned_seqs = [''.join(row) for row in _kept_arr]
    del _msa_arr, _kept_arr
    n_kept = int(keep_mask.sum())
    print(f"  MSA: {len(aligned_seqs)} seqs, {n_raw_cols} columns -> "
          f"{n_kept} after removing {n_raw_cols - n_kept} gap-dominant columns")

    scored, consensus_aligned, kept_cols = score_amplicon_positions(
        aligned_seqs, primer_len, min_amp_len, max_amp_len,
        min_gc=min_gc, max_gc=max_gc_frac,
        top_n=n_positions,
    )
    if not scored:
        print("ERROR: No valid amplicon positions found in MSA.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    print(f"  Top amplicon positions: {len(scored)}")

    # Step 2: Generate primer variants
    print("Step 2: Generating wobble-optimized primer variants...")
    primers, features = extract_primers_at_positions(
        aligned_seqs, scored, primer_len,
        primer_params["min_tm"], primer_params["max_tm"],
        primer_params["max_gc"], min_dg,
        n_variants=n_variants,
        max_pri_len=primer_params["max_pri_len"],
    )

    if not primers:
        print("ERROR: No primers passed Tm/GC/dG filters.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    n_unique_pos = len(set(p['pos_idx'] for p in primers))
    print(f"  Passed filters: {len(primers)} pairs ({n_unique_pos} positions)")

    primers.sort(key=lambda p: -p['score'])
    if len(primers) > top_pairs:
        primers = primers[:top_pairs]
    print(f"  Selected top {len(primers)} pairs")

    # Step 3: Group into batches
    print("Step 3: Grouping into batches...")
    batches = group_into_batches(primers, primer_len)
    for i, (batch, min_fwd, max_rev) in enumerate(batches):
        print(f"  Batch {i+1}: {len(batch)} pairs, MSA {min_fwd}-{max_rev}")

    # Step 4: Load ML models
    print("Step 4: Loading ML models...")
    scaler, classifier, regressor, device = load_models()
    torch.set_num_threads(args.threads)

    if args.tmpdir:
        tmpdir = str(Path(args.tmpdir) / "_intermediate")
    else:
        tmpdir = str(Path(args.out).parent)
    os.makedirs(tmpdir, exist_ok=True)

    # Step 5: Per-batch evaluation loop (same as _run_multi_region)
    all_results = []
    all_eval_re = []
    all_eval_full = []
    all_eval_cl = []
    name_to_seq = {}
    pname_to_msa_pos = {}  # pname_f -> (fwd_start_msa, rev_start_msa, fwd_len, rev_len)

    try:
        for batch_idx, (batch_primers, min_fwd, max_rev) in enumerate(batches):
            batch_start = time.time()
            print(f"\nBatch {batch_idx+1}/{len(batches)}: {len(batch_primers)} pairs, MSA {min_fwd}-{max_rev}")

            margin = 20
            win_start = max(0, min_fwd - margin)
            win_end = max_rev + margin
            all_window = extract_window_sequences(aligned_seqs, all_seq_ids, win_start, win_end)
            if not all_window:
                print("  No sequences in region. Skipping.")
                continue

            win_target_fa = Path(tmpdir) / f"batch_{batch_idx}_targets.fa"
            with open(win_target_fa, "w") as f:
                for sid, seq in all_window:
                    f.write(f">{sid}\n{seq}\n")
            n_targets = len(all_window)

            index_prefix = str(Path(tmpdir) / f"batch_{batch_idx}_bt2idx")
            print("  Building bowtie2 index...", end=" ", flush=True)
            build_index(str(win_target_fa), index_prefix, threads=args.threads)
            print("done.")

            fwd_seqs = [p['fwd_seq'] for p in batch_primers]
            rev_seqs = [p['rev_seq'] for p in batch_primers]
            fasta_path, feat_path, fwd_n2s, rev_n2s = _write_batch_fasta_and_features(
                fwd_seqs, rev_seqs, features, args.name, batch_idx, tmpdir
            )
            name_to_seq.update(fwd_n2s)
            name_to_seq.update(rev_n2s)

            # Track pname -> MSA position mapping for probe search
            for i, p in enumerate(batch_primers):
                pname_f = f"{args.name}_q{batch_idx}_{i+1}_f"
                pname_r = f"{args.name}_q{batch_idx}_{i+1}_r"
                pname_to_msa_pos[pname_f] = (
                    p['fwd_start_msa'], p['rev_start_msa'],
                    len(p['fwd_seq']), len(p['rev_seq']),
                )

            print("  Aligning...", end=" ", flush=True)
            mapped_path = _align_batch(
                fasta_path, index_prefix, n_targets, args.threads, tmpdir, batch_idx
            )
            print("done.")

            if mapped_path.stat().st_size == 0:
                print("  No alignments found. Skipping.")
                continue

            filtered_path = _apply_coverage_filter(mapped_path, n_targets, tmpdir, batch_idx, cov_frac=0.90)
            if filtered_path.stat().st_size == 0:
                print("  No primers pass coverage filter. Skipping.")
                continue

            print("  Preparing ML input...", end=" ", flush=True)
            input_path = Path(tmpdir) / f"batch_{batch_idx}.input"
            prep_args = argparse.Namespace(
                mapped=str(filtered_path),
                ml_input=str(input_path),
                reference=str(win_target_fa),
                param_file=args.param_file,
                reftype="on",
                pri_features=str(feat_path),
                prev="",
            )
            try:
                run_prepare_input(prep_args)
            except SystemExit:
                pass
            print("done.")

            if not input_path.exists() or input_path.stat().st_size == 0:
                print("  No valid pairs formed. Skipping.")
                continue

            print("  Evaluating...", end=" ", flush=True)
            res, eval_re, eval_full, eval_cl = _run_evaluate(
                input_path, str(win_target_fa), scaler, classifier, regressor, device, args.threads,
                n_targets_global=len(aligned_seqs),
            )
            print("done.")

            if res.empty:
                print("  No results. Skipping.")
                continue

            all_results.append(res)
            all_eval_re.append(eval_re)
            if not eval_full.empty:
                all_eval_full.append(eval_full)
            if not eval_cl.empty:
                all_eval_cl.append(eval_cl)
            batch_time = time.time() - batch_start

            passing = res[(res["coverage"] >= cov_min) & (res["activity"] >= act_min)]
            total_passing = sum(
                len(r[(r["coverage"] >= cov_min) & (r["activity"] >= act_min)])
                for r in all_results
            )
            top = res.iloc[0]
            print(f"  Top pair: coverage={top['coverage']:.3f}, activity={top['activity']:.3f}, score={top['score']:.3f}")
            print(f"  Pairs meeting thresholds: {len(passing)} (total: {total_passing}/{min_pairs})")
            print(f"  Batch time: {batch_time:.1f}s")

            if total_passing >= min_pairs:
                print(f"\n  SUCCESS! Found {total_passing} pairs meeting thresholds (target: {min_pairs}).")
                break

        # Step 6: Generate probes across the full MSA, then assign to primer pairs
        if all_results:
            combined = pd.concat(all_results, ignore_index=True)
            combined = combined.sort_values("activity", ascending=False).drop_duplicates(
                subset=["pname_f", "pname_r"]
            )
            combined = combined[combined["coverage"] > 0]
        else:
            combined = pd.DataFrame()

        if combined.empty:
            print("\nWARNING: No primer pairs found. Skipping probe search.")
            _write_output(args, all_results, all_eval_re, name_to_seq,
                          cov_min, act_min, len(batches), tmpdir, min_dg=min_dg,
                          features=features, all_eval_full=all_eval_full,
                          all_eval_cl=all_eval_cl)
        else:
            print(f"\nStep 6: Searching for probes across MSA...")

            # Search probes across the full MSA (single region covering everything)
            msa_len = len(aligned_seqs[0])
            full_region = [(0, msa_len, 0)]
            probes, probe_features = _generate_probes_from_msa(
                aligned_seqs, full_region, probe_params)
            print(f"  Generated {len(probes)} probe candidates (Tm/GC/homopolymer/dG filtered)")

            if not probes:
                print("  WARNING: No probes passed thermodynamic filters.")
                if args.probe_fa:
                    Path(args.probe_fa).parent.mkdir(parents=True, exist_ok=True)
                    Path(args.probe_fa).touch()
                if args.probe_feat:
                    Path(args.probe_feat).parent.mkdir(parents=True, exist_ok=True)
                    Path(args.probe_feat).touch()
                _write_output(args, all_results, all_eval_re, name_to_seq,
                              cov_min, act_min, len(batches), tmpdir, min_dg=min_dg,
                              features=features, all_eval_full=all_eval_full,
                              all_eval_cl=all_eval_cl)
            else:
                # Score probes by wobble-aware coverage (once for all probes)
                print(f"  Scoring probes by wobble-aware coverage (max_mm={max_mm}, max_indels={max_indel})...")
                _score_probes_by_coverage(probes, aligned_seqs, max_mm, max_indel,
                                         seq_ids=all_seq_ids)

                # Select top probes at distinct positions
                top_n_probes = int(params.get("PROBE_TOP_N", 50))
                selected_probes = _select_top_probes(probes, top_n=top_n_probes, min_distance=30)
                print(f"  Selected {len(selected_probes)} probes at distinct positions:")
                for i, p in enumerate(selected_probes[:10]):
                    print(f"    Probe {i+1}: {p['name']} ({p['seq'][:20]}...) "
                          f"MSA {p['msa_start']}-{p['msa_end']}, "
                          f"coverage={p['coverage']:.1%}, Tm={p['tm']}")

                selected_names = {p['name'] for p in selected_probes}
                probe_features = {k: v for k, v in probe_features.items()
                                  if k in selected_names}

                # Assign probes to primer pairs by amplicon containment
                probes_by_probe_idx = {i: [p] for i, p in enumerate(selected_probes)}
                pair_to_probe_idx = {}

                for _, row in combined.iterrows():
                    pname_f = row['pname_f']
                    if pname_f not in pname_to_msa_pos:
                        continue
                    fwd_start, rev_start, fwd_len, rev_len = pname_to_msa_pos[pname_f]
                    amp_start = fwd_start + fwd_len
                    amp_end = rev_start

                    pair_key = (pname_f, row['pname_r'])
                    probe_idxs = set()
                    for pi, p in enumerate(selected_probes):
                        if p['msa_start'] >= amp_start and p['msa_end'] <= amp_end:
                            probe_idxs.add(pi)
                    if probe_idxs:
                        pair_to_probe_idx[pair_key] = probe_idxs

                n_with_probes = sum(1 for _ in pair_to_probe_idx)
                print(f"  {n_with_probes}/{len(combined)} pairs have probes in amplicon")

                # Check probe-primer dimers
                print(f"  Checking probe-primer dimers...")
                all_probe_assignments = _check_probe_primer_dimers(
                    combined, selected_probes, name_to_seq, min_dg=min_dg)
                n_dimer_pass = sum(1 for pair_key in pair_to_probe_idx
                                   if pair_key in all_probe_assignments)
                print(f"  {n_dimer_pass}/{n_with_probes} pairs pass probe-primer dimer check")

                # Write probe mapping CSV
                if getattr(args, 'probe_csv', None):
                    _write_probe_mapping_csv(
                        selected_probes, aligned_seqs, all_seq_ids,
                        args.probe_csv, max_mm, max_indel,
                    )

                # Compute intersected coverage: primer targets ∩ probe targets
                # Primer targets come from ML eval.full; probe targets from MSA scoring
                n_total = len(all_seq_ids)
                if all_eval_full:
                    eval_full_df = pd.concat(all_eval_full, ignore_index=True)
                    if 'classifier' in eval_full_df.columns:
                        eval_full_pos = eval_full_df[eval_full_df['classifier'] >= 0.5]
                    else:
                        eval_full_pos = eval_full_df
                else:
                    eval_full_pos = pd.DataFrame()

                # Greedy probe assignment: walk pairs by score, assign best unclaimed probe
                claimed_probes = set()
                pair_assigned_probe = {}  # pair_key -> probe dict

                for _, row in combined.iterrows():
                    pair_key = (row['pname_f'], row['pname_r'])
                    # Must be in amplicon AND pass dimer check
                    position_probes = pair_to_probe_idx.get(pair_key, set())
                    dimer_probes = {p['name'] for p in all_probe_assignments.get(pair_key, [])}
                    valid_probes = []
                    for pi in position_probes:
                        p = selected_probes[pi]
                        if p['name'] in dimer_probes and p['name'] not in claimed_probes:
                            valid_probes.append(p)
                    if valid_probes:
                        # Pick probe with highest coverage
                        best = max(valid_probes, key=lambda p: p['coverage'])
                        pair_assigned_probe[pair_key] = best
                        claimed_probes.add(best['name'])

                # Compute per-pair coverage
                cov_target_map = {}
                cov_primer_only_map = {}
                assigned_probe_name_map = {}
                assigned_probe_seq_map = {}

                for idx, row in combined.iterrows():
                    pair_key = (row['pname_f'], row['pname_r'])

                    # Get primer targets from eval.full
                    primer_targets = set()
                    if not eval_full_pos.empty:
                        pair_rows = eval_full_pos[
                            (eval_full_pos['pname_f'] == row['pname_f']) &
                            (eval_full_pos['pname_r'] == row['pname_r'])
                        ]
                        for _, pr in pair_rows.iterrows():
                            targets = ast.literal_eval(str(pr['targets']))
                            primer_targets.update(targets)

                    cov_primer_only_map[idx] = f"{len(primer_targets)} / {n_total}"

                    probe = pair_assigned_probe.get(pair_key)
                    if probe:
                        probe_targets = probe.get('covered_seq_ids', set())
                        both = primer_targets & probe_targets
                        cov_target_map[idx] = f"{len(both)} / {n_total}"
                        assigned_probe_name_map[idx] = probe['name']
                        assigned_probe_seq_map[idx] = probe['seq']
                    else:
                        cov_target_map[idx] = f"{len(primer_targets)} / {n_total}"
                        assigned_probe_name_map[idx] = ''
                        assigned_probe_seq_map[idx] = ''

                combined['cov_target'] = combined.index.map(cov_target_map)
                combined['cov_primer_only'] = combined.index.map(cov_primer_only_map)
                combined['probe_name'] = combined.index.map(assigned_probe_name_map)
                combined['probe_seq'] = combined.index.map(assigned_probe_seq_map)

                n_assigned = sum(1 for v in assigned_probe_name_map.values() if v)
                print(f"  Greedy probe assignment: {n_assigned}/{len(combined)} pairs got a probe")
                for idx, row in combined.iterrows():
                    if assigned_probe_name_map.get(idx):
                        print(f"    {row['pname_f']},{row['pname_r']} → "
                              f"{assigned_probe_name_map[idx]} "
                              f"(cov_target={cov_target_map[idx]}, "
                              f"cov_primer_only={cov_primer_only_map[idx]})")

                # Build probe_data for _write_output
                probe_data = {
                    'probes': selected_probes,
                    'features': probe_features,
                    'assignments': all_probe_assignments,
                    'pair_to_probe_idx': pair_to_probe_idx,
                    'probes_by_probe_idx': probes_by_probe_idx,
                }

                _write_output(args, all_results, all_eval_re, name_to_seq,
                              cov_min, act_min, len(batches), tmpdir, min_dg=min_dg,
                              features=features, probe_data=probe_data,
                              all_eval_full=all_eval_full, all_eval_cl=all_eval_cl,
                              pair_to_probe_idx=pair_to_probe_idx,
                              combined_override=combined)

    finally:
        print(f"  Temp files kept at: {tmpdir}")

    runtime = time.time() - start_time
    print(f"\nPrimers-first quick mode completed in {runtime:.1f}s")


def _run_multi_region(args, params, primer_params, cov_min, act_min, min_pairs, start_time):
    """Position-first quick mode: score positions globally, extract primers, batch evaluate."""
    from Bio import AlignIO

    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    min_dg = primer_params["min_dg"]
    primer_len = primer_params["min_pri_len"]
    min_gc = primer_params.get("min_gc", 30) / 100.0 if primer_params.get("min_gc") else 0.30
    max_gc_frac = primer_params["max_gc"] / 100.0

    top_pairs = int(params.get("QUICK_TOP_PAIRS", TOP_PAIRS))

    n_positions = int(params.get("QUICK_N_POSITIONS", N_POSITIONS))
    n_variants = int(params.get("QUICK_N_VARIANTS", N_VARIANTS))

    print(f"Quick mode (wobble-optimized): {args.name}")
    print(f"  Thresholds: coverage >= {cov_min}, activity >= {act_min}")
    print(f"  Positions: {n_positions}, variants/direction: {n_variants}")
    print(f"  Max pairs: {n_positions} × {n_variants}² = {n_positions * n_variants**2}")

    # Step 1: Load MSA, strip gap-dominant columns, and score amplicon positions
    print("Step 1: Scoring amplicon positions across MSA...")
    aln = AlignIO.read(args.msa, "fasta-pearson")
    raw_seqs = [str(r.seq) for r in aln]
    all_seq_ids = [r.id for r in aln]
    n_raw_cols = len(raw_seqs[0])

    # Pre-filter: remove columns where gap fraction > 0.50
    gap_thresh = 0.50
    n_seqs_raw = len(raw_seqs)
    # Vectorize: convert MSA to numpy array of bytes for fast column ops
    _msa_arr = np.array([list(s) for s in raw_seqs], dtype='U1')
    _gap_counts = np.sum((_msa_arr == '-') | (_msa_arr == 'N') | (_msa_arr == 'n'), axis=0)
    keep_mask = (_gap_counts / n_seqs_raw) <= gap_thresh
    kept_indices = np.where(keep_mask)[0]
    _kept_arr = _msa_arr[:, kept_indices]
    aligned_seqs = [''.join(row) for row in _kept_arr]
    del _msa_arr, _kept_arr  # free memory
    n_kept = int(keep_mask.sum())
    print(f"  MSA: {len(aligned_seqs)} seqs, {n_raw_cols} columns → "
          f"{n_kept} after removing {n_raw_cols - n_kept} gap-dominant columns")

    scored, consensus_aligned, kept_cols = score_amplicon_positions(
        aligned_seqs, primer_len, min_amp_len, max_amp_len,
        min_gc=min_gc, max_gc=max_gc_frac,
        top_n=n_positions,
    )
    if not scored:
        print("ERROR: No valid amplicon positions found in MSA.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    print(f"  Top amplicon positions: {len(scored)}")
    for i, (fwd, rev, score) in enumerate(scored[:5]):
        print(f"    #{i+1}: fwd={fwd} rev={rev} (score={score:.4f})")

    # Step 2: Generate wobble-optimized primer variants
    print("Step 2: Generating wobble-optimized primer variants...")
    primers, features = extract_primers_at_positions(
        aligned_seqs, scored, primer_len,
        primer_params["min_tm"], primer_params["max_tm"],
        primer_params["max_gc"], min_dg,
        n_variants=n_variants,
        max_pri_len=primer_params["max_pri_len"],
    )

    if not primers:
        print("ERROR: No primers passed Tm/GC/dG filters at any position.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    n_unique_pos = len(set(p['pos_idx'] for p in primers))
    print(f"  Passed Tm/GC/dG: {len(primers)} pairs ({n_unique_pos} positions)")

    # Take top N pairs by score
    primers.sort(key=lambda p: -p['score'])
    if len(primers) > top_pairs:
        primers = primers[:top_pairs]
    n_unique_pos = len(set(p['pos_idx'] for p in primers))
    print(f"  Selected top {len(primers)} pairs ({n_unique_pos} positions)")

    # Step 4: Group into batches
    print("Step 4: Grouping into batches...")
    batches = group_into_batches(primers, primer_len)
    for i, (batch, min_fwd, max_rev) in enumerate(batches):
        print(f"  Batch {i+1}: {len(batch)} pairs, MSA {min_fwd}-{max_rev}")

    # Step 5: Load ML models once
    print("Step 5: Loading ML models...")
    scaler, classifier, regressor, device = load_models()
    torch.set_num_threads(args.threads)

    # Setup tmpdir
    if args.tmpdir:
        tmpdir = str(Path(args.tmpdir) / "_intermediate")
    else:
        tmpdir = str(Path(args.out).parent)
    os.makedirs(tmpdir, exist_ok=True)

    # Step 6: Per-batch loop
    all_results = []
    all_eval_re = []
    all_eval_full = []
    all_eval_cl = []
    name_to_seq = {}

    try:
        for batch_idx, (batch_primers, min_fwd, max_rev) in enumerate(batches):
            batch_start = time.time()
            print(f"\nBatch {batch_idx+1}/{len(batches)}: {len(batch_primers)} pairs, MSA {min_fwd}-{max_rev}")

            # Extract bowtie2 reference for this batch's region
            margin = 20
            win_start = max(0, min_fwd - margin)
            win_end = max_rev + margin
            all_window = extract_window_sequences(aligned_seqs, all_seq_ids, win_start, win_end)
            if not all_window:
                print("  No sequences in region. Skipping.")
                continue

            win_target_fa = Path(tmpdir) / f"batch_{batch_idx}_targets.fa"
            with open(win_target_fa, "w") as f:
                for sid, seq in all_window:
                    f.write(f">{sid}\n{seq}\n")
            n_targets = len(all_window)
            print(f"  Target sequences: {n_targets}")

            # Build bowtie2 index
            index_prefix = str(Path(tmpdir) / f"batch_{batch_idx}_bt2idx")
            print("  Building bowtie2 index...", end=" ", flush=True)
            build_index(str(win_target_fa), index_prefix, threads=args.threads)
            print("done.")

            # Write batch FASTA and features
            fwd_seqs = [p['fwd_seq'] for p in batch_primers]
            rev_seqs = [p['rev_seq'] for p in batch_primers]
            fasta_path, feat_path, fwd_n2s, rev_n2s = _write_batch_fasta_and_features(
                fwd_seqs, rev_seqs, features, args.name, batch_idx, tmpdir
            )

            name_to_seq.update(fwd_n2s)
            name_to_seq.update(rev_n2s)

            # Align
            print("  Aligning...", end=" ", flush=True)
            mapped_path = _align_batch(
                fasta_path, index_prefix, n_targets, args.threads, tmpdir, batch_idx
            )
            print("done.")

            if mapped_path.stat().st_size == 0:
                print("  No alignments found. Skipping.")
                continue

            # Coverage filter
            filtered_path = _apply_coverage_filter(mapped_path, n_targets, tmpdir, batch_idx, cov_frac=0.90)
            if filtered_path.stat().st_size == 0:
                print("  No primers pass coverage filter. Skipping.")
                continue

            # Prepare input
            print("  Preparing ML input...", end=" ", flush=True)
            input_path = Path(tmpdir) / f"batch_{batch_idx}.input"
            prep_args = argparse.Namespace(
                mapped=str(filtered_path),
                ml_input=str(input_path),
                reference=str(win_target_fa),
                param_file=args.param_file,
                reftype="on",
                pri_features=str(feat_path),
                prev="",
            )
            try:
                run_prepare_input(prep_args)
            except SystemExit:
                pass
            print("done.")

            if not input_path.exists() or input_path.stat().st_size == 0:
                print("  No valid pairs formed. Skipping.")
                continue

            # Evaluate
            print("  Evaluating...", end=" ", flush=True)
            res, eval_re, eval_full, eval_cl = _run_evaluate(
                input_path, str(win_target_fa), scaler, classifier, regressor, device, args.threads,
                n_targets_global=len(aligned_seqs),
            )
            print("done.")

            if res.empty:
                print("  No results. Skipping.")
                continue

            all_results.append(res)
            all_eval_re.append(eval_re)
            if not eval_full.empty:
                all_eval_full.append(eval_full)
            if not eval_cl.empty:
                all_eval_cl.append(eval_cl)
            batch_time = time.time() - batch_start

            passing = res[(res["coverage"] >= cov_min) & (res["activity"] >= act_min)]
            total_passing = sum(
                len(r[(r["coverage"] >= cov_min) & (r["activity"] >= act_min)])
                for r in all_results
            )
            top = res.iloc[0]
            print(f"  Top pair: coverage={top['coverage']:.3f}, activity={top['activity']:.3f}, score={top['score']:.3f}")
            print(f"  Pairs meeting thresholds: {len(passing)} (total: {total_passing}/{min_pairs})")
            print(f"  Batch time: {batch_time:.1f}s")

            if total_passing >= min_pairs:
                print(f"\n  SUCCESS! Found {total_passing} pairs meeting thresholds (target: {min_pairs}).")
                break

        # Output results
        _write_output(args, all_results, all_eval_re, name_to_seq,
                      cov_min, act_min, len(batches), tmpdir, min_dg=min_dg,
                      features=features, all_eval_full=all_eval_full,
                      all_eval_cl=all_eval_cl)

    finally:
        print(f"  Temp files kept at: {tmpdir}")

    runtime = time.time() - start_time
    print(f"\nQuick mode completed in {runtime:.1f}s")


def _run_single_window(args, params, primer_params, cov_min, act_min, min_pairs, start_time):
    """Legacy single-window quick mode (when --msa is not provided)."""
    if not args.selected:
        raise SystemExit("ERROR: --selected is required when --msa is not provided.")
    min_amp_len = primer_params["min_amp_len"]
    max_amp_len = primer_params["max_amp_len"]
    min_dg = primer_params["min_dg"]

    print(f"Quick mode: {args.name}")
    print(f"  Thresholds: coverage >= {cov_min}, activity >= {act_min}")
    print(f"  Max batches: {MAX_BATCHES}, primers/batch: {PRIMERS_PER_BATCH} per direction")

    # Determine which FASTA to use for bowtie2 indexing
    bt2_ref = args.target_window if args.target_window else args.target
    n_targets = sum(1 for _ in SeqIO.parse(bt2_ref, "fasta"))
    print(f"  Target sequences: {n_targets}")
    if args.target_window:
        print(f"  Using window-trimmed reference for bowtie2 index")

    # Step 1: Generate all primers from selected sequences
    print("Step 1: Generating primers...")
    selected_records = list(SeqIO.parse(args.selected, "fasta"))
    target_seqs = [sanitize_iupac(str(r.seq)) for r in selected_records]

    raw_weights = {}
    for idx, rec in enumerate(selected_records):
        frac = 1.0
        if "cluster_frac=" in rec.description:
            try:
                frac = float(rec.description.split("cluster_frac=")[1].split()[0])
            except (ValueError, IndexError):
                pass
        raw_weights[idx] = frac

    medoid_total = sum(w for w in raw_weights.values() if w > 0)
    origin_weights = {}
    for idx, w in raw_weights.items():
        if w == 0.0:
            origin_weights[idx] = CONSENSUS_QUOTA_FRAC
        else:
            origin_weights[idx] = w / medoid_total * (1.0 - CONSENSUS_QUOTA_FRAC)

    for_filt, rev_filt, features = generate_primers_multi(
        target_seqs,
        primer_params["step"],
        primer_params["min_pri_len"],
        primer_params["max_pri_len"],
        min_amp_len,
        max_amp_len,
        primer_params["max_tm"],
        primer_params["min_tm"],
        primer_params["max_gc"],
        min_dg,
    )

    if not for_filt or not rev_filt:
        print("ERROR: No primers passed filters. Cannot proceed.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    print(f"  Generated: {len(for_filt)} fwd, {len(rev_filt)} rev primers")

    # Step 2: Build origin map and compute position bins
    print("Step 2: Tracking primer origins and computing position bins...")
    origin_map = _build_origin_map(
        target_seqs, for_filt, rev_filt,
        primer_params["step"], primer_params["min_pri_len"],
        primer_params["max_pri_len"], min_amp_len,
    )
    n_origins = len(target_seqs)
    cons_idx = n_origins - 1
    cons_fwd_count = sum(1 for s in for_filt if cons_idx in origin_map.get(s, set()))
    cons_rev_count = sum(1 for s in rev_filt if cons_idx in origin_map.get(s, set()))
    print(f"  Origins: {n_origins} sequences ({n_origins-1} representatives + 1 consensus)")
    print(f"  Consensus primers: {cons_fwd_count} fwd, {cons_rev_count} rev")
    total_w = sum(origin_weights.values())
    for idx in range(n_origins):
        w = origin_weights[idx]
        q = max(1, int(round(PRIMERS_PER_BATCH * w / total_w)))
        label = "consensus" if idx == cons_idx else f"rep_{idx}"
        print(f"  Origin {label}: weight={w:.4f}, quota~{q}")

    position_bins = _compute_position_bins(for_filt, MAX_BATCHES)
    for i, (bs, be) in enumerate(position_bins):
        n_f = sum(1 for _, p in for_filt.items() if bs <= p < be)
        print(f"  Bin {i+1}: positions {int(bs)}-{int(be)}, {n_f} fwd primers")

    if not position_bins:
        print("ERROR: No valid primer positions. Try adjusting amplicon length constraints.")
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).touch()
        return

    # Step 3: Build bowtie2 index
    print("Step 3: Building bowtie2 index...")
    if args.tmpdir:
        tmpdir = str(Path(args.tmpdir) / "_intermediate")
    else:
        tmpdir = str(Path(args.out).parent)
    os.makedirs(tmpdir, exist_ok=True)
    index_prefix = str(Path(tmpdir) / "bt2idx")
    build_index(bt2_ref, index_prefix, threads=args.threads)

    # Step 4: Load ML models once
    print("Step 4: Loading ML models...")
    scaler, classifier, regressor, device = load_models()
    torch.set_num_threads(args.threads)

    # Step 5: Batch loop
    all_results = []
    all_eval_re = []
    all_eval_full = []
    all_eval_cl = []
    best_result = None
    name_to_seq = {}

    try:
        for batch_idx, (bin_start, bin_end) in enumerate(position_bins):
            batch_start_time = time.time()

            fwd_selected, rev_selected = _select_batch_primers(
                for_filt, rev_filt, features, origin_map,
                n_origins, origin_weights, bin_start, bin_end,
                min_amp_len, max_amp_len,
            )

            if not fwd_selected or not rev_selected:
                print(f"\nBatch {batch_idx+1}/{MAX_BATCHES}: No primers in bin {int(bin_start)}-{int(bin_end)}. Skipping.")
                continue

            n_cons_f = sum(1 for s in fwd_selected if cons_idx in origin_map.get(s, set()))
            n_cons_r = sum(1 for s in rev_selected if cons_idx in origin_map.get(s, set()))
            n_rep_f = len(fwd_selected) - n_cons_f
            n_rep_r = len(rev_selected) - n_cons_r
            print(f"\nBatch {batch_idx+1}/{MAX_BATCHES}: positions {int(bin_start)}-{int(bin_end)}")
            print(f"  {len(fwd_selected)} fwd ({n_cons_f} cons + {n_rep_f} rep), {len(rev_selected)} rev ({n_cons_r} cons + {n_rep_r} rep)")
            print(f"  Up to {len(fwd_selected) * len(rev_selected)} cross-pairs")

            fasta_path, feat_path, fwd_n2s, rev_n2s = _write_batch_fasta_and_features(
                fwd_selected, rev_selected, features, args.name, batch_idx, tmpdir
            )

            name_to_seq.update(fwd_n2s)
            name_to_seq.update(rev_n2s)

            print("  Aligning...", end=" ", flush=True)
            mapped_path = _align_batch(
                fasta_path, index_prefix, n_targets, args.threads, tmpdir, batch_idx
            )
            print("done.")

            if mapped_path.stat().st_size == 0:
                print("  No alignments found. Skipping batch.")
                continue

            filtered_path = _apply_coverage_filter(mapped_path, n_targets, tmpdir, batch_idx)
            if filtered_path.stat().st_size == 0:
                print("  No primers pass coverage filter. Skipping batch.")
                continue

            print("  Preparing ML input...", end=" ", flush=True)
            input_path = Path(tmpdir) / f"batch_{batch_idx}.input"
            prep_args = argparse.Namespace(
                mapped=str(filtered_path),
                ml_input=str(input_path),
                reference=bt2_ref,
                param_file=args.param_file,
                reftype="on",
                pri_features=str(feat_path),
                prev="",
            )
            try:
                run_prepare_input(prep_args)
            except SystemExit:
                pass
            print("done.")

            if not input_path.exists() or input_path.stat().st_size == 0:
                print("  No valid pairs formed. Skipping batch.")
                continue

            print("  Evaluating...", end=" ", flush=True)
            res, eval_re, eval_full, eval_cl = _run_evaluate(
                input_path, bt2_ref, scaler, classifier, regressor, device, args.threads
            )
            print("done.")

            if res.empty:
                print("  No results. Skipping batch.")
                continue

            all_results.append(res)
            all_eval_re.append(eval_re)
            if not eval_full.empty:
                all_eval_full.append(eval_full)
            if not eval_cl.empty:
                all_eval_cl.append(eval_cl)
            batch_time = time.time() - batch_start_time

            passing = res[(res["coverage"] >= cov_min) & (res["activity"] >= act_min)]
            total_passing = sum(
                len(r[(r["coverage"] >= cov_min) & (r["activity"] >= act_min)])
                for r in all_results
            )
            top = res.iloc[0]
            print(f"  Top pair: coverage={top['coverage']:.3f}, activity={top['activity']:.3f}, score={top['score']:.3f}")
            print(f"  Pairs meeting thresholds: {len(passing)} (total: {total_passing}/{min_pairs})")
            print(f"  Batch time: {batch_time:.1f}s")

            if total_passing >= min_pairs:
                print(f"\n  SUCCESS! Found {total_passing} pairs meeting thresholds (target: {min_pairs}).")
                best_result = passing
                break
            elif not passing.empty:
                best_result = passing
            else:
                if best_result is None or top["score"] > best_result.iloc[0]["score"]:
                    best_result = res

        # Output results
        _write_output(args, all_results, all_eval_re, name_to_seq,
                      cov_min, act_min, MAX_BATCHES, tmpdir, min_dg=min_dg,
                      features=features, all_eval_full=all_eval_full,
                      all_eval_cl=all_eval_cl)

    finally:
        print(f"  Temp files kept at: {tmpdir}")

    runtime = time.time() - start_time
    print(f"\nQuick mode completed in {runtime:.1f}s")


def _write_probe_outputs(args, combined, probe_data):
    """Write probe FASTA and features files for probe-first mode."""
    assignments = probe_data.get('assignments', {})
    all_probes = probe_data.get('probes', [])
    probe_features = probe_data.get('features', {})
    pair_to_probe_idx = probe_data.get('pair_to_probe_idx', {})
    probes_by_probe_idx = probe_data.get('probes_by_probe_idx', {})

    # Collect probes assigned to at least one passing pair.
    # Prefer dimer-checked assignments; fall back to position-based (pair_to_probe_idx).
    used_probe_names = set()
    for _, row in combined.iterrows():
        pair_key = (row['pname_f'], row['pname_r'])
        for probe in assignments.get(pair_key, []):
            used_probe_names.add(probe['name'])

    # If dimer assignments are empty, use position-validated assignments from probe-first design
    if not used_probe_names and pair_to_probe_idx:
        for pair_key, probe_idxs in pair_to_probe_idx.items():
            for pi in probe_idxs:
                for probe in probes_by_probe_idx.get(pi, []):
                    used_probe_names.add(probe['name'])

    # Last resort: include all probes so downstream filter can check
    if not used_probe_names:
        used_probe_names = {p['name'] for p in all_probes}

    # Write probe pair assignment CSV: records which probe is valid for which primer pair.
    # filter_primers uses this to avoid cross-region probe assignment in quick design mode.
    if getattr(args, 'probe_pair_csv', None):
        pair_probe_rows = []
        for _, row in combined.iterrows():
            pair_key = (row['pname_f'], row['pname_r'])
            # Get position-validated probes for this pair
            probe_idxs = pair_to_probe_idx.get(pair_key, set())
            valid_probe_names = set()
            for pi in probe_idxs:
                for probe in probes_by_probe_idx.get(pi, []):
                    # If dimer assignments exist, further restrict to dimer-passing probes
                    if assignments:
                        if any(p['name'] == probe['name'] for p in assignments.get(pair_key, [])):
                            valid_probe_names.add(probe['name'])
                    else:
                        valid_probe_names.add(probe['name'])
            for pname in sorted(valid_probe_names):
                pair_probe_rows.append({
                    'pname_f': row['pname_f'],
                    'pname_r': row['pname_r'],
                    'probe_name': pname,
                })
        Path(args.probe_pair_csv).parent.mkdir(parents=True, exist_ok=True)
        cols = ['pname_f', 'pname_r', 'probe_name']
        pd.DataFrame(pair_probe_rows, columns=cols).to_csv(args.probe_pair_csv, index=False)
        print(f"  Saved probe pair assignments: {args.probe_pair_csv} ({len(pair_probe_rows)} rows)")

    # Write probe FASTA
    if args.probe_fa:
        Path(args.probe_fa).parent.mkdir(parents=True, exist_ok=True)
        written = set()
        with open(args.probe_fa, "w") as f:
            for probe in all_probes:
                if probe['name'] in used_probe_names and probe['name'] not in written:
                    strand = probe.get('strand', 'fwd')
                    f.write(f">{probe['name']} strand={strand}\n{probe['seq']}\n")
                    written.add(probe['name'])
        print(f"  Saved probe FASTA: {args.probe_fa} ({len(written)} probes)")

    # Write probe features
    if args.probe_feat:
        Path(args.probe_feat).parent.mkdir(parents=True, exist_ok=True)
        feat_rows = []
        seen = set()
        for probe in all_probes:
            if probe['name'] in used_probe_names and probe['name'] not in seen:
                seen.add(probe['name'])
                feat_rows.append({
                    'pname': probe['name'],
                    'pseq': probe['seq'],
                    'len': len(probe['seq']),
                    'position': probe['msa_start'],
                    'Tm': probe['tm'],
                    'GC': probe['gc'],
                    'dG': probe['dg'],
                })
        pd.DataFrame(feat_rows).to_csv(args.probe_feat, index=False)
        print(f"  Saved probe features: {args.probe_feat} ({len(feat_rows)} probes)")


def _write_output(args, all_results, all_eval_re, name_to_seq,
                   cov_min, act_min, n_batches, tmpdir, min_dg=-7.0,
                   features=None, probe_data=None, all_eval_full=None,
                   all_eval_cl=None, pair_to_probe_idx=None,
                   combined_override=None):
    """Write combined output from all batches/windows."""
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    if all_results:
        if combined_override is not None:
            combined = combined_override
        else:
            combined = pd.concat(all_results, ignore_index=True)
            combined = combined.sort_values("activity", ascending=False).drop_duplicates(
                subset=["pname_f", "pname_r"]
            )

        combined = combined[combined["coverage"] > 0]
        if combined.empty:
            Path(args.out).touch()
            print("\nWARNING: All pairs had coverage=0. Recommend running full mode.")
        else:
            combined["fwd_seq"] = combined["pname_f"].map(name_to_seq)
            combined["rev_seq"] = combined["pname_r"].map(name_to_seq)

            # Cross-dimer dG filter for cross-position pairs
            n_before = len(combined)
            cross_pairs = list(zip(combined["fwd_seq"], combined["rev_seq"]))
            cross_dgs = compute_batch_dimer_dg(cross_pairs)
            combined["cross_dg"] = cross_dgs
            combined = combined[combined["cross_dg"] >= min_dg].drop(columns=["cross_dg"])
            n_dropped = n_before - len(combined)
            if n_dropped > 0:
                print(f"  Cross-dimer dG filter: dropped {n_dropped}/{n_before} pairs (min_dg={min_dg})")

            # Probe mode: keep all passing pairs — dedup happens in build-output
            if pair_to_probe_idx:
                combined = combined.sort_values("score", ascending=False)

            cols = ["pname_f", "pname_r", "fwd_seq", "rev_seq", "coverage", "activity", "score",
                    "cov_target", "cov_primer_only", "probe_name", "probe_seq"]
            combined = combined[[c for c in cols if c in combined.columns]]
            # Pipeline mode: eval file must have only pname_f,pname_r,coverage,activity,score
            if getattr(args, 'init_fa', None):
                combined[["pname_f", "pname_r", "coverage", "activity", "score"]].to_csv(
                    args.out, index=False)
                # Write marker so filter_primers knows to skip position validation
                Path(f"{args.out}.quick").touch()
            else:
                combined.to_csv(args.out, index=False)

            # Write pipeline-compatible init FASTA (all primers in results)
            if getattr(args, 'init_fa', None):
                Path(args.init_fa).parent.mkdir(parents=True, exist_ok=True)
                seen_seqs = set()
                with open(args.init_fa, "w") as f:
                    for _, row in combined.iterrows():
                        fname, rname = row["pname_f"], row["pname_r"]
                        fseq, rseq = row["fwd_seq"], row["rev_seq"]
                        if fname not in seen_seqs:
                            f.write(f">{fname}\n{fseq}\n")
                            seen_seqs.add(fname)
                        if rname not in seen_seqs:
                            f.write(f">{rname}\n{rseq}\n")
                            seen_seqs.add(rname)
                print(f"  Saved init FASTA: {args.init_fa}")

            # Write pipeline-compatible features file
            if getattr(args, 'init_feat', None):
                Path(args.init_feat).parent.mkdir(parents=True, exist_ok=True)
                feat_rows = []
                seen_names = set()
                for _, row in combined.iterrows():
                    for pname, seq, forrev in [
                        (row["pname_f"], row["fwd_seq"], "f"),
                        (row["pname_r"], row["rev_seq"], "r"),
                    ]:
                        if pname in seen_names:
                            continue
                        seen_names.add(pname)
                        feat = (features or {}).get(seq, {})
                        feat_rows.append({
                            "pname": pname, "pseq": seq, "forrev": forrev,
                            "len": feat.get("len", len(seq)),
                            "Tm": feat.get("Tm", 0),
                            "GC": feat.get("GC", 0),
                            "dG": feat.get("dG", 0),
                        })
                pd.DataFrame(feat_rows).to_csv(args.init_feat, index=False)
                print(f"  Saved init features: {args.init_feat}")

            # Write probe outputs (probe-first mode)
            if probe_data and getattr(args, 'probe_fa', None):
                _write_probe_outputs(args, combined, probe_data)

        if all_eval_re:
            eval_re_combined = pd.concat(all_eval_re, ignore_index=True)
            eval_re_combined = eval_re_combined.fillna(0)
            agg_dict = {c: "max" for c in eval_re_combined.columns[2:]}
            eval_re_combined = eval_re_combined.groupby(
                ["pname_f", "pname_r"]
            ).agg(agg_dict).reset_index()
            re_path = f"{args.out}.re"
            eval_re_combined.to_csv(re_path, index=False)
            print(f"  Saved per-target activity: {re_path}")

        if all_eval_full:
            eval_full_combined = pd.concat(all_eval_full, ignore_index=True)
            eval_full_path = f"{args.out}.full"
            eval_full_combined.to_csv(eval_full_path, index=False)
            print(f"  Saved eval full: {eval_full_path}")

        if all_eval_cl:
            eval_cl_combined = pd.concat(all_eval_cl, ignore_index=True)
            agg_dict = {c: "max" for c in eval_cl_combined.columns[2:]}
            eval_cl_combined = eval_cl_combined.groupby(
                ["pname_f", "pname_r"]
            ).agg(agg_dict).reset_index()
            cl_path = f"{args.out}.cl"
            eval_cl_combined.to_csv(cl_path, index=False)
            print(f"  Saved classifier table: {cl_path}")

        top = combined.iloc[0]
        passes = (top["coverage"] >= cov_min) and (top["activity"] >= act_min)

        if not passes:
            print(f"\nWARNING: No pairs met thresholds after {n_batches} batches.")
            print(f"  Best: coverage={top['coverage']:.3f}, activity={top['activity']:.3f}")
            print(f"  Recommend running full mode: snakemake --cores N")
    else:
        Path(args.out).touch()
        if getattr(args, 'init_fa', None):
            Path(f"{args.out}.quick").touch()
            Path(args.init_fa).parent.mkdir(parents=True, exist_ok=True)
            Path(args.init_fa).touch()
        if getattr(args, 'init_feat', None):
            Path(args.init_feat).parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(columns=["pname", "pseq", "forrev", "len", "Tm", "GC", "dG"]).to_csv(
                args.init_feat, index=False)
        if getattr(args, 'probe_fa', None):
            Path(args.probe_fa).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_fa).touch()
        if getattr(args, 'probe_feat', None):
            Path(args.probe_feat).parent.mkdir(parents=True, exist_ok=True)
            Path(args.probe_feat).touch()
        if getattr(args, 'probe_pair_csv', None):
            Path(args.probe_pair_csv).parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(columns=['pname_f', 'pname_r', 'probe_name']).to_csv(
                args.probe_pair_csv, index=False)
        print("\nWARNING: No results produced. Recommend running full mode.")

