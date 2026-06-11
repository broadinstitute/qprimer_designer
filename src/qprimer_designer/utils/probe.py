"""Wobble-aware probe matching utilities.

Shared between design mode (MSA-based) and evaluate mode (direct sequence matching).

Wobble weights for same-strand comparison (probe vs reference, same strand).
The probe hybridizes with complement(target).  Only G:T and T:G wobble pairs
on the hybridizing strands get a reduced penalty:

  Same-strand pair  →  Hybridization pair   →  Type
  (G, A)            →  G : complement(A)=T  →  G:T wobble
  (T, C)            →  T : complement(C)=G  →  T:G wobble

All other non-identical pairs are full mismatches (penalty 1.0).
"""

from .sequences import complement_dna, reverse_complement_dna

# Primer base-selection weights for _best_primer_base in quick_design.py.
# Maps (primer_base, template_complement_base) → hybridization stability.
# Used for direct hybridization scoring, NOT same-strand mismatch counting.
WOBBLE_W_PRIMER = {
    ('A', 'T'): 1.0, ('T', 'A'): 1.0, ('G', 'C'): 1.0, ('C', 'G'): 1.0,
    ('G', 'T'): 0.80, ('T', 'G'): 0.80,
    ('G', 'A'): 0.80, ('A', 'G'): 0.80,
}

# Wobble tolerance for probes (no ML layer, must be more stringent).
# G:T / T:G wobble on hybridizing strands = 0.50 weight → 0.50 penalty.
WOBBLE_W_PROBE = {
    ('G', 'A'): 0.50,  # G:T wobble
    ('T', 'C'): 0.50,  # T:G wobble
}

# Default alias — probe functions use WOBBLE_W_PROBE
WOBBLE_W = WOBBLE_W_PROBE

# Watson-Crick pairs for hybridization display
_WC_PAIRS = {('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')}
# DNA wobble pairs for hybridization display
_WOBBLE_PAIRS = {('G', 'T'), ('T', 'G')}


def wobble_mismatch_count(probe_seq, target_seq):
    """Count wobble-weighted mismatches between two ungapped sequences.

    Identical bases contribute 0.0, Watson-Crick pairs contribute 0.0,
    G-T/A-G wobble pairs contribute 0.20, all other mismatches contribute 1.0.

    Args:
        probe_seq: Probe sequence (ungapped, uppercase)
        target_seq: Target sequence (ungapped, uppercase, same length as probe)

    Returns:
        (effective_mismatches: float, indels: int)
        indels is always 0 for ungapped comparison.
    """
    mm = 0.0
    for pb, tb in zip(probe_seq, target_seq):
        if pb == tb:
            continue
        w = WOBBLE_W.get((pb, tb), 0.0)
        if w == 0.0:
            mm += 1.0
        else:
            mm += 1.0 - w
    return mm, 0


def wobble_mismatch_count_gapped(probe_seq, target_row, msa_start, msa_end):
    """Count wobble-weighted mismatches from an MSA window.

    Extracts target_row[msa_start:msa_end], counts gap/N positions as indels,
    and applies wobble-weighted mismatch scoring on aligned (non-gap) positions.

    The probe_seq should be in sense orientation (matching the MSA strand).

    Args:
        probe_seq: Probe consensus in sense orientation (ungapped)
        target_row: Full MSA row string for one sequence
        msa_start: Start column in MSA (inclusive)
        msa_end: End column in MSA (exclusive)

    Returns:
        (effective_mismatches: float, n_indels: int)
    """
    window = target_row[msa_start:msa_end]
    mm = 0.0
    n_indels = 0
    for pb, tb in zip(probe_seq, window):
        if tb == '-' or tb == 'N' or tb == 'n':
            n_indels += 1
            continue
        tb_upper = tb.upper()
        if pb == tb_upper:
            continue
        w = WOBBLE_W.get((pb, tb_upper), 0.0)
        if w == 0.0:
            mm += 1.0
        else:
            mm += 1.0 - w
    return mm, n_indels


def wobble_mismatch_count_cols(probe_seq, target_row, msa_cols):
    """Count wobble-weighted mismatches at specific MSA columns.

    Like wobble_mismatch_count_gapped but uses explicit column indices
    instead of a contiguous [start:end] slice. This handles probes built
    from gap-trimmed consensus where the columns may not be contiguous.

    Args:
        probe_seq: Probe consensus in sense orientation (ungapped, len = len(msa_cols))
        target_row: Full MSA row string for one sequence
        msa_cols: List/array of MSA column indices to compare

    Returns:
        (effective_mismatches: float, n_indels: int)
    """
    mm = 0.0
    n_indels = 0
    for pb, col in zip(probe_seq, msa_cols):
        tb = target_row[col]
        if tb == '-' or tb == 'N' or tb == 'n':
            n_indels += 1
            continue
        tb_upper = tb.upper()
        if pb == tb_upper:
            continue
        w = WOBBLE_W.get((pb, tb_upper), 0.0)
        if w == 0.0:
            mm += 1.0
        else:
            mm += 1.0 - w
    return mm, n_indels


def build_match_string(probe_seq, target_complement):
    """Build a pairwise match string for probe-target hybridization display.

    Compares probe sequence against the complement of the target (i.e., the
    strand the probe actually hybridizes with).

    For each position:
      '|' = Watson-Crick pair (A-T, T-A, G-C, C-G)
      '.' = wobble pair (G-T, T-G)
      ' ' = mismatch

    Args:
        probe_seq: Probe sequence (uppercase)
        target_complement: Complement of target sequence (uppercase, same length)

    Returns:
        Match string of same length as inputs.
    """
    chars = []
    for pb, tb in zip(probe_seq, target_complement):
        if (pb, tb) in _WC_PAIRS:
            chars.append('|')
        elif (pb, tb) in _WOBBLE_PAIRS:
            chars.append('.')
        else:
            chars.append(' ')
    return ''.join(chars)


def slide_probe_match(probe_seq, target_seq, max_mismatches, max_indels=0):
    """Slide a probe across a target sequence, checking both orientations.

    Uses wobble-aware mismatch counting at each position. Returns all
    positions where effective_mismatches <= max_mismatches.

    Args:
        probe_seq: Probe sequence (ungapped, uppercase)
        target_seq: Target sequence (ungapped, uppercase)
        max_mismatches: Maximum wobble-weighted mismatches allowed
        max_indels: Maximum indels allowed (always 0 for ungapped comparison)

    Returns:
        list of dicts: {start_pos, orientation, mismatches, indels}
    """
    probe_len = len(probe_seq)
    target_len = len(target_seq)
    if probe_len > target_len:
        return []

    target_upper = target_seq.upper()
    hits = []

    for orientation, seq in [('+', probe_seq.upper()),
                              ('-', reverse_complement_dna(probe_seq).upper())]:
        for i in range(target_len - probe_len + 1):
            window = target_upper[i:i + probe_len]
            mm, indels = wobble_mismatch_count(seq, window)
            if mm <= max_mismatches and indels <= max_indels:
                hits.append({
                    'start_pos': i,
                    'orientation': orientation,
                    'mismatches': round(mm, 2),
                    'indels': indels,
                })

    return hits
