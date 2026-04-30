"""Prepare input for ML evaluation from mapped alignments."""

import argparse
import sys
import time
from collections import defaultdict
from pathlib import Path

from math import log10

import numpy as np
import pandas as pd
from Bio import SeqIO
from pandas.errors import EmptyDataError

from qprimer_designer.utils import parse_params, reverse_complement_dna

# GC-based Tm with offset to approximate Tm_NN(Na=50, Mg=1.5, dNTPs=0.6)
_LOG10_NA = log10(0.05)
_TM_OFFSET = 9.6


def _gc_tm(gc_count, length):
    """Fast amplicon Tm from pre-computed GC count and length."""
    return 81.5 + 16.6 * _LOG10_NA + 41.0 * gc_count / length - 675.0 / length + _TM_OFFSET


def register(subparsers):
    """Register the prepare-input subcommand."""
    parser = subparsers.add_parser(
        "prepare-input",
        help="Prepare input for ML evaluation",
        description="""
Convert mapped primer-target alignments into a feature table suitable
for ML evaluation. Pairs forward and reverse primers into candidate
amplicons and computes amplicon-level features.
""",
    )
    parser.add_argument("--in", dest="mapped", required=True, help="Mapped alignments file")
    parser.add_argument("--out", dest="ml_input", required=True, help="Output CSV for ML input")
    parser.add_argument("--ref", dest="reference", required=True, help="Reference FASTA")
    parser.add_argument("--params", dest="param_file", required=True, help="Parameters file")
    parser.add_argument("--reftype", dest="reftype", required=True, choices=["on", "off"], help="on-target or off-target")
    parser.add_argument("--features", dest="pri_features", required=True, help="Primer features CSV")
    parser.add_argument("--prev", default="", help="Previous evaluation file (for off-target restriction)")
    parser.set_defaults(func=run)


def _ensure_list(x):
    """Ensure aggregated values are always lists."""
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    return [x]


def run(args):
    """Run the prepare-input command."""
    params = parse_params(args.param_file)
    min_amp_len = int(params.get("AMPLEN_MIN", 60))
    max_amp_len = int(params.get("AMPLEN_MAX", 200))
    min_off_len = int(params.get("OFFLEN_MIN", 50))
    max_off_len = int(params.get("OFFLEN_MAX", 5000))
    num_select = int(params.get("NUM_TOP_SENSITIVITY", 100))

    print(f"Preparing ML input from {args.mapped}...")
    start_time = time.time()
    nlines = 0

    # Load reference sequences and build GC prefix sums
    tarseqs = {s.id: str(s.seq) for s in SeqIO.parse(args.reference, "fasta")}
    gc_prefix = {}
    for sid, seq in tarseqs.items():
        arr = np.frombuffer(seq.upper().encode(), dtype='uint8')
        gc = (arr == ord('G')) | (arr == ord('C'))
        gc_prefix[sid] = np.concatenate([[0], np.cumsum(gc)])

    # Primer feature table
    feats = pd.read_csv(args.pri_features, index_col=0)

    # Previous evaluation (OFF-target restriction)
    valid_pairs = None
    if args.prev:
        try:
            teval = pd.read_csv(args.prev).iloc[:num_select]
            arr = teval[['pname_f', 'pname_r']].to_numpy()
            comb = np.vstack([
                np.column_stack([arr[:, 0], arr[:, 0]]),
                np.column_stack([arr[:, 0], arr[:, 1]]),
                np.column_stack([arr[:, 1], arr[:, 0]]),
                np.column_stack([arr[:, 1], arr[:, 1]]),
            ])
            valid_pairs = pd.DataFrame(comb, columns=['pname_f', 'pname_r']).drop_duplicates(ignore_index=True)
        except (FileNotFoundError, EmptyDataError):
            valid_pairs = None

    # Load mapped alignments
    cols = ['pname', 'orientation', 'tname', 'start', 'pseq', 'tseq', 'match']
    try:
        raw = pd.read_table(args.mapped, sep='\t', names=cols)
        if raw.empty:
            Path(args.ml_input).touch()
            sys.exit()
    except EmptyDataError:
        Path(args.ml_input).touch()
        sys.exit()

    raw['orientation'] = raw['orientation'] % 256

    maptbl = (
        raw.groupby(['pname', 'pseq', 'tseq', 'orientation'], sort=False)
        .agg(
            tnames=('tname', list),
            starts=('start', list),
            match=('match', 'first'),
        )
        .reset_index()
    )

    maptbl['indel'] = maptbl['pseq'].str.count('-') + maptbl['tseq'].str.count('-')
    maptbl['mm'] = maptbl['pseq'].str.len() - maptbl['match'].str.count(r'\|') - maptbl['indel']
    maptbl = maptbl.join(feats[['forrev', 'len', 'Tm', 'GC']], on='pname')

    drop_cols = ['orientation', 'forrev']

    if args.reftype == 'on':
        fors = maptbl[(maptbl['orientation'] == 0) & (maptbl['forrev'] == 'f')].drop(columns=drop_cols)
        revs = maptbl[(maptbl['orientation'] == 16) & (maptbl['forrev'] == 'r')].drop(columns=drop_cols)
        minl, maxl = min_amp_len, max_amp_len
        lfunc = max
    else:
        fors = maptbl[maptbl['orientation'] == 0].drop(columns=drop_cols)
        revs = maptbl[maptbl['orientation'] == 16].drop(columns=drop_cols)
        minl, maxl = min_off_len, max_off_len
        lfunc = min

    revs['pseq'] = revs['pseq'].apply(reverse_complement_dna)
    revs['tseq'] = revs['tseq'].apply(reverse_complement_dna)

    # Build sorted reverse index per target for binary search
    rev_by_target = defaultdict(list)
    for r in revs.itertuples():
        r_id = r.Index
        r_len = int(r.len)
        for t, st in zip(_ensure_list(r.tnames), _ensure_list(r.starts)):
            rev_by_target[t].append((int(st), r_id, r_len))

    # Sort by start position and convert to numpy arrays
    rev_sorted = {}
    for t, entries in rev_by_target.items():
        entries.sort()  # sort by start position
        rev_sorted[t] = (
            np.array([e[0] for e in entries]),  # starts
            np.array([e[1] for e in entries]),  # ids
            np.array([e[2] for e in entries]),  # lens
        )
    del rev_by_target

    if args.reftype == 'off' and valid_pairs is not None:
        allowed_r_by_f = valid_pairs.groupby('pname_f')['pname_r'].agg(set).to_dict()
        rname_by_id = revs['pname']

    allpairs = []
    n_fors = len(fors)
    next_log_time = start_time + 60

    for i, f in enumerate(fors.itertuples()):
        now = time.time()
        if now >= next_log_time:
            elapsed = now - start_time
            pct = (i + 1) / n_fors * 100
            print(f"  Progress: {i+1}/{n_fors} ({pct:.1f}%) — {elapsed:.0f}s elapsed", flush=True)
            next_log_time = now + 60

        f_id = f.Index
        fname = f.pname

        tnamesf = _ensure_list(f.tnames)
        startsf = _ensure_list(f.starts)
        if not tnamesf:
            continue

        targets_by_r = defaultdict(list)
        starts_by_r = defaultdict(list)
        amplens_by_r = defaultdict(set)
        tms_by_r = defaultdict(list)

        for t_f, st_f in zip(tnamesf, startsf):
            rev_data = rev_sorted.get(t_f)
            if rev_data is None:
                continue

            st_f = int(st_f)
            r_starts_arr, r_ids_arr, r_lens_arr = rev_data

            # Binary search for candidate reverse primers
            # amplicon_len = r_start + r_len - st_f + 1
            # Need: minl <= amplicon_len <= maxl
            # Use conservative bounds (min/max r_len) then exact-filter
            min_r_len = r_lens_arr[0] if len(r_lens_arr) == 1 else r_lens_arr.min()
            max_r_len = r_lens_arr[0] if len(r_lens_arr) == 1 else r_lens_arr.max()
            lo_bound = st_f - 1 + minl - max_r_len
            hi_bound = st_f - 1 + maxl - min_r_len

            lo_idx = int(np.searchsorted(r_starts_arr, lo_bound, side='left'))
            hi_idx = int(np.searchsorted(r_starts_arr, hi_bound, side='right'))

            if lo_idx >= hi_idx:
                continue

            # Vectorized length filtering on candidates
            cand_starts = r_starts_arr[lo_idx:hi_idx]
            cand_ids = r_ids_arr[lo_idx:hi_idx]
            cand_lens = r_lens_arr[lo_idx:hi_idx]
            amp_lens = cand_starts + cand_lens - st_f + 1
            valid_mask = (amp_lens >= minl) & (amp_lens <= maxl)

            if not valid_mask.any():
                continue

            v_ids = cand_ids[valid_mask]
            v_amp_lens = amp_lens[valid_mask]
            v_starts = cand_starts[valid_mask]
            v_lens = cand_lens[valid_mask]

            for j in range(len(v_ids)):
                r_id = int(v_ids[j])

                if args.reftype == 'off' and valid_pairs is not None:
                    if rname_by_id[r_id] not in allowed_r_by_f.get(fname, set()):
                        continue

                al = int(v_amp_lens[j])
                targets_by_r[r_id].append(t_f)
                starts_by_r[r_id].append(st_f)
                amplens_by_r[r_id].add(al)

                pfx = gc_prefix.get(t_f)
                if pfx is not None:
                    amp_start = st_f - 1
                    amp_end = min(int(v_starts[j]) + int(v_lens[j]), len(pfx) - 1)
                    gc_count = int(pfx[amp_end] - pfx[amp_start])
                    tms_by_r[r_id].append(_gc_tm(gc_count, al))

        if not targets_by_r:
            continue

        rev_ids = list(targets_by_r.keys())
        revsub = revs.loc[rev_ids].copy()
        revsub['targets'] = [targets_by_r[r] for r in revsub.index]
        revsub['starts'] = [starts_by_r[r] for r in revsub.index]
        revsub['prod_len'] = [lfunc(amplens_by_r[r]) for r in revsub.index]
        revsub['prod_Tm'] = [round(np.mean(tms_by_r[r]), 1) for r in revsub.index]

        forsub = fors.loc[[f_id]].copy().drop(['tnames', 'starts'], axis=1)
        pairs = forsub.merge(
            revsub.drop(['tnames'], axis=1),
            how='cross',
            suffixes=('_f', '_r'),
        )

        nlines += len(pairs)
        allpairs.append(pairs)

    if not allpairs:
        Path(args.ml_input).touch()
    else:
        pd.concat(allpairs, ignore_index=True).to_csv(args.ml_input, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {nlines} lines to {args.ml_input} ({runtime:.1f} sec)")
