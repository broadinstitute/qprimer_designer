"""Evaluate primer-target pairs using ML models."""

import argparse
import ast
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from sklearn.preprocessing import MultiLabelBinarizer
from torch.utils.data import DataLoader

from qprimer_designer.models import load_models, PcrDataset, FEATURE_COLUMNS
from qprimer_designer.utils import encode_batch_parallel


def register(subparsers):
    """Register the evaluate subcommand."""
    parser = subparsers.add_parser(
        "evaluate",
        help="Evaluate primer-target matches using ML",
        description="""
Evaluate primer-target pairs using trained ML models (classifier + regressor).
Produces raw per-pair scores, aggregated per-primer-pair scores, and final
ranked table with coverage, activity, and score.
""",
    )
    parser.add_argument("--in", dest="input", required=True, help="Input CSV from prepare-input")
    parser.add_argument("--out", dest="output", required=True, help="Output CSV for evaluation results")
    parser.add_argument("--ref", dest="reference", required=True, help="Reference FASTA")
    parser.add_argument("--reftype", dest="reftype", required=True, choices=["on", "off"], help="on-target or off-target")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for encoding")
    parser.set_defaults(func=run)


def run(args):
    """Run the evaluate command."""
    tnames = [s.id for s in SeqIO.parse(args.reference, "fasta")]

    # Load models
    scaler, classifier, regressor, device = load_models()

    print(f"Evaluating {args.input} with {device}...")
    start_time = time.time()

    if os.path.getsize(args.input) == 0:
        open(args.output, "w").close()
        sys.exit()

    header_flag = True
    mode = 'w'
    clstbl, regtbl = [], []

    for i, chunk in enumerate(pd.read_csv(args.input, chunksize=20000)):
        chunk['targets'] = chunk['targets'].apply(ast.literal_eval)

        inps_fe = chunk[FEATURE_COLUMNS]
        inps_fe = scaler.transform(inps_fe)
        inps_se = chunk[['pseq_f', 'tseq_f', 'pseq_r', 'tseq_r']]
        inps_se = encode_batch_parallel(inps_se, int(args.threads))
        dataset = PcrDataset(inps_se, inps_fe, np.array([0] * len(chunk)))
        loader = DataLoader(dataset, batch_size=64, shuffle=False)

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

        chunk.loc[:, 'classifier'] = predict_cls
        chunk.loc[:, 'regressor'] = predict_reg
        chunk.to_csv(f'{args.output}.full', mode=mode, header=header_flag, index=False)

        mlb = MultiLabelBinarizer()
        onehot = mlb.fit_transform(chunk['targets'])
        target_cols = list(mlb.classes_)

        for label, l in zip(['classifier', 'regressor'], [clstbl, regtbl]):
            targets_df = pd.DataFrame(onehot, columns=target_cols, index=chunk.index)
            targets_df = targets_df.mul(chunk[label], axis=0)
            evaltbl = pd.concat([chunk[['pname_f', 'pname_r']], targets_df], axis=1)
            agg_dict = {c: "max" for c in evaltbl.columns[2:]}
            evaltbl = evaltbl.groupby(['pname_f', 'pname_r']).agg(agg_dict).reset_index()
            l.append(evaltbl)

        header_flag = False
        mode = 'a'

    clstbl = pd.concat(clstbl, ignore_index=True)
    regtbl = pd.concat(regtbl, ignore_index=True).reindex(columns=clstbl.columns)
    agg_dict = {c: "max" for c in regtbl.columns[2:]}
    regtbl = regtbl.reset_index().groupby(['pname_f', 'pname_r']).agg(agg_dict)
    clstbl = clstbl.reset_index().groupby(['pname_f', 'pname_r']).agg(agg_dict)
    clstbl.fillna(0).round(3).to_csv(f'{args.output}.cl')
    regtbl.fillna(0).round(3).to_csv(f'{args.output}.re')

    coverage = (clstbl > .5).sum(axis=1).reset_index(name='coverage')
    if args.reftype == 'on':
        coverage['coverage'] = coverage['coverage'] / len(tnames)
        activity = (regtbl * (clstbl > .5)).mean(axis=1).reset_index(name='activity')
    else:
        activity = (regtbl * (clstbl > .5)).max(axis=1).reset_index(name='activity')

    res = coverage.merge(activity, on=['pname_f', 'pname_r'])
    res['score'] = res['coverage'] * res['activity']
    res = res.dropna().sort_values('score', ascending=False)
    res.round(4).to_csv(args.output, index=False)

    runtime = time.time() - start_time
    print(f"Wrote {len(res)} lines to {args.output} ({runtime:.1f} sec)")
