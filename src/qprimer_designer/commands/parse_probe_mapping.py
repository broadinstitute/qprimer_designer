"""Parse probe SAM alignments into mapping table."""

import argparse
import pandas as pd
from qprimer_designer.utils import parse_params


def register(subparsers):
    """Register the parse-probe-mapping subcommand."""
    parser = subparsers.add_parser(
        "parse-probe-mapping",
        help="Parse probe SAM file into mapping table",
        description="Extract probe-to-target mappings from bowtie2 SAM output."
    )
    parser.add_argument("--sam", required=True, help="Input SAM file from probe alignment")
    parser.add_argument("--out", required=True, help="Output CSV file")
    parser.add_argument("--params", required=True, help="Parameters file")
    parser.set_defaults(func=run)


def run(args):
    """Parse SAM file and create probe mapping CSV."""
    params = parse_params(args.params)
    max_mismatches = int(params.get("PROBE_MAX_MISMATCHES", 2))

    mappings = []

    with open(args.sam) as f:
        for line in f:
            if line.startswith('@'):  # Skip SAM header
                continue

            fields = line.strip().split('\t')
            probe_name = fields[0]
            flag = int(fields[1])
            target_id = fields[2]
            start_pos = int(fields[3]) - 1  # Convert to 0-based
            probe_seq = fields[9]

            # Extract NM tag (edit distance)
            nm_tag = None
            for field in fields[11:]:
                if field.startswith('NM:i:'):
                    nm_tag = int(field.split(':')[2])
                    break

            if nm_tag is None or nm_tag > max_mismatches:
                continue

            orientation = '-' if (flag & 16) else '+'

            mappings.append({
                'probe_name': probe_name,
                'probe_seq': probe_seq,
                'target_id': target_id,
                'start_pos': start_pos,
                'num_mismatches': nm_tag,
                'orientation': orientation
            })

    df = pd.DataFrame(mappings)
    df.to_csv(args.out, index=False)
    print(f"Parsed {len(df)} probe mappings to {args.out}")
