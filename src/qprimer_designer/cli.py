#!/usr/bin/env python
"""qprimer - ML-guided PCR primer design CLI."""

import argparse
import sys


def main():
    """Main entry point for the qprimer CLI."""
    parser = argparse.ArgumentParser(
        prog="qprimer",
        description="ML-guided PCR primer design with off-target minimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  qprimer generate --in targets.fa --out primers.fa --params params.txt
  qprimer evaluate --in input.csv --out output.csv --ref reference.fa --reftype on
  qprimer pick-representatives --in msa.fa --out reps.fa --params params.txt --name virus

For more information on a specific command:
  qprimer <command> --help
""",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="commands",
        description="Available subcommands",
        metavar="<command>",
    )

    # Import and register subcommands
    from qprimer_designer.commands import (
        generate,
        pick_representatives,
        prepare_input,
        evaluate,
        filter_primers,
        build_output,
        select_multiplex,
    )

    generate.register(subparsers)
    pick_representatives.register(subparsers)
    prepare_input.register(subparsers)
    evaluate.register(subparsers)
    filter_primers.register(subparsers)
    build_output.register(subparsers)
    select_multiplex.register(subparsers)

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Execute the subcommand
    args.func(args)


if __name__ == "__main__":
    main()
