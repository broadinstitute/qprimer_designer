# qprimer-designer

[![Tests](https://github.com/broadinstitute/qprimer_designer/actions/workflows/test.yml/badge.svg)](https://github.com/broadinstitute/qprimer_designer/actions/workflows/test.yml)
[![Docker Build](https://github.com/broadinstitute/qprimer_designer/actions/workflows/docker.yml/badge.svg)](https://github.com/broadinstitute/qprimer_designer/actions/workflows/docker.yml)
[![codecov](https://codecov.io/gh/broadinstitute/qprimer_designer/graph/badge.svg)](https://codecov.io/gh/broadinstitute/qprimer_designer)

ML-guided qPCR primer design with off-target minimization.

## Installation

### Using conda (recommended)

```bash
# Create environment with external tools
conda env create -f environment.yml
conda activate qprimer-designer

# Install the package
pip install .
```

### Using Docker

```bash
docker pull ghcr.io/broadinstitute/qprimer_designer:latest
docker run --rm ghcr.io/broadinstitute/qprimer_designer qprimer --help
```

## Web GUI

A Streamlit-based GUI is available as an alternative to the CLI. See the [GUI Getting Started Guide](docs/gui_getting_started.md) for detailed setup instructions.

```bash
pip install -e ".[gui]"
streamlit run gui/app.py
```

The GUI uses a multi-page workflow with sidebar navigation:

- **Home** — Entry point with getting started guide
- **Design / Evaluate / Monitor** — Step-by-step workflows: select targets, configure parameters, run, and view results
- **Past Results** — Browse and download outputs from previous runs

## Quick Start

```bash
adapt --help
adapt design --help
adapt evaluate --help
adapt fetch --help
adapt monitor --help
```

## Usage

### Directory Structure

The pipeline expects input sequences in the following structure:

```
.
└── target_seqs/
    └── original/
        ├── target1.fa
        ├── target2.fa
        └── offtarget.fa
```

FASTA files in `./target_seqs/original` should be multi-sequence, unaligned FASTAs with `.fa` extension.

### Configuration

Copy the template and fill in your parameters:

```bash
cp workflows/params.txt.template params.txt
```

Key sections in `params.txt`:
- **Primer generation** — Tm, GC%, length, tiling parameters
- **Probe generation** — Probe length, Tm, homopolymer, 5' G avoidance
- **Probe mode** — Mismatch tolerance, amplicon buffer, probes per pair
- **Amplicon** — Min/max amplicon and off-target lengths
- **Email** — Gmail sender, app password, and recipients (for `adapt monitor`)

### Design Primers

```bash
# Singleplex
adapt design --params params.txt

# With probe design
adapt design --params params.txt --probe

# Multiplex
adapt design --params params.txt --multiplex --probe

# Dry run (preview without executing)
adapt design --params params.txt --dry-run
```

Output will be in `runs/{run_id}/` (timestamped by default, or use `--runid` to specify):
- `{target}_final.csv` — Primer candidates with scores and coverage
- For probe mode: `{target}_probe.fa` and mapping results

### Evaluate Custom Primer Set

Evaluate existing primers against target sequences using the ML model.

```bash
# Option 1: Evaluate direct sequences
adapt evaluate --for ATCGATCGATCG --rev GCTAGCTAGCTA

# Option 2: With probe
adapt evaluate --for ATCGATCGATCG --rev GCTAGCTAGCTA --pro AACCGGTTAACCGG

# Option 3: Evaluate from FASTA file (multiple primer sets)
adapt evaluate --pset my_primers.fa
```

#### Primer FASTA file format

Primers must follow the `*_for`, `*_rev`, and optionally `*_pro` naming pattern:

```
>primer1_for
ATCGATCGATCGATCGATCG
>primer1_rev
GCTAGCTAGCTAGCTAGCTA
>primer1_pro
AACCGGTTAACCGGTTAACC
>primer2_for
TTTTAAAACCCCGGGG
>primer2_rev
GGGGCCCCAAAATTTT
```

Each primer set must have both a forward (`*_for`) and reverse (`*_rev`) entry. Probe (`*_pro`) is optional.

#### Output

Results will be in `evaluate/{run_id}/{pset_name}/` containing Excel reports with:
- **Summary sheet**: Primer/probe sequences, dimerization table (2x2 or 3x3 with probe), sensitivity (coverage as `covered / total`), and specificity metrics
- **Detail sheet**: Per-target alignments with classifier/regressor scores, decision/reason columns, probe match status, and mismatch counts. Unmapped sequences are included with `classifier=0, regressor=unmapped`.

Each primer set gets its own Excel file (e.g., `primer1.xlsx`, `primer2.xlsx`).

### Fetch Sequences

Download virus sequences from NCBI using a Google Sheets configuration spreadsheet.

```bash
# Fetch all queries defined in the spreadsheet
adapt fetch --params params.txt

# Fetch specific query IDs
adapt fetch --params params.txt --query-ids 7 12

# Override spreadsheet URL
adapt fetch --url "https://docs.google.com/spreadsheets/d/..." --query-ids 7
```

The spreadsheet must have columns including `query_id`, `Pathogen`, and search parameters. Configure `SPREADSHEET_URL` and optionally `QUERY_IDS` in `params.txt`, or pass them via CLI.

### Monitor Primer Performance

Automatically fetch new sequences, evaluate primers, and send email alerts with results. Designed for periodic surveillance of primer set performance against emerging variants.

```bash
# Run monitor once
adapt monitor --params params.txt

# With explicit run ID (groups results from the same spreadsheet)
adapt monitor --params params.txt --runid my_panel

# Dry run
adapt monitor --params params.txt --dry-run

# Schedule monthly cron job
adapt monitor --params params.txt --schedule

# Remove cron job
adapt monitor --unschedule
```

#### How it works

1. **Fetch** — Downloads latest sequences from NCBI via the configured Google Sheets spreadsheet
2. **Diff** — Compares accession IDs against the most recent previous fetch to identify new sequences
3. **Evaluate** — Runs the ML evaluation pipeline on new sequences using primer/probe sets from the spreadsheet (`Forward`, `Reverse`, `Probe` columns)
4. **Email** — Sends an alert with:
   - Number of new sequences detected (with length, geographic region, release date)
   - Sensitivity/specificity summary tables
   - Excel report attachments for each pathogen

#### Email configuration

In `params.txt`:

```
EMAIL_SENDER = your.alert@gmail.com
EMAIL_PASSWORD = xxxx xxxx xxxx xxxx
EMAIL_RECIPIENTS = recipient1@example.com,recipient2@example.com
```

The sender must be a Gmail account with [App Passwords](https://support.google.com/accounts/answer/185833) enabled (requires 2-Step Verification). The password is a 16-character app password, not the account password.

#### Monitor data structure

Results are organized by run ID and date. The `--runid` flag (default: derived from spreadsheet ID) groups results from the same monitoring configuration. Each run date gets a flat directory with all inputs and outputs:

```
monitor/
└── {runid}/
    └── YYYYMMDD/
        ├── mastersheet.csv                # Spreadsheet snapshot
        ├── {target}_pset.fa               # Primer set FASTA
        ├── {target}_new.fa                # New sequences (evaluate input)
        ├── {target}_accessions.txt        # Accession list (for next diff)
        ├── {target}_metadata.csv          # Sequence metadata
        └── {target}_{primer}.xlsx         # Evaluation reports
```

The full fetched FASTA is removed after extracting the new-only subset and saving the accession list, to conserve disk space. Previous accession lists are used to diff against future fetches.

The pipeline also uses internal `qprimer` subcommands via Snakemake. See [docs/qprimer_cli.md](docs/qprimer_cli.md) for details.

## GPU Support

If GPU is available, add the resource flag:

```bash
adapt design --params params.txt  # auto-detects GPU
```

Or with raw Snakemake:

```bash
snakemake -s Snakefile.example --cores all --resources gpu=1
```

CPU performance is acceptable for most use cases.

## Development

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v
```

See [CLAUDE.md](CLAUDE.md) for development guidelines.

## Pre-trained Models

Pre-trained models are bundled with the package in `src/qprimer_designer/data/`. Training scripts are available in the `training/` directory for reference (raw dataset available upon request).

## License

MIT License - see [LICENSE](LICENSE) for details.
