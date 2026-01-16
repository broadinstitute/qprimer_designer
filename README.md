# qprimer-designer

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
docker pull ghcr.io/kbryanhsu/qprimer_designer:latest
docker run --rm ghcr.io/kbryanhsu/qprimer_designer qprimer --help
```

## Quick Start

Verify installation:

```bash
qprimer --help
qprimer generate --help
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

### Singleplex qPCR

Configure `workflows/Snakefile.template` with your targets:

```python
TARGETS = ['target1']
CROSS   = []
HOST    = ['offtarget']
```

Run the pipeline:

```bash
cd workflows
snakemake -s Snakefile.example --cores all
```

Output will be in the `final/` directory as a CSV file.

### Multiplex qPCR

Configure the panel targets:

```python
PANEL = ['target1', 'target2']
HOST  = ['offtarget']
```

Run with multiplex enabled:

```bash
snakemake -s Snakefile.example --config multiplex=1 --cores all
```

Output will be `final/multiplex_output.csv` containing top candidates for each target.

## CLI Commands

The `qprimer` CLI provides the following subcommands:

| Command | Description |
|---------|-------------|
| `qprimer generate` | Generate primer candidates from target sequences |
| `qprimer pick-representatives` | Select representative sequences from MSA |
| `qprimer prepare-input` | Prepare input data for ML evaluation |
| `qprimer evaluate` | Run ML model to score primer candidates |
| `qprimer filter` | Filter primers based on evaluation scores |
| `qprimer build-output` | Build final output CSV with scores |
| `qprimer select-multiplex` | Select best multiplex primer set |

## GPU Support

If GPU is available, add the resource flag:

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
