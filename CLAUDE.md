# qprimer_designer Development Guide

## Project Structure

This is a Python package for ML-guided qPCR primer design. Key directories:

- `src/qprimer_designer/` - Main package source
  - `cli.py` - Single CLI entry point with subcommands
  - `commands/` - Individual subcommand implementations
  - `models/` - PyTorch ML model architectures
  - `external/` - Wrappers for external tools (ViennaRNA, bowtie2, MAFFT)
  - `utils/` - Shared utilities (sequence ops, encoding, params)
  - `data/` - Pre-trained ML models (bundled as package data)
- `workflows/` - Snakemake workflow templates
- `training/` - Historical model training code (not production)
- `tests/` - pytest test suite

## Development Setup

```bash
# Create conda environment with external tools
conda env create -f environment.yml
conda activate qprimer-designer

# Install package in editable mode with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/ -v
```

## CLI Usage

Single entry point with subcommands:
```bash
qprimer generate --help
qprimer evaluate --help
qprimer pick-representatives --help
qprimer prepare-input --help
qprimer filter --help
qprimer build-output --help
qprimer select-multiplex --help
```

## Key Patterns

### External Tool Wrappers
External bioinformatics tools (RNAduplex, bowtie2, mafft) are assumed to be
in PATH via conda installation. Use `shutil.which()` to verify availability.

### ML Model Loading
Models are bundled as package data. Load using `importlib.resources`:
```python
from importlib.resources import files
model_path = files('qprimer_designer.data').joinpath('combined_classifier.pth')
```

### Adding New Subcommands
1. Create module in `src/qprimer_designer/commands/`
2. Implement `register(subparsers)` function to add argparse subparser
3. Import and register in `cli.py`

## Testing

- Run all tests: `pytest tests/ -v`
- Run with coverage: `pytest tests/ -v --cov=qprimer_designer`
- Tests should not require external tools (mock them)

## Docker

Build locally:
```bash
docker build -t qprimer-designer:local .
docker run --rm qprimer-designer:local qprimer --help
```

## Snakemake Workflows

Workflows are in `workflows/`. They use the `qprimer` CLI internally:
```bash
cd workflows
snakemake -s Snakefile.example --cores all
```

Dry-run validation:
```bash
snakemake -s Snakefile.example --dry-run
```

## Environment Variables

- `QPRIMER_FONT_PATH`: Custom font directory for training plots (optional)
- `QPRIMER_TOOLPATH`: Custom tool installation path for training scripts (optional)
- `RNASTRUCTURE_DATAPATH`: Path to RNAstructure data tables (optional)
