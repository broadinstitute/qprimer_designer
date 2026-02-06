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
TARGETS = ['target']
CROSS   = []
HOST    = ['offtarget1']
```

Run the pipeline:

```bash
cd workflows
snakemake -s Snakefile.example --cores all
```

Output will be in the `final/` directory as a CSV file with the following format:

| pname_f | pname_r | cov_target | act_target | sco_target | cov_offtarget1 | act_offtarget1 | sco_offtarget1 | ... | pseq_f | pseq_r |
|---------|---------|------------|------------|------------|----------------|----------------|----------------|-----|--------|--------|
| p_19_f | p_286_r | 1.0 | 1.084 | 1.084 | 10 | 0.393 | 3.153 | ... | CATTACGTTTGGTGGACCCT | TTCCTTGCCATGTTGAGTGA |
| p_19_f | p_283_r | 1.0 | 1.083 | 1.083 | 12 | 0.412 | 4.116 | ... | CATTACGTTTGGTGGACCCT | GTCTTCCTTGCCATGTTGAG |
| p_19_f | p_285_r | 1.0 | 1.088 | 1.088 | 14 | 0.540 | 5.572 | ... | CATTACGTTTGGTGGACCCT | CTTCCTTGCCATGTTGAGTG |

**Column descriptions:**
- `pname_f/r`: Forward/reverse primer names
- `cov_target`: Coverage on target sequences (fraction)
- `act_target`: Activity score on target
- `sco_target`: Combined score on target
- `cov_offtarget`, `act_offtarget`, `sco_offtarget`: Coverage, activity, and combined scores for each off-target (can have multiple off-targets)
- `pseq_f/r`: Forward/reverse primer sequences

#### Probe Mode (Singleplex only)

For TaqMan-style qPCR assays, you can enable probe mode to generate and filter probes alongside primers:

```bash
cd workflows
snakemake -s Snakefile.example --config probe=1 --cores all
```

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

### Evaluate Custom Primer Set

You can evaluate your own primers instead of generating them. This is useful when you have existing primers and want to assess their performance.

**Option 1: Evaluate direct sequences**

```bash
snakemake -s Snakefile.template \
  --config evaluate=1 \
    for=ATCGATCGATCGATCGATCG \
    rev=GCTAGCTAGCTAGCTAGCTA \
  --cores all
```

**Option 2: Evaluate from FASTA file**

```bash
cd workflows
snakemake -s Snakefile.template \
  --config evaluate=1 pset=my_primers.fa \
  --cores all
```

#### Primer FASTA file naming convention

Primers must follow the `*_for` and `*_rev` naming pattern:

```
>primer1_for
ATCGATCGATCGATCGATCG
>primer1_rev
GCTAGCTAGCTAGCTAGCTA
>primer2_for
TTTTAAAACCCCGGGG
>primer2_rev
GGGGCCCCAAAATTTT
```

Each primer set must have both a forward (`*_for`) and reverse (`*_rev`) primer.

#### Output

Results will be in `evaluate/{pset_name}/` containing Excel reports with:
- **Summary sheet**: Dimerization, sensitivity, and specificity metrics
- **Detail sheet**: Per-target alignments with scores and coverage

Each primer set gets its own Excel file (e.g., `primer1.xlsx`, `primer2.xlsx`).

**Summary sheet format:**

| Dimerization | | | | | |
|--------------|---------|---------|-----|-----|-----|
| | forward | reverse | | | |
| forward | -4.9 | -4.4 | | | |
| reverse | -4.4 | -4.1 | | | |

| Sensitivity | Coverage | Act_mean | Act_median | Act_min | Act_max |
|-------------|----------|----------|------------|---------|---------|
| target | 1 | 1.084 | 1.084 | 1.084 | 1.084 |

| Specificity | Coverage | Act_mean | Act_median | Act_min | Act_max |
|-------------|----------|----------|------------|---------|---------|
| offtarget1 | 10 | 0.248 | 0.212 | 0.138 | 0.393 |

- **Dimerization**: Dimer ΔG values for forward and reverse primers
- **Sensitivity**: On-target coverage and activity statistics
- **Specificity**: Off-target coverage and activity statistics

**Detail sheet format:**

| seq_id | target | eval_type | classifier | regressor | pname_f | pname_r | prod_len |
|--------|--------|-----------|------------|-----------|---------|---------|----------|
| SC2N | SC2N | on | 1 | 1.084 | 1_for | 1_rev | 146 |
| ENST00000580258.1 | HUMAN | off | 1 | 0.393 | 1_rev | 1_rev | 413 |

Additional columns (not shown above): `align_f`, `align_r`, `prod_Tm`, `mm_f`, `indel_f`, `len_f`, `Tm_f`, `GC_f`, `mm_r`, `indel_r`, `len_r`, `Tm_r`, `GC_r`

**Example alignment strings** (from `align_f` and `align_r` columns):

```
align_f (on-target):
1    CATTACGTTTGGTGGACCCT   20
     ||||||||||||||||||||
42   GTAATGCAAACCACCTGGGA   61

align_r (on-target):
1    TTCCTTGCCATGTTGAGTGA   20
     ||||||||||||||||||||
168  AAGGAACGGTACAACTCACT  187

align_f (off-target with mismatches):
1    TTCCTTGCCATGTTGAGTGA   20
      || |  || || |||||||
33   GAGTATGGGCACGACTCACT   52
```

**Column descriptions:**
- **seq_id**: Target sequence identifier
- **target**: Target name (on-target or off-target)
- **eval_type**: Evaluation type ("on" or "off")
- **classifier**: Binary classification score (0-1)
- **regressor**: Activity prediction score
- **pname_f/r**: Forward/reverse primer names
- **align_f/r**: Multi-line alignment strings showing primer-target matches with positions
- **prod_len**: Amplicon length (bp)
- **mm_f/r**: Number of mismatches in forward/reverse primer
- **indel_f/r**: Number of indels in forward/reverse primer

## CLI Commands

The `qprimer` CLI provides the following subcommands:

### Core Commands

| Command | Description |
|---------|-------------|
| `qprimer generate` | Generate primer candidates from target sequences |
| `qprimer prepare-features` | Compute features (Tm, GC%, dG) for existing primers |
| `qprimer pick-representatives` | Select representative sequences from MSA |
| `qprimer prepare-input` | Prepare input data for ML evaluation |
| `qprimer evaluate` | Run ML model to score primer candidates |
| `qprimer filter` | Filter primers based on evaluation scores |
| `qprimer build-output` | Build final output CSV with scores |
| `qprimer select-multiplex` | Select best multiplex primer set |
| `qprimer export-report` | Export evaluation results to Excel reports |

### Probe Mode Commands

| Command | Description |
|---------|-------------|
| `qprimer generate-probe` | Generate probe candidates from target sequences (for TaqMan assays) |
| `qprimer parse-probe-mapping` | Parse probe SAM alignments into mapping table |

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
