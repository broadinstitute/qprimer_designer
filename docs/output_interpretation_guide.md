# qPrimer Designer Output Interpretation Guide

This guide describes the columns and values in the CSV and XLSX files produced by the qPrimer Designer pipeline.

---

## Final Output CSV

**Location:** `final/{target_name}.csv`

This is the main output file. Each row represents a candidate primer pair, sorted by `offtarget_score_sum` (ascending — lower off-target scores are better).

### Index columns

| Column | Description |
|--------|-------------|
| `pname_f` | Forward primer name |
| `pname_r` | Reverse primer name |

### On-target columns

| Column | Description |
|--------|-------------|
| `cov_target` | **Coverage** — fraction of on-target sequences predicted to amplify (0–1). Higher is better. |
| `act_target` | **Activity** — mean predicted amplification activity across covered on-target sequences. Higher is better. Normalized to 1.0 being performance of CDC SARS-CoV-2 N1 primer set on perfect matching sequence.|
| `sco_target` | **Score** — `cov_target × act_target`. Combined metric; higher is better. |

### Off-target columns (one set per off-target panel)

For each off-target panel named `{offname}`:

| Column | Description |
|--------|-------------|
| `cov_{offname}` | Number of off-target sequences predicted to amplify (sum of positives). Lower is better. |
| `act_{offname}` | Maximum predicted activity across off-target sequences. Lower is better. |
| `sco_{offname}` | `cov_{offname} × act_{offname}`. Combined off-target metric; lower is better. |

### Sequence columns

| Column | Description |
|--------|-------------|
| `pseq_f` | Forward primer sequence (5'→3') |
| `pseq_r` | Reverse primer sequence (5'→3') |

### Aggregate off-target column

| Column | Description |
|--------|-------------|
| `offtarget_score_sum` | Sum of all `sco_{offname}` values. The table is sorted by this column ascending. Lower is better. |

### Probe-mode columns (when probes are enabled)

| Column | Description |
|--------|-------------|
| `valid_probes` | Comma-separated names of probes compatible with this primer pair (filtered to avoid off-target amplicon overlap). |
| `valid_probe_sequences` | Comma-separated sequences for the valid probes. |
| `amplicon_seq` | Representative amplicon sequence extracted from the reference representative sequence. |

---

## Evaluation Report XLSX

**Location:** `evaluate/{primer_set_id}.xlsx`

Each XLSX report describes a single primer pair in detail. It contains two sheets.

### Sheet: "summary"

Contains three tables:

#### Dimerization table

A 2×2 matrix showing the predicted ΔG (kcal/mol) for all primer-primer dimer combinations:

|  | forward | reverse |
|--|---------|---------|
| **forward** | F–F dimer ΔG | F–R dimer ΔG |
| **reverse** | R–F dimer ΔG | R–R dimer ΔG |

More negative values indicate stronger (worse) dimer formation. Values above approximately −6 kcal/mol are generally acceptable.

#### Sensitivity table (on-target)

| Column | Description |
|--------|-------------|
| `Coverage` | Number of on-target sequences predicted to amplify (classifier > 0.5). |
| `Act_mean` | Mean predicted activity across on-target sequences. |
| `Act_median` | Median predicted activity. |
| `Act_min` | Minimum predicted activity. |
| `Act_max` | Maximum predicted activity. |

#### Specificity table (off-target)

Same columns as the sensitivity table, but grouped by off-target panel. Lower values indicate better specificity (less off-target amplification).

### Sheet: "detail"

Per-target sequence alignment and scoring data. One row per target sequence.

| Column | Description |
|--------|-------------|
| `classifier` | Binary classification (1 = predicted to amplify, 0 = not). Binarized at a 0.5 threshold from the raw ML output. |
| `regressor` | Predicted amplification activity (continuous score, higher = more active). |
| `align_f` | Visual alignment of the forward primer against the target sequence (monospaced). |
| `align_r` | Visual alignment of the reverse primer against the target sequence (monospaced). |
| `prod_len` | Predicted amplicon product length (bp). |
| `prod_Tm` | Predicted amplicon melting temperature (°C). |
| `mm_f` / `mm_r` | Number of mismatches in the forward / reverse primer alignment. |
| `indel_f` / `indel_r` | Number of indels in the forward / reverse primer alignment. |
| `len_f` / `len_r` | Forward / reverse primer length (nt). |
| `Tm_f` / `Tm_r` | Forward / reverse primer melting temperature (°C). |
| `GC_f` / `GC_r` | Forward / reverse primer GC content (fraction, 0–1). |

---

## Intermediate Files

These files are produced during the pipeline and may be useful for debugging or further analysis.

### Evaluation summary CSV

**Location:** `evaluate/{target_name}.{panel}.csv`

| Column | Description |
|--------|-------------|
| `pname_f` | Forward primer name |
| `pname_r` | Reverse primer name |
| `coverage` | For on-target: fraction of sequences amplified (0–1). For off-target: count of sequences amplified. |
| `activity` | For on-target: mean activity across covered sequences. For off-target: max activity. |
| `score` | `coverage × activity` |

### Filtered primer pairs CSV

**Location:** `filter/{target_name}_pairs.csv`

| Column | Description |
|--------|-------------|
| `pname_f` / `pname_r` | Primer pair names |
| `coverage` | Coverage metric from evaluation |
| `activity` | Activity metric from evaluation |
| `score` | Combined score |
| `pseq_f` / `pseq_r` | Primer sequences |
| `primer_dimer_dg` | Predicted ΔG of primer dimer formation (kcal/mol). More negative = stronger dimer. |
| `num_probes` | (Probe mode) Number of compatible probes for this pair. |
| `probe_names` | (Probe mode) Comma-separated compatible probe names. |

### Primer features file

**Location:** `*.feat`

| Column | Description |
|--------|-------------|
| `pname` | Primer or probe name |
| `pseq` | Sequence |
| `forrev` | Orientation: `f` (forward) or `r` (reverse). Absent for probes. |
| `len` | Length in nucleotides |
| `Tm` | Melting temperature (°C) |
| `GC` | GC content (fraction, 0–1) |
| `dG` | Self-dimer free energy (kcal/mol) |

---

## Key Concepts

- **Classifier**: An ML model that predicts whether a primer pair will successfully amplify a given target (binary outcome, threshold 0.5). Values > 0.5 are predicted to amplify, values < 0.5 are predicted to not amplify.
- **Regressor**: An ML model that predicts the amplification activity/efficiency (continuous value; higher is better for on-target, lower is better for off-target).
- **Coverage**: The proportion (on-target) or count (off-target) of sequences in a panel that a primer pair is predicted to amplify.
- **Activity**: A measure of predicted amplification strength. For on-target panels, the mean across covered targets is used. For off-target panels, the max is used (worst-case).
- **Score**: `coverage × activity`. A single metric combining breadth and strength of amplification.
- **ΔG (dG)**: Gibbs free energy of dimer formation. More negative values indicate more stable (and thus more problematic) dimers. A typical cutoff is −6 kcal/mol.
