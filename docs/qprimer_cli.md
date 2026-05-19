# `qprimer` CLI Reference

Low-level subcommands used internally by the Snakemake workflow. These are not intended to be called directly by users — use the `adapt` CLI instead.

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
| `qprimer generate-probe` | Generate probe candidates for TaqMan assays |
| `qprimer parse-probe-mapping` | Parse probe SAM alignments into mapping table |
| `qprimer quick-design` | Batch primer generation with early stopping for multi-sequence targets |

For usage details on any subcommand:

```bash
qprimer <subcommand> --help
```
