# Transcriptional sensitivity

| Script | Analysis |
|---|---|
| `ir_expression_diversity.R` | Expression diversity and the IR-versus-Ctrl response |
| `dmso_expression_diversity.R` | Expression diversity and the DMSO-versus-Ctrl response |
| `nfeature_dsb_score.R` | Detected genes and the DNA double-strand-break expression score |
| `meld.py` | Archived MELD analysis and parameter exploration |

The R scripts use objects from the original analysis session. The MELD script requires `1_exp_matrix.csv`, `1_metadata.csv` and the saved parameter-search table. Set `MELD_PARAMETER_SEARCH_FILE` and `MELD_OUTPUT_FILE`, then run it from the input directory. Its original `Lin-con` control label must match the input metadata.

The supplied Python source is separate from the damaged historical MELD R plot. See the [source notes](../../docs/provenance.md) for file history and input availability.
