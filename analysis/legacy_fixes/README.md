# Historical source repairs and unresolved damage

Two pre-existing public R files do not parse. Their original root copies remain unchanged for provenance.

- `Figure2_3_ATAC_Homer_Plot_fixed.R` replaces invalid single-backslash Windows paths with the distributed HOMER table paths. It also corrects the y-axis label to **−ln(P)**, consistent with the stored HOMER LogP column and current manuscript. Numeric values and selection logic are unchanged. It remains a plotting snippet requiring the original shared plotting packages and `rdwhbu` palette function.
- The root `Figure3_4_IR_vs_Ctrl_Transcriptional_Diversity_Correlates_with_Radiation_Sensitivity(Meld).R` contains binary-corrupted bytes starting at byte 4096 and cannot be reconstructed from its sole Git version. It is an archival damaged file, not a supported execution entry point.
- `run_meld_archived.py` provides the original, separately saved Python MELD analysis for reference. IPython reload commands were commented and two personal paths became environment variables; calculations are otherwise unchanged. It requires original `1_exp_matrix.csv`, `1_metadata.csv` and a parameter-search table. It is not a replacement reconstruction of the lost R plotting commands, and includes exploratory PHATE/vertex-frequency steps.

Set `MELD_PARAMETER_SEARCH_FILE` to the original trusted parameter-search pickle and `MELD_OUTPUT_FILE` to a new output filename before running the archived Python source from its input directory. The historical `Lin-con` control label is preserved in this source; changing it without the corresponding original metadata would change the input subset.
