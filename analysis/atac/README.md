# ATAC motif enrichment

`motif_enrichment.R` plots the archived HOMER enrichment results for regions with increased or decreased accessibility after irradiation (Fig. 2). Inputs are in [`resources/atac/`](../../resources/atac/).

Run from the repository root with the original plotting packages and `rdwhbu` palette function loaded. The script uses the corrected −ln(P) axis label and replaces the earlier copy with invalid Windows path escapes. Peak selection and plotted values are unchanged.
