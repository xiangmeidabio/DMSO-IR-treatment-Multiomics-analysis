Figure 6h descriptive-expression verification
Prepared: 25 November 2026

These are audit records for the existing plot, not a rerun of the statistical analysis or an original plotting script. No models were trained, predictions generated, or significance tests performed.

figure6h_box_alignment.csv records 180 comparisons: nine genes x four conditions x five box-plot statistics (first quartile, median, third quartile, lower whisker and upper whisker). The PDF values were recovered from vector positions calibrated against labelled y-axis ticks. The archived RNA values were calculated from the original normalized expression matrices. The maximum absolute discrepancy is 0.0007388472 on the plotted expression scale, consistent with PDF coordinate precision.

archived_expression_quantiles.csv records the cell count and descriptive summaries for each of the 36 gene-condition combinations. Quartiles use linear interpolation at (n-1)*p, equivalent to the default type-7 quantile convention. Whiskers in the alignment table are the most extreme observed values within 1.5 interquartile ranges of the box.

figure6h_verification_summary.json states the comparison scope and limits. source_manifest.json identifies the original expression CSV files, their SHA-256 checksums and cell counts. Source paths are relative to the original Radiation project root. The source matrices and cell-level expression values are not redistributed in this folder.

This numerical correspondence supports identification of the plotted values as single-cell log(1 + counts per 10,000) RNA expression. It does not recover the original plotting command or establish the method used to generate the original significance annotations. The original significance tests and any multiple-testing adjustment were not verified. The label Integrated Dataset in the former plot title does not establish use of a Seurat integrated assay.

No inferential conclusion about independent animals is implied by the cell counts. This folder supports descriptive plot provenance only.
