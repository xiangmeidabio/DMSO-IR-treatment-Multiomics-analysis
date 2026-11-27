Figure 6h plotting code and expression verification
Prepared: 25 November 2026; author-supplied plotting code added 26 November 2026

Author-supplied code
plot_Fig6h_author_supplied.R preserves the full R code supplied by the author on 26 November 2026. Comments and whitespace were organized for archiving; analysis choices, omitted defaults, plotting settings, title and output filename were retained. The script was not executed during this audit, and no significance tests were rerun.

The script requires the original Seurat object merge_inter. It sets the default assay to RNA, fetches the data layer for nine genes and orig.ident, and retains Ctrl, IR, DMSO and DMSO-IR cells in that order. The nine facets follow Uqcrb, Ptma, Rps2, Gapdh, Impa2, Prtn3, Ramp1, Snrpf and Cenpw. The expression scale is supported by the archived RNA checks described below; the existing plot title does not indicate use of a Seurat integrated assay.

The four comparisons for each gene are Ctrl versus IR, Ctrl versus DMSO, Ctrl versus DMSO-IR, and IR versus DMSO-IR. stat_compare_means uses method="wilcox.test", label="p.signif" and vjust=0.5, with other statistical arguments left at their defaults. The comparisons branch calls ggsignif directly: paired=FALSE and the default two-sided alternative give unpaired Wilcoxon rank-sum tests on cell-level expression. The displayed P values are unadjusted; the compare_means default Holm procedure is not used by this branch. With these group sizes, the default wilcox.test uses a normal approximation with continuity correction. No multiple-testing adjustment was added to the archived code.

The four-star mapping used in ggpubr v0.6.0 and v0.6.1 is consistent with the figure: * P<0.05, ** P<0.01, *** P<0.001, **** P<0.0001, and ns P>=0.05. Source inspection confirms the code path and defaults, not the version originally installed. The author-supplied code does not report an original package session. Each facet's four brackets match the specified comparisons; this visual check does not recalculate the 36 P values.

Original numerical audit
figure6h_box_alignment.csv records 180 comparisons: nine genes x four conditions x five box-plot statistics (first quartile, median, third quartile, lower whisker and upper whisker). The PDF values were recovered from vector positions calibrated against labelled y-axis ticks. The archived RNA values were calculated from the original normalized expression matrices. The maximum absolute discrepancy is 0.0007388472 on the plotted expression scale, consistent with PDF coordinate precision.

archived_expression_quantiles.csv records cell counts and descriptive summaries for the 36 gene-condition combinations. Quartiles use linear interpolation at (n-1)*p, equivalent to the default type-7 quantile convention. Whiskers in the alignment table are the most extreme observed values within 1.5 interquartile ranges of the box.

figure6h_verification_summary.json is the unchanged record of that earlier numerical audit. Its statement that statistical methods were not recovered describes the audit before the author's later code submission. The newly supplied method is documented here and in source_manifest.json; no numerical audit values were replaced.

source_manifest.json identifies the four original expression matrices and their checksums, and records the later author-supplied script and its checksum. Matrix paths are relative to the original Radiation project root. The source matrices, original Seurat object and cell-level expression values are not redistributed here. These records support single-cell log(1 + counts per 10,000) RNA expression and the stated plotting method; cell-level comparisons do not quantify variation among independent animals.

Official implementation references
https://github.com/kassambara/ggpubr/blob/v0.6.0/R/stat_compare_means.R
https://github.com/kassambara/ggpubr/blob/v0.6.1/R/stat_compare_means.R
https://github.com/const-ae/ggsignif/blob/0.6.4/R/geom_signif.R
https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html
