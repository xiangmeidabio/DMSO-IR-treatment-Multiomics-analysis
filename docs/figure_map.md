# Figure guide

The table links the current manuscript figures to their analysis code and supporting resources. R scripts retain the original computations; their filenames have been changed to describe the analysis rather than earlier figure numbering. Input objects and working-directory requirements are described in the [execution instructions](reproduction.md).

| Figure | Analysis | Code and resources |
|---|---|---|
| Fig. 1 | Cell populations, composition and overall treatment responses | [Cell atlas](../analysis/cell_atlas/): preprocessing, annotation, marker expression, composition and sample correlation |
| Fig. 2 | Radiation-related expression changes and ATAC motif enrichment | [Radiation-response scripts](../analysis/radiation_response/); [program-scoring resources](../resources/reanalysis/Result2/); [motif plotting](../analysis/atac/motif_enrichment.R) and [ATAC inputs](../resources/atac/) |
| Fig. 3 | Treatment comparisons, expression preservation and inferred communication | [Augur](../analysis/radiation_response/augur.R); [program-scoring resources](../resources/reanalysis/Result2/); [CellChat database, settings and saved-object audit](../resources/reanalysis/integration_new/cellchat/) |
| Fig. 4 | HSC states, stemness, lineage expression and coordinated gene modules | [HSC scripts](../analysis/hsc/); [130-gene HSC module](../resources/HSC_module3_genes_130.txt) |
| Fig. 5 | Expression diversity, transcriptional sensitivity and response-associated gene selection | [Sensitivity scripts](../analysis/transcriptional_sensitivity/); [100-gene panel](../resources/response_genes_100.txt) |
| Fig. 6a–d | Treatment classification, within-cell-type evaluation and SHAP | [Classification workflow](../analysis/classification/); [training protocols](../resources/classifier_training_protocols.csv); [cell partitions](../resources/classifier_cell_split.csv) |
| Fig. 6e,f | Predicted and measured DMSO-IR expression | [Squidiff workflow](../analysis/squidiff/): transfer of the Ctrl-to-IR semantic shift to DMSO cells, joint PCA and expression comparison |
| Fig. 6g | Candidate-set intersection | [Nine-gene panel](../resources/classifier_genes_9.txt); see [source notes](provenance.md) for availability of the intermediate gene sets |
| Fig. 6h | Candidate-gene expression across treatment conditions | [Archived expression summaries and plot-value checks](../resources/fig6h/) |
| Fig. 6i | Nine-gene treatment classifier | [Nine-gene classification scripts](../analysis/classification/7-best_model_9genes/) |
| Supplementary Fig. S7e | Squidiff training-loss curve | [Recorded training configuration](../analysis/squidiff/); the complete original training log is not included |

The CellChat audit operates on saved objects. RNA-velocity commands and some final figure-layout steps are not included in this code collection. See [source notes](provenance.md) for the scope of the archived material.
