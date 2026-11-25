# Current manuscript figure-to-code map

Historical `Figure*.R` filenames have been preserved. The mapping below is at the analysis level: a matching historical filename is not a claim that it regenerates every panel of the current composite figure without its original inputs and plotting session.

| Current figure | Analysis | Relevant sources and resources |
|---|---|---|
| Fig. 1 | Cell populations, composition and overall treatment responses | Root `Figure1_1_Seurat.R` through `Figure1_5_Sample_Correlation.R` |
| Fig. 2 | Radiation-related expression changes, fixed programs and ATAC motif enrichment | Root `Figure2_*.R`; `resources/reanalysis/Result2/`; `resources/atac/` |
| Fig. 3 | Augur treatment comparisons, expression preservation and inferred communication | Root `Figure3_1_Augur.R`; fixed scoring resources under `resources/reanalysis/Result2/`; archived CellChat database/settings and audit under `resources/reanalysis/integration_new/cellchat/` |
| Fig. 4 | HSC states, stemness/lineage expression and coordinated gene modules | Root `Figure4_*.R` and `Figure5_*.R`; `resources/HSC_module3_genes_130.txt` |
| Fig. 5 | Baseline expression diversity, transcriptional sensitivity and response-associated gene selection | Root `Figure3_2_*` through `Figure3_5_*`; fixed `resources/response_genes_100.txt`; source locations/hashes in `resources/source_inventory.txt` |
| Fig. 6a–d | Four-condition classification, evaluation within cell types and SHAP interpretation | `analysis/classification/`; `resources/classifier_training_protocols.csv`; `resources/classifier_cell_split.csv` |
| Fig. 6e,f | Joint semantic PCA and measured-versus-virtual expression | `analysis/squidiff/`; 100-gene panel, Ctrl-to-IR latent shift applied to DMSO |
| Fig. 6g | Historical intermediate-set intersection | Final nine-gene panel and HSC module membership are distributed; the original complete 84-gene and intermediate 14-gene sets are not available |
| Fig. 6h | Descriptive expression of the nine candidates | `resources/fig6h/` contains archived expression summaries and plot-value verification, not a regenerated significance test |
| Fig. 6i | Final nine-gene condition classifier | `analysis/classification/7-best_model_9genes/`; `resources/classifier_genes_9.txt` |
| Supplementary Fig. S7e | Training-loss curve retained in the manuscript | Recorded training configuration in `analysis/squidiff/`; a complete original training log is not included |

RNA-velocity and some final figure-layout commands are not present in this code snapshot. CellChat audit code reads original saved objects; it is not a complete reconstruction of the inference workflow. These differences are documented rather than filled with new, unverified procedures.
