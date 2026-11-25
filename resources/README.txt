Reproducibility materials for the DMSO irradiation study
Prepared: 25 November 2026

This directory contains gene lists, cell-split metadata, fixed scoring definitions, the database extracted from saved CellChat objects, archived ATAC-seq motif-analysis resources, Figure 6h descriptive-expression verification, and selected original scripts. It does not contain complete expression matrices, trained model weights or a complete executable analysis environment. No models were retrained to prepare these files, and no reported performance values or model rankings were changed.

FILES

response_genes_100.txt
The ordered 100-gene panel used by the saved four-condition classifiers. One gene symbol per line, without a header. The list matches the saved training, validation and test feature columns and the independently reconstructed feature-selection list. It is also the reference panel for the 100-gene Squidiff analysis. Original gene symbols are retained.

HSC_module3_genes_130.txt
The 130 genes assigned to module 3 in the saved 387-gene HSC analysis. One gene symbol per line, without a header. Membership was extracted from the saved module assignments without reclustering or applying new thresholds. This is an expression module, not a functionally validated HSC-protection signature.

classifier_genes_9.txt
The ordered final nine-gene classifier panel used in Figure 6h,i: Rps2, Uqcrb, Gapdh, Impa2, Ptma, Prtn3, Ramp1, Snrpf and Cenpw. This list is distinct from an earlier nine-gene panel retained in the original analysis directory.

classifier_cell_split.csv
The 23,006 cell identifiers and condition labels in the saved classifier splits. Columns are cell_id, split, row_in_saved_split_1based, condition_code and condition. The row number preserves the original order within each split. No expression values or predictions are included. Codes are 0 = Ctrl, 1 = IR, 2 = DMSO and 3 = DMSO-IR.

source_inventory.txt
Locations and SHA-256 checksums of the source lists, saved split files and analysis scripts. All paths are relative to the Radiation project root. These paths locate the original files in the project and are not download links. Only the resources and scripts described below are copied into this directory.

SHA256SUMS.txt
SHA-256 checksums for the distributed files. File paths are relative to this directory.

fig6h/
Descriptive verification of Figure 6h against archived RNA expression. The alignment table compares 180 values recovered from the PDF with expression-derived summaries: nine genes x four conditions x five statistics (first quartile, median, third quartile, lower whisker and upper whisker). The maximum absolute discrepancy was 0.0007388472 on the plotted scale. The correspondence supports identification of the plot values as single-cell log(1 + counts per 10,000) RNA expression. Quartiles use linear interpolation, and whiskers identify the most extreme observed values within 1.5 interquartile ranges of the box. The folder contains only summary statistics and provenance records, not the full expression matrices or cell-level expression values. The original plotting command, significance tests and multiple-testing adjustment were not recovered by this verification.

atac/
Archived differential-peak tables, the exact increased/decreased region selections and background used for Figure 2e,f, HOMER known-motif results and parameters, and library-processing evidence. The accompanying README and validation summary distinguish the verified motif inputs from the unverified differential-testing design. Large BAM files and the complete per-library peak tables are not included.

reanalysis/Result2/resources/
fixed_program_gene_lists.json contains the 12 fixed program gene sets. score_definitions_and_matched_controls.json records the exact expressed target genes, matched control genes and random seeds used for the scores. download_manifest.json records upstream resource URLs and checksums, including the MSigDB mouse 2025.1 collections. The full GMT collections are not included.

reanalysis/integration_new/cellchat/resources/
original_saved_CellChatDB.mouse.rds is the database extracted from the original saved CellChat objects. It has 3,379 interactions and 337 complexes and was identical in the four saved condition objects. It differs from the separately downloaded database snapshot. saved_object_options.json records the parameters stored in those objects. These files preserve the database and settings used for the saved analysis; the package does not infer a software version from the database.

Selected source scripts retain their original relative paths:
reanalysis/Result2/analyze_programs.py
reanalysis/Result2/check_depth_and_complexity.py
reanalysis/integration_new/cellchat/audit_saved_cellchat.R
They are unmodified copies. The Python scripts implement program scoring and UMI-depth checks. The R script extracts and verifies saved CellChat outputs; it does not rerun CellChat inference. These scripts require the original project inputs, including expression data and saved analysis objects, which are not included. They may also require intermediate tables and output directories from earlier analysis steps. They have not been executed as a standalone package.

classifier_training_protocols.csv
Model-specific fitting pools, tuning procedures, validation use and common-test evaluation. The table retains all 22 source configurations in their original order. Twenty were included in Figure 6b. The additional GBM configuration and LightGBM are identified explicitly. Table row order is not a performance ranking.

CLASSIFIER DESIGN

The saved training, validation and test splits contain 16,104, 4,601 and 2,301 cells, respectively, with no overlapping cell identifiers. Their condition counts are:

Split       Ctrl    IR      DMSO    DMSO-IR
Training    4068    3754    3676    4606
Validation  1162    1073    1050    1316
Test         581     537     525     658

The preprocessing code uses stratified random splitting with random_state = 42. StandardScaler was fitted to the 16,104 training cells and then applied to all three splits. Some model-specific routines combined training and validation cells for fitting and cross-validation; others used the training cells alone or used validation AUC for model selection, as detailed in Additional_file_3.csv. Scaling preceded model-specific cross-validation and was not refitted separately in each cross-validation fold.

The 100-gene panel was selected before the classifier split. The upstream selection involved 4,890 Ctrl cells, including 490 cells later assigned to the classifier test split. Thus, the classifier comparison evaluates a fixed, previously selected panel within this dataset; it is not an evaluation of a fully nested feature-selection pipeline. The separate subcluster-grouped regression analysis repeated feature selection within each fold and addresses a different prediction task.

All 20 models in Figure 6b were evaluated on the same 2,301 test cells. The additional saved GBM configuration used the same test cells. No saved LightGBM test result was available. The final nine-gene classifier also used the same test cells. Comparisons are internal and cell-level: cells originate from one pooled scRNA-seq sample per condition, and the split does not establish performance in independent animals or cohorts. The available records do not establish a model ranking selected solely on validation data before test-set comparison.

For the neural-network classifiers, the source code uses validation AUC for epoch selection but stores a state dictionary without a deep copy. Preservation of the best-epoch weights was not verified. The protocol table retains this limitation and does not relabel saved results as newly trained or externally validated results.

SQUIDIFF PROVENANCE

The source inventory identifies the preparation notebook, inference notebook and local training entry point for the 100-gene branch. The classifier panel above is the reference gene order. Saved notebook outputs support the 100-gene branch and its relation to the manuscript data, but they do not themselves contain the complete 100-gene list. This package does not claim to include the original trained checkpoint or a successful training log. The presence of the recorded training command alone does not establish that execution completed.

The inference analysis uses a 60-dimensional semantic representation, a Ctrl-to-IR mean latent shift applied to DMSO cells, and a comparison of 5,251 virtual cells with 6,580 measured DMSO-IR cells. The PCA comparison is in joint semantic space. Generated-expression differences are model-relative observations and do not alone establish biological non-additivity.

REMAINING REPRODUCIBILITY GAPS

The reported 84-gene measured-versus-virtual selection is retained as an author-confirmed historical result. The earlier manuscript records positive log2 fold change and adjusted P < 0.01 as the selection criteria. On 25 November 2026, the author confirmed that the complete original result table was no longer available following disk formatting. This package therefore does not include that table or claim an independent reconstruction of all 84 gene identities. The original 14-gene SHAP selection list and its exact selection rule also remain unavailable and have not been confirmed. The supplied final nine-gene list documents the saved classifier panel; it does not by itself reconstruct the intermediate-set intersection.

The complete original Squidiff training log and the corresponding trained checkpoint are not included. Recorded commands and notebook outputs do not supply a complete training history.

For ATAC-seq, the original differential-testing sample sheet, peak-by-library count matrix, normalization, design and contrast, backend settings, and software session record have not been recovered. In particular, the archived differential table does not establish how the two technical libraries per pooled sample were handled in the model. Its reported FDR values should not be interpreted as estimates of reproducibility across independent animals or biological pools.

These materials preserve verified inputs, selected source code and descriptive checks. They do not constitute a complete executable pipeline or independent biological validation.

USE

Preserve gene order when matching saved model inputs. Join cell identifiers exactly, including their suffixes, and retain the supplied condition-code mapping. To reconstruct the original analysis, obtain the corresponding expression objects, trained models and software environment separately, then follow the source scripts with their original inputs. The gene lists and split file alone are insufficient to rerun the complete workflow.

qc/
Sequencing quality-control tables: author-supplied Cell Ranger summary (four scRNA-seq conditions), eight ATAC library duplicate/insert-size reports, and two archived vendor raw-read reports. See QC_table_notes.txt for units, provenance and distinctions between estimated cells and the final retained analysis cells.
