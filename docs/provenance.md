# Source notes

The code combines the original analysis scripts and supporting resources for the current manuscript. R scripts are grouped by analysis and named by function; their calculations are unchanged. The [file history](provenance/file_moves.json) maps the earlier figure-based names to the current paths.

## Source records

- [Source manifest](provenance/source_files.json): original source locations, checksums and preparation notes.
- [File checksums](provenance/SHA256SUMS.txt): SHA-256 values for the distributed files.
- [Resource inventory](../resources/source_inventory.txt): saved analysis objects and inputs.
- [Resource descriptions](../resources/README.txt): gene lists, cell partitions, scoring definitions and archived results.

Classification scripts retain the recorded fitting procedures. Notebook excerpts contain the study's 100-gene Squidiff workflow, with outputs removed and personal paths replaced. The command-line input and inference scripts add explicit paths and checks for cell identity, gene order and output locations.

## Script maintenance

The ATAC motif plot uses the repaired script at `analysis/atac/motif_enrichment.R`. Invalid Windows path escapes were replaced with paths to the supplied HOMER tables, and the axis label was corrected to −ln(P). The plotted values and selection criteria are unchanged.

The earlier MELD R file contained binary damage and has been removed from the working tree. The separately saved Python analysis is provided at `analysis/transcriptional_sensitivity/meld.py`; it is not a reconstruction of the damaged R plotting code. Both removed R originals remain available in [the preceding commit](https://github.com/xiangmeidabio/DMSO-IR-treatment-Multiomics-analysis/tree/fad90ec8a11590fc62b6b3d6dcbaa7cb2dc48a31).

## Intermediate files and analysis scope

The 100-gene panel, final nine candidates, 130-gene HSC module and exact classifier cell partitions are included. The [SHAP selection records](../resources/shap_selection/) preserve all 100 local values for the example cell and the 14 positive contributors among its 19 individually displayed genes. The selection rule follows the author clarification and archived values; no SHAP values were recomputed. The complete historical 84-gene measured-versus-predicted table remains unavailable. Its reported selection used positive log2 fold change and adjusted P < 0.01; that table and the full intermediate-set intersection have not been independently reconstructed.

The [author-supplied Figure 6h plotting script](../resources/fig6h/plot_Fig6h_author_supplied.R), received on 26 November 2026, is included with its [method and provenance notes](../resources/fig6h/README.txt). It requires the original Seurat object. Its statistical and plotting settings were preserved, and the script was not executed to generate new results. Earlier numerical audit records are retained unchanged; statements about an unrecovered plotting method in those records describe the earlier audit.

The original 100-gene Squidiff checkpoint and complete training log are not included. The recorded training schedule is provided, and Supplementary Fig. S7e remains part of the manuscript. Raw sequencing files and processed sequencing data are deposited in GEO; full analysis matrices and model weights are not bundled with this code.

Each condition has one pooled biological sample per assay. Its two ATAC libraries were prepared from the same pool. The available ATAC peak and HOMER tables are included; the original differential-testing design matrix and session are unavailable. Classifier evaluation uses a cell-level split, with the gene panel selected before partitioning. Fitting and validation procedures differ among the original classifiers and are specified in the [protocol table](../resources/classifier_training_protocols.csv).

## Validation

File hashes, gene panels, cell partitions, Python syntax and stripped notebooks were checked. A synthetic input check verified cell-order alignment, training-set standardization and protection against overwriting outputs. All retained R sources parse. These checks do not rerun model training, GPU inference or the original R analyses.
