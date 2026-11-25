# Source provenance and reproducibility scope

This release was assembled from the archived analysis project and the supporting resources accompanying the current manuscript on 25 November 2026. Root `Figure*.R` scripts remain unchanged. Added original classification scripts retain their existing computations, hyperparameters and file conventions. Original notebook excerpts retain only the manuscript's 100-gene branch, remove embedded outputs and execution counts, and replace personal absolute paths with configurable/relative paths. Their source-cell indices and original SHA-256 values are recorded.

- `provenance/source_files.json`: source paths relative to the original project, hashes before copying, destination paths, and transformation notes.
- `resources/source_inventory.txt`: additional source objects and verified inputs; entries are provenance records, not download links.
- `provenance/SHA256SUMS.txt`: checksums of distributed files, excluding itself and Git metadata.
- `resources/README.txt`: detailed interpretation of gene lists, cell splits, scoring and archived analyses.

The portable entry points are convenience transcriptions of recorded procedures. They add input checks, explicit output directories and stable cell/gene ordering; they do not add new study results. The provided final 100- and nine-gene panels document saved model inputs.

## Available and unavailable intermediates

- The final nine genes and 130-gene HSC module are available, together with exact 23,006-cell classifier partitions.
- The complete 84-gene measured-versus-virtual table was lost after disk formatting. The manuscript preserves the author-confirmed count and recorded positive log2FC/adjusted-P < 0.01 criteria. No substitute gene list is generated here.
- The exact intermediate 14-gene SHAP selection rule/list was not recovered. Global top-14 importance and an individual-cell waterfall are not interchangeable selection procedures; the repository does not assume either rule to rebuild the intersection.
- The original 100-gene Squidiff checkpoint and complete training log are not included. The recorded 100,000-step schedule is a configuration, not proof that all steps completed. The manuscript's existing Supplementary Fig. S7e is retained.
- Eight ATAC libraries were made from four condition-specific biological pools, two libraries per pool. The original differential-testing design/count matrix/session was not recovered; available differential-peak and HOMER tables are preserved with their documented origin.
- The classifier split is an internal cell-level comparison within the pooled dataset. It does not represent independent animal-level validation. The feature panel was selected before this split; some classifier routines additionally combine training and validation cells for tuning. See the protocol table for exact differences.

No raw sequence reads, full expression matrices, checkpoint files, private reviewer-access tokens or unpublished manuscript text are included. The repository does not claim a complete end-to-end reconstruction when required original inputs are absent.
