# DMSO pretreatment and hematopoietic responses to irradiation

Analysis code and supporting resources for **“Single-cell analysis identifies expression programs and candidate genes for DMSO-mediated hematopoietic radioprotection.”**

The study profiles lineage-negative-enriched mouse bone marrow under Ctrl, IR, DMSO and DMSO-IR conditions. Single-cell expression analysis characterizes radiation responses and expression programs retained after DMSO pretreatment. Bulk ATAC-seq provides complementary regulatory context. Treatment classification, SHAP and a Squidiff-generated expression reference prioritize candidate genes.

## Start here

- [Current manuscript figure-to-code map](docs/figure_map.md)
- [Inputs, environments and execution](docs/reproduction.md)
- [Scope and source provenance](docs/provenance.md)
- [Sequencing QC tables and source notes](resources/qc/QC_table_notes.txt)
- [Historical R source repairs](analysis/legacy_fixes/README.md)
- [Archived classifier training protocols](resources/classifier_training_protocols.csv)
- [Squidiff 100-gene workflow](analysis/squidiff/README.md)

The `Figure*.R` files at the repository root are retained historical analysis scripts. **Their filenames follow the earlier manuscript numbering.** They are analysis snippets with shared R objects and historical input paths, not a new standalone pipeline. Use the figure map for the current manuscript.

The November 2026 update adds the previously missing classification/SHAP code, the manuscript's **100-gene** Squidiff branch, fixed gene lists and cell partitions, scoring resources, and ATAC motif-analysis inputs. New command-line entry points validate input identities and gene order and write new runs to separate directories. No models or reported figure values were regenerated for this release.

## Repository layout

| Path | Contents |
|---|---|
| `Figure*.R` | Original R analysis snippets, preserved under their historical names |
| `analysis/classification/` | Original preprocessing, 22 classifier configurations, SHAP and final nine-gene scripts |
| `analysis/squidiff/` | Portable input preparation/inference, recorded notebook excerpts and local Squidiff implementation |
| `resources/` | Ordered 100-, 9- and 130-gene lists; exact cell splits; fixed scoring definitions; saved CellChat database; ATAC and expression-summary resources |
| `scripts/` | Fixed-split classifier input preparation and dependency-free repository checks |
| `provenance/` | Source hashes and a checksum manifest for this release |

## Data

Raw and processed sequencing data were deposited at NCBI GEO:

- [GSE353967](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353967): single-cell RNA-seq.
- [GSE353968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353968): bulk ATAC-seq.

The records are private pending public release. Reviewer access is provided through the manuscript/submission materials; access tokens are not included in this repository. Raw sequencing files, full expression matrices, model weights and unpublished manuscript files are not bundled here.

The final scRNA-seq dataset contains 23,006 cells (Ctrl 5,811; IR 5,364; DMSO 5,251; DMSO-IR 6,580). The classifier split is internal and cell-level: 16,104 training, 4,601 validation and 2,301 test cells. Each condition has one pooled biological sample per assay. The two ATAC libraries per condition were prepared from the same pool and are technical libraries.

## Validate the distributed files

From the repository root, with Python 3.10 or later:

```bash
python scripts/verify_repository.py
```

This checks gene lists, cell identities and ordering, condition labels, classifier protocol counts, Python syntax, stripped notebook outputs and SHA-256 hashes. It does not train models or claim independent biological replication.

## Interpretation and availability

The 100-gene panel was selected before classifier partitioning. Classification and SHAP identify expression features that distinguish the four conditions; they do not establish a protective effect of a gene. Squidiff transfers the Ctrl-to-IR shift in learned semantic space to DMSO cells and generates a model-based expression reference. Its output is not a simple addition of gene-expression values.

The historical 84-gene measured-versus-virtual selection is retained in the manuscript as an author-confirmed result (positive log2 fold change and adjusted P < 0.01). The complete original table is unavailable. The exact selection rule/list for the intermediate 14-gene SHAP set and the original 100-gene Squidiff checkpoint/training log are also not distributed. The final nine-gene classifier panel is provided; missing intermediate lists have not been recreated from that panel. See [provenance](docs/provenance.md).

## Citation and third-party software

Please cite the associated study when using its analysis resources; publication details will be added when available. Squidiff is described by [He et al., *Nature Methods*, DOI: 10.1038/s41592-025-02877-y](https://doi.org/10.1038/s41592-025-02877-y). The vendored Squidiff source retains its upstream BSD license and attribution; see [the vendor notice](analysis/squidiff/vendor/NOTICE.md). Other dependencies retain their own licenses. This repository update does not assign a new license to the authors' existing code.
