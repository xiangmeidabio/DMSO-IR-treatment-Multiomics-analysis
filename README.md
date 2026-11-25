# DMSO and hematopoietic radioprotection

Code and analysis resources for **Single-cell analysis identifies expression programs and candidate genes for DMSO-mediated hematopoietic radioprotection**.

The study examines lineage-negative-enriched mouse bone marrow under four conditions: Ctrl, irradiation (IR), DMSO, and DMSO pretreatment followed by irradiation (DMSO-IR). Single-cell RNA-seq, supported by bulk ATAC-seq, characterizes radiation responses and expression programs retained after pretreatment. Classification, SHAP and Squidiff analyses identify candidate genes for further functional investigation.

## Analyses

| Workflow | Code |
|---|---|
| Cell populations and composition | [Cell atlas](analysis/cell_atlas/) |
| Radiation responses and DMSO-associated expression changes | [Expression analysis](analysis/radiation_response/) · [ATAC motif enrichment](analysis/atac/) |
| HSC stemness, lineage expression and gene modules | [HSC analysis](analysis/hsc/) |
| Expression diversity and transcriptional sensitivity | [Sensitivity analysis](analysis/transcriptional_sensitivity/) |
| Treatment classification and SHAP | [Classification](analysis/classification/) |
| Predicted and measured combined-treatment expression | [Squidiff](analysis/squidiff/) |

See the [figure guide](docs/figure_map.md) for links between manuscript figures, scripts and supporting resources.

## Data and use

Sequencing data: [GSE353967 — scRNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353967) and [GSE353968 — bulk ATAC-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353968). The GEO records are private until release; reviewer access is provided with the manuscript. Gene panels, cell partitions and [sequencing QC summaries](resources/qc/) are provided in `resources/`.

[Execution instructions](docs/reproduction.md) describe inputs and software requirements. [Source notes](docs/provenance.md) document the archived analyses and available intermediate files. To check the distributed files:

```bash
python scripts/verify_repository.py
```

Please cite the associated study when using these resources. Squidiff attribution and licensing are provided in the [vendor notice](analysis/squidiff/vendor/NOTICE.md).
