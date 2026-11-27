# Inputs, environments and execution

## Software records

The analyses use separate R and Python environments. The table distinguishes recorded package versions from dependencies identified in the source code.

| Component | Available record | Meaning |
|---|---|---|
| Gradient boosting | scikit-learn 1.2.2 in the archived analysis evidence | Version used for the manuscript's saved gradient-boosting model |
| Original classifier scripts | NumPy, pandas, SciPy, scikit-learn, scikit-optimize, matplotlib, seaborn; optional torch/torchvision, xgboost, lightgbm | Imports in the preserved source; not all original package versions were recovered |
| SHAP | `KernelExplainer`, `kmeans` and first 100 validation cells | Exact package version was not established; historical code expects the pre-0.45 list-per-class SHAP output convention |
| Local Squidiff directory | NumPy 1.24.4, pandas 2.0.3, Scanpy 1.9.8, AnnData 0.9.2, SciPy 1.10.1, torch 2.4.1, Squidiff 1.0.8 listed in its requirements file | A local software record, not proof of the exact original training environment; the local source is vendored with hashes |
| ATAC alignment/peak calling | Bowtie2 2.3.3.1 and MACS2 2.2.7.1 | Recovered from BAM headers and peak-call headers; see `resources/atac/` |

Use Python 3.10 for an environment compatible with these older packages. A fresh compatible environment is a new execution environment and cannot guarantee byte-identical stochastic results. Install only dependencies needed for the chosen analysis. The classifier protocol table documents differences in fitting pool and selection among models.

## 1. Fixed classifier inputs

Required input: the original cells-by-genes expression CSV `All_Groups_Combined_transposed.csv`, with cell IDs in its first column, gene-symbol columns, and the original numeric `label` (0=Ctrl, 1=IR, 2=DMSO, 3=DMSO-IR). The code does not guess the expression scale or substitute raw counts for the original classifier matrix.

From the repository root:

```bash
python scripts/prepare_classifier_inputs.py \
  --expression /path/to/All_Groups_Combined_transposed.csv \
  --genes resources/response_genes_100.txt \
  --out analysis/classification/results
```

This uses the **archived cell identities and within-split ordering**, validates all labels and fits `StandardScaler` to training cells only. It fails if any panel gene is missing. No feature selection or new random partitioning is performed. The original `1-data/05_data_preprocessing.py` is also preserved, but new users should use the fixed-split entry point to avoid input-row-order differences.

To run the saved gradient-boosting configuration (a new computational run):

```bash
MPLBACKEND=Agg python 'analysis/classification/2-model codes/06_gradient_boosting_model.py'
```

The classifier performs 30 Bayesian search iterations with three-fold cross-validation and accuracy scoring. The original training-set standardization precedes those folds. This preserves the recorded analysis; it is not a fully nested feature-selection evaluation. Other model scripts use different fitting/selection procedures: consult `resources/classifier_training_protocols.csv` before running them. Running all scripts can be expensive and was not required to prepare this repository.

For the final nine-gene model, use the same expression input, changing `--genes` to `resources/classifier_genes_9.txt` and `--out` to `analysis/classification/7-best_model_9genes`, then run `2_gradient_boosting_model_key_genes.py` there. The input-preparation script allows existing source files in that directory but refuses to replace existing data/output files. Do not overwrite original saved analyses.

## 2. SHAP

The preserved `analysis/classification/4-SHAP/SHAP.py` loads the 100-gene gradient-boosting model and standardized data from the sibling `results` directory. It uses 50 k-means background centers and the first 100 validation cells; its displayed waterfall is for the first cell. For a fresh run, use a SHAP version with the list-per-class return convention expected by this archived script, or adapt the returned dimensions explicitly. The [archived selection records](../resources/shap_selection/) contain the example cell's 100 local SHAP values and the 14 positive contributors among the 19 individually displayed genes for the DMSO-IR class. This selection uses the local waterfall display and sign rule, separately from global mean-absolute SHAP importance. These records can be inspected without rerunning SHAP.

Only load locally trusted pickle/model files. Model weights and pickle data are not published in this repository.

## 3. Squidiff

See [the dedicated instructions](../analysis/squidiff/README.md). Count normalization for Squidiff is documented explicitly and is separate from the classifier's archived standardized feature matrices. Both analyses use the same ordered 100-gene panel.

## 4. R analyses and supporting resources

R scripts are grouped under `analysis/cell_atlas/`, `analysis/radiation_response/`, `analysis/hsc/` and `analysis/transcriptional_sensitivity/`. They use Seurat, Monocle and Augur objects from the original analysis session. Set the input paths and working directory for the corresponding local data before use. The repaired ATAC motif plot is in `analysis/atac/`. `resources/README.txt` lists the saved inputs available in this release and those obtained separately. The program-scoring and CellChat audit scripts there retain relative paths from the original project, including intermediate objects not distributed here.

## Verification performed for this release

The dependency-free repository checker validates resource counts, labels, ordering, hashes and Python ASTs. CLI help was checked for the portable entry points. A synthetic 23,006-row, nine-feature input check recovered the archived order after shuffling input rows, confirmed training-only zero-mean/unit-variance scaling and confirmed that existing outputs cannot be overwritten. This test did not use or regenerate study expression results. Source review compared the portable preprocessing/inference steps with the archived notebook and model code. All retained R sources parse. The repaired ATAC plot replaces the old copy; a damaged MELD R file is retained only in Git history, with the separately saved Python analysis provided for reference. Expensive model training, GPU inference, R analysis execution and exact historical-output regeneration were not performed for this release.
