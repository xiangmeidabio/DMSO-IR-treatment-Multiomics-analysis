# Squidiff: manuscript 100-gene branch

This directory documents the analysis used as the model-based reference in Fig. 6e,f. It uses the same ordered 100 response-associated genes as the treatment classifier. The 2,000-gene example branch present elsewhere in the old notebook is deliberately not included in the curated excerpts.

The recorded procedure is:

1. Normalize filtered raw counts across all genes to 10,000 counts per cell and log1p-transform; then subset to the ordered 100-gene panel.
2. Train on Ctrl (5,811 cells) and IR (5,364 cells). The recorded command configures batch size 128, learning rate 1e-4 and checkpoint saving every 5,000 steps. Local defaults configure AdamW, weight decay 0 and linear learning-rate decay over 100,000 steps.
3. Encode all 23,006 cells into 60-dimensional semantic states. Add the mean Ctrl-to-IR semantic shift to each of the 5,251 DMSO states.
4. Generate virtual expression using DDIM with 1,000 diffusion steps, in batches of 64. Compare it with measured DMSO-IR expression from 6,580 cells using Scanpy's Wilcoxon test.
5. Construct the PCA jointly from measured and virtual **semantic states**, not directly from the expression matrix.

The shift is applied in learned semantic space and is decoded through the diffusion model. It is not a simple expression-space additive baseline, and a measured-versus-predicted difference does not by itself establish biological non-additivity.

## Input preparation

From the repository root:

```bash
python analysis/squidiff/prepare_inputs.py \
  --matrix /path/to/matrix.mtx \
  --genes-tsv /path/to/genes.tsv \
  --barcodes-tsv /path/to/barcodes.tsv \
  --metadata /path/to/metadata.csv \
  --out runs/squidiff100/input
```

`matrix.mtx` must contain genes-by-cells **raw counts from the retained cells**. `genes.tsv` contains symbols in the first column, not a different Ensembl-ID mapping; barcodes match matrix columns. `metadata.csv` is indexed by those cell IDs and has a `Group` column with the exact four condition names. The script verifies all group counts and all 100 panel genes. It writes `train_targeted.h5ad` and `predict_targeted.h5ad`.

## Recorded training command

The original notebook command was:

```bash
python train_squidiff.py \
  --logger_path logs/mouse_targeted_genes \
  --data_path train_targeted.h5ad \
  --resume_checkpoint mouse_targeted_100 \
  --gene_size 100 --output_dim 100 \
  --batch_size 128 --save_interval 5000 --lr 1e-4
```

The local training entry point and library are preserved under `vendor/`. In this local `TrainLoop`, initial checkpoint loading is commented out and the `resume_checkpoint` argument is used as a directory for saving checkpoints. Thus the directory name in the historical command must not be interpreted as proof of a pretrained model or a resumed run. The source requires the corresponding CUDA/PyTorch environment; a fresh run is not an archived model artifact. The local configuration/log-plot code is preserved unchanged, including historical plot labels. The recorded command does not establish a completed 100,000-step training run.

To reproduce that configuration using the vendored source, from the repository root, provide absolute local paths:

```bash
python analysis/squidiff/vendor/train_squidiff.py \
  --logger_path /path/to/new-run/logs \
  --data_path /path/to/train_targeted.h5ad \
  --resume_checkpoint /path/to/new-run/checkpoints \
  --gene_size 100 --output_dim 100 \
  --batch_size 128 --save_interval 5000 --lr 1e-4
```

The original checkpoint is not bundled. Running this command is new training and may be computationally expensive.

## Inference and comparison

With a compatible, locally trusted 100-gene checkpoint:

```bash
python analysis/squidiff/infer_and_compare.py \
  --input runs/squidiff100/input/predict_targeted.h5ad \
  --checkpoint /path/to/original/model.pt \
  --out runs/squidiff100/new-comparison
```

This entry point follows the selected notebook branch but avoids duplicated interactive generation cells and validates the 60-dimensional semantic representation. It writes expression comparisons and PCA objects to a new directory. Filenames explicitly mark them as **new outputs**, not recovered historical results. No seed was recorded in the original inference notebook, so a fresh generation is stochastic. The filtered output follows positive log2FC and adjusted P < 0.01; it is not asserted to reproduce the unavailable historical 84-gene table or the intermediate SHAP intersection.

## Notebook excerpts and upstream source

`notebooks/1-process.ipynb` and `2-data.ipynb` are archival, output-free excerpts. They preserve historical diagnostics and selected source-cell indices for review, rather than presenting an automatically executable top-to-bottom workflow. The portable scripts above are the entry points for fresh execution.

The local dependency record listed Scanpy 1.9.8, AnnData 0.9.2, NumPy 1.24.4, pandas 2.0.3, SciPy 1.10.1 and torch 2.4.1; it does not establish the full historical training environment. See the [vendor notice](vendor/NOTICE.md) and [study-wide provenance](../../docs/provenance.md).
