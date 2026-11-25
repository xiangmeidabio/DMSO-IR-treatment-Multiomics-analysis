#!/usr/bin/env python3
"""Prepare the manuscript's 100-gene Squidiff inputs from filtered raw counts.

All-gene normalization precedes panel selection, as in 1-process.ipynb.
"""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="Genes-by-cells raw-count Matrix Market file")
    parser.add_argument("--genes-tsv", type=Path, required=True, help="One symbol per row, first column, no header")
    parser.add_argument("--barcodes-tsv", type=Path, required=True, help="One cell ID per row, no header")
    parser.add_argument("--metadata", type=Path, required=True, help="CSV indexed by cell ID, with Group column")
    parser.add_argument("--panel", type=Path, default=Path("resources/response_genes_100.txt"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    import numpy as np
    import pandas as pd
    import scanpy as sc

    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Output directory must be empty.")
    panel = args.panel.read_text().splitlines()
    if len(panel) != 100 or len(set(panel)) != 100:
        raise ValueError("The ordered panel must contain 100 unique symbols.")
    adata = sc.read_mtx(args.matrix).T
    symbols = pd.read_csv(args.genes_tsv, sep="\t", header=None)[0].astype(str)
    cells = pd.read_csv(args.barcodes_tsv, sep="\t", header=None)[0].astype(str)
    if adata.shape != (len(cells), len(symbols)) or cells.duplicated().any():
        raise ValueError("Count matrix dimensions or cell IDs are inconsistent.")
    adata.var_names = symbols.to_numpy()
    adata.obs_names = cells.to_numpy()
    adata.var_names_make_unique()
    metadata = pd.read_csv(args.metadata, index_col=0)
    metadata.index = metadata.index.astype(str)
    if not metadata.index.is_unique or set(metadata.index) != set(cells):
        raise ValueError("Metadata and matrix cell IDs must match exactly.")
    adata.obs = metadata.loc[adata.obs_names].copy()
    expected = {"Ctrl": 5811, "IR": 5364, "DMSO": 5251, "DMSO-IR": 6580}
    if adata.obs["Group"].value_counts().to_dict() != expected:
        raise ValueError("Group labels/counts differ from the manuscript's filtered data.")
    values = adata.X.data
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.round(values)):
        raise ValueError("Expected finite, nonnegative raw integer counts.")
    if (np.asarray(adata.X.sum(axis=1)).ravel() <= 0).any():
        raise ValueError("Zero-count cells cannot be normalized.")
    missing = set(panel) - set(adata.var_names)
    if missing:
        raise ValueError(f"Missing panel genes: {sorted(missing)}")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    subset = adata[:, panel].copy()
    args.out.mkdir(parents=True, exist_ok=True)
    subset[subset.obs.Group.isin(["Ctrl", "IR"])].copy().write_h5ad(args.out / "train_targeted.h5ad")
    subset.write_h5ad(args.out / "predict_targeted.h5ad")
    (args.out / "preparation.json").write_text(json.dumps({
        "genes": panel, "group_counts": expected, "training_groups": ["Ctrl", "IR"],
        "normalization": "Normalize total counts to 10,000 across all genes; log1p; then select panel"
    }, indent=2) + "\n")
    print("Wrote 11,175 training cells and 23,006 prediction-input cells, each with 100 genes.")


if __name__ == "__main__":
    main()
