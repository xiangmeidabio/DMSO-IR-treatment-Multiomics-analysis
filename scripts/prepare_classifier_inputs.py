#!/usr/bin/env python3
"""Prepare original classifier-compatible inputs with the archived cell split.

Input is the original, unstandardized cells-by-genes expression CSV; this script
does not infer a normalization method or derive gene panels from treatment labels.
"""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", type=Path, required=True,
                        help="CSV with SampleID, label (0,1,2,3), and gene columns")
    parser.add_argument("--genes", type=Path, required=True,
                        help="Ordered one-symbol-per-line gene panel")
    parser.add_argument("--split", type=Path, default=Path("resources/classifier_cell_split.csv"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    import pickle
    import numpy as np
    import pandas as pd
    from sklearn.preprocessing import StandardScaler

    output_names = ["metadata.pkl", "preparation.json"] + [
        f"{split}_data.{suffix}" for split in ("train", "val", "test")
        for suffix in ("pkl", "csv")]
    if any((args.out / name).exists() for name in output_names):
        raise SystemExit("Output files already exist; choose a new directory to preserve the previous run.")
    genes = [g.strip() for g in args.genes.read_text().splitlines() if g.strip()]
    if len(genes) not in (9, 100) or len(set(genes)) != len(genes):
        raise ValueError("Expected a unique ordered 9- or 100-gene panel.")
    columns = pd.read_csv(args.expression, nrows=0).columns.tolist()
    id_column = columns[0]
    needed = [id_column, "label", *genes]
    absent = set(needed) - set(columns)
    if absent:
        raise ValueError(f"Missing required columns: {sorted(absent)}")
    data = pd.read_csv(args.expression, usecols=needed, index_col=id_column)
    data.index = data.index.astype(str)
    if not data.index.is_unique:
        raise ValueError("Duplicate input cell IDs.")
    missing = set(genes) - set(data.columns)
    if missing:
        raise ValueError(f"Missing panel genes: {sorted(missing)}")
    partition = pd.read_csv(args.split, dtype={"cell_id": str})
    if not partition.cell_id.is_unique or len(partition) != 23006:
        raise ValueError("Expected 23,006 unique cells in the archived split.")
    if set(partition.cell_id) != set(data.index):
        raise ValueError("Input cell IDs do not exactly match the archived split.")
    if "label" not in data:
        raise ValueError("Input CSV must include the original numeric label column.")
    labels = pd.to_numeric(data["label"], errors="raise")
    aligned_labels = labels.loc[partition.cell_id].to_numpy()
    if not np.array_equal(aligned_labels, partition.condition_code.to_numpy()):
        raise ValueError("Condition labels disagree with the archived split.")
    expected = {"train": 16104, "val": 4601, "test": 2301}
    split_names = {"training": "train", "validation": "val", "test": "test",
                   "train": "train", "val": "val"}
    partition["short_split"] = partition["split"].map(split_names)
    if partition.short_split.isna().any():
        raise ValueError("Unexpected split labels.")
    prepared = {}
    for name, count in expected.items():
        rows = partition.loc[partition.short_split == name].sort_values("row_in_saved_split_1based")
        if len(rows) != count or rows.row_in_saved_split_1based.tolist() != list(range(1, count + 1)):
            raise ValueError(f"Invalid saved order/count for {name}.")
        names = rows.cell_id.to_numpy()
        X = data.loc[names, genes].to_numpy(dtype=float)
        if not np.isfinite(X).all():
            raise ValueError("Expression must contain finite numeric values.")
        y = pd.Series(rows.condition_code.astype(str).to_numpy())
        prepared[name] = (X, y, names)
    scaler = StandardScaler().fit(prepared["train"][0])
    args.out.mkdir(parents=True, exist_ok=True)
    for name, (X, y, names) in prepared.items():
        scaled = scaler.transform(X)
        obj = {f"X_{name}": scaled, f"y_{name}": y, f"names_{name}": names}
        with (args.out / f"{name}_data.pkl").open("wb") as handle:
            pickle.dump(obj, handle)
        frame = pd.DataFrame(scaled, columns=genes)
        frame.insert(0, "SampleID", names)
        frame["label"] = y.to_numpy()
        frame.to_csv(args.out / f"{name}_data.csv", index=False)
    metadata = {"feature_names": genes, "available_genes": genes, "scaler": scaler,
                "class_distribution": labels.astype(str).value_counts().to_dict()}
    with (args.out / "metadata.pkl").open("wb") as handle:
        pickle.dump(metadata, handle)
    (args.out / "preparation.json").write_text(json.dumps({
        "panel_size": len(genes), "split_cells": expected,
        "standardization": "StandardScaler fitted on archived training cells only",
        "feature_selection": "Fixed panel selected before classifier splitting",
        "input_scale": "Supplied original expression scale; no new normalization applied"
    }, indent=2) + "\n")
    print(f"Prepared {len(genes)} features with archived train/validation/test order.")


if __name__ == "__main__":
    main()
