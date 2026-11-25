#!/usr/bin/env python3
"""Portable transcription of the 100-gene inference branch, not archived results.

A compatible trusted checkpoint is required and is not distributed here.
Newly generated results are not replacements for the historical 84-gene table.
"""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="predict_targeted.h5ad")
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--panel", type=Path, default=Path("resources/response_genes_100.txt"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=64)
    args = parser.parse_args()
    if not args.checkpoint.is_file():
        raise SystemExit("Provide the original compatible 100-gene checkpoint; it is not included.")
    if args.out.exists() and any(args.out.iterdir()):
        raise SystemExit("Output directory must be empty.")
    if args.batch_size < 1:
        raise ValueError("Batch size must be positive.")
    sys.path.insert(0, str(Path(__file__).parent / "vendor"))
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    from sample_squidiff import sampler

    adata = sc.read_h5ad(args.input)
    panel = args.panel.read_text().splitlines()
    if len(panel) != 100 or adata.var_names.tolist() != panel:
        raise ValueError("Input gene order must exactly match the archived 100-gene panel.")
    expected = {"Ctrl": 5811, "IR": 5364, "DMSO": 5251, "DMSO-IR": 6580}
    if adata.obs.Group.value_counts().to_dict() != expected:
        raise ValueError("Input groups do not match the manuscript dataset.")
    generator = sampler(str(args.checkpoint), gene_size=100, output_dim=100,
                        use_drug_structure=False)
    device = next(generator.model.parameters()).device
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    if not np.isfinite(matrix).all():
        raise ValueError("Input expression must be finite.")
    chunks = []
    with torch.no_grad():
        for start in range(0, len(adata), args.batch_size):
            out = generator.model.encoder(torch.as_tensor(
                matrix[start:start + args.batch_size], dtype=torch.float32, device=device))
            if isinstance(out, (tuple, list)):
                out = out[0]
            chunks.append(out.detach().cpu().numpy())
    semantic = np.concatenate(chunks)
    if semantic.shape != (23006, 60):
        raise ValueError(f"Expected 23,006 x 60 semantic matrix; received {semantic.shape}.")
    group = adata.obs.Group.to_numpy()
    delta = semantic[group == "IR"].mean(axis=0) - semantic[group == "Ctrl"].mean(axis=0)
    virtual_semantic = semantic[group == "DMSO"] + delta
    generated = []
    with torch.no_grad():
        for start in range(0, len(virtual_semantic), args.batch_size):
            z = torch.as_tensor(virtual_semantic[start:start + args.batch_size],
                                dtype=torch.float32, device=device)
            generated.append(generator.pred(z, gene_size=100).detach().cpu().numpy())
    virtual_expression = np.concatenate(generated)
    measured = matrix[group == "DMSO-IR"]
    combined = sc.AnnData(np.vstack([measured, virtual_expression]))
    combined.var_names = panel
    combined.obs_names = (["measured_" + c for c in adata.obs_names[group == "DMSO-IR"]]
                          + ["virtual_" + c for c in adata.obs_names[group == "DMSO"]])
    combined.obs["Group"] = (["Measured DMSO-IR"] * len(measured)
                              + ["Virtual DMSO-IR"] * len(virtual_expression))
    sc.tl.rank_genes_groups(combined, "Group", groups=["Measured DMSO-IR"],
                            reference="Virtual DMSO-IR", method="wilcoxon")
    table = sc.get.rank_genes_groups_df(combined, group="Measured DMSO-IR")
    selected = table.loc[(table.logfoldchanges > 0) & (table.pvals_adj < 0.01)]
    pca = sc.AnnData(np.vstack([semantic[group == "DMSO-IR"], virtual_semantic]))
    pca.obs_names = combined.obs_names
    pca.obs["Group"] = combined.obs["Group"].to_numpy()
    sc.tl.pca(pca)
    args.out.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out / "new_measured_vs_virtual_gene_statistics.csv", index=False)
    selected.to_csv(args.out / "new_positive_log2FC_adjustedP_lt_0p01.csv", index=False)
    combined.write_h5ad(args.out / "new_measured_and_virtual_expression.h5ad")
    pca.write_h5ad(args.out / "new_joint_semantic_pca.h5ad")
    np.savez_compressed(args.out / "new_semantic_states.npz", measured=semantic,
                        virtual=virtual_semantic, radiation_shift=delta)
    (args.out / "run_description.json").write_text(json.dumps({
        "run_type": "New execution of archived inference logic; not historical recovered results",
        "panel_size": 100, "semantic_dimensions": 60, "batch_size": args.batch_size,
        "measured_DMSO_IR_cells": len(measured), "virtual_DMSO_IR_cells": len(virtual_expression),
        "diffusion_steps": generator.arg["diffusion_steps"], "use_ddim": generator.arg["use_ddim"],
        "random_seed": "Not set in archived notebook; fresh generation is stochastic",
        "comparison": "Measured versus virtual, Wilcoxon, positive log2FC and adjusted P < 0.01",
        "SHAP_intersection": "Not reconstructed: original 14-gene selection rule is not available"
    }, indent=2) + "\n")
    print(f"New run selected {len(selected)} genes. This count is not asserted to reproduce 84.")


if __name__ == "__main__":
    main()
