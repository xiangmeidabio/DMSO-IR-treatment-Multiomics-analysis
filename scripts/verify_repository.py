#!/usr/bin/env python3
"""Validate distributed resource integrity and source syntax without model fitting."""
import ast
import collections
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main():
    resources = ROOT / "resources"
    panels = {}
    for filename, count in [("response_genes_100.txt", 100),
                            ("classifier_genes_9.txt", 9),
                            ("HSC_module3_genes_130.txt", 130)]:
        panel = (resources / filename).read_text().splitlines()
        assert len(panel) == len(set(panel)) == count, filename
        assert all(g and g == g.strip() for g in panel), filename
        panels[filename] = panel
    assert set(panels["classifier_genes_9.txt"]).issubset(panels["response_genes_100.txt"])
    assert {"Uqcrb", "Ptma"}.issubset(panels["HSC_module3_genes_130.txt"])
    assert (ROOT / "analysis/classification/data/key_gene.csv").read_text().splitlines() == panels["response_genes_100.txt"]
    with (resources / "classifier_cell_split.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 23006 == len({r["cell_id"] for r in rows})
    counts = collections.Counter(r["split"] for r in rows)
    assert dict(counts) == {"train": 16104, "val": 4601, "test": 2301}, counts
    names = {"0": "Ctrl", "1": "IR", "2": "DMSO", "3": "DMSO-IR"}
    for split, count in counts.items():
        selected = [r for r in rows if r["split"] == split]
        assert [int(r["row_in_saved_split_1based"]) for r in selected] == list(range(1, count + 1))
    assert all(names[r["condition_code"]] == r["condition"] for r in rows)
    protocol = list(csv.DictReader((resources / "classifier_training_protocols.csv").open(encoding="utf-8-sig")))
    assert len(protocol) == 22
    assert sum(r["Included in Figure 6b"] == "Yes" for r in protocol) == 20
    python_files = [p for p in ROOT.rglob("*.py") if ".git" not in p.parts]
    for path in python_files:
        ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path.relative_to(ROOT)))
    notebooks = list((ROOT / "analysis/squidiff/notebooks").glob("*.ipynb"))
    for path in notebooks:
        notebook = json.loads(path.read_text())
        assert notebook["nbformat"] == 4
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                assert cell["execution_count"] is None and cell["outputs"] == []
                source = "".join(cell["source"])
                assert "/home/" not in source and "/media/" not in source
    checksums = ROOT / "docs/provenance/SHA256SUMS.txt"
    verified = 0
    for line in checksums.read_text().splitlines():
        digest, name = line.split("  ", 1)
        target = ROOT / name
        assert target.is_file(), name
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, name
        verified += 1
    print(json.dumps({"gene_counts": {k: len(v) for k, v in panels.items()},
                      "cell_splits": dict(counts), "classifier_protocols": len(protocol),
                      "Python_sources_parsed": len(python_files), "stripped_notebooks": len(notebooks),
                      "verified_SHA256_files": verified,
                      "scope": "Static/resource validation; no model fitting or scientific result regeneration"}, indent=2))


if __name__ == "__main__":
    main()
