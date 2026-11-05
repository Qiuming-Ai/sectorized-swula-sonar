"""Freeze only the legacy Fig. 2 curve slices and Fig. 7 diagram geometry.

These optional archival input files are not required to replot the distributed
data. Pass their paths explicitly if you hold the originals.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def export_fig2(source, output):
    import h5py

    with h5py.File(source, "r") as mat:
        angle = mat["grid/tgt_deg"][()].ravel()
        rx, tx, product = (mat[key][()] for key in ("PSF_Rx", "AF_Tx", "G2"))
    # HDF5 exposes the transpose of the MATLAB matrices. These are exactly
    # the columns corresponding to the rows selected in legacy PSF.m.
    rows = np.array([850, 790, 910])
    np.testing.assert_allclose(angle[rows], [0, -6, 6], atol=1e-10)
    assert rx.shape == product.shape == (len(angle), len(angle))
    assert tx.shape == (len(angle), 1)
    data = np.column_stack((
        angle, rx[:, rows] / rx.max(), tx[:, 0] / tx.max(),
        product[:, rows] / product.max(),
    ))
    np.savetxt(
        output / "figure2_legacy_psf_slices.csv", data, delimiter=",",
        header="angle_deg,rx_focus_0,rx_focus_minus6,rx_focus_plus6,tx_focus_0,"
               "product_focus_0,product_focus_minus6,product_focus_plus6",
        comments="", fmt="%.17g",
    )
    return {
        "archive_file": source.name, "sha256": sha256(source),
        "rows": int(len(angle)), "source_recipe": "PSF.m",
        "selected_matlab_rows": [851, 791, 911],
        "normalization": "Linear magnitudes; each family divided by its full-array maximum.",
    }


def export_fig7(source, output):
    import fitz

    with fitz.open(source) as pdf:
        page = pdf[0]
        drawings = page.get_drawings()
        spans = [s for b in page.get_text("dict")["blocks"] if "lines" in b
                 for line in b["lines"] for s in line["spans"] if s["text"].strip()]
    nodes = sorted(
        (d for d in drawings if d.get("fill") is not None
         and d["rect"].width > 30 and 20 < d["rect"].height < 40),
        key=lambda d: d["rect"].x0,
    )
    links = [d for d in drawings if d not in nodes and d["rect"].height < 10]
    assert len(nodes) == 8 and len(links) == 14
    assert sum(d.get("fill") is not None for d in links) == 7
    ink = fitz.Rect(nodes[0]["rect"])
    for node in nodes[1:]:
        ink |= node["rect"]
    pad = 3.5
    width, height = ink.width + 2 * pad, ink.height + 2 * pad

    def point(p):
        return [p.x - ink.x0 + pad, height - (p.y - ink.y0 + pad)]

    labels = []
    for node in nodes:
        selected = [s for s in spans if node["rect"].contains(fitz.Rect(s["bbox"]))]
        selected.sort(key=lambda s: (s["bbox"][1], s["bbox"][0]))
        center = (node["rect"].tl + node["rect"].br) / 2
        labels.append({"center_pt": point(center), "lines": [s["text"] for s in selected]})
    assert [" ".join(x["lines"]) for x in labels] == [
        "Signal transmission", "Channel propagation", "Matched filtering", "Shift and padding",
        "FFT", "Phase compensation", "Weighted sum", "IFFT",
    ]
    paths = []
    for d in drawings:
        if d not in nodes and d not in links:
            continue
        assert all(item[0] in ("l", "c") for item in d["items"])
        paths.append({
            "items": [[item[0], *[point(p) for p in item[1:]]] for item in d["items"]],
            "fill_rgb": d.get("fill"), "stroke_rgb": d.get("color"),
            "width_pt": d.get("width") or .75, "closed": bool(d.get("closePath")),
        })
    payload = {"width_pt": width, "height_pt": height, "paths": paths, "nodes": labels,
               "operation": "Preserve original vector geometry; crop external whitespace only."}
    (output / "figure7_pipeline_geometry.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8",
    )
    return {"archive_file": source.name, "sha256": sha256(source),
            "stages": len(labels), "arrows": len(links) // 2, "vector_paths": len(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-mat", type=Path, required=True)
    parser.add_argument("--pipeline-pdf", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    record = {
        "figure_2": export_fig2(args.legacy_mat, args.output),
        "figure_7": export_fig7(args.pipeline_pdf, args.output),
    }
    (args.output / "archive_extraction_provenance.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8",
    )


if __name__ == "__main__":
    main()
