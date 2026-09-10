#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.tree import nj

def sanitize(labels):
    return [str(x).strip().replace(" ", "_").replace("(", "_").replace(")", "_").replace(",", "_")
            for x in labels]

def read_matrix_xlsx(path, sheet=None, meta_cols=None, drop_unnamed=True):
    xls = pd.ExcelFile(path)
    target_sheet = sheet if sheet is not None else xls.sheet_names[0]
    df = pd.read_excel(xls, sheet_name=target_sheet)

    df.columns = [str(c).strip() for c in df.columns]

    if drop_unnamed:
        df = df.loc[:, ~df.columns.str.match(r"^Unnamed", na=False)]

    default_meta = ["Protist", "Mode", "genus", "Phylum", "Class", "Order", "Family"]
    meta_set = set([c.lower() for c in (meta_cols if meta_cols else default_meta)])

    sample_cols = [c for c in df.columns if c.lower() not in meta_set]

    if len(sample_cols) < 2:
        sys.exit(f"Insufficient sample columns (less than 2 after removing metadata). sheet={target_sheet}")

    Xdf = df[sample_cols].copy()

    for c in Xdf.columns:
        Xdf[c] = pd.to_numeric(Xdf[c], errors="coerce")

    all_nan_cols = [c for c in Xdf.columns if Xdf[c].notna().sum() == 0]
    if all_nan_cols:
        Xdf.drop(columns=all_nan_cols, inplace=True)

    if Xdf.shape[1] < 2:
        sys.exit(f"Insufficient numeric sample columns (less than 2 after removing non-numeric columns). sheet={target_sheet}")

    Xdf = Xdf.fillna(0.0)

    Xdf = Xdf.loc[:, Xdf.sum(axis=0) > 0]

    if Xdf.shape[1] < 2:
        sys.exit(f"Insufficient valid sample columns (less than 2 after removing all-zero columns). sheet={target_sheet}")

    Xdf = Xdf.loc[Xdf.sum(axis=1) > 0, :]

    if Xdf.shape[0] < 2:
        sys.exit(f"Insufficient valid feature rows (less than 2 after removing all-zero rows). sheet={target_sheet}")

    return Xdf, target_sheet

def main():
    ap = argparse.ArgumentParser(description="Build Neighbor-Joining tree using any distance metric (default Bray-Curtis) and output Newick; automatically ignores metadata and Unnamed columns.")
    ap.add_argument("-i", "--input", required=True, help="Input Excel file")
    ap.add_argument("-s", "--sheet", default=None, help="Sheet name (defaults to first sheet)")
    ap.add_argument("-m", "--metric", default="braycurtis", help="Distance metric (default braycurtis, accepts SciPy metrics like euclidean)")
    ap.add_argument("-o", "--out", default="genus_NJ.newick", help="Output Newick filename")
    ap.add_argument("--meta-cols", default="Protist,Life,genus,Phylum,Class,Order,Family",
                    help="Comma-separated metadata column names (case-insensitive) to be excluded (default: Protist,Life,genus,Phylum,Class,Order,Family)")
    ap.add_argument("--no-drop-unnamed", action="store_true", help="Do not automatically drop columns starting with 'Unnamed:'")
    ap.add_argument("--save-samples", default=None, help="Optional: Save the final sample column names to this file")
    args = ap.parse_args()

    meta_cols = [c.strip() for c in args.meta_cols.split(",") if c.strip()]
    Xdf, used_sheet = read_matrix_xlsx(
        args.input,
        sheet=args.sheet,
        meta_cols=meta_cols,
        drop_unnamed=(not args.no_drop_unnamed)
    )

    X = Xdf.T.values
    labels = sanitize(Xdf.columns.tolist())

    D = squareform(pdist(X, metric=args.metric))
    if not np.isfinite(D).all():
        sys.exit("Distance matrix contains non-finite values. Please check for remaining NaNs or all-zero columns/rows.")

    dm = DistanceMatrix(D, labels)
    tree = nj(dm)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(str(tree) + "\n")

    if args.save_samples:
        with open(args.save_samples, "w", encoding="utf-8") as fo:
            for s in Xdf.columns.tolist():
                fo.write(f"{s}\n")

    print(f"[OK] NJ tree saved -> {args.out}  (sheet={used_sheet}, samples={len(labels)}, metric={args.metric})")
    if args.save_samples:
        print(f"[OK] Samples used in calculation saved to -> {args.save_samples}")

if __name__ == "__main__":
    main()
