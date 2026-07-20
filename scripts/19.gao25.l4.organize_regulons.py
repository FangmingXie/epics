#!/usr/bin/env python
"""Stage 7 (gao25 · l4) — Organize SCENIC+ eRegulons into a clean regulon -> gene table (mouse L4/5).

Fork of scripts/19.gao25.organize_regulons.py for the L4 layer (plan/06); logic unchanged, repointed
at the l4 work dir. From the Stage 5 direct/extended eRegulon tables, build one long table where each
row is a (regulon, gene) pair, after:
  - keeping only R2G-positive regulons  (name signs `++` and `-+`; drops all `*/-`)
  - deduplicating direct vs extended by TF, favoring the direct version
  - collapsing the multiple region rows per gene into aggregated scores

Run: conda run --no-capture-output -n epics python -u scripts/19.gao25.l4.organize_regulons.py
"""
import argparse
import os
import sys

import pandas as pd

# --- paths (capitalized per CLAUDE.md) ---
PROJ = "/data/qlyu/project/epics"
WORK = f"{PROJ}/data/gao25/l4/work"

EREGULON_DIRECT = f"{WORK}/eRegulon.tsv"             # Stage 5c (is_extended=False)
EREGULON_EXTENDED = f"{WORK}/eRegulon_extended.tsv"  # Stage 5c extended (is_extended=True)
OUT_TABLE = f"{WORK}/regulon_gene_table.tsv"


def load_eregulon(path):
    """Load an eRegulon TSV; fail loudly if missing/empty."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        sys.exit(f"ERROR: required input missing/empty: {path}")
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        sys.exit(f"ERROR: no rows in {path}")
    return df


def add_signs(df):
    """Parse TF2G/R2G signs from eRegulon_name suffix `..._{s1}/{s2}`."""
    signs = df["eRegulon_name"].str.rsplit("_", n=1).str[-1]  # e.g. '+/+'
    parts = signs.str.split("/", expand=True)
    df = df.copy()
    df["TF2G_sign"] = parts[0]
    df["R2G_sign"] = parts[1]
    if not df["R2G_sign"].isin(["+", "-"]).all():
        sys.exit(f"ERROR: unexpected sign format in eRegulon_name: {signs.unique()}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--direct", default=EREGULON_DIRECT, help="direct eRegulon TSV")
    ap.add_argument("--extended", default=EREGULON_EXTENDED, help="extended eRegulon TSV")
    ap.add_argument("--out", default=OUT_TABLE, help="output regulon->gene TSV")
    args = ap.parse_args()

    direct = add_signs(load_eregulon(args.direct))
    extended = add_signs(load_eregulon(args.extended))

    # keep only R2G-positive regulons (++ and -+)
    direct = direct[direct["R2G_sign"] == "+"]
    extended = extended[extended["R2G_sign"] == "+"]

    # dedup by TF: a TF present in direct drops its extended rows entirely
    direct_tfs = set(direct["TF"])
    extended_keep = extended[~extended["TF"].isin(direct_tfs)]
    combined = pd.concat([direct, extended_keep], ignore_index=True)

    # collapse the multiple region rows per gene into one (regulon, gene) row
    grouped = (
        combined.groupby(
            ["TF", "eRegulon_name", "is_extended", "TF2G_sign", "R2G_sign", "Gene"],
            as_index=False,
        )
        .agg(
            n_regions=("Region", "size"),
            max_importance_x_rho=("importance_x_rho", "max"),
            best_triplet_rank=("triplet_rank", "min"),
        )
    )

    out = grouped[
        ["TF", "eRegulon_name", "is_extended", "TF2G_sign", "R2G_sign", "Gene",
         "n_regions", "max_importance_x_rho", "best_triplet_rank"]
    ].sort_values(["TF", "best_triplet_rank"]).reset_index(drop=True)

    out.to_csv(args.out, sep="\t", index=False)

    n_reg = out["eRegulon_name"].nunique()
    n_direct = out.loc[~out["is_extended"], "eRegulon_name"].nunique()
    n_ext = out.loc[out["is_extended"], "eRegulon_name"].nunique()
    print(f"Wrote {args.out}")
    print(f"  regulons: {n_reg}  (direct={n_direct}, extended={n_ext})")
    print(f"  TFs:      {out['TF'].nunique()}")
    print(f"  (regulon, gene) rows: {len(out)}")


if __name__ == "__main__":
    main()
