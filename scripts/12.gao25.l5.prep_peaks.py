#!/usr/bin/env python
"""Stage 3a (gao25 · l5) — L5-accessible peak matrix from the precomputed ATAC count matrix.

Fork of scripts/12.gao25.prep_peaks_matrix.py for the L5 layer (plan/06). Method and the 5% peak
filter are unchanged; only the barcode/BED/npz paths point at the l5 namespace. The DevVIS ATAC h5ad
ships a peaks×cells count matrix (882,075 fixed-width peaks), so there is no peak calling. We:
  1. read the ATAC h5ad (backed) and subset rows to the final multiome L5 barcodes (ATAC _index),
  2. keep peaks accessible in >= MIN_CELL cells (5% of cells) — cell-type-relevant, lighter DB,
  3. convert peak names 'chr1_3079841_3080341' -> 'chr1:3079841-3080341' (colon form for cisTopic),
  4. write the L5 peak BED (for the cisTarget DB build) and a regions×cells CSR npz (for cisTopic).

Run in the epics env:
  conda run --no-capture-output -n epics python -u scripts/12.gao25.l5.prep_peaks.py
"""
import numpy as np
import scipy.sparse as sp
import anndata as ad

# --- paths (capitalized per CLAUDE.md) ---
ORIG = "/home/qlyu/mydata/data/gao25/orig"
ATAC_H5AD = f"{ORIG}/DevVIS_multiome_snATAC_processed.h5ad"    # 45G, peaks×cells matrix (X = cells×peaks CSR)
DATA = "/data/qlyu/project/epics/data/gao25/l5"
L5_BARCODES = f"{DATA}/vis_l5_barcodes.txt"                   # Stage 1 (common key == ATAC _index)
WORK = f"{DATA}/work"
L5_PEAKS_BED = f"{WORK}/vis_l5_peaks.bed"                     # output: chr\tstart\tend (for DB build)
L5_MATRIX_NPZ = f"{WORK}/vis_l5_matrix.npz"                  # output: regions×cells CSR + names (cisTopic)

MIN_CELL_FRAC = 0.05                                          # keep peaks accessible in >= 5% of L5 cells


def main():
    valid_bc = set(open(L5_BARCODES).read().split())
    print(f"[Stage 3a gao25·l5] {len(valid_bc):,} target L5 barcodes")

    print(f"[Stage 3a gao25·l5] opening {ATAC_H5AD} (backed)")
    adata = ad.read_h5ad(ATAC_H5AD, backed="r")
    print(f"  full ATAC: {adata.n_obs:,} cells x {adata.n_vars:,} peaks")

    mask = adata.obs_names.isin(valid_bc)
    n_sel = int(mask.sum())
    print(f"  matched barcodes in ATAC: {n_sel:,}")
    if n_sel != len(valid_bc):
        raise ValueError(f"barcode mismatch: {n_sel} matched of {len(valid_bc)} — expected all present.")

    sub = adata[mask].to_memory()             # cells×peaks CSR, int counts (float64 on disk)
    X = sub.X.tocsr()
    cell_names = sub.obs_names.to_numpy().astype(str)     # ATAC _index (common key)
    n_cells = X.shape[0]
    print(f"  L5 submatrix: {X.shape[0]:,} cells x {X.shape[1]:,} peaks; nnz={X.nnz:,}")

    # per-peak nonzero cell count -> keep peaks accessible in >= MIN_CELL cells
    min_cell = int(round(MIN_CELL_FRAC * n_cells))
    per_peak = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = per_peak >= min_cell
    n_keep = int(keep.sum())
    print(f"  peak filter: >= {min_cell} cells ({MIN_CELL_FRAC:.0%}); kept {n_keep:,} of {X.shape[1]:,} peaks")
    if n_keep < 10000:
        raise ValueError(f"only {n_keep} peaks passed — filter too strict or matrix wrong.")

    keep_idx = np.where(keep)[0]
    Xk = X[:, keep_idx]                        # cells×peaks(kept)
    peak_names_raw = sub.var_names.to_numpy().astype(str)[keep_idx]   # 'chr1_3079841_3080341'

    # convert 'chr1_3079841_3080341' -> colon/BED forms (rsplit: chrom may itself contain no '_')
    region_names, bed_lines = [], []
    for r in peak_names_raw:
        chrom, start, end = r.rsplit("_", 2)
        region_names.append(f"{chrom}:{start}-{end}")
        bed_lines.append(f"{chrom}\t{start}\t{end}\n")
    region_names = np.asarray(region_names, dtype=object)

    with open(L5_PEAKS_BED, "w") as fh:
        fh.writelines(bed_lines)
    print(f"  wrote {L5_PEAKS_BED} ({n_keep:,} peaks)")

    # regions×cells CSR (transpose) for cisTopic create_cistopic_object.
    # counts are integer-valued but float64 on disk; pycisTopic's LDA needs an integer matrix -> int32.
    if not np.allclose(Xk.data, np.rint(Xk.data)):
        raise ValueError("ATAC matrix has non-integer values — not raw counts.")
    RxC = Xk.T.tocsr().astype(np.int32)
    density = RxC.nnz / (RxC.shape[0] * RxC.shape[1])
    print(f"  regions×cells matrix: {RxC.shape[0]:,} x {RxC.shape[1]:,}; nnz={RxC.nnz:,}; density={density:.4f}")
    np.savez(
        L5_MATRIX_NPZ,
        data=RxC.data, indices=RxC.indices, indptr=RxC.indptr, shape=np.asarray(RxC.shape),
        cell_names=cell_names, region_names=region_names,
    )
    print(f"  wrote {L5_MATRIX_NPZ}")
    print("[Stage 3a gao25·l5] complete.")


if __name__ == "__main__":
    main()
