# Enhancer-driven GRN (SCENIC+) for V1 L2/3 IT neurons — Wang et al. 2025

**Date:** 2026-07-16 · **Env:** `epics` (SCENIC+ 1.0a2, pycisTopic 2.0a0) + `chrombpnet` (htslib/bedtools) + `macs3`
**Plan:** [`plan/03.wang25_l23_scenicplus_executable.md`](../plan/03.wang25_l23_scenicplus_executable.md) · **Scripts:** `scripts/09`–`scripts/18`

---

## 1. Summary

We inferred an enhancer-driven gene-regulatory network (eGRN / eRegulons) for **V1 layer-2/3 intratelencephalic
excitatory neurons** (`EN-L2_3-IT`, primary visual cortex) from the paired 10x snMultiome atlas of
**Wang L. et al., *Nature* 647:169–178 (2025)**, "Molecular and cellular dynamics of the developing human
neocortex" (hg38).

**Headline result:** a biologically coherent **39-eRegulon network (37 TFs)** whose activity traces a
developmental cascade across the four stages V1 L2/3 spans (2nd trimester → adolescence), anchored by
canonical cortical upper-layer regulators (POU3F2/BRN2, POU3F3/BRN1, MEF2C, FOXP1, TCF4, ZBTB18) and an
activity-dependent Jun/JunB AP-1 module.

---

## 2. Data & scope

| | |
| --- | --- |
| Assay | Paired 10x snMultiome (RNA + ATAC, same nuclei), hg38 |
| Source | CELLxGENE (RNA h5ad + ATAC fragments) — chosen over the 32.8 GB Dryad Seurat `.rds` |
| Cell type | `Type == "EN-L2_3-IT"` **and** `Region == "V1"` |
| Cells | **5,558** V1 L2/3 nuclei (~94% postnatal: Infancy 3,672 + Adolescence 1,740; prenatal 146) |
| Stages present | Second_trimester (38), Third_trimester (108), Infancy (3,672), Adolescence (1,740) |

Scope decision: **V1 only** (PFC and the tiny fetal "General" pool excluded) to avoid the PFC↔V1
regional-regulatory confound that is a central theme of the paper.

---

## 3. Pipeline (Stages 0–6)

| Stage | Script | What | Key result |
| --- | --- | --- | --- |
| 0 | `09_download.sh` | Fetch CELLxGENE RNA + ATAC fragments | ~45 GB (RNA 2.6 GB, fragments 40 GB) |
| 1 | `10_select_l23.py` | Subset V1 L2/3; raw counts → `.X`; **gene symbols** as var_names | 5,558 × 35,468 genes |
| 2 | `11_refs.sh` | chrom sizes, ENCODE blacklist, biomart gene annotation → **UCSC + primary chroms** | annot 103,810 rows, blacklist 636 |
| 2b | `11b_fragments_to_ucsc.sh` | Drop Ensembl-named scaffolds from fragments (primary chroms already UCSC) | UCSC fragments (40 GB, `chr1..chrY`) |
| 3a | `12_call_peaks.sh` + `12b_consensus_peaks.py` | MACS3 `-g hs` on L2/3 fragments → **consensus peaks** | 209,551 raw → **141,169** non-overlapping (500 bp) |
| 3b | `13_build_ctx_db.sh` | Custom cisTarget motif DB on the consensus peaks | **141,169 regions × 10,249 motifs** (rankings 3.8 GB) |
| 3c | `14a_extract_l23_fragments.sh` + `14_run_pycistopic.py` | cisTopic object + LDA topic model | 141,169 × 5,558; **40 topics** → 40 region-sets (120,000 regions) |
| 4 | `15_prep_tf_names.py`, `16_prepare_menr.sh` | MuData, search space, motif enrichment, cistromes | MuData 5,558 × (35,468 + 141,169); **756 cistrome TFs** |
| 5 | `17_run_grn.sh` (+ `17b_egrn_extended.sh`) | TF→gene, region→gene, eGRN, AUCell | **39 eRegulons / 37 TFs / 9,766 triplets** |
| 6 | `18_qc.py` | eRegulon summary + developmental-activity QC | positive controls + developmental cascade |

### Stage detail & data-derived decisions

- **Peak set — our own, not the authors'.** Wang25 publishes no standalone peak BED; their peaks live only
  in the 32.8 GB atlas `.rds` (R/Signac) or a 3.4 GB Cellranger h5, and are a **global** atlas set (all
  cell types / regions / ages). For a single-cell-type L2/3 eGRN, cell-type pseudobulk MACS peaks are the
  SCENIC+-recommended input, so we called our own.
- **Consensus peaks.** MACS3 `--call-summits` produced 208,641 heavily overlapping 501 bp windows; replaced
  with `pycisTopic.get_consensus_peaks` (Corces-style iterative filtering + blacklist) → **141,169**
  non-overlapping 500 bp peaks (also ~40% cheaper for the Stage 3b DB build).
- **LDA topics.** Sweep `[10,20,30,40]`, n_iter=150; auto-selected **40** (the top of the range — optimum
  at the boundary; accepted since topics are only candidate region-sets and the eGRN re-filters).
- **TF list.** Stage 5 uses the **756 cistrome TFs** from `prepare_menr` (only motif-backed, data-derived
  TFs can form eRegulons), superseding a broader pre-generated 1,292-TF candidate list.

---

## 4. Key results

### 4.1 eGRN — direct vs extended cistromes

| | Direct (`eRegulon.tsv`) | Extended (`eRegulon_extended.tsv`) |
| --- | --- | --- |
| eRegulons | **39** | 31 |
| unique TFs | 37 | 27 |
| triplets | 9,766 | 6,345 |
| target genes / eRegulon | median 125 (10–844) | — |

The two are **complementary, not nested** (orthology annotation reassigns region↔TF links): the extended run
**recovers SATB2** (+ RBPJ, TCF12, HES4, MYEF2) but loses POU3F2 / OLIG1 / BCL11A / ZBTB18. Their **union
covers 4/6 positive controls**.

### 4.2 Positive-control L2/3 TFs

| TF | Direct | Extended | Note |
| --- | --- | --- | --- |
| POU3F2 (BRN2) | ✅ | – | upper-layer marker |
| SATB2 | – | ✅ | had cistrome; dropped in direct thresholding |
| MEF2C | ✅ | ✅ | |
| FOXP1 | ✅ | ✅ | |
| CUX2 | ✗ | ✗ | **no cistrome** (motif not enriched in L2/3 topics) |
| RORB | ✗ | ✗ | **no cistrome** |

The eRegulon TF set is otherwise cortical-neuron-appropriate: POU3F3 (BRN1), TCF4, ZBTB18 (RP58), MEF2A,
RFX3, BCL11A, PBX1, MEIS3, OLIG1, NR3C1, plus expected ubiquitous factors (YY1, JUN/JUNB, ATF4, NRF1).

### 4.3 Developmental dynamics (AUCell, z-scored across stages)

Single cell type → grouped by developmental stage rather than by cell type. A coherent cascade emerged
(figure: `data/wang25/work/qc/eRegulon_activity_by_stage.png`):

```
2nd trimester  → OLIG1, BCL11A, PBX1              (early / progenitor, deep-layer)
3rd trimester  → POU3F2, POU3F3, MEF2C, ZBTB18    (upper-layer maturation)
Infancy        → PBX1(+/+), RFX3, MEIS3, EGR3, YY1
Adolescence    → FOXP1, MEF2A, NFAT5, NR3C1, TCF4, RORA, SMAD3  (mature / activity-dependent)
```

*Caveat:* prenatal stages are low-N (38 / 108 cells) vs Infancy / Adolescence (3,672 / 1,740), so the
2nd/3rd-trimester columns are noisier.

### 4.4 AP-1 family — a Jun-anchored, activity-dependent module

Three AP-1 members surfaced (both runs): **JUN** (588 target genes), **JUNB** (185), **JDP2** (144).
The **FOS half of the canonical Jun/Fos dimer is entirely absent** — not for lack of the (shared) AP-1
motif, but because of **expression**: JUN/JUNB/JDP2 are detected in 21–42% of nuclei, whereas FOS/FOSB/FOSL1
are detected in ≤15% (FOSB/FOSL1 ≈ 0). The eGRN's co-expression step can't build a robust regulon from such
sparse signal, so well-expressed **JUN "claims" the shared AP-1 motif** — compounded by snRNA's poor capture
of transient immediate-early FOS transcripts. Consistently, **JUNB's activating regulon is a textbook
activity-dependent program** — targets include **ARC, HOMER1, NR4A2/3, EGR2/3, DUSP1, VGF, SCG2, and FOS
itself** — i.e., FOS appears as a *target*, not a driver.

---

## 5. Runtime issues resolved (lessons)

The skeleton→executable transition surfaced many real-data gotchas; all fixed and committed:

| Issue | Resolution |
| --- | --- |
| Fragments mixed UCSC (primary) + Ensembl (scaffolds); gene annotation fully Ensembl; fasta/blacklist UCSC | Standardize on **UCSC + primary chroms** everywhere (fragments, gene annotation → chr-prefixed) |
| `bgzip`/`tabix`/`bedtools`/`samtools`/`macs3` not in `epics` | Use `chrombpnet` env (htslib/bedtools) and `macs3` env by explicit path |
| pycisTopic reads a fragment file as text unless it ends in `.gz` | Rename UCSC fragments `.bgz` → `.gz` |
| `regions.join(fragments)` (ray, `n_cpu>1`) returns NaN at scale → `IntCastingNaNError` | Build cisTopic with **`n_cpu=1`** (serial join); LDA stays parallel. Also pre-extract an L2/3-only 1.6 GB fragment file for fast reads |
| RNA `var_names` were Ensembl IDs; gene annotation + tf_names use symbols → `search_spance KeyError 'Gene'` | Remap RNA to **gene symbols** (Stage 1) |
| `search_spance` chromsizes needs a `Chromosome/Start/End` header file | Build it inline (not the 2-col `chrom.sizes`, which bedtools needs) |
| `prepare_GEX_ACC` metacell flags are non-multiome-only | Run **per-cell multiome**; match barcodes with `--bc_transform_func "x+'___l23'"`, `--do_not_use_raw` |
| `motif_enrichment_cistarget --region_set_folder` iterates **subfolders** | Point at the parent `region_sets/`, not the `topics_top_3k` leaf |
| eGRN's `--cistromes_fname` needs `prepare_menr` output (missing from skeleton) | Added Stage **4d prepare_menr** → `cistromes_{direct,extended}.h5ad` |

---

## 6. Reproducibility

- **Scripts** `scripts/09_download.sh` … `scripts/18_qc.py` (+ `17b_egrn_extended.sh`); all paths defined at
  the top of each script per repo convention.
- **Envs:** `epics` (SCENIC+/pycisTopic/pycistarget), `chrombpnet` (htslib 1.22.1, bedtools, samtools),
  `macs3` (MACS3 3.0.4). Nothing was installed into `epics`.
- **Commit trail** (main): `cad95c0` (Stage 0/1) → `647c2d1`/`6bf7562` (chrom naming) → `285cd1a`/`307b5a1`
  (Stage 3) → `43c8863`/`eb444f6` (Stage 3c fixes) → `908be52` (Stage 4) → `5c8ab6e` (Stage 5) →
  `ba54ba8` (Stage 6) → `4491fcc` (extended variant).

### Output manifest (`data/wang25/`, gitignored)

| File | Contents |
| --- | --- |
| `work/cistopic_obj_model.pkl` | cisTopic object + 40-topic LDA (607 MB) |
| `db/wang25_l23.regions_vs_motifs.{rankings,scores}.feather` | custom cisTarget DB (3.8 / 2.6 GB) |
| `work/ACC_GEX.h5mu` | multiome MuData 5,558 × (35,468 + 141,169) |
| `work/cistromes_{direct,extended}.h5ad` | motif→TF cistromes |
| `work/eRegulon.tsv` / `work/eRegulon_extended.tsv` | **final eGRNs** (2.2 / 1.4 MB) |
| `work/AUCell.h5mu` / `work/AUCell_extended.h5mu` | per-cell eRegulon activity (Gene + Region based) |
| `work/qc/` | summary CSV, activity-by-stage CSV, 2 figures, Jun-side target lists |

---

## 7. Suggested next steps

- **Merge direct + extended** eRegulons into one non-redundant network (recovers both POU3F2 and SATB2).
- Add **DEM** motif enrichment and merge with cisTarget cistromes (full-pipeline default).
- Validate JUN/JUNB targets against known experience-dependent visual-cortex genes; test whether CUX2/RORB
  targets can be recovered with a relaxed motif-enrichment threshold or an alternative motif annotation.
