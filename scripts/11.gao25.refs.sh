#!/usr/bin/env bash
# Stage 2 (gao25) — Mouse references for the DevVIS L2/3 SCENIC+ run (mm10 / GRCm38, UCSC chrom style).
# mm10 fasta/chromsizes/blacklist are REUSED in place (no copy). Only the gene annotation is produced,
# and it MUST be GRCm38 to match the mm10 peaks (verified: mm10.chrom.sizes chr1 = 195,471,971 = GRCm38).
#
# HIGH-RISK NOTE (plan open risk #1): the biomart route is unusable here — SCENIC+
# `download_genome_annotations` against the nov2020 GRCm38 archive host returns malformed XML, and the
# v1astro precedent annotation is actually GRCm39 (its Xkr4 = chr1:3,276,124-3,741,721 matches GRCm39,
# not GRCm38 3,214,482-3,671,498 — a latent assembly mismatch). We therefore derive the SCENIC+
# annotation TSV directly from GENCODE vM25 — the FINAL GRCm38 GENCODE mouse release — which ships
# UCSC-style `chr` names natively. Output columns replicate the biomart schema exactly:
#   Chromosome, Start(gene), End(gene), Strand, Gene(symbol), Transcription_Start_Site(per transcript),
#   Transcript_type   — protein_coding transcripts only.
# The GENCODE GTF lives with the shared mm10 genome under the common ref dir (reusable).
# Run: bash scripts/11.gao25.refs.sh
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/gao25
DB=${DATA}/db
MM10_DIR=/home/qlyu/mydata/data/common/mm10

MM10_CHROMSIZES=${MM10_DIR}/mm10.chrom.sizes
GTF_GZ=${MM10_DIR}/gencode.vM25.annotation.gtf.gz       # GENCODE vM25 = final GRCm38 mouse release (shared ref)
GENE_ANNOT=${DB}/mm10_gene_annotation.tsv               # final: UCSC chr names, GRCm38 coords
MM10_CHROMSIZES_SP=${DB}/mm10_chromsizes_scenicplus.tsv # search_spance format: Chromosome/Start/End header

mkdir -p "${DB}"
[[ -s "${MM10_CHROMSIZES}" ]] || { echo "ERROR: mm10 chrom.sizes missing: ${MM10_CHROMSIZES}" >&2; exit 1; }
[[ -s "${GTF_GZ}" ]] || { echo "ERROR: GENCODE GTF missing: ${GTF_GZ} (download it first)" >&2; exit 1; }
gzip -t "${GTF_GZ}" || { echo "ERROR: GTF is truncated/corrupt: ${GTF_GZ}" >&2; exit 1; }

# sanity: confirm mm10 == GRCm38 (chr1 length) before trusting GRCm38 coords downstream
CHR1=$(awk -F'\t' '$1=="chr1"{print $2}' "${MM10_CHROMSIZES}")
[[ "${CHR1}" == "195471971" ]] || { echo "ERROR: mm10 chr1=${CHR1} != GRCm38 195471971 — assembly mismatch." >&2; exit 1; }
echo "[$(date)] mm10 == GRCm38 confirmed (chr1=${CHR1})"

# --- derive the SCENIC+ gene annotation TSV from the GENCODE GTF (GRCm38, UCSC chr names) ---
echo "[$(date)] building SCENIC+ gene annotation -> ${GENE_ANNOT}"
conda run --no-capture-output -n epics python -u - "${GTF_GZ}" "${GENE_ANNOT}" <<'PY'
import gzip, re, sys
GTF, OUT = sys.argv[1], sys.argv[2]
# GENCODE mouse already uses UCSC chr names. Keep the mm10 primary contigs; drop scaffolds.
PRIMARY = {f"chr{i}" for i in range(1, 20)} | {"chrX", "chrY", "chrM"}
def attr(s, key):
    m = re.search(key + r' "([^"]+)"', s)
    return m.group(1) if m else None

genes = {}          # gene_id -> (chrom, start, end, strand, name)
rows = []           # per protein_coding transcript: (gene_id, tss)
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        chrom, feat, start, end, strand, attrs = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        if chrom not in PRIMARY:
            continue
        if feat == "gene":
            genes[attr(attrs, "gene_id")] = (chrom, start, end, strand, attr(attrs, "gene_name"))
        elif feat == "transcript":
            if attr(attrs, "transcript_type") != "protein_coding":
                continue
            tss = start if strand == "+" else end        # + -> transcript start, - -> transcript end
            rows.append((attr(attrs, "gene_id"), tss))

n_skipped = 0
with open(OUT, "w") as out:
    out.write("Chromosome\tStart\tEnd\tStrand\tGene\tTranscription_Start_Site\tTranscript_type\n")
    for gid, tss in rows:
        g = genes.get(gid)
        if g is None or g[4] is None:
            n_skipped += 1
            continue
        chrom, gstart, gend, strand, name = g
        out.write(f"{chrom}\t{gstart}\t{gend}\t{strand}\t{name}\t{tss}\tprotein_coding\n")
print(f"  protein_coding transcript rows: {len(rows) - n_skipped:,} (skipped {n_skipped}); "
      f"protein_coding genes with symbol: "
      f"{len({genes[g][4] for g in genes if genes[g][4]}):,}")
PY
[[ -s "${GENE_ANNOT}" ]] || { echo "ERROR: gene annotation not written: ${GENE_ANNOT}" >&2; exit 1; }
echo "  gene annotation: $(($(wc -l < "${GENE_ANNOT}") - 1)) rows -> ${GENE_ANNOT}"
echo "  chroms kept: $(cut -f1 "${GENE_ANNOT}" | tail -n +2 | sort -u | tr '\n' ' ')"

# --- search_spance chromsizes (Chromosome/Start/End header) from mm10 primary chroms (chr1-19,X,Y,M) ---
printf 'Chromosome\tStart\tEnd\n' > "${MM10_CHROMSIZES_SP}"
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|[XYM])$/ {print $1, 0, $2}' "${MM10_CHROMSIZES}" >> "${MM10_CHROMSIZES_SP}"
echo "  search-space chromsizes: $(($(wc -l < "${MM10_CHROMSIZES_SP}") - 1)) chroms -> ${MM10_CHROMSIZES_SP}"

# --- spot-check GRCm38 assembly on a known gene (Xkr4 GRCm38 gene ≈ chr1:3,214,482-3,671,498) ---
echo "  Xkr4 (expect GRCm38 ~chr1:3,214,482-3,671,498):"
awk -F'\t' '$5=="Xkr4"{print "   ", $1, $2, $3, $4; exit}' "${GENE_ANNOT}"

echo "[$(date)] done. references:"
ls -lh "${GENE_ANNOT}" "${MM10_CHROMSIZES_SP}"
