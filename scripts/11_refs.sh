#!/usr/bin/env bash
# Stage 2 — Human references for the Wang25 L2/3 SCENIC+ run (hg38, UCSC chrom style).
# Produces: chrom sizes, ENCODE hg38 blacklist v2, and SCENIC+ gene annotation.
set -euo pipefail

# --- paths (capitalized per CLAUDE.md) ---
PROJ=/data/qlyu/project/epics
DATA=${PROJ}/data/wang25
DB=${DATA}/db

HG38_FASTA=/home/qlyu/mydata/data/v1_multiome_ai/chrombpnet_tutorial/data/downloads/hg38.fa
HG38_FAI=${HG38_FASTA}.fai
HG38_CHROMSIZES=${DB}/hg38.chrom.sizes
HG38_BLACKLIST=${DB}/hg38-blacklist.v2.bed
GENE_ANNOT=${DB}/hg38_gene_annotation.tsv            # final: UCSC style, primary chroms
GENE_ANNOT_RAW=${DB}/hg38_gene_annotation.raw.tsv    # raw biomart output (Ensembl style)
CHROMSIZES_UCSC=${DB}/hg38_chromsizes_ucsc.tsv

BLACKLIST_URL=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz

mkdir -p "${DB}"

# --- chrom sizes from the genome fasta index ---
# samtools faidx writes ${HG38_FASTA}.fai next to the fasta; reuse it if it already exists.
if [[ ! -f "${HG38_FAI}" ]]; then
    echo "[$(date)] indexing genome fasta -> ${HG38_FAI}"
    samtools faidx "${HG38_FASTA}"
else
    echo "[$(date)] reusing existing fasta index ${HG38_FAI}"
fi
cut -f1,2 "${HG38_FAI}" > "${HG38_CHROMSIZES}"
echo "  chrom sizes: $(wc -l < "${HG38_CHROMSIZES}") contigs -> ${HG38_CHROMSIZES}"

# --- ENCODE hg38 blacklist v2 (Boyle-Lab) ---
echo "[$(date)] downloading ENCODE hg38 blacklist v2 -> ${HG38_BLACKLIST}"
wget -qO- "${BLACKLIST_URL}" | gunzip -c > "${HG38_BLACKLIST}"
echo "  blacklist: $(wc -l < "${HG38_BLACKLIST}") regions"

# --- gene annotation via SCENIC+ (biomart) ---
# NOTE: biomart returns Ensembl-style chrom names (1, MT, GL...) and does not convert to
# UCSC when NCBI is unreachable. We normalize to UCSC (chr-prefixed) + primary chromosomes
# below so the annotation matches the chr-prefixed fragments/peaks/fasta used downstream.
echo "[$(date)] fetching SCENIC+ gene annotation (hsapiens) -> ${GENE_ANNOT_RAW}"
conda run --no-capture-output -n epics scenicplus prepare_data download_genome_annotations \
    --species hsapiens \
    --genome_annotation_out_fname "${GENE_ANNOT_RAW}" \
    --chromsizes_out_fname "${CHROMSIZES_UCSC}" || true   # chromsizes step may skip (NCBI); annot still written

# Normalize Chromosome col: 1..22/X/Y -> chr*, MT -> chrM, drop scaffolds/patches/alts.
echo "[$(date)] normalizing gene annotation to UCSC + primary chroms -> ${GENE_ANNOT}"
awk -F'\t' 'BEGIN{OFS="\t";
             split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y",a," ");
             for(i in a) keep[a[i]]=1}
            NR==1{print; next}
            { c=$1;
              if(c=="MT") c="chrM";
              else if(c in keep) c="chr" c;
              else next;
              $1=c; print }' "${GENE_ANNOT_RAW}" > "${GENE_ANNOT}"
echo "  gene annotation: $(($(wc -l < "${GENE_ANNOT}") - 1)) rows (from $(($(wc -l < "${GENE_ANNOT_RAW}") - 1)) raw)"
echo "  chroms kept: $(cut -f1 "${GENE_ANNOT}" | tail -n +2 | sort -u | tr '\n' ' ')"

echo "[$(date)] done. references:"
ls -lh "${HG38_CHROMSIZES}" "${HG38_BLACKLIST}" "${GENE_ANNOT}"
