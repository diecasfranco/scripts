#!/bin/bash
set -e

# === Usage ===
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <mt_ref.fasta> <output_dir> <threads>"
    exit 1
fi

# === Input Arguments ===
R1="$1"
R2="$2"
REF="$3"
OUTDIR="$4"
THREADS="$5"

# === Derived paths ===
TMPDIR="${OUTDIR}/tmp"
mkdir -p "$OUTDIR" "$TMPDIR"

# === Step 1: Index the reference ===
bwa index "$REF"
samtools faidx "$REF"
lofreq faidx "$REF"

# === Step 2: Align reads to mitochondrial reference ===
bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | \
    samtools view -bS - | \
    samtools sort -@ "$THREADS" -o "${TMPDIR}/aligned.sorted.bam"

samtools index "${TMPDIR}/aligned.sorted.bam"

# === Step 3: Add indel qualities (required by LoFreq) ===
lofreq indelqual --dindel -f "$REF" -o "${TMPDIR}/aligned.indelqual.bam" "${TMPDIR}/aligned.sorted.bam"
samtools index "${TMPDIR}/aligned.indelqual.bam"

# === Step 4: Variant calling with LoFreq ===
lofreq call-parallel --pp-threads "$THREADS" -f "$REF" -o "${OUTDIR}/raw_variants.vcf" "${TMPDIR}/aligned.indelqual.bam"

# === Step 5: Filter variants for heteroplasmy (1% < AF < 50%) ===
bcftools filter -i 'INFO/AF>0.01 && INFO/AF<0.5 && INFO/DP>100' "${OUTDIR}/raw_variants.vcf" -o "${OUTDIR}/filtered_heteroplasmy.vcf"

# === Optional: Summary ===
echo "Heteroplasmic variant count:"
grep -v "^#" "${OUTDIR}/filtered_heteroplasmy.vcf" | wc -l
