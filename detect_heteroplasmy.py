#!/usr/bin/env python3

import os
import subprocess
import argparse
from pathlib import Path

def run_cmd(cmd, log_msg):
    print(f"[INFO] {log_msg}")
    print(f"[CMD]  {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main(r1, r2, ref, outdir, threads):
    outdir = Path(outdir)
    tmpdir = outdir / "tmp"
    outdir.mkdir(parents=True, exist_ok=True)
    tmpdir.mkdir(parents=True, exist_ok=True)

    # === Step 1: Index reference ===
    run_cmd(f"bwa index {ref}", "Indexing reference for BWA")
    run_cmd(f"samtools faidx {ref}", "Indexing reference for samtools")
    run_cmd(f"lofreq faidx {ref}", "Indexing reference for LoFreq")

    # === Step 2: Align reads to mitochondrial reference ===
    bam_sorted = tmpdir / "aligned.sorted.bam"
    run_cmd(
        f"bwa mem -t {threads} {ref} {r1} {r2} | "
        f"samtools view -bS - | "
        f"samtools sort -@ {threads} -o {bam_sorted}",
        "Aligning reads and sorting BAM"
    )
    run_cmd(f"samtools index {bam_sorted}", "Indexing sorted BAM")

    # === Step 3: Add indel qualities ===
    bam_indel = tmpdir / "aligned.indelqual.bam"
    run_cmd(
        f"lofreq indelqual --dindel -f {ref} -o {bam_indel} {bam_sorted}",
        "Adding indel quality scores"
    )
    run_cmd(f"samtools index {bam_indel}", "Indexing indel-qualified BAM")

    # === Step 4: Variant calling ===
    raw_vcf = outdir / "raw_variants.vcf"
    run_cmd(
        f"lofreq call-parallel --pp-threads {threads} -f {ref} "
        f"-o {raw_vcf} {bam_indel}",
        "Calling variants with LoFreq"
    )

    # === Step 5: Filter for heteroplasmic variants ===
    filtered_vcf = outdir / "filtered_heteroplasmy.vcf"
    run_cmd(
        f"bcftools filter -i 'INFO/AF>0.01 && INFO/AF<0.5 && INFO/DP>100' "
        f"{raw_vcf} -o {filtered_vcf}",
        "Filtering variants for heteroplasmy"
    )

    # === Summary ===
    print("[INFO] Heteroplasmic variant count:")
    run_cmd(f"grep -v '^#' {filtered_vcf} | wc -l", "Counting variants")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect heteroplasmy in Illumina data.")
    parser.add_argument("r1", help="Forward reads FASTQ.gz")
    parser.add_argument("r2", help="Reverse reads FASTQ.gz")
    parser.add_argument("ref", help="Mitochondrial reference genome (FASTA)")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads (default: 8)")
    args = parser.parse_args()

    main(args.r1, args.r2, args.ref, args.outdir, args.threads)
