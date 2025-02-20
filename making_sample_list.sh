#!/bin/bash

usage() {
  echo "Usage: $0 [-h] [-lsd|-split]"
  echo "  -h                  Display help message"
  echo "  -lsd                Generate sample list for files with pattern _R1_*.fastq and _R2_*.fastq"
  echo "  -split              Generate sample list for files with pattern _R1.fastq and _R2.fastq"
}

generate_sample_list_lsd() {
  touch sample_list.txt
  for r1_file in *_R1_*.fastq; do
    r2_file=${r1_file/_R1_/_R2_}
    sample_name=${r1_file%%_R1_*}
    sample_name=${sample_name//-/_}
    echo "${sample_name} ${r1_file} ${r2_file}" >> sample_list.txt
  done
}

generate_sample_list_split() {
  touch sample_list.txt
  for r1_file in *_R1.fastq; do
    r2_file=${r1_file/_R1.fastq/_R2.fastq}
    sample_name=${r1_file%_R1.fastq}
    sample_name=${sample_name//-/_}
    echo "${sample_name} ${r1_file} ${r2_file}" >> sample_list.txt
  done
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -lsd)
      generate_sample_list_lsd
      exit 0
      ;;
    -split)
      generate_sample_list_split
      exit 0
      ;;
    *)
      echo "Invalid option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

echo "No option provided. Use -h or --help for usage instructions." >&2
exit 1

