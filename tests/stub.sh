#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/data"
touch "$tmp/data/a.bam" "$tmp/data/b.bam" "$tmp/data/c.bam" "$tmp/data/r1.fq.gz" "$tmp/data/r2.fq.gz" "$tmp/data/short.cram" "$tmp/data/reference.fa" "$tmp/data/transcripts.fa"
printf 'sample,bam,short_read_1,short_read_2,short_read_cram,cram_reference\na,a.bam,,,short.cram,reference.fa\na,b.bam,,,short.cram,reference.fa\nb,c.bam,,,,\nc,a.bam,r1.fq.gz,r2.fq.gz,,\n' > "$tmp/data/samples.csv"
printf 'apptainer.enabled = false\n' > "$tmp/stub.config"
for mode in full skip asmgene; do
  extra=()
  if [[ "$mode" == skip ]]; then extra+=(--skip_flagger --skip_nucflag); fi
  if [[ "$mode" == asmgene ]]; then extra+=(--asmgene_reference "$tmp/data/reference.fa" --asmgene_transcripts "$tmp/data/transcripts.fa"); fi
  nextflow run "$root/main.nf" -c "$tmp/stub.config" -profile test --input "$tmp/data/samples.csv" --outdir "$tmp/$mode" -stub-run -work-dir "$tmp/work" "${extra[@]}"
  test -s "$tmp/$mode/a/summary.tsv"
  test -s "$tmp/$mode/b/summary.tsv"
  test -s "$tmp/$mode/c/summary.tsv"
  test -f "$tmp/$mode/c/qc/hap1.yak.qv.txt"
  test -f "$tmp/$mode/a/qc/hap1.yak.qv.txt"
  test -f "$tmp/$mode/a/hifiasm_extra/a.ec.fa.gz"
  test -f "$tmp/$mode/a/hifiasm_extra/a.ovlp.paf.gz"
  test ! -e "$tmp/$mode/b/qc/hap1.yak.qv.txt"
  test -f "$tmp/$mode/a/qc/hap1.gfastats.txt"
  if [[ "$mode" == asmgene ]]; then test -f "$tmp/$mode/a/asmgene/hap1.asmgene.tsv"; fi
  if [[ "$mode" != skip ]]; then test -f "$tmp/$mode/a/nucflag/nucflag.bed"; fi
  if [[ "$mode" == full ]]; then test -f "$tmp/$mode/a/flagger/flagger_results/prediction_summary_final.tsv"; fi
  if [[ "$mode" == skip ]]; then test ! -e "$tmp/$mode/a/flagger"; fi
done
