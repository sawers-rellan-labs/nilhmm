#!/usr/bin/env bash
# Build a per-TAXON subset of the SNP50K joint cohort VCF for nilhmm.
#
# Mirrors the RTIGER subset scheme (fit_rtiger_by_taxa.R): sample -> donor_id via
# data/sample_metadata_master.csv (project=="bzea", not a check, donor_id!=NA),
# taxon = substr(donor_id,1,2) in {Zx,Zv,Zd,Zl,Zh}. Subsets the cohort VCF to the
# samples of ONE taxon so nilhmm can start on a small, single-divergence group.
#
# Usage:  bash agent/make_nilhmm_taxon_subset.sh <TAXON>   # e.g. Zh
set -euo pipefail

TAXON=${1:?usage: make_nilhmm_taxon_subset.sh <TAXON e.g. Zh>}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
META=$ROOT/data/sample_metadata_master.csv
# cohort.vcf.gz is the only joint VCF carrying FORMAT/AD (counts) — required for the
# nilHMM count path (the *_cohort.vcf.gz files are GT-only, 73% missing, unusable at 0.4x).
VCF=${VCF:-/Volumes/BZea/bzeaseq/50K/results/joint/cohort.vcf.gz}
OUTDIR=$ROOT/data/nilhmm
mkdir -p "$OUTDIR"

[ -r "$META" ] || { echo "ERROR: missing $META" >&2; exit 1; }
[ -r "$VCF" ]  || { echo "ERROR: cannot read cohort VCF (mount down?): $VCF" >&2; exit 1; }

SAMPLES=$OUTDIR/${TAXON}_samples.txt
# Sample list from metadata. Header:
#   field_row,project,sample,is_check,sample_group,founder_group,donor_id,pedigree
# Teosinte taxa (Zx/Zv/Zd/Zl/Zh): non-check NILs keyed by substr(donor_id,1,2).
# Check groups (B73/Purple): checks have donor_id=NA, so key on sample_group instead
# (negative controls; mirrors RTIGER fit_rtiger_checks.R). B73 control = PN10_SID893
# (see memory b73-control-impostors).
case "$TAXON" in
  B73|Purple)
    awk -F',' -v g="$TAXON" 'NR>1 && $2=="bzea" && toupper($4)=="TRUE" && $5==g {print $3}' \
      "$META" | sort -u > "$SAMPLES" ;;
  *)
    awk -F',' -v tx="$TAXON" 'NR>1 && $2=="bzea" && toupper($4)=="FALSE" && $7!="NA" && $7!="" && substr($7,1,2)==tx {print $3}' \
      "$META" | sort -u > "$SAMPLES" ;;
esac

# keep only samples actually present in the VCF (intersection), preserving VCF presence
bcftools query -l "$VCF" | sort -u > "$OUTDIR/.vcf_samples.tmp"
comm -12 "$SAMPLES" "$OUTDIR/.vcf_samples.tmp" > "$OUTDIR/.kept.tmp"
mv "$OUTDIR/.kept.tmp" "$SAMPLES"; rm -f "$OUTDIR/.vcf_samples.tmp"

N=$(wc -l < "$SAMPLES" | tr -d ' ')
echo "Taxon $TAXON: $N samples present in cohort VCF -> $SAMPLES"
[ "$N" -gt 0 ] || { echo "ERROR: no samples for taxon $TAXON" >&2; exit 1; }

OUTVCF=$OUTDIR/${TAXON}_counts.vcf.gz
# biallelic SNPs only, keep GT+AD (drop the heavy INFO/PL annotations)
bcftools view -S "$SAMPLES" --force-samples -m2 -M2 -v snps "$VCF" 2>/dev/null \
  | bcftools annotate -x '^FORMAT/GT,FORMAT/AD' -Oz -o "$OUTVCF"
bcftools index -t "$OUTVCF"
echo "Wrote $OUTVCF ($(bcftools query -l "$OUTVCF" | wc -l | tr -d ' ') samples, $(bcftools index -n "$OUTVCF") records)"
