#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=qc-vcf
#SBATCH --output=slurm-qc-vcf-%j.out
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 04a_qc-vcf.sh --vcf <vcf_file> --bed <bed_file> --outdir <outdir>

Tabulate genotypes and depth for a single-sample VCF file (.vcf.gz).

This is the only step that reads genotypes out of the VCFs: 04b and 04c both work
from its output, so that every downstream statistic is counted the same way and
the VCFs are parsed once instead of once per consumer.

Everything comes from one pass over the VCF. Multi-allelic records are merged back
(bcftools norm -m+) and de-duplicated (-d all) first, so that a site is counted
once, with its combined genotype.

Output files, all prefixed with the sample ID:
  _sites.tsv.gz          One row per site whose genotype is NOT homozygous
                         reference: chrom, pos, locus, code (het / alt / miss).
                         Hom-ref is the implicit default and is left out, which is
                         what keeps this file small -- it is ~5% of the full
                         sample-by-site matrix. Consumers reconstruct what they
                         need: at a site, n_het is the number of samples with
                         'het', and n_called is the sample count minus the number
                         with 'miss'.
  _counts-per-locus.tsv  Per-locus genotype counts: locus, n_het, n_alt, n_ref,
                         n_miss, and mean FORMAT/DP over the called genotypes
  _depth.tsv             Mean FORMAT/DP over all non-missing genotypes
  _het_per-locus.tsv     Per-locus count of heterozygous genotypes, in BED layout
  _hom-alt_per-locus.tsv Per-locus count of homozygous-alternate genotypes, ditto
  _hom-ref_per-locus.tsv Per-locus count of homozygous-reference genotypes, ditto

Required arguments:
  --vcf             Input .vcf.gz file
  --bed             BED file with loci/intervals
  -o, --outdir      Output directory

Optional arguments:
  -h, --help        Show this help and exit
EOF
}

# Parse options
vcf=
bed=
outdir=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)          vcf="$2"; shift 2 ;;
        --bed)          bed="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$vcf" || -z "$bed" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -f "$vcf" ]] && { echo "Error: VCF file not found: $vcf" >&2; exit 1; }
[[ ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

# Infer sample ID
sample_id=$(basename "$vcf" .vcf.gz)

# Report
echo "# Starting script 04a_qc-vcf.sh"
date
echo "# Input VCF:          $vcf"
echo "# Sample ID:          $sample_id"
echo "# BED file:           $bed"
echo "# Output directory:   $outdir"

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains bcftools and bedtools

# Create output directory
mkdir -p "$outdir"/logs

sites_file="$outdir/${sample_id}_sites.tsv.gz"
counts_file="$outdir/${sample_id}_counts-per-locus.tsv"
dp_file="$outdir/${sample_id}_depth.tsv"

# One pass over the VCF. Sites are assigned to loci through a coarse position
# index over the BED, so that each site is only compared against the few loci that
# can possibly contain it; sites outside every locus are dropped, which applies the
# BED restriction at the same time.
echo -e "\n# Tabulating genotypes and depth..."
bcftools norm -m+ "$vcf" 2>/dev/null |
    bcftools norm -d all 2>/dev/null |
    bcftools query -f '%CHROM\t%POS[\t%GT\t%DP]\n' 2>/dev/null |
    awk -F'\t' -v OFS='\t' \
        -v sample="$sample_id" -v counts="$counts_file" -v dpa="$dp_file" '
    FNR == NR {
        for (b = int($2 / 5000); b <= int($3 / 5000); b++) {
            idx[$1 ":" b] = idx[$1 ":" b] " " $2 "," $3 "," $4
        }
        loci[$4] = 1; lchrom[$4] = $1; lstart[$4] = $2; lend[$4] = $3
        next
    }
    # A position inside a locus that has no VCF record at all is as uncalled as one
    # with a missing genotype, and consumers divide by the number of called samples,
    # so both have to reach the sites file. Records arrive sorted and loci do not
    # overlap, so all records of a locus are contiguous: the positions seen in the
    # current locus are tracked, and the rest are emitted when the locus ends.
    function flush_locus(l,   p) {
        if (l == "") return
        for (p = lstart[l] + 1; p <= lend[l]; p++)
            if (!((l SUBSEP p) in seen)) { print lchrom[l], p, l, "miss"; n_miss[l]++ }
        delete seen
    }
    {
        bucket = $1 ":" int($2 / 5000)
        if (!(bucket in idx)) next
        locus = ""
        split(idx[bucket], cand, " ")
        for (i in cand) {
            if (cand[i] == "") continue
            split(cand[i], c, ",")
            if ($2 > c[1] + 0 && $2 <= c[2] + 0) { locus = c[3]; break }
        }
        if (locus == "") next

        if (locus != cur) { flush_locus(cur); cur = locus }
        seen[locus, $2] = 1

        if ($3 ~ /\./) {
            code = "miss"
            n_miss[locus]++
        } else {
            split($3, g, /[\/|]/)
            if (g[1] != g[2])      { code = "het"; n_het[locus]++ }
            else if (g[1] + 0 > 0) { code = "alt"; n_alt[locus]++ }
            else                   { code = "ref"; n_ref[locus]++ }
            # Depth is only meaningful where a genotype was called
            if ($4 != ".") { dp_sum[locus] += $4; dp_n[locus]++; dp_all += $4; dp_all_n++ }
        }
        # Hom-ref is the implicit default and is not written out
        if (code != "ref") print $1, $2, locus, code
    }
    END {
        flush_locus(cur)
        # Loci with no record at all are uncalled over their whole length
        for (l in loci) if (l != cur && !(l in n_het) && !(l in n_ref) && !(l in n_alt) && !(l in n_miss))
            for (p = lstart[l] + 1; p <= lend[l]; p++) { print lchrom[l], p, l, "miss"; n_miss[l]++ }

        print "locus", "n_het", "n_alt", "n_ref", "n_miss", "mean_dp", "n_dp" > counts
        for (l in loci) {
            print l, n_het[l] + 0, n_alt[l] + 0, n_ref[l] + 0, n_miss[l] + 0,
                  (dp_n[l] ? sprintf("%.2f", dp_sum[l] / dp_n[l]) : "NA"), dp_n[l] + 0 > counts
        }
        if (dp_all_n > 0) printf "%s\t%.2f\n", sample, dp_all / dp_all_n > dpa
        else              printf "%s\tNA\n", sample > dpa
    }' "$bed" - |
    gzip > "$sites_file"

# The per-locus count files keep the BED layout (chrom, start, end, locus, count)
# that 05a reads, so the interval coordinates are joined back on here.
echo -e "\n# Writing per-locus genotype counts in BED layout..."
for pair in "het:2" "hom-alt:3" "hom-ref:4"; do   # columns of the counts file
    name=${pair%%:*} col=${pair##*:}
    awk -F'\t' -v OFS='\t' -v col="$col" '
        FNR == NR { if (FNR > 1) cnt[$1] = $col; next }
        { print $1, $2, $3, $4, cnt[$4] + 0 }' "$counts_file" "$bed" \
        > "$outdir/${sample_id}_${name}_per-locus.tsv"
done

# Report
echo -e "\n# Sites with a non-hom-ref genotype:"
zcat "$sites_file" | awk -F'\t' '{n[$4]++} END {for (c in n) printf "  %-6s %d\n", c, n[c]}'
echo -e "\n# Mean counts per locus:"
for f in "$outdir/${sample_id}"_*_per-locus.tsv; do
    mean=$(awk '{sum += $NF; n++} END {if (n) printf "%.2f", sum/n}' "$f")
    echo "  $(basename "$f"): $mean"
done
echo -e "\n# Mean depth over non-missing genotypes:"
cat "$dp_file"

echo -e "\n# Output files:"
ls -lh "$outdir/${sample_id}"_*

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 04a_qc-vcf.sh"
date
