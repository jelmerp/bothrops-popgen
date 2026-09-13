#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=120
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vcf-stats
#SBATCH --output=slurm-vcf-stats-%j.out
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 04f_vcf-stats.sh --vcf <merged_vcf> --outdir <outdir> [--bed <bed_file>]

Run 'bcftools stats' on the merged multi-sample VCF from 04d and pull out the
per-sample counts.

The point of this step is Ts/Tv. Real SNPs in a vertebrate genome sit around a
transition:transversion ratio of 2, and a sample whose ratio drifts towards 0.5
-- the ratio expected if the calls were random -- is calling noise, whether from
contamination, low depth, or degraded DNA. It is the cheapest single number for
spotting a bad sample, and nothing else in this pipeline measures it, because
04a records only whether a genotype is het/alt/missing, not which bases are
involved.

Singleton counts come out of the same pass and are worth reading next to Ts/Tv:
a sample carrying far more private alleles than its peers is the other signature
of contamination.

Output:
  <outdir>/per_sample_vcf_stats.tsv   One row per sample: n_ref_hom, n_non_ref_hom,
                                      n_het, n_transitions, n_transversions, ts_tv,
                                      n_indels, mean_dp, n_singletons, n_missing
  <outdir>/bcftools_stats.txt         The full 'bcftools stats' output

Required arguments:
  --vcf             Input merged .vcf.gz file
  -o, --outdir      Output directory

Optional arguments:
  --bed             Only count records inside these intervals (BED, 0-based).
                    Needed when the merged VCF combines samples restricted to
                    different locus sets: outside the shared set, the samples
                    that do have records there would look like singletons.
  -h, --help        Show this help and exit
USAGE
}

# Parse options
vcf=
outdir=
bed=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)          vcf="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        --bed)          bed="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$vcf" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -f "$vcf" ]] && { echo "Error: VCF file not found: $vcf" >&2; exit 1; }
[[ -n "$bed" && ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

# Report
echo "# Starting script 04f_vcf-stats.sh"
date
echo "# Input VCF:          $vcf"
echo "# Output directory:   $outdir"
echo "# Restrict to BED:    ${bed:-no}"

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools

mkdir -p "$outdir"/logs
stats_file="$outdir"/bcftools_stats.txt
outfile="$outdir"/per_sample_vcf_stats.tsv
threads=${SLURM_CPUS_PER_TASK:-4}

echo -e "\n# Running bcftools stats..."
# -T streams through the file rather than jumping via the index, and a file ending
# in .bed is read as 0-based, as 03c and 04c also rely on
bcftools stats --threads "$threads" -s - ${bed:+-T "$bed"} "$vcf" > "$stats_file"

echo -e "\n# Extracting the per-sample counts..."
# PSC lines hold the per-sample counts; the column order is fixed by bcftools and
# is repeated in the '# PSC' header line of the stats file
awk -F'\t' -v OFS='\t' '
BEGIN {
    print "sample", "n_ref_hom", "n_non_ref_hom", "n_het", "n_transitions",
          "n_transversions", "ts_tv", "n_indels", "mean_dp", "n_singletons", "n_missing"
}
$1 == "PSC" {
    ts = $7; tv = $8
    print $3, $4, $5, $6, ts, tv, (tv > 0 ? sprintf("%.4f", ts / tv) : "NA"),
          $9, $10, $11, $14
}' "$stats_file" > "$outfile"

# Report
echo -e "\n# First rows of the output table:"
head -3 "$outfile" | column -t
echo -e "\n# Summary across samples:"
awk -F'\t' 'NR > 1 && $7 != "NA" {
        n++; tstv += $7; single += $10
        if (min == "" || $7 < min) { min = $7; min_smp = $1 }
        if (max == "" || $7 > max) { max = $7; max_smp = $1 }
        if (max_s == "" || $10 > max_s) { max_s = $10; max_s_smp = $1 }
    }
    END {
        printf "  Samples:            %d\n", n
        printf "  Mean Ts/Tv:         %.3f\n", tstv / n
        printf "  Lowest Ts/Tv:       %.3f (%s)\n", min, min_smp
        printf "  Highest Ts/Tv:      %.3f (%s)\n", max, max_smp
        printf "  Mean singletons:    %.0f\n", single / n
        printf "  Most singletons:    %d (%s)\n", max_s, max_s_smp
    }' "$outfile"
echo -e "\n# The 10 samples with the lowest Ts/Tv:"
awk -F'\t' 'NR > 1' "$outfile" | sort -t$'\t' -k7,7g | head -10 |
    awk -F'\t' -v OFS='\t' 'BEGIN {print "sample", "ts_tv", "n_singletons", "mean_dp"}
                            {print $1, $7, $10, $9}' | column -t
echo -e "\n# Output files:"
ls -lh "$outfile" "$stats_file"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 04f_vcf-stats.sh"
date
