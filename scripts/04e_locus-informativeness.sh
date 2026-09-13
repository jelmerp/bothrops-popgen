#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=informativeness
#SBATCH --output=slurm-informativeness-%j.out
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 04e_locus-informativeness.sh --qc_dir <dir> --bed <bed_file> --outdir <outdir>

Count the variable, parsimony-informative and singleton sites of each locus.

Everything else in the pipeline counts how far each sample sits from the
reference, which says nothing about whether a locus separates the samples from
each other. A locus can carry hundreds of differences from the reference and
still be constant across all 125 samples, which makes it useless for anything
downstream. This is the step that measures that.

Like 04b and 04c, this reads the per-site output of 04a rather than the VCFs, so
that genotypes are counted the same way everywhere.

Per site, over the samples that pass the call-rate threshold:
  carriers_alt = n_het + n_alt      samples carrying at least one ALT allele
  carriers_ref = n_called - n_alt   samples carrying at least one REF allele
  variable     both carrier counts > 0     (segregating among the samples)
  informative  both carrier counts >= 2    (parsimony-informative: each of the two
                                            alleles is seen in at least two samples,
                                            so the site can support a grouping)
  singleton    the rarer allele is in exactly one sample

Sites absent from the 04a output are homozygous reference in every sample, so
they are counted as scored and invariant without being read.

Caveat: 04a merges multi-allelic records and records only 'het' / 'alt' / 'miss',
so a site with two ALT alleles is treated here as biallelic. Those are a small
minority and this is a QC measure, not an input to a downstream analysis.

Output:
  <outdir>/per_locus_informativeness.tsv   One row per locus

Required arguments:
  --qc_dir          Directory with the '*_sites.tsv.gz' files from 04a
  --bed             BED file with the loci to score (use the final locus set,
                    results/stats/excess_het/retained_loci_final.bed)
  -o, --outdir      Output directory

Optional arguments:
  --min_call_frac   Min. proportion of samples with a called genotype for a site
                    to be scored                                    [default: 0.25]
  -h, --help        Show this help and exit
USAGE
}

# Parse options
qc_dir=
bed=
outdir=
min_call_frac=0.25
while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc_dir)           qc_dir="$2"; shift 2 ;;
        --bed)              bed="$2"; shift 2 ;;
        -o|--outdir)        outdir="$2"; shift 2 ;;
        --min_call_frac)    min_call_frac="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *)                  echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$qc_dir" || -z "$bed" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -d "$qc_dir" ]] && { echo "Error: QC dir not found: $qc_dir" >&2; exit 1; }
[[ ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

shopt -s nullglob
site_files=("$qc_dir"/*_sites.tsv.gz)
n_samples=${#site_files[@]}
[[ $n_samples -eq 0 ]] && { echo "Error: no '*_sites.tsv.gz' files in $qc_dir" >&2; exit 1; }

# Report
echo "# Starting script 04e_locus-informativeness.sh"
date
echo "# Input QC dir:                $qc_dir"
echo "# Input BED file:              $bed"
echo "# Output directory:            $outdir"
echo "# Number of samples:           $n_samples"
echo "# Min. call fraction per site: $min_call_frac"

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools

mkdir -p "$outdir"/logs
outfile="$outdir"/per_locus_informativeness.tsv
threads=${SLURM_CPUS_PER_TASK:-4}

echo -e "\n# Aggregating per-site genotype counts and scoring loci..."
# One row per non-hom-ref genotype comes in; sorting by locus and position brings
# all samples of a site together, so the counts can be built in a stream.
zcat "${site_files[@]}" |
    sort -S 4G --parallel="$threads" -k3,3 -k2,2n |
    awk -F'\t' -v OFS='\t' -v ns="$n_samples" -v mcf="$min_call_frac" -v bed="$bed" '
    BEGIN {
        while ((getline line < bed) > 0) {
            split(line, f, "\t")
            locus = f[4]; sub(/^Bothrops_I9999_/, "", locus)
            keep[locus] = 1
            chrom[locus] = f[1]; start[locus] = f[2]; end[locus] = f[3]
            len[locus] = f[3] - f[2]
        }
        close(bed)
    }
    # Score the site that just ended, then reset for the next one
    function flush_site(   n_called, c_alt, c_ref, minor) {
        if (key == "") return
        n_called = ns - n_miss
        if (n_called < mcf * ns) { unscored[locus]++; return }
        c_alt = n_het + n_alt
        c_ref = n_called - n_alt
        if (c_alt > 0 && c_ref > 0) {
            n_var[locus]++
            if (c_alt >= 2 && c_ref >= 2) n_pi[locus]++
            minor = (c_alt < c_ref ? c_alt : c_ref)
            if (minor == 1) n_single[locus]++
        }
    }
    {
        if (!($3 in keep)) next
        k = $3 SUBSEP $2
        if (k != key) {
            flush_site()
            key = k; locus = $3; n_het = 0; n_alt = 0; n_miss = 0
        }
        if ($4 == "het")      n_het++
        else if ($4 == "alt") n_alt++
        else                  n_miss++
    }
    END {
        flush_site()
        print "locus", "chrom", "start", "end", "length_bp", "n_samples",
              "n_sites_scored", "n_variable", "n_pi", "n_singleton",
              "prop_variable", "prop_pi"
        for (l in keep) {
            # Positions never seen above are hom-ref in every sample: scored and invariant
            scored = len[l] - unscored[l]
            print l, chrom[l], start[l], end[l], len[l], ns,
                  scored, n_var[l] + 0, n_pi[l] + 0, n_single[l] + 0,
                  (scored > 0 ? sprintf("%.6f", n_var[l] / scored) : "NA"),
                  (scored > 0 ? sprintf("%.6f", n_pi[l] / scored) : "NA")
        }
    }' |
    { IFS= read -r header; echo "$header"; sort -t$'\t' -k1,1; } > "$outfile"

# Report
echo -e "\n# First rows of the output table:"
head -3 "$outfile" | column -t
echo -e "\n# Summary across loci:"
awk -F'\t' 'NR > 1 {
        n++; scored += $7; var += $8; pi += $9; single += $10
        if ($8 == 0) invariant++
        if ($9 == 0) uninformative++
    }
    END {
        printf "  Loci:                          %d\n", n
        printf "  Scored sites:                  %d\n", scored
        printf "  Variable sites:                %d (%.2f%% of scored)\n", var, 100 * var / scored
        printf "  Parsimony-informative sites:   %d (%.2f%% of scored)\n", pi, 100 * pi / scored
        printf "  Singleton sites:               %d\n", single
        printf "  Loci with no variable site:    %d\n", invariant + 0
        printf "  Loci with no informative site: %d\n", uninformative + 0
        printf "  Mean informative sites/locus:  %.1f\n", pi / n
    }' "$outfile"
echo -e "\n# The 10 least informative loci:"
awk -F'\t' 'NR > 1' "$outfile" | sort -t$'\t' -k9,9n | head -10 |
    awk -F'\t' -v OFS='\t' 'BEGIN {print "locus", "length_bp", "n_sites_scored", "n_variable", "n_pi"}
                            {print $1, $5, $7, $8, $9}' | column -t
echo -e "\n# Output file:"
ls -lh "$outfile"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 04e_locus-informativeness.sh"
date
