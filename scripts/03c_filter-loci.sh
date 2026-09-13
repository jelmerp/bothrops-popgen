#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=filter-loci
#SBATCH --output=slurm-filter-loci-%j.out

set -euo pipefail

# Help text
usage() {
    cat <<'EOF'
Usage: filter-loci.sh --fasta_dir <dir> --vcf_dir <dir> --fasta_out <dir> --vcf_out <dir> --stats_dir <dir> --max_missing_sample <float> --max_missing_locus <float> --min_locus_length <int> [--remove_samples <file>]

Filter samples and loci based on missing data (N) and minimum length in FASTA files.
First removes samples exceeding a per-sample missingness threshold, then filters
loci that are too short or exceed a per-locus missingness threshold.
Final VCFs exclude both removed samples and removed loci.

Arguments:
  --fasta_dir_in        Input directory with FASTA files (.fa)
  --vcf_dir_in          Input directory with VCF files (.vcf.gz)
  --fasta_dir_out       Output directory for filtered FASTA files
  --vcf_dir_out         Output directory for filtered VCF files
  --stats_dir           Output directory for missingness stats
  --max_missing_sample  Maximum mean proportion of N per sample (0-1, e.g., 0.5 for 50%)
  --max_missing_locus   Maximum mean proportion of N per locus (0-1, e.g., 0.25 for 25%)
  --min_locus_length    Minimum alignment length in bp to retain a locus (e.g., 100)
  --remove_samples      Optional: file with sample IDs to remove a priori (one per line)
EOF
}

# Options
fasta_dir=
vcf_dir=
fasta_out=
vcf_out=
stats_dir=
max_missing_sample=
max_missing_locus=
min_locus_length=
remove_samples=

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta_dir_in)         fasta_dir="$2"; shift 2 ;;
        --vcf_dir_in)           vcf_dir="$2"; shift 2 ;;
        --fasta_dir_out)        fasta_out="$2"; shift 2 ;;
        --vcf_dir_out)          vcf_out="$2"; shift 2 ;;
        --stats_dir)            stats_dir="$2"; shift 2 ;;
        --max_missing_sample)   max_missing_sample="$2"; shift 2 ;;
        --max_missing_locus)    max_missing_locus="$2"; shift 2 ;;
        --min_locus_length)     min_locus_length="$2"; shift 2 ;;
        --remove_samples)       remove_samples="$2"; shift 2 ;;
        -h|--help)              usage; exit 0 ;;
        *)                      echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$fasta_dir" || -z "$vcf_dir" || -z "$fasta_out" || -z "$vcf_out" || -z "$stats_dir" || -z "$max_missing_sample" || -z "$max_missing_locus" || -z "$min_locus_length" ]]; then
    usage >&2; exit 1
fi

# Report
echo "# Starting script filter-loci.sh"
date
echo "# Input FASTA directory:       $fasta_dir"
echo "# Input VCF directory:         $vcf_dir"
echo "# Output FASTA directory:      $fasta_out"
echo "# Output VCF directory:        $vcf_out"
echo "# Stats directory:             $stats_dir"
echo "# Max sample missing threshold: $max_missing_sample"
echo "# Max locus missing threshold:  $max_missing_locus"
echo "# Min locus length:             $min_locus_length bp"
echo "# A priori sample removal file:  ${remove_samples:-none}"
echo

# Check inputs
[[ ! -d "$fasta_dir" ]] && { echo "Error: FASTA directory not found: $fasta_dir" >&2; exit 1; }
orig_fasta_dir="$fasta_dir"
[[ ! -d "$vcf_dir" ]] && { echo "Error: VCF directory not found: $vcf_dir" >&2; exit 1; }
[[ -n "$remove_samples" && ! -f "$remove_samples" ]] && { echo "Error: Sample removal file not found: $remove_samples" >&2; exit 1; }

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains seqkit and GNU-parallel too

# Convert thresholds to percentages
max_sample_pct=$(awk -v x="$max_missing_sample" 'BEGIN{print x*100}')
max_locus_pct=$(awk -v x="$max_missing_locus" 'BEGIN{print x*100}')

# Enable nullglob
shopt -s nullglob

# Create output directories
mkdir -p "$fasta_out" "$vcf_out" "$stats_dir"/logs

# Output stats files
sample_locus_stats="$stats_dir"/per_sample_locus_missingness.tsv
sample_mean_stats="$stats_dir"/per_sample_mean_missingness.tsv
removed_samples_file="$stats_dir"/removed_samples.txt
locus_stats="$stats_dir"/per_locus_missingness.tsv
locus_lengths="$stats_dir"/per_locus_length.tsv
short_loci_file="$stats_dir"/removed_short_loci.txt
retained_locus_stats="$stats_dir"/per_locus_missingness_filtered.tsv
retained_bed="$stats_dir"/retained_loci.bed


# ===============================================================================
echo "# Step 1: Removing a priori blacklisted samples..."
# ===============================================================================
first_fa=("$fasta_dir"/*.fa)
n_samples_init=$(seqkit seq -n "${first_fa[0]}" 2>/dev/null | wc -l)
echo "# Initial number of samples: $n_samples_init"

if [[ -n "$remove_samples" ]]; then
    n_remove=$(grep -c . "$remove_samples")
    echo "# Samples to remove a priori: $n_remove (from $remove_samples)"
    echo "# A priori removed samples:"
    cat "$remove_samples"

    fasta_apriori="$stats_dir"/fasta_apriori_filtered
    mkdir -p "$fasta_apriori"

    # Exact match on the sample ID. This used to be a '-r' regex anchored at the
    # end, to make the short IDs in the metadata file match the prefixed IDs that
    # Sarek produced ('I1674' vs 'P0021_FG_I1674'). Step 8 filters VCFs by an exact
    # match on the file name, so the prefixed samples were dropped from the FASTA
    # set but kept in the VCF set. Sample IDs are the short form throughout now
    # (see run/1_main.md section 2) and both steps match the same way.
    for fa in "$fasta_dir"/*.fa; do
        out_fa="$fasta_apriori"/$(basename "$fa")
        seqkit grep -v -f "$remove_samples" "$fa" -o "$out_fa" 2>/dev/null
    done

    first_fa_after=("$fasta_apriori"/*.fa)
    n_samples_after=$(seqkit seq -n "${first_fa_after[0]}" 2>/dev/null | wc -l)
    echo "# Samples after a priori removal: $n_samples_after (removed $((n_samples_init - n_samples_after)))"
    fasta_dir="$fasta_apriori"
else
    echo "# No a priori sample removal file provided, skipping"
fi


# ===============================================================================
echo -e "\n# Step 2: Computing per-sample missingness across all loci..."
# ===============================================================================
echo -e "locus\tsample\tpct_N" > "$sample_locus_stats"

for fa in "$fasta_dir"/*.fa; do
    locus=$(basename "$fa" .fa | sed 's/_Binsu.*//')
    seqkit fx2tab -n -i -l -C N "$fa" 2>/dev/null |
        awk -F'\t' -v locus="$locus" -v outfile="$sample_locus_stats" '
        {
            pct = ($2 > 0 ? 100 * $3 / $2 : 100)
            printf "%s\t%s\t%.2f\n", locus, $1, pct >> outfile
        }'
done

# Compute mean missingness per sample across all loci
echo -e "sample\tmean_pct_N" > "$sample_mean_stats"
tail -n +2 "$sample_locus_stats" |
    awk -F'\t' '{sum[$2] += $3; n[$2]++} END {for (s in sum) printf "%s\t%.2f\n", s, sum[s]/n[s]}' |
    sort -t$'\t' -k2,2nr >> "$sample_mean_stats"

echo "# Per-sample-locus stats: $sample_locus_stats"
echo "# Per-sample mean stats:  $sample_mean_stats"


# ===============================================================================
echo -e "\n# Step 3: Filtering samples with mean missingness > ${max_sample_pct}%..."
# ===============================================================================
if [[ -n "$remove_samples" ]]; then
    cp "$remove_samples" "$removed_samples_file"
else
    : > "$removed_samples_file"
fi
kept_samples=0 && removed_samples=0

while IFS=$'\t' read -r sample mean_pct_N; do
    if awk -v m="$mean_pct_N" -v t="$max_sample_pct" 'BEGIN{exit !(m<=t)}'; then
        ((kept_samples += 1))
    else
        echo "$sample" >> "$removed_samples_file"
        ((removed_samples += 1))
    fi
done < <(tail -n +2 "$sample_mean_stats")

echo "# Kept samples: $kept_samples // removed samples: $removed_samples"
if [[ $removed_samples -gt 0 ]]; then
    echo "# Removed samples listed in: $removed_samples_file"
    cat "$removed_samples_file"
fi


# ===============================================================================
echo -e "\n# Step 4: Creating FASTA files without removed samples..."
# ===============================================================================
fasta_intermediate="$stats_dir"/fasta_sample_filtered
mkdir -p "$fasta_intermediate"

if [[ $removed_samples -gt 0 ]]; then
    for fa in "$fasta_dir"/*.fa; do
        out_fa="$fasta_intermediate"/$(basename "$fa")
        seqkit grep -v -f "$removed_samples_file" "$fa" -o "$out_fa" 2>/dev/null
    done
else
    for fa in "$fasta_dir"/*.fa; do
        cp "$fa" "$fasta_intermediate"/
    done
fi

fa_count=$(find "$fasta_intermediate" -maxdepth 1 -type f -name "*.fa" | wc -l)
echo "# Created $fa_count sample-filtered FASTA files in: $fasta_intermediate"


# ===============================================================================
echo -e "\n# Step 5: Computing per-locus missingness from retained samples..."
# ===============================================================================
echo -e "locus\tmean_pct_N" > "$locus_stats"

for fa in "$fasta_intermediate"/*.fa; do
    locus=$(basename "$fa" .fa | sed 's/_Binsu.*//')
    seqkit fx2tab -n -i -l -C N "$fa" 2>/dev/null |
        awk -F'\t' -v locus="$locus" -v locus_file="$locus_stats" '
        {
            pct = ($2 > 0 ? 100 * $3 / $2 : 100)
            sum += pct
            n++
        }
        END {
            if (n > 0) {
                mean = sum / n
                printf "%s\t%.2f\n", locus, mean >> locus_file
            }
        }'
done

echo "# Per-locus stats (from retained samples): $locus_stats"


# ===============================================================================
echo -e "\n# Step 6: Computing locus lengths and blacklisting short loci (< ${min_locus_length} bp)..."
# ===============================================================================
echo -e "locus\tlength_bp" > "$locus_lengths"
: > "$short_loci_file"
short_loci=0

for fa in "$fasta_intermediate"/*.fa; do
    locus=$(basename "$fa" .fa | sed 's/_Binsu.*//')
    len=$(seqkit stats -T "$fa" 2>/dev/null | awk -F'\t' 'NR==2 {print $8}')
    echo -e "$locus\t$len" >> "$locus_lengths"
    if [[ "$len" -lt "$min_locus_length" ]]; then
        echo "$locus" >> "$short_loci_file"
        ((short_loci += 1))
    fi
done

echo "# Locus lengths file: $locus_lengths"
tail -n +2 "$locus_lengths" | awk -F'\t' '
    {a[NR] = $2; sum += $2}
    END {
        n = asort(a)
        median = (n % 2) ? a[(n+1)/2] : (a[n/2] + a[n/2+1]) / 2
        printf "# Locus length (bp): min=%d, mean=%.1f, median=%.1f, max=%d\n", a[1], sum/n, median, a[n]
    }'
echo "# Loci shorter than ${min_locus_length} bp: $short_loci"
if [[ $short_loci -gt 0 ]]; then
    echo "# Short loci listed in: $short_loci_file"
fi


# ===============================================================================
echo -e "\n# Step 7: Filtering loci with mean missingness > ${max_locus_pct}% or length < ${min_locus_length} bp..."
# ===============================================================================
echo -e "locus\tmean_pct_N" > "$retained_locus_stats"
kept_loci=0 && removed_loci=0

while IFS=$'\t' read -r locus mean_pct_N; do
    if grep -qxF "$locus" "$short_loci_file"; then
        ((removed_loci += 1))
        continue
    fi
    if awk -v m="$mean_pct_N" -v t="$max_locus_pct" 'BEGIN{exit !(m<=t)}'; then
        for fa in "$fasta_intermediate"/*"$locus"_*.fa; do
            [[ ! -f "$fa" ]] && break
            cp "$fa" "$fasta_out"/
            echo -e "$locus\t$mean_pct_N" >> "$retained_locus_stats"
            ((kept_loci += 1))
            break
        done
    else
        ((removed_loci += 1))
    fi
done < <(tail -n +2 "$locus_stats")

mean_n=$(tail -n +2 "$retained_locus_stats" | awk -F'\t' '{sum+=$2} END {if (NR>0) printf "%.2f%%", sum/NR; else print "NA"}')
echo "# Retained per-locus stats: $retained_locus_stats"
echo "# Kept loci: $kept_loci // removed loci: $removed_loci"
echo "# Mean % Ns across retained loci: $mean_n"


# ===============================================================================
echo -e "\n# Step 8: Creating filtered VCFs (excluding removed samples and loci)..."
# ===============================================================================

# First create BED file with retained loci coordinates for filtering VCFs
: > "$retained_bed"
for fa in "$fasta_out"/*.fa; do
    base=$(basename "$fa")
    if [[ "$base" =~ _(L[0-9]+)_((Binsu_[^_]+))_([0-9]+)-([0-9]+)\.fa$ ]]; then
        locus_id="${BASH_REMATCH[1]}"
        chrom="${BASH_REMATCH[2]}"
        start="${BASH_REMATCH[4]}"
        end="${BASH_REMATCH[5]}"

        # Adjust coordinates using ClipKit log (accounts for end-trimming)
        log_file="$orig_fasta_dir/${base}.log"
        if [[ -f "$log_file" ]]; then
            trim_from_start=$(awk '$2 == "keep" {print NR-1; exit}' "$log_file")
            trim_from_end=$(awk '$2 == "keep" {last=NR} END {print NR-last}' "$log_file")
            start=$((start + trim_from_start))
            end=$((end - trim_from_end))
        else
            echo "Warning: No ClipKit log found for $base, using original coordinates" >&2
        fi

        printf "%s\t%s\t%s\t%s\n" "$chrom" "$((start - 1))" "$end" "$locus_id" >> "$retained_bed"
    else
        echo "Warning: Could not parse coordinates from FASTA filename: $base" >&2
    fi
done

echo "# First lines of retained loci BED file: $retained_bed"
head -n 5 "$retained_bed"
[[ ! -s "$retained_bed" ]] && { echo "Error: No retained intervals parsed; cannot filter VCFs" >&2; exit 1; }

# Then filter the VCFs
rm -f "$vcf_out"/*.vcf.gz "$vcf_out"/*.vcf.gz.csi "$vcf_out"/*.vcf.gz.tbi

filter_vcf() {
    local vcf="$1" retained_bed="$2" vcf_out="$3" removed_samples_file="$4"
    local sample_id out_vcf
    sample_id=$(basename "$vcf" .vcf.gz)
    if [[ -s "$removed_samples_file" ]] && grep -qxF "$sample_id" "$removed_samples_file"; then
        echo "# Skipping removed sample: $sample_id"
        return
    fi
    out_vcf="$vcf_out/${sample_id}.vcf.gz"
    bcftools view -T "$retained_bed" -Oz -o "$out_vcf" "$vcf" 2>/dev/null
    bcftools index -f -c "$out_vcf"
}
export -f filter_vcf

parallel -j "$SLURM_CPUS_ON_NODE" \
    filter_vcf {} "$retained_bed" "$vcf_out" "$removed_samples_file" \
    ::: "$vcf_dir"/*.vcf.gz

vcf_count=$(find "$vcf_out" -maxdepth 1 -type f -name "*.vcf.gz" | wc -l)
echo "# Created $vcf_count filtered VCF file(s)"


# ===============================================================================
# Clean up intermediate files
rm -rf "$fasta_intermediate"
[[ -n "${fasta_apriori:-}" ]] && rm -rf "$fasta_apriori"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$stats_dir"/logs/ 2>/dev/null || true
    fi
fi

# Report
echo -e "\n# Done with script filter-loci-by-missing.sh"
date
