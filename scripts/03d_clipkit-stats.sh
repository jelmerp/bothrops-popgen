#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=30
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=clipkit-stats
#SBATCH --output=slurm-clipkit-stats-%j.out
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 03d_clipkit-stats.sh --clipped_dir <dir> --outdir <outdir>

Tabulate how much of each locus the ClipKIT ends-trimming step (03b) removed.

ClipKIT writes one '<locus>.fa.log' per alignment with a line per alignment
column: position, 'keep' or 'trim', a site class, and the gap fraction. That is
all the input needed here -- the alignments themselves are not re-read.

A locus trimmed back by 80-99% started with a buffer that carried almost no
sequence, so the retained part is short and often poorly covered. Those are the
loci to look at alongside the length and missingness filters in 03c.

The site classes come with an important caveat and are therefore prefixed 'pre_':
ClipKIT computed them on the *untrimmed* alignment, over *all* samples that came
out of step 2, before the contaminated and high-missingness samples were removed
in 03c and before the sex-linked and excess-het loci were removed in 04b/04c. Use
them to judge the trimming, not as informativeness measures for the final data
set -- script 04e computes those on the final sample and locus set.

Output:
  <outdir>/per_locus_clipkit.tsv   One row per locus:
                                   locus, n_sites_pre, n_kept, n_trimmed,
                                   pct_trimmed, pre_n_constant, pre_n_singleton,
                                   pre_n_pi, pre_n_other, trim_failed

Loci where ClipKIT failed are reported with trim_failed=TRUE. It signals failure
by writing an empty '.fa.log' next to an empty '<locus>.fa.untrimmed' marker and
copying the alignment through untrimmed, so for those loci n_trimmed is 0 and
n_sites_pre is read off the alignment instead of the log.

Required arguments:
  --clipped_dir     Directory with the ClipKIT output and '.fa.log' files
  -o, --outdir      Output directory

Optional arguments:
  -h, --help        Show this help and exit
USAGE
}

# Parse options
clipped_dir=
outdir=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --clipped_dir)  clipped_dir="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$clipped_dir" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -d "$clipped_dir" ]] && { echo "Error: dir not found: $clipped_dir" >&2; exit 1; }

# Report
echo "# Starting script 03d_clipkit-stats.sh"
date
echo "# ClipKIT directory:  $clipped_dir"
echo "# Output directory:   $outdir"

shopt -s nullglob
mkdir -p "$outdir"/logs
outfile="$outdir"/per_locus_clipkit.tsv

log_files=("$clipped_dir"/*.fa.log)
[[ ${#log_files[@]} -eq 0 ]] && { echo "Error: no '*.fa.log' files in $clipped_dir" >&2; exit 1; }
echo "# ClipKIT logs found: ${#log_files[@]}"

# ClipKIT writes an empty log for a locus it failed to trim, and an empty log has
# no records for awk to key a new locus on, so the two cases are handled apart
nonempty_logs=() empty_logs=()
for log_file in "${log_files[@]}"; do
    if [[ -s "$log_file" ]]; then nonempty_logs+=("$log_file"); else empty_logs+=("$log_file"); fi
done
echo "# Loci ClipKIT failed to trim: ${#empty_logs[@]}"

echo -e "\n# Tabulating per-locus trimming..."
# The locus ID is the 'L<n>' field of the FASTA file name, matching the rest of
# the pipeline. FNR == 1 starts a new locus; the counts are flushed per file.
awk -F' ' -v OFS='\t' '
BEGIN {
    print "locus", "n_sites_pre", "n_kept", "n_trimmed", "pct_trimmed",
          "pre_n_constant", "pre_n_singleton", "pre_n_pi", "pre_n_other",
          "trim_failed"
}
function flush() {
    if (locus == "") return
    printf "%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%s\n",
        locus, n, keep, trim, (n ? 100 * trim / n : 0),
        constant, singleton, pi, other, "FALSE"
}
FNR == 1 {
    flush()
    # e.g. Bothrops_I9999_L1000_Binsu_ma-1_125612849-125613648.fa.log -> L1000
    locus = FILENAME
    sub(/.*\//, "", locus)
    sub(/^Bothrops_I9999_/, "", locus)
    sub(/_Binsu.*/, "", locus)
    n = keep = trim = constant = singleton = pi = other = 0
}
{
    n++
    if ($2 == "keep") keep++; else trim++
    if ($3 == "constant")                  constant++
    else if ($3 == "singleton")            singleton++
    else if ($3 == "parsimony-informative") pi++
    else                                   other++
}
END { flush() }' "${nonempty_logs[@]}" > "$outfile"

# The failed loci pass through untrimmed, so their alignment length is the length
# of the first sequence in the clipped FASTA; ClipKIT logged no site classes
for log_file in "${empty_logs[@]:-}"; do
    [[ -z "$log_file" ]] && continue
    fasta=${log_file%.log}
    locus=$(basename "$fasta" | sed -E 's/^Bothrops_I9999_//; s/_Binsu.*//')
    aln_len=$(awk '/^>/ {if (seen) exit; seen = 1; next} {len += length($0)} END {print len + 0}' "$fasta")
    printf '%s\t%d\t%d\t0\t0.00\tNA\tNA\tNA\tNA\tTRUE\n' "$locus" "$aln_len" "$aln_len"
done >> "$outfile"

# Keep the table in locus order rather than log-then-failure order
{ head -1 "$outfile"; tail -n +2 "$outfile" | sort -t$'\t' -k1,1; } > "$outfile".tmp &&
    mv "$outfile".tmp "$outfile"

# Report
echo -e "\n# First rows of the output table:"
head -3 "$outfile" | column -t
echo -e "\n# Summary across loci:"
awk -F'\t' 'NR > 1 {
        n++; tot_pre += $2; tot_kept += $3; pct += $5
        if ($5 > 75) heavy++
        if ($10 == "TRUE") failed++
    }
    END {
        printf "  Loci:                       %d\n", n
        printf "  Mean %% trimmed:             %.2f\n", pct / n
        printf "  Total bp before / after:    %d / %d\n", tot_pre, tot_kept
        printf "  Loci trimmed by >75%%:       %d\n", heavy + 0
        printf "  Loci where ClipKIT failed:  %d\n", failed + 0
    }' "$outfile"
echo -e "\n# The 10 most heavily trimmed loci:"
awk -F'\t' 'NR > 1' "$outfile" | sort -t$'\t' -k5,5gr | head -10 |
    awk -F'\t' -v OFS='\t' 'BEGIN {print "locus", "n_sites_pre", "n_kept", "pct_trimmed"}
                            {print $1, $2, $3, $5}' | column -t
echo -e "\n# Output file:"
ls -lh "$outfile"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 03d_clipkit-stats.sh"
date
