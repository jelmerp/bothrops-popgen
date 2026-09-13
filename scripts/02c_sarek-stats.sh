#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sarek-stats
#SBATCH --output=slurm-sarek-stats-%j.out
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 02c_sarek-stats.sh --sarek_dir <dir> --bed <bed_file> --outdir <outdir>

Collect the per-sample mapping and coverage metrics that nf-core Sarek already
computed, into one table. These explain why a sample has low depth or high
missingness downstream: a poor mapping rate, a high duplicate rate, or simply
few reads.

Reads four report files per sample from '<sarek_dir>/reports':
  samtools/*/*.md.cram.stats            read counts, mapping and pairing rates,
                                        MQ0 rate, mismatch error rate, insert size
  markduplicates/*/*.md.cram.metrics    Picard PERCENT_DUPLICATION
  mosdepth/*/*.md.mosdepth.summary.txt  genome length and total sequenced bases,
                                        from the 'total' row
  mosdepth/*/*.md.regions.bed.gz        per-window depth, intersected with --bed
                                        to get the on-target figures

ON-TARGET FIGURES
mosdepth was run over the genome in 500 bp windows rather than over the probe
intervals, so the on-target numbers are obtained here by intersecting those
windows with the probe BED:

  pct_on_target   percentage of all sequenced bases that fall inside the probe
                  intervals -- the capture's hit rate
  target_mean_dp  mean depth across the probe intervals, counting uncovered
                  positions as zero
  enrichment      target_mean_dp / genome_mean_dp, i.e. how far the capture
                  concentrated the reads onto the targets

Depth is assumed uniform within a 500 bp window, so a window straddling the edge
of a target contributes its depth in proportion to the overlap. On this data the
window-summed total agrees with mosdepth's own exact total to within 0.001%.

'target_mean_dp' is not the same quantity as the 'mean_dp' that 04a reports, and
the two are expected to differ: this one averages over every position of every
target including the uncovered ones, while 04a averages FORMAT/DP over the sites
where a genotype was actually called.

The probe intervals must not overlap each other; the script checks and refuses to
run if they do, since overlapping targets would be double-counted.

Sample IDs are reduced to the short form (I1603) to match the rest of the
pipeline, by stripping everything up to and including the last underscore.

Output:
  <outdir>/sarek_stats.tsv   One row per sample

Required arguments:
  --sarek_dir       Sarek output dir (the one containing 'reports')
  --bed             BED file with the probe/capture target intervals
  -o, --outdir      Output directory

Optional arguments:
  -h, --help        Show this help and exit
USAGE
}

# Parse options
sarek_dir=
bed=
outdir=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sarek_dir)    sarek_dir="$2"; shift 2 ;;
        --bed)          bed="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$sarek_dir" || -z "$bed" || -z "$outdir" ]] && { usage >&2; exit 1; }
reports_dir="$sarek_dir"/reports
[[ ! -d "$reports_dir" ]] && { echo "Error: reports dir not found: $reports_dir" >&2; exit 1; }
[[ ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

# Overlapping targets would be counted twice in the on-target sum
n_overlap=$(sort -k1,1 -k2,2n "$bed" |
    awk '$1 == pc && $2 < pe {n++} {pc = $1; ps = $2; pe = $3} END {print n + 0}')
[[ "$n_overlap" -gt 0 ]] && {
    echo "Error: $n_overlap overlapping intervals in $bed; merge them first" >&2; exit 1
}

# Report
echo "# Starting script 02c_sarek-stats.sh"
date
echo "# Sarek directory:    $sarek_dir"
echo "# Target BED file:    $bed"
echo "# Target intervals:   $(wc -l < "$bed") ($(awk '{s += $3 - $2} END {print s}' "$bed") bp)"
echo "# Output directory:   $outdir"

shopt -s nullglob
mkdir -p "$outdir"/logs
outfile="$outdir"/sarek_stats.tsv

stats_files=("$reports_dir"/samtools/*/*.md.cram.stats)
[[ ${#stats_files[@]} -eq 0 ]] && { echo "Error: no samtools stats files in $reports_dir/samtools" >&2; exit 1; }
echo "# Samples found:      ${#stats_files[@]}"

echo -e "\n# Collecting per-sample metrics..."
{
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        sample n_reads pct_mapped pct_properly_paired pct_mq0 \
        error_rate insert_size_avg pct_duplication genome_mean_dp \
        pct_on_target target_mean_dp enrichment

    for stats_file in "${stats_files[@]}"; do
        sarek_id=$(basename "$stats_file" .md.cram.stats)
        sample=${sarek_id##*_}   # Long Sarek IDs (P0021_FG_I1603) become I1603

        dup_file="$reports_dir/markduplicates/$sarek_id/$sarek_id.md.cram.metrics"
        dp_file="$reports_dir/mosdepth/$sarek_id/$sarek_id.md.mosdepth.summary.txt"
        regions_file="$reports_dir/mosdepth/$sarek_id/$sarek_id.md.regions.bed.gz"

        # PERCENT_DUPLICATION is the 9th field of the single metrics row, which is
        # the line after the header line starting with 'LIBRARY'
        pct_dup=NA
        [[ -f "$dup_file" ]] && pct_dup=$(
            awk -F'\t' '/^LIBRARY\t/ {getline; printf "%.4f", $9 * 100; exit}' "$dup_file"
        )

        # The 'total' row is the whole-genome summary; 'total_region' repeats it
        # when mosdepth was not given a BED, so match 'total' exactly. Its own
        # mean column is rounded to 2 decimals, so divide the exact counts instead.
        genome_dp=NA
        [[ -f "$dp_file" ]] && genome_dp=$(
            awk -F'\t' '$1 == "total" {printf "%.6f", $3 / $2; exit}' "$dp_file"
        )

        # Intersect the per-window depths with the target intervals. The targets
        # are indexed into 100 kb buckets so each window is only tested against
        # the few that could contain it, rather than against all of them.
        on_target="NA\tNA\tNA"
        [[ -f "$regions_file" ]] && on_target=$(
            zcat "$regions_file" |
            awk -F'\t' -v OFS='\t' -v bed="$bed" -v genome_dp="$genome_dp" '
            BEGIN {
                while ((getline line < bed) > 0) {
                    split(line, f, "\t")
                    n = ++n_targets
                    tchrom[n] = f[1]; tstart[n] = f[2]; tend[n] = f[3]
                    target_bp += f[3] - f[2]
                    for (b = int(f[2] / 100000); b <= int((f[3] - 1) / 100000); b++)
                        idx[f[1] ":" b] = idx[f[1] ":" b] " " n
                }
                close(bed)
            }
            {
                total += $4 * ($3 - $2)
                # A target spanning a bucket boundary is listed in both, so keep
                # track of which ones this window has already been credited with
                seen = ""
                for (b = int($2 / 100000); b <= int(($3 - 1) / 100000); b++) {
                    if (!(($1 ":" b) in idx)) continue
                    split(idx[$1 ":" b], cand, " ")
                    for (i in cand) {
                        if (cand[i] == "") continue
                        k = cand[i]
                        if (index(seen, "," k ",")) continue
                        seen = seen "," k ","
                        lo = ($2 > tstart[k] ? $2 : tstart[k])
                        hi = ($3 < tend[k]   ? $3 : tend[k])
                        if (hi > lo) on += $4 * (hi - lo)
                    }
                }
            }
            END {
                printf "%.3f\t%.2f\t%.1f",
                    (total > 0 ? 100 * on / total : 0),
                    (target_bp > 0 ? on / target_bp : 0),
                    (genome_dp > 0 && target_bp > 0 ? (on / target_bp) / genome_dp : 0)
            }'
        )

        awk -F'\t' -v OFS='\t' -v sample="$sample" \
            -v pct_dup="${pct_dup:-NA}" -v genome_dp="${genome_dp:-NA}" \
            -v on_target="$on_target" '
        /^SN/ {
            sub(/:$/, "", $2)
            v[$2] = $3
        }
        END {
            n = v["raw total sequences"]
            printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\n",
                sample, n,
                (n ? 100 * v["reads mapped"] / n : 0),
                (n ? 100 * v["reads properly paired"] / n : 0),
                (n ? 100 * v["reads MQ0"] / n : 0),
                v["error rate"], v["insert size average"], pct_dup, genome_dp, on_target
        }' "$stats_file"
    done
} > "$outfile"

# Report
echo -e "\n# First rows of the output table:"
head -3 "$outfile" | column -t
echo -e "\n# Summary across samples:"
awk -F'\t' 'NR > 1 {
        n++; reads += $2; map += $3; dup += $8; ont += $10; tdp += $11; enr += $12
        if (min_map == "" || $3 < min_map) { min_map = $3; min_map_smp = $1 }
        if (min_ont == "" || $10 < min_ont) { min_ont = $10; min_ont_smp = $1 }
    }
    END {
        printf "  Samples:                %d\n", n
        printf "  Mean read count:        %.0f\n", reads / n
        printf "  Mean %% mapped:          %.2f\n", map / n
        printf "  Mean %% duplicates:      %.2f\n", dup / n
        printf "  Lowest %% mapped:        %.2f (%s)\n", min_map, min_map_smp
        printf "  Mean %% on target:       %.2f\n", ont / n
        printf "  Lowest %% on target:     %.2f (%s)\n", min_ont, min_ont_smp
        printf "  Mean on-target depth:   %.2f\n", tdp / n
        printf "  Mean enrichment:        %.1fx\n", enr / n
    }' "$outfile"
echo -e "\n# Output file:"
ls -lh "$outfile"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 02c_sarek-stats.sh"
date
