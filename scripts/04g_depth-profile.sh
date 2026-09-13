#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=depth-profile
#SBATCH --output=slurm-depth-profile-%j.out
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 04g_depth-profile.sh --vcf_dir <dir> --bed <probe_bed> --outdir <outdir>
                            [--buffer <int>] [--binsize <int>] [--n_samples <int>]

Profile FORMAT/DP against distance from the nearest capture-probe edge.

Locus intervals are a probe's BLAST hit plus a fixed buffer on each side, so the
share of a locus that sits in the low-coverage hybridisation shoulder depends on
the probe's length. This script measures that shoulder directly: every called
position is binned by its distance from the nearest edge of the probe hit that
claims it, with negative distances for positions out in the buffer.

--bed must be the UNBUFFERED probe BED (probes_no-buffer.bed); the buffer is added
here, so that probe and buffer positions can be told apart.

Probes are split into a 'short' and a 'long' class at 600 bp, which separates the
two length classes of the bait set (~200 bp and ~1,950 bp).

Output files:
  depth_profile.tsv   probe_class, dist_bin, n_pos, n_called, mean_dp
                      dist_bin is the left edge of the bin, in bp from the nearest
                      probe edge; 0 is the first on-probe base, negative is buffer
  depth_by_region.tsv probe_class, region (probe / buffer), n_pos, n_called,
                      call_frac, mean_dp

Required arguments:
  --vcf_dir         Directory with single-sample .vcf.gz files
  --bed             Unbuffered probe BED file
  -o, --outdir      Output directory

Optional arguments:
  --buffer          Buffer size per side, in bp, as used to build the locus
                    intervals [default: 300]
  --binsize         Distance bin width in bp [default: 50]
  --n_samples       Use only the first N VCFs, for a quick profile
                    [default: all of them]
  -h, --help        Show this help and exit
EOF
}

# Parse options
vcf_dir=
bed=
outdir=
buffer=300
binsize=50
n_samples=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf_dir)      vcf_dir="$2"; shift 2 ;;
        --bed)          bed="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        --buffer)       buffer="$2"; shift 2 ;;
        --binsize)      binsize="$2"; shift 2 ;;
        --n_samples)    n_samples="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$vcf_dir" || -z "$bed" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -d "$vcf_dir" ]] && { echo "Error: VCF dir not found: $vcf_dir" >&2; exit 1; }
[[ ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

mapfile -t vcfs < <(find "$vcf_dir" -maxdepth 1 -name "*.vcf.gz" | sort)
[[ ${#vcfs[@]} -eq 0 ]] && { echo "Error: no .vcf.gz files in $vcf_dir" >&2; exit 1; }
[[ -n "$n_samples" ]] && vcfs=("${vcfs[@]:0:$n_samples}")

# Report
echo "# Starting script 04g_depth-profile.sh"
date
echo "# VCF directory:      $vcf_dir"
echo "# Number of VCFs:     ${#vcfs[@]}"
echo "# Probe BED file:     $bed"
echo "# Buffer per side:    $buffer"
echo "# Bin size:           $binsize"
echo "# Output directory:   $outdir"

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains bcftools and bedtools

mkdir -p "$outdir"/logs

profile_file="$outdir/depth_profile.tsv"
region_file="$outdir/depth_by_region.tsv"

# One pass over all VCFs. Positions are assigned to probes through the same coarse
# 5-kb position index that 04a uses, widened by the buffer so that buffer positions
# find their probe too.
echo -e "\n# Profiling depth against distance from the probe edge..."
for vcf in "${vcfs[@]}"; do
    bcftools query -f '%CHROM\t%POS[\t%GT\t%DP]\n' "$vcf" 2>/dev/null
done |
    awk -F'\t' -v OFS='\t' -v buffer="$buffer" -v binsize="$binsize" \
        -v profile="$profile_file" -v region_out="$region_file" '
    FNR == NR {
        n++; pstart[n] = $2; pend[n] = $3; plen[n] = $3 - $2
        for (b = int(($2 - buffer) / 5000); b <= int(($3 + buffer) / 5000); b++)
            idx[$1 ":" b] = idx[$1 ":" b] " " n
        next
    }
    {
        # Find the probe whose buffered interval contains this position
        hit = 0
        for (b = int(($2 - buffer) / 5000); b <= int(($2 + buffer) / 5000); b++) {
            if (!(($1 ":" b) in idx)) continue
            split(idx[$1 ":" b], cand, " ")
            for (i in cand) {
                if (cand[i] == "") continue
                j = cand[i] + 0
                if ($2 > pstart[j] - buffer && $2 <= pend[j] + buffer) { hit = j; break }
            }
            if (hit) break
        }
        if (!hit) next

        # Distance from the nearest probe edge: 1 at the first and last on-probe
        # base, <= 0 out in the buffer
        dist_left = $2 - pstart[hit]
        dist_right = pend[hit] - $2 + 1
        dist = (dist_left < dist_right ? dist_left : dist_right)

        cls = (plen[hit] < 600 ? "short" : "long")
        reg = (dist > 0 ? "probe" : "buffer")
        # Bin on the on-probe side starts at 0 for the first base; the buffer side
        # is binned away from the edge, so that -1 falls in the [-binsize, -1) bin
        bin = (dist > 0 ? int((dist - 1) / binsize) * binsize \
                        : -int(-dist / binsize) * binsize - binsize)

        pos_n[cls SUBSEP bin]++
        reg_n[cls SUBSEP reg]++
        # A genotype that was not called carries no usable depth
        if ($3 !~ /\./ && $4 != ".") {
            dp_sum[cls SUBSEP bin] += $4; dp_n[cls SUBSEP bin]++
            rdp_sum[cls SUBSEP reg] += $4; rdp_n[cls SUBSEP reg]++
        }
    }
    END {
        print "probe_class", "dist_bin", "n_pos", "n_called", "mean_dp" > profile
        for (k in pos_n) {
            split(k, a, SUBSEP)
            printf "%s\t%d\t%d\t%d\t%s\n", a[1], a[2], pos_n[k], dp_n[k] + 0,
                (dp_n[k] ? sprintf("%.2f", dp_sum[k] / dp_n[k]) : "NA") > profile
        }
        print "probe_class", "region", "n_pos", "n_called", "call_frac", "mean_dp" > region_out
        for (k in reg_n) {
            split(k, a, SUBSEP)
            printf "%s\t%s\t%d\t%d\t%.4f\t%s\n", a[1], a[2], reg_n[k], rdp_n[k] + 0,
                rdp_n[k] / reg_n[k],
                (rdp_n[k] ? sprintf("%.2f", rdp_sum[k] / rdp_n[k]) : "NA") > region_out
        }
    }' "$bed" -

# Sort the profile by class and distance, which awk's hash order does not preserve
for f in "$profile_file"; do
    { head -n 1 "$f"; tail -n +2 "$f" | sort -k1,1 -k2,2n; } > "$f.tmp" && mv "$f.tmp" "$f"
done
{ head -n 1 "$region_file"; tail -n +2 "$region_file" | sort -k1,1 -k2,2; } \
    > "$region_file.tmp" && mv "$region_file.tmp" "$region_file"

# Report
echo -e "\n# Depth by region:"
column -t "$region_file"
echo -e "\n# Listing output files:"
ls -lh "$profile_file" "$region_file"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 04g_depth-profile.sh"
date
