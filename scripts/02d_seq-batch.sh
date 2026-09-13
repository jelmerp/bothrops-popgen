#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=seq-batch
#SBATCH --output=slurm-seq-batch-%j.out
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: 02d_seq-batch.sh --fastq_dir <dir> --outdir <outdir>

Record which sequencing run and flowcell each sample came from, by reading the
read name of the first read of its R1 FASTQ file.

This matters because the capture efficiency differs sharply between the two runs
in this data set, so the sequencing batch is a confounder for every depth- and
missingness-based metric downstream. Without it, plots of depth against read
count show two separate trends with no visible explanation.

Illumina read names (Casava >= 1.8) are colon-separated:
  @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:...
so the first three fields give the instrument, the run number, and the flowcell
ID. Only the first read is inspected: one FASTQ file here holds one flowcell
(checked over the leading 400k reads of several files), since the files are
per-sample demultiplexed output of a single run.

The 'seq_batch' column combines run and flowcell ('605_HCTGLBCXX'). The run
number alone repeats across instruments and the flowcell alone can be re-used,
whereas the pair identifies a run unambiguously; it is also the column meant for
grouping and plotting.

Sample IDs are taken from the FASTQ file name and reduced to the short form
(I1603) to match the rest of the pipeline, by stripping everything up to and
including the last underscore.

Output:
  <outdir>/seq_batch.tsv   One row per sample

Required arguments:
  --fastq_dir       Dir with the per-sample FASTQ files (searched recursively)
  -o, --outdir      Output directory

Optional arguments:
  --r1_suffix       Suffix of the R1 files             [default: _R1.fastq.gz]
  -h, --help        Show this help and exit
USAGE
}

# Parse options
fastq_dir=
outdir=
r1_suffix=_R1.fastq.gz
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq_dir)    fastq_dir="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        --r1_suffix)    r1_suffix="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$fastq_dir" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -d "$fastq_dir" ]] && { echo "Error: FASTQ dir not found: $fastq_dir" >&2; exit 1; }

# Report
echo "# Starting script 02d_seq-batch.sh"
date
echo "# FASTQ directory:    $fastq_dir"
echo "# R1 suffix:          $r1_suffix"
echo "# Output directory:   $outdir"

mkdir -p "$outdir"/logs
outfile="$outdir"/seq_batch.tsv

mapfile -t r1_files < <(find -L "$fastq_dir" -name "*$r1_suffix" | sort)
[[ ${#r1_files[@]} -eq 0 ]] &&
    { echo "Error: no *$r1_suffix files in $fastq_dir" >&2; exit 1; }
echo "# Samples found:      ${#r1_files[@]}"

echo -e "\n# Reading the first read name of each R1 file..."
{
    printf '%s\t%s\t%s\t%s\t%s\n' sample instrument run flowcell seq_batch

    for r1 in "${r1_files[@]}"; do
        fastq_id=$(basename "$r1" "$r1_suffix")
        sample=${fastq_id##*_}   # Long IDs (P0021_FG_I1603) become I1603

        # zcat exits non-zero when head closes the pipe early, which `set -o
        # pipefail` would otherwise turn into a script-wide failure
        header=$( { zcat "$r1" || true; } | head -1)
        [[ -z "$header" ]] && { echo "Error: $r1 is empty" >&2; exit 1; }

        # Read names from other pipelines (SRA, Casava < 1.8) have too few
        # fields to locate the flowcell, so fail loudly rather than guess
        n_field=$(awk -F: '{print NF}' <<< "${header%% *}")
        [[ "$n_field" -lt 7 ]] && {
            echo "Error: unexpected read name in $r1 ($n_field colon-separated" \
                 "fields, expected >= 7):" >&2
            echo "  $header" >&2
            exit 1
        }

        awk -F: -v smp="$sample" 'NR == 1 {
            sub(/^@/, "", $1)
            printf "%s\t%s\t%s\t%s\t%s_%s\n", smp, $1, $2, $3, $2, $3
        }' <<< "$header"
    done
} > "$outfile"

# Report
echo -e "\n# First rows of the output table:"
head -3 "$outfile" | column -t
echo -e "\n# Samples per sequencing batch:"
awk -F'\t' 'NR > 1 {n[$5]++} END {for (b in n) printf "  %-20s %d\n", b, n[b]}' \
    "$outfile" | sort
echo -e "\n# Output file:"
ls -lh "$outfile"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 02d_seq-batch.sh"
date
