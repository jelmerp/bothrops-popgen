#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=vcf-depth
#SBATCH --output=slurm-vcf-depth-%j.out

set -euo pipefail

# Help text
usage() {
    cat <<'EOF'
Usage: vcf-depth.sh -i <vcf_dir> -o <output_stats_file>

Compute mean FORMAT/DP for non-missing genotypes in each single-sample VCF
(.vcf.gz) in a directory.
EOF
}

# Options
vcf_dir=
stats_file=
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i)
            vcf_dir="$2"
            shift 2
            ;;
        -o)
            stats_file="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$vcf_dir" || -z "$stats_file" ]]; then
    usage >&2
    exit 1
fi

# Report
echo "# Starting script vcf-depth.sh"
date
echo "# Input VCF directory: $vcf_dir"
echo "# Output stats file: $stats_file"
echo

# Load Conda env
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools

if [[ ! -d "$vcf_dir" ]]; then
    echo "Error: VCF directory not found: $vcf_dir" >&2
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "Error: bcftools not found in PATH" >&2
    exit 1
fi

mkdir -p "$(dirname "$stats_file")"
: > "$stats_file"

shopt -s nullglob
vcf_files=("$vcf_dir"/*.vcf.gz)
shopt -u nullglob

if [[ ${#vcf_files[@]} -eq 0 ]]; then
    echo "Error: No .vcf.gz files found in: $vcf_dir" >&2
    exit 1
fi

for vcf in "${vcf_files[@]}"; do
    sample=$(basename "$vcf" .vcf.gz)
    bcftools query -f '[%GT\t%DP\n]' "$vcf" 2>/dev/null |
        awk -v sample="$sample" '$1 !~ /\./ && $2 != "." {sum += $2; n++} END {if (n) printf "%s\t%.2f\n", sample, sum / n; else printf "%s\tNA\n", sample}' |
        tee -a "$stats_file"
done

# Report
echo -e "\n# Depth statistics written to: $stats_file"
ls -lh "$stats_file"
echo "# Done with script vcf-depth.sh"
date
