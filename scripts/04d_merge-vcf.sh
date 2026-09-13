#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=merge-vcf
#SBATCH --output=slurm-merge-vcf-%j.out
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: merge-vcf.sh --vcf <vcf_file> --bed <bed_file> --outdir <outdir>

Merge multiple single-sample VCF files (.vcf.gz) into a combined multi-sample VCF file.

Required arguments:
  -i/--indir             Input directory with .vcf.gz files
  -o, --outfile          Output merged VCF file

Optional arguments:
  -h, --help        Show this help and exit
EOF
}

# Parse options
indir=
outfile=
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--indir)     indir="$2"; shift 2 ;;
        -o|--outfile)   outfile="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$indir" || -z "$outfile" ]] && { usage >&2; exit 1; }
[[ ! -d "$indir" ]] && { echo "Error: VCF directory not found: $indir" >&2; exit 1; }

# Report
echo "# Starting script merge-vcf.sh"
date
echo
echo "# Input VCF directory: $indir"
echo "# Output merged VCF:   $outfile"
echo

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains bcftools and bedtools

# Create output directory
outdir=$(dirname "$outfile")
mkdir -p "$outdir"/logs

# Merge the VCF file
echo "# Merging VCF files with bcftools merge..."
bcftools merge --threads 8 -Oz -o "$outfile" "$indir"/*.vcf.gz
bcftools index "$outfile"

# Report
echo "# Number of samples in merged VCF:"
bcftools query -l "$outfile" | wc -l
echo -e "\n# Output file:"
ls -lh "$outfile"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script merge-vcf.sh"
date
