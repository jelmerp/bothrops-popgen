#!/bin/bash
#SBATCH --account=PAS0471
#SBATCH --time=15
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=qc-consensus
#SBATCH --output=slurm-qc-consensus-%j.out

# Process command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --vcf) vcf="$2"; shift ;;
        --bed) bed="$2"; shift ;;
        -o|--outdir) outdir="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$vcf" || -z "$bed" || -z "$outdir" ]]; then
    echo "Usage: $0 --vcf <path_to_vcf> --bed <path_to_bed> -o|--outdir <path_to_outdir>"
    exit 1
fi
# Check that input files exist
if [[ ! -f "$vcf" ]]; then
    echo "Error: VCF file '$vcf' not found!"
    exit 1
fi
if [[ ! -f "$bed" ]]; then
    echo "Error: BED file '$bed' not found!"
    exit 1
fi

# Load the conda environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools

# Create output directory if it doesn't exist
mkdir -p "$outdir"/logs

# Get sample ID from VCF filename (assuming format: sample_id.vcf.gz)
sample_id=$(basename "$vcf" .vcf.gz)

# Report
date
echo "Running script qc-consensus.sh"
echo "VCF file:         $vcf"
echo "BED file:         $bed"
echo "Sample ID:        $sample_id"
echo "Output directory: $outdir"
echo

# Get the stats
bedtools intersect -c -a "$bed" -b <(bcftools view -g het "$vcf") \
    > "$outdir/${sample_id}_het_per-locus.txt"
bedtools intersect -c -a "$bed" -b <(bcftools view -i 'GT="AA"' "$vcf") \
    > "$outdir/${sample_id}_hom-alt_per-locus.txt"
bedtools intersect -c -a "$bed" -b <(bcftools view -i 'GT="RR"' "$vcf") \
    > "$outdir/${sample_id}_hom-ref_per-locus.txt"

# Report:
echo
date
echo "Done! Listing resulting files:"
ls -lh "$outdir"/"${sample_id}"_*.txt
