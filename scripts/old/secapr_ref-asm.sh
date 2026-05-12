#!/bin/bash
#SBATCH --account=PAS1533
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --job-name=secapr_ref-asm
#SBATCH --output=slurm-secapr_ref-asm-%j.out
set -euo pipefail

# Load the Conda env
module load miniconda3/24.1.2-py310 
conda activate /fs/ess/PAS1533/users/jelmer/software/conda/secapr-2.2.8

# Process arguments
indir=$1            # Input dir with trimmed FASTQ files
ref=$2              # Reference genome FASTA file
outdir=$3           # Output dir for results

# Other variables/constants
threads=8
REF_TYPE="user-ref-lib"  # Type of reference library (see secapr docs for options)

# Report
echo "Starting script secapr_ref-asm.sh"
date
echo "Input dir:            $indir"
echo "Output dir:           $outdir"
echo "Working dir:          $PWD"
echo -e "======================\n"

# Create the output dir
if [[ -d "$outdir" ]]; then
    rmdir "$outdir" || echo "Output dir already exists and is not empty, exiting..." && exit 1
fi

# Run clean_reads
secapr reference_assembly \
    --reads "$indir" \
    --output "$outdir" \
    --reference "$ref" \
    --reference_type "$REF_TYPE" \
    --cores "$threads"

#? --reference_type     "user-ref-lib" enables to input one single
#?                      fasta file created by the user which will be used as a reference library for all samples.

#? --keep_duplicates    Use this flag if you do not want to discard all
#?                      duplicate reads with Picard.

# Report
echo "Successfully finished script secapr_ref-asm.sh"
date
