#!/bin/bash
#SBATCH --account=PAS1533
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --job-name=secapr_clean-reads
#SBATCH --output=slurm-secapr_clean-reads-%j.out
set -euo pipefail

# Load the Conda env
module load miniconda3/24.1.2-py310 
conda activate /fs/ess/PAS1533/users/jelmer/software/conda/secapr-2.2.8

# Process arguments
indir=$1
outdir=$2

# Other variables
threads=8

# Report
echo "Starting script secapr_clean-reads.sh"
date
echo "Input dir:            $indir"
echo "Output dir:           $outdir"
echo "Working dir:          $PWD"
echo -e "======================\n"

# Create the output dir
mkdir -p "$outdir"

# Run clean_reads
secapr clean_reads \
    --input "$indir" \
    --output "$outdir" \
    --cores "$threads"

#? --read_min              Min. per-sample read count, default 200,000
#? --index {single,double} Specify if single- or double-indexed adapters were
#?                         used for the library preparation (essential
#?                         information in order to interpret the control-file correctly).

# Report
echo "Successfully finished script secapr_clean-reads.sh"
date
