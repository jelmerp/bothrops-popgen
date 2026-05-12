--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## Read QC

```bash
for fq in data/fastq/*/*fastq.gz; do
    sbatch mcic-scripts/fq/fastqc.sh -i "$fq" -o results/fastqc/raw
done
sbatch mcic-scripts/fq/fastqc.sh -i results/fastqc/raw -o results/multiqc/raw
```

## Read trimming

```bash
for R1 in data/fastq/*/*R1.fastq.gz; do
    sbatch mcic-scripts/fq/trimgalore.sh -i "$R1" -o results/trimgalore
done
```

## Map reads with bwa

```bash
# Create the BWA index for the reference genome
sbatch mcic-scripts/map/bwa_index.sh -i "$ref" -o results/bwa/index
# Align to the index
for R1 in results/trimgalore/trimmed/*_R1.fastq.gz; do
    sample=$(basename "$R1" | sed 's/_R1.fastq.gz//')
    sbatch mcic-scripts/map/bwa_mem2.sh \
        -i "$R1" \
        -o results/bwa \
        --index_dir results/bwa/index \
        --readgroup "$sample"
done
#> grep -A1 "primary mapped" slurm-bwa-mem2-6897*
```

- TODO: - DUPLICATE MARKING

## Subset BAM files to only keep probe-regions

- TODO: Add buffer around probe regions, e.g. 500 bp on either side

```bash
micromamba activate /fs/ess/PAS0471/conda/bedtools_2.31.1

smp=P0021_FG_I12277

bam_in=results/geno/bwa/"$smp".bam && ls -lh "$bam_in"
bed=results/blast_probes/probes_sorted.bed && ls -lh "$bed"
bam_regionfilt=results/geno/filt-bam/"$smp".bam
bedtools intersect -a "$bam_in" -b "$bed" > "$bam_regionfilt"
samtools index "$bam_regionfilt"
```

## Check depths

```bash
samtools depth "$bam_regionfilt" | awk '{sum+=$3} END {print "Global Average Depth: ", sum/NR}' # 30.1693
samtools coverage "$bam_regionfilt"

# Per-region depths:
samtools bedcov "$bed" "$bam_regionfilt"
# Add column with mean coverage:
samtools bedcov "$bed" "$bam_regionfilt" | awk '{print $0, $5/($3-$2)}' > "$smp"_coverage.tsv

#? Try mosdepth?
#? Try this for coverage stats: <https://joss.theoj.org/papers/10.21105/joss.09774>
```

## Call variants with FreeBayes

- TODO: Do separately for each locus? This will go faster and then the consensus
  sequences are automatically generated per locus.

```bash
bam_regionfilt=results/geno/filt-bam/"$smp".bam && ls -lh "$bam_regionfilt"

sbatch mcic-scripts/popgenom/freebayes.sh \
    --ref_fasta "$ref" \
    --bam "$bam_regionfilt" \
    -o results/geno/freebayes/"$smp".vcf.gz \
    --allsites
```
