
## 2026-03-23 -- Coverage drop-off around probes

```bash
ref=data/ref/binsu_all.fna && ls -lh "$ref"
bed=results/blast_probes/probes_sorted.bed && ls -lh "$bed"

# Convert CRAM to BAM
cram=results/geno/nfc-sarek/preprocessing/markduplicates/P0021_FG_I1485/P0021_FG_I1485.md.cram
samtools view -b -T "$ref" -o P0021_FG_I1485.bam "$cram"

# Subtract 1000 from column 2 and add 1000 to column 3, if column 2 is below 0 after subtracting 1000, change it to 0:
awk -F'\t' -v OFS='\t' '{print $1, ($2-1000<0?0:$2-1000), $3+1000, $4}' "$bed" > buffer.bed

# Quantify depth
bedtools coverage -d -a buffer.bed -b P0021_FG_I1485.bam > depth.tsv
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-23 -- Checking the VCF filtering steps manually

```bash
# Manually
bcftools filter -e 'FORMAT/DP < 5 || FORMAT/GQ < 20' --set-GTs . "$vcf" -O u |
    bcftools filter --SnpGap 5 -O u |
    bcftools view -v snps -O z -o TEST.vcf.gz

bcftools filter -e 'FORMAT/DP < 5 || FORMAT/GQ < 20' --set-GTs . -O z -o filt1.vcf.gz "$vcf"
bcftools filter --SnpGap 10 -O z -o filt2.vcf.gz filt1.vcf.gz
bcftools norm -m -any --atomize -O u filt2.vcf.gz | bcftools view -V indels -O z -o filt3.vcf.gz

#? -a, --atomize                   Decompose complex variants (e.g. MNVs become consecutive SNVs)

zgrep -vc "^#" "$vcf" filt1.vcf.gz filt2.vcf.gz filt3.vcf.gz

bcftools index TEST.vcf.gz
bcftools consensus --haplotype I --missing N -f "$ref" TEST.vcf.gz > consensus.fasta
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-22 -- SNPs near indels with BCFtools call

I was calling SNPs with BCFtools call, but when I checked the distribution
of SNPS around indels, there were losts of SNPs near indels.

```bash
# Run BCFtools call
ref=data/ref/binsu_all.fna && ls -lh "$ref"
vcf=results/geno/bcftools/"$smp".vcf.gz && ls -lh "$vcf"
sbatch mcic-scripts/popgenom/bcftools-call.sh \
    --ref_fasta "$ref" \
    --bam_dir results/geno/filt-bam \
    --vcf "$vcf" \
    --allsites

# Check if many SNPs occur near indels - section `IDD` has indel distance
bcftools stats "$vcf" > vcf_stats.txt

grep IDD vcf_stats.txt
```

- Many more SNPs near indels, only trailing off at around 15 bp distance.
  Therefore it's better to use FreeBayes (or GATK), which will perform local
  assembly around indels and thus should be more robust to this issue.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-01

- Andrew's Python script was run like so:

```bash
python ~/python_scripts/Anchored_mapping_consensus.py \
    -r ../../BothropsjararacaGenomeAssemblyFinalFinal.normalized.fasta \
    -bt ../../probes_blast_tab \
    -b ../../02_mapping/P0021_FG_I1698/P0021_FG_I1698_mark_dups.bam \
    -o P0021_FG_I1698.fasta
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-02-23 - SECAPR => Nextflow workflow

- Made a Nextflow workflow based on `SECAPR` `reference_assembly.py`.
  But later decided to not use this, because it doesn't really do any variant calling.

```bash
reads=data/test2
outdir=results/nf-refasm
workdir=/fs/scratch/PAS1533/jelmer/nf-refasm
sbatch workflows/ref_assembly/main.sh \
    --local_wf workflows/nf-probes \
    --reads "$reads" \
    --reference "$ref" \
    --outdir "$outdir" \
    --workdir "$workdir"
```
