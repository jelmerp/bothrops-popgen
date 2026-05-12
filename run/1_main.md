

## Settings

```bash
buffer_size=1000 # In bp, on each side of the BLAST hit. 
snp_gap=10 # Minimum distance of SNPs to indels
max_missing_at_end=0.20 # Maximum proportion of sequences with a base (not N) at the end of the alignment to consider it for trimming
max_missing=0.20 # Maximum proportion of N in the final consensus sequences to keep a locus
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## BLAST the probe sequences against the reference genome

- Check probe lengths:
 
```bash
seqkit fx2tab -l -n -i data/probes/probe_seqs.fasta > results/blast_probes/probe_lengths.txt
```

- BLAST the probes against the reference genome:

```bash
sbatch mcic-scripts/align/blastdb_make.sh -i $ref -o results/blast_probes/blastdb
sbatch mcic-scripts/align/blast.sh \
    --local_db results/blast_probes/blastdb/binsu_all \
    -i data/probes/probe_seqs.fasta \
    -o results/blast_probes \
    --top_n_query 5 --pct_id 90 --pct_qcov 80 \
    --qcov_metric hsp --no_taxinfo
```

- Create a BED file with the probe locations, with 500 bp buffer on either side
  (this script includes BLAST output filtering and BED-file creation):

```bash
micromamba activate /fs/ess/PAS0471/jelmer/conda/R
blast_file=results/blast_probes/blast_out_final.tsv && ls -lh "$blast_file"
Rscript scripts/filter-loci.R --blast-file "$blast_file" \
    --outdir results/blast_probes --buffer "$buffer_size"
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## Generating consensus sequences

```bash
# Files
fqdir=data/fastq && ls -lh "$fqdir"
ref=data/ref/binsu_all.fna && ls -lh "$ref"
bed=results/blast_probes/probes_"$buffer_size"bp-buffer.bed && ls -lh "$bed"
samplesheet=config/nfc-sarek_samplesheet.csv
cdir=results/geno/consensus
```

- Create a sample sheet for nf-core Sarek:

```bash
echo "patient,sample,lane,fastq_1,fastq_2" > "$samplesheet"
find "$fqdir" -name "*fastq.gz" | sort | paste -d, - - |
    sed -E 's@data/fastq/([^/]+)/.*fastq.gz@\1,\1,1,&@' >> "$samplesheet"
```

- Run nf-core Sarek:
  - NOTE: `config/nfc-sarek.yml` contains a reference to the BED file with the probe locations,
    so make sure to update that if you change the buffer size or file name.

```bash
sbatch -t 10:00:00 mcic-scripts/popgenom/nfc-sarek.sh \
    --samplesheet "$samplesheet" \
    --ref_fasta "$ref" \
    --outdir results/geno/nfc-sarek \
    --config config/nfc-sarek.config \
    --params config/nfc-sarek.yml
```

- Variant filtering and consensus sequence generation:

```bash
vcfdir=results/geno/nfc-sarek/variant_calling/freebayes
csdir="$cdir"/sample_init
for vcf in "$vcfdir"/*/*.freebayes.vcf.gz; do
    smp=$(basename "$vcf" .freebayes.vcf.gz)
    sbatch mcic-scripts/popgenom/bcftools-consensus.sh -o "$csdir" \
        --sample_id "$smp" --vcf "$vcf" --ref "$ref" --bed "$bed" --snp_gap "$snp_gap"
done

# Check number of variable sites in the consensus sequences:
grep "Variable sites" "$csdir"/logs/slurm*
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## Pivot and trim consensus sequences

- Pivot to locus-centric FASTA files:

```bash
sdir="$cdir"/sample_init/fa && ls -lh "$sdir"
ldir="$cdir"/locus_init
sbatch scripts/pivot-fasta-by-locus.sh "$sdir" "$ldir"
```

- Remove terminal Ns (due to buffers):

```bash
echo "$max_missing_at_end"
clipdir="$cdir"/locus_clipped
for fa in "$ldir"/*.fa; do
    # Make sure it's not submitting too many jobs at once (SLURM limit is 1000):
    while (( $(squeue -u "$USER" -h | wc -l) >= 950 )); do sleep 60; done
    sbatch -t5 mcic-scripts/align/clipkit.sh -i "$fa" -o "$clipdir"/$(basename "$fa") \
        --mode gappy --more_opts "--ends_only --gaps "$max_missing_at_end" -gc N"
done

# Check
grep "Percentage of alignment trimmed" "$clipdir"/logs/slurm-clipkit* |
    cut -d: -f1,3 | sort -k2,2nr | head -n20
```

- Filter loci with too much missing data (N):

```bash
clipdir="$cdir"/locus_clipped && ls -lh "$clipdir"
odir="$cdir"/locus_filtered && mkdir -p "$odir"
max_missing_pct=$(awk -v x="$max_missing" 'BEGIN{print x*100}') && echo "$max_missing_pct"
kept=0 && removed=0
for fa in "$clipdir"/*.fa; do
    locus=$(basename "$fa" .fa | sed 's/_Binsu.*//')
    mean_pct_N=$(seqkit fx2tab -n -i -l -C N "$fa" |
        awk -F'\t' '{pct=($2>0?100*$3/$2:0); s+=pct; n++} END{if(n>0) print s/n; else print 100}')
    echo "$locus: mean % N = $mean_pct_N"
    if awk -v m="$mean_pct_N" -v t="$max_missing_pct" 'BEGIN{exit !(m<=t)}'; then
        cp "$fa" "$odir"/
        ((kept++))
    else
        ((removed++))
    fi
done
echo "Kept loci: $kept // Removed loci: $removed"
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## QC consensus sequences

- For each sample and locus, tabulate the VCF genotype counts (`0/0`, `0/1`, `1/1`, etc.):

```bash
vcfdir="$csdir"/vcf && ls -lh "$vcfdir"
countdir="$cdir"/stats/sample-geno-counts
for vcf in "$vcfdir"/*.vcf.gz; do
    sbatch scripts/qc-consensus.sh --vcf "$vcf" --bed "$bed" --outdir "$countdir"
done
```

- Percentage of missing data (N) in the final consensus sequences:

```bash
micromamba activate /fs/project/PAS0471/jelmer/conda/seqkit
clipdir="$bdir"/locus_clipped && ls -lh "$clipdir"
missing_file="$bdir"/stats/missing_data.tsv
missing_locus_file="$bdir"/stats/missing_data_locus.tsv

> "$missing_file"
for file in "$clipdir"/*.fa; do
    locus=$(basename "$file" .fa) && echo "Processing locus: $locus"
    seqkit fx2tab -n -i -l -C N "$file" |
        awk -F'\t' '{pct=($2>0?100*$3/$2:0); printf "%s\t%d\t%d\t%.2f\n",$1,$2,$3,pct}' |
        awk -v fname="$locus" 'BEGIN{OFS="\t"} {print fname, $0}' \
        >> "$missing_file"
done

# Mean of column 5 (percentage of Ns) for each locus (column 1):
cat "$missing_file" |
    awk -F'\t' -v OFS="\t" '{g[$1]+=$5; count[$1]++} END {for (loc in g) print loc, g[loc]/count[loc]}' \
    > "$missing_locus_file"
awk -F'\t' '{sum+=$2} END {print "Mean % Ns across loci:", sum/NR}' "$missing_locus_file"
sort -k2,2nr "$missing_locus_file" | head
```

- Get depth statistics from VCF:

```bash
vcfdir="$bdir"/sample_init/vcf && ls -lh "$vcfdir"
stats_file="$bdir"/stats/vcf-depth.tsv
sbatch scripts/vcf-depth.sh -i "$vcfdir" -o "$stats_file"
```


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## Check nr of reads mapping to mtDNA

```bash
# Inputs
mt_ref=data/ref/Binsularis_mitogenome.fasta && ls -lh "$mt_ref"
fqdir=data/fastq
idxdir=results/mtDNA/index
bamdir=results/mtDNA/bam && mkdir -p "$bamdir"

# Build BWA-MEM2 index for mtDNA
sbatch mcic-scripts/map/bwa_index.sh -i "$mt_ref" -o "$idxdir"

# Map all samples to mtDNA only
for r1 in "$fqdir"/*/*_R1.fastq.gz; do
    sbatch mcic-scripts/map/bwa_mem2.sh -i "$r1" --index_dir "$idxdir" -o "$bamdir"
done

# Get stats
out=results/mtDNA/mt_mapped_counts.tsv
echo -e "sample\tprimary_mapped" > "$out"
for f in "$bamdir"/flagstat/*.flagstat; do
    s=$(basename "$f" .flagstat)
    n=$(awk '/primary mapped/ {print $1; exit}' "$f")
    echo -e "$s\t$n" >> "$out"
done
# =30x: 161 samples
# =50x: 141 samples
# =100x: 89 samples
```
