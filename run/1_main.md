## 0. Define settings and files

```bash
# Settings
buffer_size=300             # Size around BLAST hit in bp to keep on each side
snp_gap=5                   # SNPs within this distance of an indel are set to missing (N)
max_missing_at_end=0.25     # Maximum proportion of sequences with a base (not N) at the end of the alignment to consider it for trimming
max_missing_sample=0.1      # Maximum proportion of missing data per sample
max_missing_locus=0.2       # Maximum proportion of N in the final consensus sequences to keep a locus
min_locus_length=150        # Minimum length of the final consensus sequences to keep a locus
min_qual=30                 # Minimum QUAL at ALT-bearing sites (monomorphic sites are exempt)
min_ab=0.25                 # Minimum allele balance for heterozygous genotypes
max_ab=0.75                 # Maximum allele balance for heterozygous genotypes
ab_tol=2                    # Max. reads by which AO + RO may differ from DP for a het to be judged on allele balance
min_call_frac=0.25          # Min. proportion of samples with a called genotype for a site to be scored for excess het
max_het_frac=0.6            # A site with a higher proportion of het genotypes is an 'excess-het site' (HWE maximum is 0.5)
max_site_frac=0.005         # A locus with a higher proportion of excess-het sites is dropped as a likely paralog
min_het_sites=3             # ...and it must have at least this many excess-het sites
min_ab_hets=5               # Min. evaluated hets for a sample to count toward a locus's allele-balance test
ab_skew_cutoff=0.4          # A locus is 'skewed' in a sample when the median allele balance there is below this
min_ab_samples=20           # Min. number of counted samples for a locus to be judged on allele balance
max_ab_skew_frac=0.5        # A locus is dropped when it is skewed in a larger fraction of counted samples
sex_chrom=Binsu_Z           # Name of the sex chromosome in the reference
z_policy=all                # Which sex-linked loci to drop: 'all' / 'paralogs' / 'gametologs' / 'none'

# Files
fqdir=data/fastq
ref=data/ref/binsu_all.fna
bed=results/blast_probes/probes_"$buffer_size"bp-buffer.bed
samplesheet=config/nfc-sarek_samplesheet.csv
kept_bed=results/stats/locus_filtering/retained_loci.bed
```

--------------------------------------------------------------------------------

## 1. Create a BED file with probe locations

- Check probe lengths:
 
```bash
seqkit fx2tab -l -n -i data/probes/probe_seqs.fasta > results/blast_probes/probe_lengths.txt
```

### A) and B) BLAST the probes against the reference genome

```bash
sbatch scripts/01a_blastdb.sh -i $ref -o results/blast_probes/blastdb
sbatch scripts/01b_blast.sh \
    --local_db results/blast_probes/blastdb/binsu_all \
    -i data/probes/probe_seqs.fasta \
    -o results/blast_probes \
    --top_n_query 5 \
    --pct_id 90 \
    --pct_qcov 80 \
    --qcov_metric hsp \
    --no_taxinfo
```

### C) Check the BLAST results

`scripts/01c_explore-blast.qmd`

### D) Create a BED file with the probe locations, with 300 bp buffer on either side

(this script includes BLAST output filtering and BED-file creation)

```bash
blast_file=results/blast_probes/blast_out_final.tsv && ls -lh "$blast_file"
conda activate /fs/ess/PAS0471/jelmer/conda/R
Rscript scripts/01d_filter-loci.R \
    --blast-file "$blast_file" \
    --buffer "$buffer_size" \
    --outdir results/blast_probes 
```

--------------------------------------------------------------------------------

## 2. Generate consensus sequences

- Create a sample sheet for nf-core Sarek:
  - Sample IDs are the short form (`I1603`), matching the `sample` column of
    `metadata/metadata-final.tsv`. 
  - For the sample name inside the VCF, Sarek sets it to `<patient>_<sample>`,
    i.e. `I1603_I1603`. Script 02b then reheaders its output VCF to `--sample_id`.

```bash
echo "patient,sample,lane,fastq_1,fastq_2" > "$samplesheet"
find "$fqdir" -name "*fastq.gz" | sort | paste -d, - - |
    sed -E 's@data/fastq/([^/]*_)?([^/_]+)/.*fastq.gz@\2,\2,1,&@' >> "$samplesheet"

# The short IDs must be unique and must all occur in the metadata file:
awk -F, 'NR > 1 {print $2}' "$samplesheet" | sort | uniq -d          # must print nothing
comm -23 <(awk -F, 'NR > 1 {print $2}' "$samplesheet" | sort) \
         <(tail -n +2 metadata/metadata-final.tsv | cut -f1 | sort)  # must print nothing
```

### A) Run nf-core Sarek

- `config/nfc-sarek.yml` points to the BED file with the probe locations.
- `config/nfc-sarek.config` passes `--genotype-qualities` to FreeBayes
  so FreeBayes emits `FORMAT/GQ`.

```bash
sbatch scripts/02a_nfc-sarek.sh \
    --samplesheet "$samplesheet" \
    --ref_fasta "$ref" \
    --outdir results/nfc-sarek \
    --config config/nfc-sarek.config \
    --params config/nfc-sarek.yml

# Check that there are GQ values in the VCFs (otherwise the consensus script will fail):
vcf=results/nfc-sarek/variant_calling/freebayes/I1603/I1603.freebayes.vcf.gz
bcftools query -f '[%GQ\n]' "$vcf" | head # must show numbers, not `.`
```

### B) Collect the Sarek mapping and coverage stats into one table

These are what explain a low-depth sample later on: read count, mapping
rate, duplicate rate, and how much of the sequencing landed on the probes.

- mosdepth was run over the genome in 500 bp windows rather than over the
  probe intervals, so the script intersects those windows with `--bed` to get
  `pct_on_target`, `target_mean_dp` and `enrichment`.
- `target_mean_dp` is not the same as 04a's `mean_dp`: this one averages over
  every target position including the uncovered ones, 04a averages FORMAT/DP
  over the sites where a genotype was called.
- This runs on all 191 samples. The 52 flagged as contaminated are removed by
  03c before any locus is judged; their other per-sample metrics come from
  step 4I.

```bash
sbatch scripts/02c_sarek-stats.sh \
    --sarek_dir results/nfc-sarek \
    --bed results/blast_probes/probes_"$buffer_size"bp-buffer.bed \
    -o results/stats/sarek

# Mean 3.0M reads, 99.2% mapped, 4.6% duplicates; lowest mapping rate 96.3% (I1594)
# Mean 14.2% of bases on target (range 6.2-25.4%), mean on-target depth 34.2x,
# mean enrichment 100x over the genome-wide depth
column -t results/stats/sarek/sarek_stats.tsv | head -5
```

### C) Variant filtering and consensus sequence generation

- This uses `*.freebayes.vcf.gz` and not `*.freebayes.filtered.vcf.gz`.
  The latter is produced by `vcffilter -f 'QUAL > 30'`, but that would remove
  all monomorphic sites. A QUAL filter restricted to records with ALT alleles
  is applied inside script 02b instead.

```bash
for vcf in results/nfc-sarek/variant_calling/freebayes/*/*.freebayes.vcf.gz; do
    smp=$(basename "$vcf" .freebayes.vcf.gz)
    smp=${smp##*_}  # No-op for short IDs; strips the prefix off Sarek output
    sbatch scripts/02b_bcftools-consensus.sh \
        -o results/consensus_init \
        --sample_id "$smp" \
        --vcf "$vcf" \
        --ref "$ref" \
        --bed "$bed" \
        --snp_gap "$snp_gap" \
        --min_qual "$min_qual" \
        --min_ab "$min_ab" \
        --max_ab "$max_ab" \
        --ab_tol "$ab_tol"
done
```

### D) Record the sequencing run/flowcell of each sample

The samples came off two runs, and the capture worked about twice as well in
one of them, so the batch is a confounder for every depth- and
missingness-based metric downstream.
This reads the run and flowcell out of the first read name of each R1 FASTQ file.

```bash
sbatch scripts/02d_seq-batch.sh --fastq_dir data/fastq -o results/stats/seq_batch
#> 95 samples on run 605 (flowcell HCTGLBCXX), 96 on run 705 (H2HLLBCXY)
#> Of the 125 retained: 68 on 605, 57 on 705
#> 605 averages 18.8% of bases on target vs 10.2% for 705, i.e. ~2x the depth
#> per read; 705 was sequenced ~30% deeper, which partly offsets this
```

--------------------------------------------------------------------------------

## 3. Pivot, trim, and filter consensus sequences

### A) Pivot to locus-centric FASTA files

```bash
sbatch scripts/03_pivot-fasta.sh results/consensus_init/fa results/seq/init
```

### B) Remove terminal Ns, which are present due to buffers

```bash
for fa in results/seq/init/*.fa; do
    # Make sure the max. number of jobs is not exceeded (Slurm limit is 1,000):
    while (( $(squeue -u "$USER" -h | wc -l) >= 950 )); do sleep 60; done
    sbatch -t5 scripts/03b_clipkit.sh \
        -i "$fa" \
        -o results/seq/clipped/$(basename "$fa") \
        --mode gappy  \
        --more_opts "--ends_only --gaps "$max_missing_at_end" -gc N -l"
    #? --'ends_only': only trim ends of alignment
    #? - 'gc N': inform clipkit that gap character is N
    #? -l: create alignment log 
done

# List loci where ClipKIT failed to trim:
clipkit_failed=results/seq/clipped/clipkit_failed_loci.txt
find results/seq/clipped -maxdepth 1 -name '*.fa.untrimmed' -printf '%f\n' |
    sed 's/\.untrimmed$//' | sort > "$clipkit_failed"
echo "# Loci that failed ClipKit: $(grep -c . "$clipkit_failed")" # 34
```

### C) Filter samples and loci with too much missing data (N)

Produces filtered FASTA and VCF files,
and removes the samples flagged as possible contaminants by Felipe.
The flagged samples are removed here so they cannot later influence a locus filter
or any statistic pooled over samples.
Their per-sample stats are computed separately in step 4I.

```bash
# Make a list of the samples (52/191) flagged as contaminated in the metadata file:
awk -F'\t' '$6 == "TRUE" {print $1}' metadata/metadata-final.tsv > metadata/contam_samples.txt

# Run the filter-loci script
sbatch scripts/03c_filter-loci.sh \
    --fasta_dir_in results/seq/clipped \
    --vcf_dir_in results/consensus_init/vcf \
    --fasta_dir_out results/seq/filt \
    --vcf_dir_out results/vcf/filt \
    --stats_dir results/stats/locus_filtering \
    --max_missing_sample "$max_missing_sample" \
    --max_missing_locus "$max_missing_locus" \
    --min_locus_length "$min_locus_length" \
    --remove_samples metadata/contam_samples.txt

#> Filtering samples with mean missingness >10%...: Kept samples: 125 // removed samples: 14
#> Filtering loci with mean missingness >20% or length <150 bp...:  Kept loci: 1762 // removed loci: 57
```

### D) Tabulate the ClipKIT trimming per locus

```bash
sbatch scripts/03d_clipkit-stats.sh --clipped_dir results/seq/clipped -o results/stats/locus_filtering
#> Mean 29.1% of each alignment trimmed; 38 loci trimmed by >75%
```

--------------------------------------------------------------------------------

## 4. VCF QC and merging

### A) VCF QC

Tabulate genotypes and depth per sample and locus.
Besides the per-locus count tables that step 5 reads, it writes
`*_sites.tsv.gz`: one row per site whose genotype is not homozygous reference
(`het` / `alt` / `miss`, the last also covering positions with no VCF record).
Hom-ref is the implicit default and is left out to keep the file size manageable.

```bash
for vcf in results/vcf/filt/*.vcf.gz; do
    sbatch scripts/04a_qc-vcf.sh --vcf "$vcf" --bed "$kept_bed" --outdir results/stats/vcf_qc
done
```

### B) Assign sample sex and classify sex-linked (Z) loci

Snakes are ZW, so females are hemizygous over the differentiated part of the Z.
This shows up in depth (females have half the male depth there) and
heterozygosity (female het collapses to ~0, except where a second copy still
maps onto the locus and samples look near-fixed heterozygous).
To distinguish `gametolog` (only the heterogametic group carries excess-het sites)
or `paralog` (both groups do), excess-heterozygosity is used.

```bash
sbatch scripts/04b_sex-and-z-class.sh \
    --qc_dir results/stats/vcf_qc \
    --bed "$kept_bed" \
    --outdir results/stats/sex_z \
    --sex_chrom "$sex_chrom" \
    --z_policy "$z_policy" \
    --min_call_frac "$min_call_frac" \
    --max_het_frac "$max_het_frac" \
    --paralog_site_frac "$max_site_frac" \
    --min_het_sites "$min_het_sites"

#> Classes: ~50 hemizygous, ~5 PAR-like, 4 gametologs (L1970, L489, L918, L879), 1 paralog (L1839)
#> Z:autosome depth ratio is bimodal, ~0.44 (female) vs ~0.81 (male)
#> Sex: 89 female, 36 male; 5 calls are not confident (see the `flags` column)
```

### C) Allele balance of het genotypes

02b blanks hets with allele balance outside `min_ab`-`max_ab`,
but a blanked genotype looks like any other missing one in the filtered VCF.
This reruns the 02b filtering on the raw VCF and records every het that reaches the
allele-balance step, with whether it was kept, blanked, or exempt
(AO + RO more than `ab_tol` reads from DP).
The filter arguments must match those given to 02b.
This is used in the next step to drop loci whose allele balance is skewed in most samples.

```bash
for vcf in results/nfc-sarek/variant_calling/freebayes/*/*.freebayes.vcf.gz; do
    smp=$(basename "$vcf" .freebayes.vcf.gz)
    smp=${smp##*_}  # As in 2C
    [[ -f results/vcf/filt/"$smp".vcf.gz ]] || continue
    sbatch -t 30 scripts/04h_allele-balance.sh \
        --vcf "$vcf" \
        --sample_id "$smp" \
        --bed "$kept_bed" \
        --min_ab "$min_ab" \
        --max_ab "$max_ab" \
        --ab_tol "$ab_tol" \
        --min_qual "$min_qual" \
        --snp_gap "$snp_gap" \
        -o results/stats/allele_balance
done
```

### D) Remove loci with excess heterozygosity or skewed allele balance

`f` is the proportion of samples heterozygous at a site, and a locus is dropped
when too many of its sites have `f` above `$max_het_frac` (0.6). Under HWE the
maximum is 0.5, so `f > 0.6` is an excess, and a collapsed paralog gives a run
of such sites where every sample carries the fixed difference between the two copies.

A locus is also dropped when its allele balance is skewed in most samples: the
median balance of a sample's hets there is below `$ab_skew_cutoff` (0.4) in more
than `$max_ab_skew_frac` (half) of the samples with enough hets. 

```bash
sbatch scripts/04c_filter-excess-het.sh \
    --qc_dir results/stats/vcf_qc \
    --bed "$kept_bed" \
    --vcf_in results/vcf/filt \
    --fasta_in results/seq/filt \
    --fasta_out results/seq/final \
    --vcf_out results/vcf/final \
    --stats_dir results/stats/excess_het \
    --exclude_loci results/stats/sex_z/loci_to_exclude.txt \
    --min_call_frac "$min_call_frac" \
    --max_het_frac "$max_het_frac" \
    --max_site_frac "$max_site_frac" \
    --min_het_sites "$min_het_sites" \
    --ab_dir results/stats/allele_balance \
    --min_ab_hets "$min_ab_hets" \
    --ab_skew_cutoff "$ab_skew_cutoff" \
    --min_ab_samples "$min_ab_samples" \
    --max_ab_skew_frac "$max_ab_skew_frac"

# Check which loci were dropped, and reconsider the thresholds if needed:
cat results/stats/excess_het/excess_het_loci.txt # L1306, L1813, L473, L663
cat results/stats/excess_het/ab_skew_loci.txt
```

### E) Create a combined multi-sample VCF by merging single-sample VCFs

```bash
merged_vcf=results/vcf/merged/vcf_merged.vcf.gz
sbatch scripts/04d_merge-vcf.sh -i results/vcf/final -o "$merged_vcf"
```

### F) Count variable and parsimony-informative sites per locus

Parsimony-informative = each allele is present in at least two samples.

```bash
sbatch scripts/04e_locus-informativeness.sh \
    --qc_dir results/stats/vcf_qc \
    --bed results/stats/excess_het/retained_loci_final.bed \
    -o results/stats/informativeness \
    --min_call_frac "$min_call_frac"

#> 1.62M scored sites: 6.47% variable, 4.28% parsimony-informative
#> Mean 41 informative sites per locus (across all samples); only 1 locus (L379) has none
```

### G) Per-sample Ts/Tv and singleton counts from the merged VCF

Ts/Tv to check for noise: real SNPs should be about 2, random calls about 0.5.

```bash
sbatch scripts/04f_vcf-stats.sh --vcf "$merged_vcf" -o results/stats/vcf_tstv
#> Ts/Tv 1.90-2.11 across the 125 samples (mean 2.012)
#> Singletons average 289/sample; I20937 has 3,587
```

### H) Depth against distance from the capture-probe edge

Locus intervals are a probe hit plus a fixed 300-bp buffer, so how much of a locus
falls in the low-coverage hybridisation shoulder depends on the probe's length.
This bins every called position on its distance from the nearest probe edge
to separate the above from a real difference in capture efficiency.
Takes the unbuffered probe BED, since it adds the buffer itself.

```bash
sbatch scripts/04g_depth-profile.sh \
    --vcf_dir results/vcf/filt \
    --bed results/blast_probes/probes_no-buffer.bed \
    -o results/stats/depth_profile
```

### I) Per-sample stats for the flagged contaminant samples

03c removes the 52 flagged samples before any locus is judged, which leaves them
without most per-sample stat -- so those stats are computed here.

```bash
flagged=metadata/contam_samples.txt   # Written in 3C

# Genotypes, depth and missingness per locus (04a)
while read -r smp; do
    sbatch scripts/04a_qc-vcf.sh \
        --vcf results/consensus_init/vcf/"$smp".vcf.gz \
        --bed "$kept_bed" \
        --outdir results/stats/vcf_qc_flagged
done < "$flagged"

# Allele balance of het genotypes (04h), with the same filter arguments as 02b
for vcf in results/nfc-sarek/variant_calling/freebayes/*/*.freebayes.vcf.gz; do
    smp=$(basename "$vcf" .freebayes.vcf.gz)
    smp=${smp##*_}
    grep -qxF "$smp" "$flagged" || continue
    sbatch -t 30 scripts/04h_allele-balance.sh \
        --vcf "$vcf" \
        --sample_id "$smp" \
        --bed "$kept_bed" \
        --min_ab "$min_ab" \
        --max_ab "$max_ab" \
        --ab_tol "$ab_tol" \
        --min_qual "$min_qual" \
        --snp_gap "$snp_gap" \
        -o results/stats/allele_balance_flagged
done

# Check that all flagged samples have tables before continuing:
ls results/stats/vcf_qc_flagged/*_counts-per-locus.tsv | wc -l # 52
ls results/stats/allele_balance_flagged/*_ab-summary.tsv | wc -l # 52

# Sex calls (04b), scored together with the retained samples
qc_with_flagged=results/stats/vcf_qc_with-flagged
mkdir -p "$qc_with_flagged"
ln -sf "$PWD"/results/stats/vcf_qc/*_{sites.tsv.gz,counts-per-locus.tsv} "$qc_with_flagged"/
ln -sf "$PWD"/results/stats/vcf_qc_flagged/*_{sites.tsv.gz,counts-per-locus.tsv} "$qc_with_flagged"/

sbatch scripts/04b_sex-and-z-class.sh \
    --qc_dir "$qc_with_flagged" \
    --bed "$kept_bed" \
    --outdir results/stats/sex_z_flagged \
    --sex_chrom "$sex_chrom" \
    --z_policy "$z_policy" \
    --min_call_frac "$min_call_frac" \
    --max_het_frac "$max_het_frac" \
    --paralog_site_frac "$max_site_frac" \
    --min_het_sites "$min_het_sites"

# Ts/Tv and singletons (04d + 04f), from a merged VCF of the retained samples'
# final VCFs and the flagged samples' unfiltered ones. The flagged VCFs cover
# every locus, so 04f is restricted to the final loci; outside them, the flagged
# samples would be the only ones with records and everything would be a singleton.
vcf_with_flagged=results/vcf/with-flagged
mkdir -p "$vcf_with_flagged"
ln -sf "$PWD"/results/vcf/final/*.vcf.gz{,.csi} "$vcf_with_flagged"/
while read -r smp; do
    ln -sf "$PWD"/results/consensus_init/vcf/"$smp".vcf.gz{,.csi} "$vcf_with_flagged"/
done < "$flagged"

merged_vcf_flagged=results/vcf/merged/vcf_merged_with-flagged.vcf.gz
sbatch scripts/04d_merge-vcf.sh -i "$vcf_with_flagged" -o "$merged_vcf_flagged"

# After the merge has finished:
sbatch scripts/04f_vcf-stats.sh \
    --vcf "$merged_vcf_flagged" \
    --bed results/stats/excess_het/retained_loci_final.bed \
    -o results/stats/vcf_tstv_flagged
```

--------------------------------------------------------------------------------

## 5. Final QC and stats in R

### A) Build the QC summary tables

Reads per-locus counts and depth from 04a, sex and Z locus classes from 04b,
excess-het and allele-balance locus stats from 04c, per-sample allele balance from
04h, informativeness counts from 04e, Ts/Tv counts from 04f, mapping stats from
02c, and filtering and trimming stats from 03c/03d, and writes four tables to
`results/stats/summaries`. The flagged samples' stats from 4I fill in their own
rows of the per-sample table and nothing else.

```bash
conda activate /fs/ess/PAS0471/jelmer/conda/R
Rscript scripts/05a_prep-qc-tables.R
```

### B) Locus and sample QC report

```bash
quarto render scripts/05b_qc-loci.qmd
```

### C) Genotype count plots by species, land type, state, and sex

```bash
quarto render scripts/05c_genotype-counts.qmd
```
