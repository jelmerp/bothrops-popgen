#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=filter-excess-het
#SBATCH --output=slurm-filter-excess-het-%j.out

set -euo pipefail

# Help text
usage() {
    cat <<'EOS'
Usage: 04c_filter-excess-het.sh --qc_dir <dir> --bed <file> --vcf_in <dir> --fasta_in <dir> --fasta_out <dir> --vcf_out <dir> --stats_dir <dir> [--exclude_loci <file>] [--ab_dir <dir>] [thresholds]

Remove loci with excess heterozygosity, which is the standard signature of
paralogous copies collapsing onto a single reference locus in target-capture data.

The statistic is computed PER SITE, not per locus. For each site, 'f' is the
fraction of called samples that are heterozygous. Under Hardy-Weinberg the
maximum is 0.5 (at allele frequency 0.5), so a site with f > 0.6 is already in
excess, and a collapsed paralog produces a run of such sites where every sample
carries the fixed difference between the two copies. A locus is flagged when a
large enough fraction of its sites are excess-het sites.

This replaces the earlier per-locus het-density cutoff, which took the top 1% of
the density distribution. That distribution is smooth and unimodal, so the
percentile removed a fixed 1% of loci whether or not any were anomalous, and it
could not tell a genuinely diverse locus from a collapsed paralog. The per-site
statistic can: across all loci its median and 95th percentile are both zero.
See doc/diary.md, 2026-09-11.

Genotypes are read from the tables written by 04a_qc-vcf.sh, which is also what
04b_sex-and-z-class.sh reads, so both filters score excess heterozygosity the same
way. The VCFs are only opened to write the filtered output.

Het density is still computed and reported per locus, but no longer filters.

A second, independent test uses allele balance (optional, see --ab_dir). Where
only one of two collapsed copies carries a variant, a heterozygous genotype shows
it on about a quarter to a third of the reads rather than half, and it does so in
most samples that are heterozygous there. For each sample with enough evaluated
hets at a locus, the locus is 'skewed' in that sample when the median allele
balance of those hets is below --ab_skew_cutoff; a locus is flagged when too large
a fraction of samples are skewed. Scoring per sample keeps a few noisy samples
from dominating a locus, as they would in a pooled share of skewed hets. The
excess-het test needs many samples to be heterozygous at the same sites and the
allele-balance test does not, so it catches paralogs that differ at fewer, or less
widely shared, sites. A locus is removed when it fails either test. The het tables
come from 04h_allele-balance.sh, which also records which hets 02b judged.

Arguments:
  --qc_dir            Input dir with the per-sample tables from 04a_qc-vcf.sh
                      ('*_sites.tsv.gz' and '*_counts-per-locus.tsv')
  --bed               BED file with the retained loci (from 03c_filter-loci.sh)
  --vcf_in            Input dir with per-sample VCF files (from 03c_filter-loci.sh)
  --fasta_in          Input dir with per-locus FASTA files (from 03c_filter-loci.sh)
  --fasta_out         Output dir for FASTA files of retained loci
  --vcf_out           Output dir for VCF files restricted to retained loci
  --stats_dir         Output dir for the per-site/per-locus stats and the blacklist
  --exclude_loci      File with locus IDs (one per line) to drop outright, e.g. the
                      sex-chromosome loci written by 04b_sex-and-z-class.sh. These
                      are removed from the FASTA and VCF sets and are also left out
                      of the stats table.
  --min_call_frac     A site is only scored if at least this fraction of samples
                      has a non-missing genotype there (default: 0.25). Kept low on
                      purpose: the most diverged paralogs map poorly in most samples,
                      so their excess-het sites sit in low-call-rate windows and a
                      stricter threshold discards exactly the loci worth catching.
  --max_het_frac      A site is an 'excess-het site' when f exceeds this (default: 0.6)
  --max_site_frac     A locus is flagged when the fraction of its scored sites that
                      are excess-het sites exceeds this (default: 0.005)
  --min_het_sites     ...and when it has at least this many excess-het sites
                      (default: 3). Guards against short loci where a couple of
                      sites are enough to clear --max_site_frac.
  --ab_dir            Input dir with the per-sample het tables from
                      04h_allele-balance.sh ('*_het-ab.tsv.gz'). Without it, loci
                      are judged on excess heterozygosity alone. Only samples that
                      also have tables in --qc_dir are counted, and hets that 02b
                      did not judge on allele balance ('exempt') are left out.
  --min_ab_hets       A sample counts toward a locus's allele-balance test when it
                      has at least this many evaluated hets there (default: 5)
  --ab_skew_cutoff    ...and it is skewed there when the median allele balance of
                      those hets is below this (default: 0.4)
  --min_ab_samples    A locus is only judged on allele balance when at least this
                      many samples count toward it (default: 20)
  --max_ab_skew_frac  A locus is flagged when the fraction of counted samples in
                      which it is skewed exceeds this (default: 0.5)
EOS
}

# Options
qc_dir=
bed=
vcf_in=
fasta_in=
fasta_out=
vcf_out=
stats_dir=
exclude_loci=
min_call_frac=0.25
max_het_frac=0.6
max_site_frac=0.005
min_het_sites=3
ab_dir=
min_ab_hets=5
ab_skew_cutoff=0.4
min_ab_samples=20
max_ab_skew_frac=0.5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc_dir)           qc_dir="$2"; shift 2 ;;
        --bed)              bed="$2"; shift 2 ;;
        --vcf_in)           vcf_in="$2"; shift 2 ;;
        --fasta_in)         fasta_in="$2"; shift 2 ;;
        --fasta_out)        fasta_out="$2"; shift 2 ;;
        --vcf_out)          vcf_out="$2"; shift 2 ;;
        --stats_dir)        stats_dir="$2"; shift 2 ;;
        --exclude_loci)     exclude_loci="$2"; shift 2 ;;
        --min_call_frac)    min_call_frac="$2"; shift 2 ;;
        --max_het_frac)     max_het_frac="$2"; shift 2 ;;
        --max_site_frac)    max_site_frac="$2"; shift 2 ;;
        --min_het_sites)    min_het_sites="$2"; shift 2 ;;
        --ab_dir)           ab_dir="$2"; shift 2 ;;
        --min_ab_hets)      min_ab_hets="$2"; shift 2 ;;
        --ab_skew_cutoff)   ab_skew_cutoff="$2"; shift 2 ;;
        --min_ab_samples)   min_ab_samples="$2"; shift 2 ;;
        --max_ab_skew_frac) max_ab_skew_frac="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *)                  echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$qc_dir" || -z "$bed" || -z "$vcf_in" || -z "$fasta_in" || -z "$fasta_out" || -z "$vcf_out" || -z "$stats_dir" ]]; then
    usage >&2; exit 1
fi

# Report
echo "# Starting script 04c_filter-excess-het.sh"
date
echo "# Input QC dir:                $qc_dir"
echo "# Input BED file:              $bed"
echo "# Input VCF dir:               $vcf_in"
echo "# Input FASTA dir:             $fasta_in"
echo "# Output FASTA dir:            $fasta_out"
echo "# Output VCF dir:              $vcf_out"
echo "# Stats dir:                   $stats_dir"
echo "# Locus exclusion list:        ${exclude_loci:-none}"
echo "# Min. call fraction per site: $min_call_frac"
echo "# Excess-het site threshold:   f > $max_het_frac"
echo "# Locus flagged when:          >${max_site_frac} of sites are excess-het sites"
echo "#                              AND >= $min_het_sites such sites"
echo "# Allele-balance het tables:   ${ab_dir:-none (allele-balance test skipped)}"
if [[ -n "$ab_dir" ]]; then
    echo "# Sample counts at a locus:    >= $min_ab_hets evaluated hets there"
    echo "# ...and is skewed there when: median allele balance < $ab_skew_cutoff"
    echo "# Locus flagged when:          >${max_ab_skew_frac} of >= $min_ab_samples counted samples are skewed"
fi
echo

# Check inputs
[[ ! -d "$qc_dir" ]]   && { echo "Error: QC dir not found: $qc_dir" >&2; exit 1; }
[[ ! -f "$bed" ]]      && { echo "Error: BED file not found: $bed" >&2; exit 1; }
[[ ! -d "$vcf_in" ]]   && { echo "Error: VCF dir not found: $vcf_in" >&2; exit 1; }
[[ ! -d "$fasta_in" ]] && { echo "Error: FASTA dir not found: $fasta_in" >&2; exit 1; }
[[ -n "$exclude_loci" && ! -f "$exclude_loci" ]] && { echo "Error: exclusion file not found: $exclude_loci" >&2; exit 1; }
[[ -n "$ab_dir" && ! -d "$ab_dir" ]] && { echo "Error: allele-balance dir not found: $ab_dir" >&2; exit 1; }
# Step 6 clears the output VCF dir, so refuse to run if it is also the input dir
[[ "$(realpath -m "$vcf_out")" == "$(realpath -m "$vcf_in")" ]] && { echo "Error: --vcf_out must differ from --vcf_in" >&2; exit 1; }
[[ "$(realpath -m "$fasta_out")" == "$(realpath -m "$fasta_in")" ]] && { echo "Error: --fasta_out must differ from --fasta_in" >&2; exit 1; }

# Load software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains seqkit and GNU-parallel too

shopt -s nullglob
mkdir -p "$fasta_out" "$vcf_out" "$stats_dir"/logs

site_het_file="$stats_dir"/per_site_het.tsv.gz
locus_stats_file="$stats_dir"/per_locus_het_stats.tsv
excess_het_list="$stats_dir"/excess_het_loci.txt
ab_skew_list="$stats_dir"/ab_skew_loci.txt
blacklist="$stats_dir"/removed_loci.txt
retained_bed="$stats_dir"/retained_loci_final.bed

threads=${SLURM_CPUS_PER_TASK:-4}
workdir=$(mktemp -d) && trap 'rm -rf "$workdir"' EXIT


# ===============================================================================
echo "# Step 1: Aggregating per-site heterozygote counts across samples..."
# ===============================================================================
site_files=("$qc_dir"/*_sites.tsv.gz)
n_samples=${#site_files[@]}
[[ $n_samples -eq 0 ]] && { echo "Error: no '*_sites.tsv.gz' files in $qc_dir" >&2; exit 1; }
echo "# Number of samples: $n_samples"

# 04a writes one row per non-hom-ref genotype, so the records are few enough to
# sort and aggregate in a stream instead of holding every site in memory.
# Hom-alt rows say nothing about heterozygosity and are dropped here.
zcat "${site_files[@]}" |
    awk -F'\t' '$4 != "alt"' |
    sort -S 4G --parallel="$threads" -k3,3 -k2,2n |
    awk -F'\t' -v OFS='\t' -v ns="$n_samples" '
    function flush() {
        if (key == "") return
        called = ns - n_miss
        print locus, chrom, pos, called, n_het, (called > 0 ? n_het / called : 0)
    }
    {
        k = $3 SUBSEP $2
        if (k != key) { flush(); key = k; locus = $3; chrom = $1; pos = $2; n_het = 0; n_miss = 0 }
        if ($4 == "het") n_het++; else n_miss++
    }
    END { flush() }' |
    gzip > "$site_het_file"
echo "# Sites with at least one non-reference or missing genotype: $(zcat "$site_het_file" | wc -l)"
echo "# Per-site counts: $site_het_file"


# ===============================================================================
echo -e "\n# Step 2: Aggregating per-site counts into per-locus statistics..."
# ===============================================================================
# Only sites that appear above can be excess-het sites or fall below the call-rate
# threshold; every other position in a locus is called in all samples and
# homozygous reference, so the scored-site count is the locus length minus the
# positions that failed the call-rate threshold.
zcat "$site_het_file" |
    awk -F'\t' -v OFS='\t' \
        -v ns="$n_samples" -v mcf="$min_call_frac" -v mhf="$max_het_frac" '
    FNR == NR { chrom[$4] = $1; start[$4] = $2; end[$4] = $3; len[$4] = $3 - $2; next }
    {
        # Het genotypes feed the reported het density regardless of call rate
        tot_het[$1] += $5
        if ($4 < mcf * ns) { unscored[$1]++; next }
        if ($6 + 0 > mhf)  n_excess[$1]++
        if ($6 + 0 > max_f[$1]) max_f[$1] = $6 + 0
    }
    END {
        print "locus", "chrom", "start", "end", "length_bp", "n_samples",
              "n_sites", "n_excess_sites", "site_frac", "max_het_frac",
              "total_het", "het_density"
        for (l in len) {
            if (len[l] <= 0) continue
            ns_l = len[l] - (unscored[l] + 0)
            ne_l = n_excess[l] + 0
            print l, chrom[l], start[l], end[l], len[l], ns,
                  ns_l, ne_l, (ns_l > 0 ? ne_l / ns_l : 0), max_f[l] + 0,
                  tot_het[l] + 0, (tot_het[l] + 0) / (ns * len[l])
        }
    }' "$bed" - |
    { read -r header; echo "$header"; sort -k9,9gr -k8,8nr; } > "$locus_stats_file"

# Allele-balance skew per locus. The median needs every value of a locus-sample
# pair at once, so the hets are sorted by locus, sample, and balance, and each
# pair is resolved when the next one starts.
ab_locus_file="$workdir"/ab_per_locus.tsv
if [[ -n "$ab_dir" ]]; then
    ab_files=("$ab_dir"/*_het-ab.tsv.gz)
    [[ ${#ab_files[@]} -eq 0 ]] && { echo "Error: no '*_het-ab.tsv.gz' files in $ab_dir" >&2; exit 1; }
    n_ab_used=0
    for f in "${ab_files[@]}"; do
        [[ -f "$qc_dir/$(basename "$f" _het-ab.tsv.gz)_sites.tsv.gz" ]] && ((n_ab_used += 1))
    done
    echo "# Allele-balance tables: ${#ab_files[@]}, of which $n_ab_used for samples in --qc_dir"
    [[ $n_ab_used -eq 0 ]] && { echo "Error: no allele-balance table matches a sample in $qc_dir" >&2; exit 1; }

    for f in "${ab_files[@]}"; do
        smp=$(basename "$f" _het-ab.tsv.gz)
        [[ -f "$qc_dir/${smp}_sites.tsv.gz" ]] || continue
        zcat "$f" | awk -F'\t' -v OFS='\t' -v s="$smp" 'NR > 1 && $8 != "exempt" { print $3, s, $7 }'
    done |
        sort -S 4G --parallel="$threads" -k1,1 -k2,2 -k3,3g |
        awk -F'\t' -v OFS='\t' -v mh="$min_ab_hets" -v cut="$ab_skew_cutoff" '
        function flush(   med) {
            if (key == "") return
            if (n >= mh) {
                med = (n % 2) ? v[(n + 1) / 2] : (v[n / 2] + v[n / 2 + 1]) / 2
                n_samp[locus]++
                if (med < cut) n_skew[locus]++
            }
            split("", v); n = 0
        }
        {
            k = $1 SUBSEP $2
            if (k != key) { flush(); key = k; locus = $1 }
            v[++n] = $3 + 0
        }
        END {
            flush()
            for (l in n_samp) print l, n_samp[l], n_skew[l] + 0, (n_skew[l] + 0) / n_samp[l]
        }' > "$ab_locus_file"
    echo "# Loci with at least one counted sample for the allele-balance test: $(wc -l < "$ab_locus_file")"
fi

# Append the allele-balance columns: NA throughout without --ab_dir, and 0 for a
# locus where no sample had enough hets to count
tmp="$workdir"/stats_ab.tsv
awk -F'\t' -v OFS='\t' -v use_ab="$([[ -n "$ab_dir" ]] && echo 1 || echo 0)" -v abf="$ab_locus_file" '
    BEGIN {
        if (use_ab) while ((getline line < abf) > 0) {
            split(line, a, "\t"); ab[a[1]] = a[2] "\t" a[3] "\t" a[4]
        }
    }
    NR == 1 { print $0, "n_ab_samples", "n_ab_skewed", "ab_skew_frac"; next }
    {
        if (!use_ab)       cols = "NA\tNA\tNA"
        else if ($1 in ab) cols = ab[$1]
        else               cols = "0\t0\t0"
        print $0, cols
    }' "$locus_stats_file" > "$tmp"
mv "$tmp" "$locus_stats_file"

if [[ -n "$exclude_loci" ]]; then
    n_excl=$(grep -c . "$exclude_loci" || true)
    tmp="$workdir"/stats_tmp.tsv
    awk -F'\t' -v E="$exclude_loci" '
        BEGIN { while ((getline x < E) > 0) excl[x] = 1 }
        NR == 1 || !($1 in excl)' "$locus_stats_file" > "$tmp"
    n_dropped=$(( $(wc -l < "$locus_stats_file") - $(wc -l < "$tmp") ))
    mv "$tmp" "$locus_stats_file"
    echo "# Excluded up front: $n_dropped of the $n_excl loci in $exclude_loci"
fi

n_loci=$(( $(wc -l < "$locus_stats_file") - 1 ))
echo "# Loci with statistics: $n_loci"
echo "# Per-locus stats file: $locus_stats_file"

echo "# Distribution of the fraction of excess-het sites per locus:"
tail -n +2 "$locus_stats_file" | awk -F'\t' '{print $9}' | sort -g |
    awk '{a[NR]=$1} END {printf "#   median=%.5f  p90=%.5f  p95=%.5f  p99=%.5f  max=%.5f\n",
          a[int(NR*0.5)], a[int(NR*0.9)], a[int(NR*0.95)], a[int(NR*0.99)], a[NR]}'
echo "# Loci with no excess-het site at all: $(tail -n +2 "$locus_stats_file" | awk -F'\t' '$8 == 0' | wc -l) of $n_loci"


# ===============================================================================
echo -e "\n# Step 3: Flagging loci with excess heterozygosity..."
# ===============================================================================
tail -n +2 "$locus_stats_file" |
    awk -F'\t' -v sf="$max_site_frac" -v mn="$min_het_sites" \
        '$9 + 0 > sf && $8 + 0 >= mn {print $1}' | sort > "$excess_het_list"
n_excess=$(grep -c . "$excess_het_list" || true)
echo "# Loci flagged as excess-het: $n_excess"

if [[ -n "$ab_dir" ]]; then
    tail -n +2 "$locus_stats_file" |
        awk -F'\t' -v mf="$max_ab_skew_frac" -v ms="$min_ab_samples" \
            '$13 != "NA" && $13 + 0 >= ms && $15 + 0 > mf {print $1}' | sort > "$ab_skew_list"
else
    : > "$ab_skew_list"
fi
n_ab_skew=$(grep -c . "$ab_skew_list" || true)
echo "# Loci flagged for skewed allele balance: $n_ab_skew"
echo "#   ...of which also flagged as excess-het: $(comm -12 "$excess_het_list" "$ab_skew_list" | grep -c . || true)"
if [[ $n_ab_skew -gt 0 ]]; then
    echo "# Allele-balance locus list: $ab_skew_list"
    echo "# Flagged loci (locus, length, n_ab_samples, n_ab_skewed, ab_skew_frac, n_excess_sites):"
    awk -F'\t' -v mf="$max_ab_skew_frac" -v ms="$min_ab_samples" \
        'NR > 1 && $13 != "NA" && $13 + 0 >= ms && $15 + 0 > mf {
            printf "#   %-8s %6s %6s %6s %7.3f %6s\n", $1, $5, $13, $14, $15, $8
        }' "$locus_stats_file"
fi

# The blacklist is what gets removed: loci failing either test, plus any up-front
# exclusions
sort -u "$excess_het_list" "$ab_skew_list" ${exclude_loci:+"$exclude_loci"} > "$blacklist"
n_black=$(grep -c . "$blacklist" || true)
echo "# Loci on the removal blacklist (excess-het + allele balance + excluded): $n_black"
if [[ $n_excess -gt 0 ]]; then
    echo "# Excess-het locus list: $excess_het_list"
    echo "# Full removal blacklist: $blacklist"
    echo "# Flagged loci (locus, length, n_sites, n_excess_sites, site_frac, max_f, het_density):"
    awk -F'\t' -v sf="$max_site_frac" -v mn="$min_het_sites" \
        'NR > 1 && $9 + 0 > sf && $8 + 0 >= mn {
            printf "#   %-8s %6s %7s %6s %9.5f %7.3f %10.5f\n", $1, $5, $7, $8, $9, $10, $12
        }' "$locus_stats_file"
fi


# ===============================================================================
echo -e "\n# Step 4: Copying FASTA files of retained loci..."
# ===============================================================================
kept=0 && removed=0
for fa in "$fasta_in"/*.fa; do
    base=$(basename "$fa")
    # FASTA names look like: Bothrops_I9999_L1000_Binsu_ma-1_125612149-125614348.fa
    locus=$(echo "$base" | sed -E 's/.*_(L[0-9]+)_Binsu.*/\1/')
    # A failed match returns the whole filename, which would silently retain the locus
    [[ ! "$locus" =~ ^L[0-9]+$ ]] && { echo "Error: cannot parse a locus ID from $base" >&2; exit 1; }
    if grep -qxF "$locus" "$blacklist"; then
        ((removed += 1))
    else
        cp "$fa" "$fasta_out"/
        ((kept += 1))
    fi
done
echo "# Kept loci: $kept // removed loci: $removed"


# ===============================================================================
echo -e "\n# Step 5: Writing BED file of retained loci..."
# ===============================================================================
tail -n +2 "$locus_stats_file" |
    awk -F'\t' -v OFS='\t' -v bl="$blacklist" '
    BEGIN { while ((getline l < bl) > 0) black[l] = 1 }
    !($1 in black) { print $2, $3, $4, $1 }' |
    sort -k1,1 -k2,2n > "$retained_bed"
echo "# Retained loci in BED: $(wc -l < "$retained_bed")"
[[ ! -s "$retained_bed" ]] && { echo "Error: no retained intervals; refusing to filter VCFs" >&2; exit 1; }


# ===============================================================================
echo -e "\n# Step 6: Filtering VCFs to retained loci..."
# ===============================================================================
rm -f "$vcf_out"/*.vcf.gz "$vcf_out"/*.vcf.gz.csi "$vcf_out"/*.vcf.gz.tbi

filter_vcf() {
    local vcf="$1" retained_bed="$2" vcf_out="$3"
    local sample_id out_vcf
    sample_id=$(basename "$vcf" .vcf.gz)
    out_vcf="$vcf_out/${sample_id}.vcf.gz"
    bcftools view -T "$retained_bed" -Oz -o "$out_vcf" "$vcf"
    bcftools index -f -c "$out_vcf"
}
export -f filter_vcf

parallel -j "$threads" \
    filter_vcf {} "$retained_bed" "$vcf_out" \
    ::: "$vcf_in"/*.vcf.gz

echo "# Created $(find "$vcf_out" -maxdepth 1 -name '*.vcf.gz' | wc -l) filtered VCF file(s)"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$stats_dir"/logs/ 2>/dev/null || true
    fi
fi

# Report
echo -e "\n# Done with script 04c_filter-excess-het.sh"
date
