#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=sex-and-z-class
#SBATCH --output=slurm-sex-and-z-class-%j.out

set -euo pipefail

# Help text
usage() {
    cat <<'EOS'
Usage: 04b_sex-and-z-class.sh --qc_dir <dir> --bed <file> --outdir <dir> [options]

Assign the sex of each sample and classify sex-chromosome (Z) loci, then write a
locus exclusion list for the downstream filter.

Reads the genotype and depth tables written by 04a_qc-vcf.sh; it does not open the
VCFs. The excess-heterozygosity statistic here is the same per-site one that
04c_filter-excess-het.sh uses, applied within each sex group.

Snakes are ZW: females (ZW) are hemizygous over the differentiated part of the Z,
while males (ZZ) are diploid there. This produces two signals that are used here:

  1. Sequencing depth. Females have ~half the depth of males over the
     differentiated Z, but equal depth on autosomes. The Z:autosome depth ratio
     is therefore bimodal, with the two modes a factor of two apart.
  2. Heterozygosity. Females are hemizygous over the differentiated Z, so the
     diploid caller emits (pseudo-)homozygous calls and their het counts collapse
     to ~0. Where a second copy still maps onto the locus, the opposite happens:
     samples show near-fixed heterozygosity at the sites that tell the copies
     apart.

Z loci are classified as:
  gametolog     excess-het sites in the heterogametic group only -- the W copy
                collapsed onto the Z. True paralogs, and female-specific.
  paralog       excess-het sites in both groups -- a collapsed duplication that
                happens to sit on the Z. A female:male het ratio cannot see these,
                because the excess is present in both sexes; they were previously
                misclassified as PAR-like.
  hemizygous    female:male het density at or below --hemi_ratio; differentiated
                Z, females haploid
  PAR-like      none of the above; undifferentiated, behaves autosomally

Sex is called from the Z:autosome depth ratio by splitting the sorted ratios at
their largest internal gap, so no absolute threshold has to be hardcoded. Because
the ratio is cleanest when restricted to hemizygous loci, and hemizygous loci can
only be identified once sex is known, the two steps are iterated: sex is first
bootstrapped from all Z loci, loci are classified, and then sex is recalled using
the hemizygous loci only and the loci reclassified.

Arguments:
  --qc_dir            Input dir with the per-sample tables from 04a_qc-vcf.sh
                      ('*_sites.tsv.gz' and '*_counts-per-locus.tsv')
  --bed               BED file of retained loci (from 03c_filter-loci.sh)
  --outdir            Output dir for sex calls, locus classes and exclusion lists
  --sex_chrom         Sex chromosome name in the reference        [default: Binsu_Z]
  --z_policy          Which loci to write to the exclusion list:  [default: all]
                        'all'        - every locus on the sex chromosome
                        'paralogs'   - the gametolog and paralog loci
                        'gametologs' - only the gametolog loci
                        'none'       - write an empty list
  --min_call_frac     Min. proportion of a group with a called genotype for a site
                      to be scored                                [default: 0.25]
  --max_het_frac      A site above this proportion of het genotypes within a group
                      is an 'excess-het site' (HWE maximum is 0.5)  [default: 0.6]
  --paralog_site_frac Min. proportion of excess-het sites for a locus to count as
                      a paralog in that group                     [default: 0.005]
  --min_het_sites     ...and the minimum number of such sites     [default: 3]
  --hemi_ratio        Max female:male het ratio for hemizygous     [default: 0.6]
  --min_auto_dp       Min autosomal depth to trust a sex call      [default: 15]
  --max_group_z       Max robust (MAD-based) z-score within the assigned sex
                      group for the depth ratio, hemizygous het density and
                      gametolog het count; a sample exceeding it on any of the
                      three is flagged as a non-confident call  [default: 3]
EOS
}

# Options
qc_dir=
bed=
outdir=
sex_chrom=Binsu_Z
z_policy=all
min_call_frac=0.25
max_het_frac=0.6
paralog_site_frac=0.005
min_het_sites=3
hemi_ratio=0.6
min_auto_dp=15
max_group_z=3

while [[ $# -gt 0 ]]; do
    case "$1" in
        --qc_dir)            qc_dir="$2"; shift 2 ;;
        --bed)               bed="$2"; shift 2 ;;
        --outdir)            outdir="$2"; shift 2 ;;
        --sex_chrom)         sex_chrom="$2"; shift 2 ;;
        --z_policy)          z_policy="$2"; shift 2 ;;
        --min_call_frac)     min_call_frac="$2"; shift 2 ;;
        --max_het_frac)      max_het_frac="$2"; shift 2 ;;
        --paralog_site_frac) paralog_site_frac="$2"; shift 2 ;;
        --min_het_sites)     min_het_sites="$2"; shift 2 ;;
        --hemi_ratio)        hemi_ratio="$2"; shift 2 ;;
        --min_auto_dp)       min_auto_dp="$2"; shift 2 ;;
        --max_group_z)       max_group_z="$2"; shift 2 ;;
        -h|--help)           usage; exit 0 ;;
        *)                   echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -z "$qc_dir" || -z "$bed" || -z "$outdir" ]]; then
    usage >&2; exit 1
fi
case "$z_policy" in
    all|paralogs|gametologs|none) ;;
    *) echo "Error: --z_policy must be 'all', 'paralogs', 'gametologs' or 'none'" >&2; exit 1 ;;
esac

# Report
echo "# Starting script 04b_sex-and-z-class.sh"
date
echo "# Input QC dir:                $qc_dir"
echo "# Input BED file:              $bed"
echo "# Output dir:                  $outdir"
echo "# Sex chromosome:              $sex_chrom"
echo "# Exclusion-list policy:       $z_policy"
echo "# Min. call fraction per site: $min_call_frac"
echo "# Excess-het site threshold:   f > $max_het_frac"
echo "# Paralog when:                >${paralog_site_frac} of sites are excess-het sites"
echo "#                              AND >= $min_het_sites such sites"
echo "# Hemizygous het-ratio cutoff: <= $hemi_ratio"
echo "# Min autosomal depth:         $min_auto_dp"
echo "# Max within-group z-score:    $max_group_z"
echo

# Check inputs
[[ ! -d "$qc_dir" ]] && { echo "Error: QC dir not found: $qc_dir" >&2; exit 1; }
[[ ! -f "$bed" ]]    && { echo "Error: BED file not found: $bed" >&2; exit 1; }

# Load software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools

shopt -s nullglob
mkdir -p "$outdir"/logs

z_bed="$outdir"/z_loci.bed
depth_file="$outdir"/depth_per_sample.tsv
class_file="$outdir"/z_locus_classes.tsv
sex_file="$outdir"/sex_calls.tsv
exclude_file="$outdir"/loci_to_exclude.txt
het_long="$outdir"/logs/z_het_long.tsv
z_sites="$outdir"/logs/z_sites_long.tsv

awk -v c="$sex_chrom" '$1 == c' "$bed" > "$z_bed"
n_z=$(wc -l < "$z_bed")
[[ $n_z -eq 0 ]] && { echo "Error: no loci on '$sex_chrom' in $bed" >&2; exit 1; }
echo "# Loci on $sex_chrom: $n_z (of $(wc -l < "$bed") total)"


# ===============================================================================
echo -e "\n# Step 1: Reading per-locus depth and het counts..."
# ===============================================================================
count_files=("$qc_dir"/*_counts-per-locus.tsv)
n_samples=${#count_files[@]}
[[ $n_samples -eq 0 ]] && { echo "Error: no '*_counts-per-locus.tsv' files in $qc_dir" >&2; exit 1; }
echo "# Number of samples: $n_samples"

# Depth in long format: sample, locus ('__AUTOSOME__' for the autosomal aggregate),
# mean FORMAT/DP, and the number of sites it was averaged over.
# Het in long format: sample, locus, n_het, locus length -- Z loci only.
: > "$depth_file"
: > "$het_long"
for f in "${count_files[@]}"; do
    s=$(basename "$f" _counts-per-locus.tsv)
    awk -F'\t' -v OFS='\t' -v s="$s" -v c="$sex_chrom" -v B="$bed" \
        -v DP="$depth_file" -v HET="$het_long" '
        BEGIN {
            while ((getline l < B) > 0) { split(l, f, "\t"); chrom[f[4]] = f[1]; len[f[4]] = f[3] - f[2] }
        }
        FNR == 1 { next }
        {
            if (chrom[$1] == c) {
                if ($7 + 0 > 0) print s, $1, $6, $7 >> DP
                print s, $1, $2, len[$1] >> HET
            } else if ($7 + 0 > 0) {
                as += $6 * $7; an += $7
            }
        }
        END { if (an > 0) print s, "__AUTOSOME__", as / an, an >> DP }
    ' "$f"
done
echo "# Depth records: $(wc -l < "$depth_file")"
echo "# Het records on $sex_chrom: $(wc -l < "$het_long")"


# ===============================================================================
echo -e "\n# Step 2: Reading per-site genotypes on $sex_chrom..."
# ===============================================================================
# Long format: sample, locus, pos, code. Restricted to the sex chromosome, which is
# a small share of the loci, so this stays well under the size of the full tables.
z_loci_file="$outdir"/logs/z_loci_list.txt
cut -f4 "$z_bed" | sort > "$z_loci_file"

read_z_sites() {
    local f="$1" z_loci_file="$2"
    local s
    s=$(basename "$f" _sites.tsv.gz)
    zcat "$f" | awk -F'\t' -v OFS='\t' -v s="$s" -v L="$z_loci_file" '
        BEGIN { while ((getline x < L) > 0) keep[x] = 1 }
        $3 in keep { print s, $3, $2, $4 }'
}
export -f read_z_sites

site_files=("$qc_dir"/*_sites.tsv.gz)
[[ ${#site_files[@]} -eq 0 ]] && { echo "Error: no '*_sites.tsv.gz' files in $qc_dir" >&2; exit 1; }
if [[ ${#site_files[@]} -ne $n_samples ]]; then
    echo "Error: $n_samples count files but ${#site_files[@]} site files - the two sets must match" >&2
    exit 1
fi
parallel -j "${SLURM_CPUS_PER_TASK:-4}" read_z_sites {} "$z_loci_file" \
    ::: "${site_files[@]}" > "$z_sites"
echo "# Per-site records on $sex_chrom: $(wc -l < "$z_sites")"


# ===============================================================================
echo -e "\n# Step 3: Calling sex and classifying loci (2 iterations)..."
# ===============================================================================
# Split a sorted list of ratios at its largest internal gap (searching the middle
# 60% only, so a single outlier at either end cannot define the split).
split_at_gap() {
    sort -k2,2g | awk -v OFS='\t' '
        { s[NR] = $1; r[NR] = $2 }
        END {
            lo = int(NR * 0.2); if (lo < 1) lo = 1
            hi = int(NR * 0.8); if (hi >= NR) hi = NR - 1
            best = -1
            for (i = lo; i <= hi; i++) if (r[i+1] - r[i] > best) { best = r[i+1] - r[i]; bi = i }
            printf "#   split between %.4f and %.4f (gap %.4f): %d low, %d high\n", r[bi], r[bi+1], best, bi, NR - bi > "/dev/stderr"
            for (i = 1; i <= NR; i++) print s[i], (i <= bi ? "hetero" : "homo"), r[i]
        }'
}

# Depth ratio over a given set of loci ('-' = all Z loci)
depth_ratio() {
    local loci="$1"
    awk -F'\t' -v OFS='\t' -v L="$loci" '
        BEGIN { all = (L == "-"); if (!all) while ((getline x < L) > 0) keep[x] = 1 }
        $2 == "__AUTOSOME__" { auto[$1] = $3; next }
        (all || $2 in keep) { zs[$1] += $3 * $4; zn[$1] += $4 }
        END { for (s in zn) if (auto[s] > 0 && zn[s] > 0) print s, (zs[s] / zn[s]) / auto[s] }' "$depth_file"
}

# Classify Z loci given a sex assignment. Within each group, 'f' at a site is the
# proportion of called samples that are heterozygous; a locus carrying too many
# sites above --max_het_frac in a group has a collapsed second copy in that group.
classify_loci() {
    local sex="$1"
    awk -F'\t' -v OFS='\t' -v S="$sex" -v H="$het_long" \
        -v mcf="$min_call_frac" -v mhf="$max_het_frac" \
        -v psf="$paralog_site_frac" -v mns="$min_het_sites" -v hr="$hemi_ratio" '
    BEGIN {
        while ((getline l < S) > 0) { split(l, f, "\t"); grp[f[1]] = f[2]; gn[f[2]]++ }
        # Per-locus het density per group, and locus lengths
        while ((getline l < H) > 0) {
            split(l, f, "\t")
            len[f[2]] = f[4]
            hs[grp[f[1]], f[2]] += f[3]
            hn[grp[f[1]], f[2]]++
        }
    }
    # sample, locus, pos, code
    {
        g = grp[$1]
        if (g == "") next
        if ($4 == "het")       het[g, $2, $3]++
        else if ($4 == "miss") miss[g, $2, $3]++
        if ($4 != "ref") touched[g, $2, $3] = 1
    }
    END {
        for (k in touched) {
            split(k, p, SUBSEP)
            g = p[1]; l = p[2]
            called = gn[g] - miss[k]
            if (called < mcf * gn[g]) { unscored[g, l]++; continue }
            if (het[k] / called > mhf) excess[g, l]++
        }
        for (l in len) {
            split("hetero homo", grps, " ")
            for (i = 1; i <= 2; i++) {
                g = grps[i]
                ns[i] = len[l] - (unscored[g, l] + 0)
                ne[i] = excess[g, l] + 0
                sf[i] = (ns[i] > 0 ? ne[i] / ns[i] : 0)
                hd[i] = (hn[g, l] ? hs[g, l] / hn[g, l] / len[l] : 0)
                para[i] = (sf[i] > psf && ne[i] >= mns)
            }
            ratio = (hd[2] > 0 ? hd[1] / hd[2] : (hd[1] > 0 ? 999 : 1))
            if (para[1] && para[2])      cls = "paralog"
            else if (para[1])            cls = "gametolog"
            else if (ratio <= hr)        cls = "hemizygous"
            else                         cls = "PAR-like"
            printf "%s\t%d\t%d\t%d\t%.5f\t%d\t%d\t%.5f\t%.6f\t%.6f\t%.3f\t%s\n",
                   l, len[l], ns[1], ne[1], sf[1], ns[2], ne[2], sf[2], hd[1], hd[2], ratio, cls
        }
    }' "$z_sites" | sort -k5,5gr
}

# --- Iteration 1: bootstrap sex from all Z loci, then classify loci
echo "# Iteration 1 - sex from all $sex_chrom loci:"
depth_ratio - | split_at_gap > "$outdir"/logs/sex_iter1.tsv
classify_loci "$outdir"/logs/sex_iter1.tsv > "$outdir"/logs/z_classes_iter1.tsv
awk -F'\t' '$12 == "hemizygous" {print $1}' "$outdir"/logs/z_classes_iter1.tsv | sort > "$outdir"/logs/hemi_loci_iter1.txt
n_hemi1=$(wc -l < "$outdir"/logs/hemi_loci_iter1.txt)
echo "#   Hemizygous loci found: $n_hemi1"

# --- Iteration 2: recall sex using hemizygous loci only, then reclassify
if [[ $n_hemi1 -eq 0 ]]; then
    echo "# Warning: no hemizygous loci found; keeping the iteration-1 sex calls" >&2
    cp "$outdir"/logs/sex_iter1.tsv "$outdir"/logs/sex_iter2.tsv
else
    echo "# Iteration 2 - sex from the $n_hemi1 hemizygous loci:"
    depth_ratio "$outdir"/logs/hemi_loci_iter1.txt | split_at_gap > "$outdir"/logs/sex_iter2.tsv
fi
classify_loci "$outdir"/logs/sex_iter2.tsv > "$outdir"/logs/z_classes_iter2.tsv

n_switched=$(join -1 1 -2 1 <(sort -k1,1 "$outdir"/logs/sex_iter1.tsv) <(sort -k1,1 "$outdir"/logs/sex_iter2.tsv) |
    awk '$2 != $4' | wc -l)
echo "# Samples whose sex call changed between iterations: $n_switched"


# ===============================================================================
echo -e "\n# Step 4: Writing locus classes..."
# ===============================================================================
{ echo -e "locus\tlength_bp\tn_sites_hetero\tn_excess_hetero\tsite_frac_hetero\tn_sites_homo\tn_excess_homo\tsite_frac_homo\thet_density_hetero\thet_density_homo\thet_ratio\tclass"
  cat "$outdir"/logs/z_classes_iter2.tsv
} > "$class_file"

echo "# Locus classes on $sex_chrom:"
tail -n +2 "$class_file" | awk -F'\t' '{n[$12]++} END {for (c in n) printf "#   %-12s %d\n", c, n[c]}'
echo "# Loci with a collapsed second copy (gametolog = heterogametic group only, paralog = both):"
awk -F'\t' '$12 == "gametolog" || $12 == "paralog" {
    printf "#   %-8s %-10s len=%-6s excess_sites: hetero=%-4s homo=%-4s  het_ratio=%.2f\n", $1, $12, $2, $4, $7, $11
}' "$class_file"
echo "# Locus class file: $class_file"

awk -F'\t' '$12 == "gametolog"  {print $1}' "$class_file" | sort > "$outdir"/z_gametolog_loci.txt
awk -F'\t' '$12 == "paralog"    {print $1}' "$class_file" | sort > "$outdir"/z_paralog_loci.txt
awk -F'\t' '$12 == "hemizygous" {print $1}' "$class_file" | sort > "$outdir"/z_hemizygous_loci.txt
awk -F'\t' '$12 == "PAR-like"   {print $1}' "$class_file" | sort > "$outdir"/z_par_loci.txt
cut -f4 "$z_bed" | sort > "$outdir"/z_all_loci.txt


# ===============================================================================
echo -e "\n# Step 5: Writing sex calls..."
# ===============================================================================
# The heterogametic (low-depth-ratio) group is female in a ZW system.
#
# A call is flagged as not confident when the sample is an outlier relative to
# the group it was assigned to, rather than against absolute cutoffs: the sex
# groups are defined by a gap split, so the only thing that can be judged is how
# well a sample sits inside its own group. Three statistics are z-scored within
# the assigned group using the median and a MAD-based SD estimate (robust, so a
# handful of aberrant samples cannot inflate the spread they are tested against),
# and any one of them exceeding --max_group_z flags the sample:
#   depth_ratio        catches samples sitting between the two modes, which the
#                      gap split has to assign to one side regardless
#   het_density_hemi   catches a sample whose heterozygosity over the
#                      differentiated Z does not match its assigned sex
#   n_het_gametolog    the female-specific signal (paralog loci are heterozygous
#                      in both sexes and do not discriminate)
# A statistic whose group MAD is zero is skipped, since there is no spread to
# score against; female het_density_hemi is usually such a case.
echo -e "sample\tauto_depth\tz_depth_hemi\tdepth_ratio\thet_density_hemi\tn_het_gametolog\tsex\tconfident\tflags" > "$sex_file"
awk -F'\t' -v OFS='\t' \
    -v D="$depth_file" -v H="$het_long" \
    -v HEMI="$outdir"/z_hemizygous_loci.txt -v GAM="$outdir"/z_gametolog_loci.txt \
    -v mindp="$min_auto_dp" -v maxz="$max_group_z" '
    # Median of the values in an array (the array itself is left untouched)
    function med(src,   cp, n, i) {
        n = asort(src, cp)
        if (n == 0) return 0
        if (n % 2) return cp[(n + 1) / 2]
        return (cp[n / 2] + cp[n / 2 + 1]) / 2
    }
    # MAD scaled to a normal SD, so the z-scores are on the usual scale
    function mad_sd(src,   m, cp, n, i, d) {
        m = med(src)
        n = asort(src, cp)
        for (i = 1; i <= n; i++) d[i] = (cp[i] < m ? m - cp[i] : cp[i] - m)
        return 1.4826 * med(d)
    }
    # Zero when the group has no spread in this statistic, which skips the check
    function rz(x, m, s) {
        if (s <= 0) return 0
        return (x < m ? m - x : x - m) / s
    }
    BEGIN {
        while ((getline x < HEMI) > 0) hemi[x] = 1
        while ((getline x < GAM)  > 0) gam[x] = 1
        while ((getline l < D) > 0) {
            split(l, f, "\t")
            if (f[2] == "__AUTOSOME__") { auto[f[1]] = f[3]; continue }
            if (f[2] in hemi) { zs[f[1]] += f[3] * f[4]; zn[f[1]] += f[4] }
        }
        while ((getline l < H) > 0) {
            split(l, f, "\t")
            if (f[2] in hemi) { hs[f[1]] += f[3]; hl[f[1]] += f[4] }
            if (f[2] in gam)  { gs[f[1]] += f[3] }
        }
    }
    { sex[$1] = ($2 == "hetero" ? "female" : "male"); ratio[$1] = $3 }
    END {
        # Per-sample statistics, collected per sex group for the group summaries
        for (s in sex) {
            zd[s] = (zn[s] ? zs[s] / zn[s] : 0)
            hd[s] = (hl[s] ? hs[s] / hl[s] : 0)
            g = sex[s]
            n[g]++
            gr[g, n[g]] = ratio[s]; gh[g, n[g]] = hd[s]; gg[g, n[g]] = gs[s] + 0
        }
        for (g in n) {
            for (i = 1; i <= n[g]; i++) { vr[i] = gr[g, i]; vh[i] = gh[g, i]; vg[i] = gg[g, i] }
            mr[g] = med(vr); sr[g] = mad_sd(vr)
            mh[g] = med(vh); sh[g] = mad_sd(vh)
            mg[g] = med(vg); sg[g] = mad_sd(vg)
            printf "#   %-8s n=%-4d ratio: median=%.4f sd=%.4f | het_hemi: median=%.5f sd=%.5f | gametolog hets: median=%.0f sd=%.1f\n", \
                g, n[g], mr[g], sr[g], mh[g], sh[g], mg[g], sg[g] > "/dev/stderr"
            delete vr; delete vh; delete vg
        }
        for (s in sex) {
            g = sex[s]
            zr = rz(ratio[s], mr[g], sr[g])
            zh = rz(hd[s], mh[g], sh[g])
            zg = rz(gs[s] + 0, mg[g], sg[g])
            flags = ""
            if (auto[s] < mindp)  flags = flags (flags ? "," : "") sprintf("low_depth(%.1f)", auto[s])
            if (zr > maxz)        flags = flags (flags ? "," : "") sprintf("depth_ratio(z=%.1f)", zr)
            if (zh > maxz)        flags = flags (flags ? "," : "") sprintf("het_density(z=%.1f)", zh)
            if (zg > maxz)        flags = flags (flags ? "," : "") sprintf("gametolog_hets(z=%.1f)", zg)
            printf "%s\t%.1f\t%.1f\t%.4f\t%.6f\t%d\t%s\t%s\t%s\n", \
                s, auto[s], zd[s], ratio[s], hd[s], gs[s], g, (flags ? "no" : "yes"), (flags ? flags : "-")
        }
    }' "$outdir"/logs/sex_iter2.tsv | sort -k4,4g >> "$sex_file"

echo "# Sex calls:"
tail -n +2 "$sex_file" | awk -F'\t' '{n[$7]++; if ($8 == "no") u[$7]++} END {for (s in n) printf "#   %-8s %d (%d not confident)\n", s, n[s], u[s] + 0}'
echo "# Mean depth ratio by sex:"
tail -n +2 "$sex_file" | awk -F'\t' '{s[$7] += $4; n[$7]++} END {for (k in s) printf "#   %-8s %.3f\n", k, s[k] / n[k]}'
echo "# Samples with a non-confident call:"
tail -n +2 "$sex_file" | awk -F'\t' '$8 == "no" {printf "#   %-10s auto_dp=%-6.1f ratio=%-7.3f het_hemi=%-9.5f gametolog_hets=%-4s sex=%-7s flags=%s\n", $1, $2, $4, $5, $6, $7, $9}'
echo "# Sex call file: $sex_file"


# ===============================================================================
echo -e "\n# Step 6: Writing the locus exclusion list..."
# ===============================================================================
case "$z_policy" in
    all)        cp "$outdir"/z_all_loci.txt "$exclude_file" ;;
    paralogs)   sort -u "$outdir"/z_gametolog_loci.txt "$outdir"/z_paralog_loci.txt > "$exclude_file" ;;
    gametologs) cp "$outdir"/z_gametolog_loci.txt "$exclude_file" ;;
    none)       : > "$exclude_file" ;;
esac
echo "# Loci to exclude ($z_policy): $(grep -c . "$exclude_file" || true)"
echo "# Exclusion list: $exclude_file"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

# Report
echo -e "\n# Done with script 04b_sex-and-z-class.sh"
date
