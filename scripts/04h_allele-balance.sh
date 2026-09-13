#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=60
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=allele-balance
#SBATCH --output=slurm-allele-balance-%j.out
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: 04h_allele-balance.sh --vcf <raw_vcf> --sample_id <id> --bed <bed_file> --outdir <outdir>
                             [--min_ab <num>] [--max_ab <num>] [--ab_tol <int>]
                             [--min_qual <num>] [--snp_gap <int>] [--filter_expr <str>]

Allele balance of heterozygous genotypes, as seen by the allele-balance filter in
02b, for a single sample.

02b sets het genotypes whose allele balance, AO / (AO + RO), falls outside
[min_ab, max_ab] to missing. The filtered VCF therefore only shows the hets that
passed, and it cannot say how many were removed: a blanked genotype looks like any
other missing one. This script replays the 02b filter chain on the raw FreeBayes
VCF up to the allele-balance step, and records every het that reaches that step
together with what the filter did to it.

The replayed steps, which must stay in sync with 02b:
  0) split multi-bp monomorphic block records
  1) low-quality genotypes to missing (--filter_expr)
  2) genotypes at low-QUAL ALT-bearing sites to missing (--min_qual)
  3) SNPs within --snp_gap bp of an indel to missing
  4) split multiallelics and MNVs (norm -m -any --atomize)
  5) remove indels
  7) drop positions left with more than one record
Step 6, the allele-balance filter itself, is evaluated rather than applied.

Only sites inside --bed are reported, and the BED restriction is applied after the
SnpGap step, since an indel just outside a locus can still blank a SNP inside it.

Output files, prefixed with the sample ID:
  _het-ab.tsv.gz      One row per het genotype reaching the allele-balance step:
                      chrom, pos, locus, dp, ro, ao, ab, status
                      ao is the first ALT's count, which is what 02b uses
                      status is 'kept', 'blanked' (outside [min_ab, max_ab]), or
                      'exempt' (AO + RO more than --ab_tol reads from DP, which
                      02b does not evaluate; ab is NA for these)
  _ab-summary.tsv     sample, n_het, n_kept, n_blanked, n_exempt

Required arguments:
  --vcf             Raw FreeBayes .vcf.gz file (the 02b input)
  --sample_id       Sample ID, used for the output file names
  --bed             BED file with loci (4th column: locus name)
  -o, --outdir      Output directory

Optional arguments (defaults match run/1_main.md):
  --min_ab          Minimum allele balance for het genotypes    [default: 0.25]
  --max_ab          Maximum allele balance for het genotypes    [default: 0.75]
  --ab_tol          Max. |AO + RO - DP| for a het to be judged  [default: 2]
  --min_qual        Minimum QUAL at ALT-bearing sites           [default: 30]
  --snp_gap         Minimum distance of SNPs to indels, in bp   [default: 5]
  --filter_expr     Genotype filter expression   [default: 'FORMAT/DP < 5 || FORMAT/GQ < 20']
  -h, --help        Show this help and exit
EOF
}

# Parse options
vcf=
sample_id=
bed=
outdir=
min_ab=0.25
max_ab=0.75
ab_tol=2
min_qual=30
snp_gap=5
filter_expr='FORMAT/DP < 5 || FORMAT/GQ < 20'
while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)          vcf="$2"; shift 2 ;;
        --sample_id)    sample_id="$2"; shift 2 ;;
        --bed)          bed="$2"; shift 2 ;;
        -o|--outdir)    outdir="$2"; shift 2 ;;
        --min_ab)       min_ab="$2"; shift 2 ;;
        --max_ab)       max_ab="$2"; shift 2 ;;
        --ab_tol)       ab_tol="$2"; shift 2 ;;
        --min_qual)     min_qual="$2"; shift 2 ;;
        --snp_gap)      snp_gap="$2"; shift 2 ;;
        --filter_expr)  filter_expr="$2"; shift 2 ;;
        -h|--help)      usage; exit 0 ;;
        *)              echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

# Check inputs
[[ -z "$vcf" || -z "$sample_id" || -z "$bed" || -z "$outdir" ]] && { usage >&2; exit 1; }
[[ ! -f "$vcf" ]] && { echo "Error: VCF file not found: $vcf" >&2; exit 1; }
[[ ! -f "$bed" ]] && { echo "Error: BED file not found: $bed" >&2; exit 1; }

# Report
echo "# Starting script 04h_allele-balance.sh"
date
echo "# Input VCF:          $vcf"
echo "# Sample ID:          $sample_id"
echo "# BED file:           $bed"
echo "# Output directory:   $outdir"
echo "# Allele balance:     [$min_ab, $max_ab]"
echo "# AB tolerance:       |AO + RO - DP| <= $ab_tol"
echo "# Min. QUAL:          $min_qual"
echo "# SnpGap:             $snp_gap"
echo "# Genotype filter:    $filter_expr"

# Load software environment
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/bcftools # Contains bcftools and bedtools

mkdir -p "$outdir"/logs

het_file="$outdir/${sample_id}_het-ab.tsv.gz"
summary_file="$outdir/${sample_id}_ab-summary.tsv"

# Replay 02b steps 0-5, then restrict to the loci
echo -e "\n# Replaying the 02b filter chain and tabulating het allele balance..."
bcftools view "$vcf" 2>/dev/null |
    awk -F'\t' -v OFS='\t' '
        /^#/ { print; next }
        $5 == "." && length($4) > 1 {
            ref = $4; pos = $2
            for (i = 0; i < length(ref); i++) { $2 = pos + i; $4 = substr(ref, i+1, 1); print }
            next
        }
        { print }' |
    bcftools view -Ou 2>/dev/null |
    bcftools filter -e "$filter_expr" --set-GTs . -Ou 2>/dev/null |
    bcftools filter -e 'ALT!="." && QUAL < '"$min_qual" --set-GTs . -Ou 2>/dev/null |
    bcftools filter --SnpGap "$snp_gap" -s SnpGap -m + -Ou 2>/dev/null |
    bcftools filter -e 'FILTER~"SnpGap"' --set-GTs . -Ou 2>/dev/null |
    bcftools norm --multiallelics -any --atomize -Ou 2>/dev/null |
    bcftools view -V indels -Ou 2>/dev/null |
    bcftools view -T "$bed" -Ou 2>/dev/null |
    bcftools query -f '%CHROM\t%POS[\t%GT\t%DP\t%RO\t%AO]\n' 2>/dev/null |
    awk -F'\t' -v OFS='\t' -v min_ab="$min_ab" -v max_ab="$max_ab" -v tol="$ab_tol" \
        -v sample="$sample_id" -v summary="$summary_file" '
    FNR == NR {
        for (b = int($2 / 5000); b <= int($3 / 5000); b++)
            idx[$1 ":" b] = idx[$1 ":" b] " " $2 "," $3 "," $4
        next
    }
    # Step 7 of 02b drops every position that carries more than one record. Records
    # arrive sorted, so duplicates are contiguous: each record is held back until
    # the next one shows whether its position repeats.
    function emit(line,   f, g, locus, bucket, cand, c, i, dp, ro, ao, ao_vals, ab, status) {
        split(line, f, "\t")
        if (f[3] ~ /\./) return
        split(f[3], g, /[\/|]/)
        if (g[1] == g[2]) return

        bucket = f[1] ":" int(f[2] / 5000)
        if (!(bucket in idx)) return
        locus = ""
        split(idx[bucket], cand, " ")
        for (i in cand) {
            if (cand[i] == "") continue
            split(cand[i], c, ",")
            if (f[2] > c[1] + 0 && f[2] <= c[2] + 0) { locus = c[3]; break }
        }
        if (locus == "") return

        # AO has one value per ALT. A record can keep a second ALT that is not a SNP
        # (the spanning-deletion allele '*' survives -V indels), and 02b judges
        # balance on the first value (FMT/AO[0:0]), so this does too.
        dp = f[4]; ro = f[5]; split(f[6], ao_vals, ","); ao = ao_vals[1]
        # The same condition 02b applies before judging balance: AO + RO within
        # tol reads of DP. A zero denominator cannot be judged either way.
        if (ao ~ /^[0-9]+$/ && ro ~ /^[0-9]+$/ && dp ~ /^[0-9]+$/ &&
            ao + ro - dp <= tol && dp - ao - ro <= tol && ao + ro > 0) {
            ab = ao / (ao + ro)
            status = (ab < min_ab || ab > max_ab) ? "blanked" : "kept"
            ab = sprintf("%.4f", ab)
        } else {
            ab = "NA"; status = "exempt"
        }
        n[status]++; n_het++
        print f[1], f[2], locus, dp, ro, ao, ab, status
    }
    {
        key = $1 ":" $2
        if (key == prev_key) { dup = 1; next }
        if (prev_line != "" && !dup) emit(prev_line)
        prev_line = $0; prev_key = key; dup = 0
    }
    END {
        if (prev_line != "" && !dup) emit(prev_line)
        print "sample", "n_het", "n_kept", "n_blanked", "n_exempt" > summary
        print sample, n_het + 0, n["kept"] + 0, n["blanked"] + 0, n["exempt"] + 0 > summary
    }' "$bed" - |
    { echo -e "chrom\tpos\tlocus\tdp\tro\tao\tab\tstatus"; cat; } |
    gzip > "$het_file"

# Report
echo -e "\n# Summary:"
column -t "$summary_file"
echo -e "\n# Listing output files:"
ls -lh "$het_file" "$summary_file"

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$outdir"/logs/ 2>/dev/null || true
    fi
fi

echo -e "\n# Done with script 04h_allele-balance.sh"
date
