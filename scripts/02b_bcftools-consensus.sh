#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=bcftools-consensus
#SBATCH --output=slurm-%x-%j.out

# Strict Bash settings
set -euo pipefail

# ==============================================================================
#                          CONSTANTS AND DEFAULTS
# ==============================================================================
# Constants - generic
DESCRIPTION="Filter a VCF and then produce a per-region consensus FASTA sequence with bcftools.
Assumes that the VCF only contains a single sample. The VCF is filtered in eight steps:
0) Multi-bp monomorphic block records are split into one record per base.
   NOTE: FreeBayes, called with --report-monomorphic, occasionally merges a run of
   invariant positions with identical stats into a single record (REF several bp long,
   ALT '.'). Left intact, such a record can span past the end of the small per-region
   reference slice 'bcftools consensus' is given below and crash it with
   '[E::bcf_get_variant_type] Requested allele outside valid range'.
1) Low-quality genotypes are soft-filtered to missing (./.) based on a user-provided expression
2) Genotypes at variant sites with QUAL below a threshold are set to missing.
   NOTE: this is deliberately restricted to ALT-bearing records. The VCF is called with
   --report-monomorphic, and FreeBayes gives invariant sites QUAL ~ 0, so an unrestricted
   QUAL filter would wipe out every reference call.
3) SNPs within a user-defined distance of indels are set to missing (./.).
   NOTE: this is done as a two-pass soft filter. 'bcftools filter --SnpGap' on its own
   hard-filters, i.e. drops the records, after which 'bcftools consensus' falls back to
   the reference base -- turning a call we don't trust into a positive hom-ref assertion.
   '--set-GTs' does not apply to '--SnpGap' failures, so the records are tagged first
   (FILTER=SnpGap) and blanked in a second pass, which makes them 'N' in the consensus.
4) The VCF is normalized: multiallelics and MNVs are split, indels left-aligned
5) Indels are removed
6) Heterozygous genotypes with skewed allele balance are set to missing.
   NOTE: restricted to records where AO+RO is within --ab_tol reads of DP. After
   'bcftools norm -m -any', a true 1/2 genotype is split across two records whose
   AO+RO falls well short of DP, and its allele balance looks extreme, so applying
   this filter unrestricted would discard legitimate multiallelic hets. A small
   tolerance still lets through the many records that are off by a read or two
   (a read with a third allele, or one freebayes left out of both counts).
7) Positions left with more than one record after normalization are dropped
   (optional, see --drop_multiallelic), because bcftools consensus cannot render them.
The resulting filtered VCF is then used to generate the consensus sequence(s) with bcftools consensus.
Genotypes that were set to missing in the steps above will be set to N in the consensus sequence.
Heterozygous sites will be represented with IUPAC codes.
"
SCRIPT_VERSION="2026-09-12"
SCRIPT_AUTHOR="Jelmer Poelstra"
REPO_URL=https://github.com/mcic-osu/mcic-scripts
FUNCTION_SCRIPT_URL=https://raw.githubusercontent.com/mcic-osu/mcic-scripts/main/dev/bash_functions.sh
TOOL_BINARY=
TOOL_NAME=bcftools
TOOL_DOCS=https://samtools.github.io/bcftools/bcftools.html
VERSION_COMMAND="bcftools --version; samtools --version; bedtools --version; seqkit version"

# Defaults - generics
env_type=conda                                 # Use a 'conda' env or a Singularity 'container'
conda_path=/fs/ess/PAS0471/jelmer/conda/bcftools  # Must also contain samtools, bedtools, and seqkit
container_url=
container_dir="$HOME/containers"
container_path=

# Defaults - tool parameters
filter_expr='FORMAT/DP < 5 || FORMAT/GQ < 20'  # Filter expression
snp_gap=5                                      # Minimum distance of SNPs to indels (bp)
min_qual=30                                    # Minimum QUAL at ALT-bearing sites
min_ab=0.25                                    # Minimum allele balance for het genotypes
max_ab=0.75                                    # Maximum allele balance for het genotypes
ab_tol=2                                       # Max. reads by which AO+RO may differ from DP for the AB filter to apply
drop_multiallelic=true                         # Remove positions left with >1 record after normalization

# ==============================================================================
#                                   FUNCTIONS
# ==============================================================================
script_help() {
    echo -e "
                        $0
    v. $SCRIPT_VERSION by $SCRIPT_AUTHOR, $REPO_URL
            =================================================

DESCRIPTION:
$DESCRIPTION
    
USAGE / EXAMPLE COMMANDS:
  - Basic usage example:
      sbatch $0 -i results/my.vcf --sample_id my_sample --bed my.bed --ref data/ref.fa -o results/consensus
    
REQUIRED OPTIONS:
  -i/--vcf            <file>  Input VCF file
  --sample_id         <str>   Sample ID to use in output file names
  --ref               <file>  Input reference FASTA file
  --bed               <file>  Input BED file with regions to include in the consensus sequence
  -o/--outdir         <dir>   Output dir (will be created if needed)
    
OTHER KEY OPTIONS:
  --filter_expr       <str>   Expression to filter genotypes with, e.g.:        [default: $filter_expr]
                                'FORMAT/DP < 5 || FORMAT/GQ < 20'
                                Genotypes matching this expression are set to
                                missing (./.) in the output VCF, and thus to 'N' 
                                in the consensus sequence.
  --snp_gap           <int>   Minimum distance of SNPs to indels (bp)           [default: $snp_gap]
                                Genotypes at SNPs within this distance of an
                                indel are set to missing, and thus to 'N' in the
                                consensus sequence. The records are kept and
                                tagged with FILTER=SnpGap. Use 0 to disable.
  --min_qual          <num>   Minimum QUAL at ALT-bearing sites                 [default: $min_qual]
                                Genotypes at variant sites below this QUAL are
                                set to missing. Monomorphic records are exempt.
                                Use 0 to disable.
  --min_ab            <num>   Minimum allele balance for het genotypes          [default: $min_ab]
  --max_ab            <num>   Maximum allele balance for het genotypes          [default: $max_ab]
                                Het genotypes outside [min_ab, max_ab] are set
                                to missing, but only where AO+RO is within
                                --ab_tol reads of DP. Set min_ab to 0 and max_ab
                                to 1 to disable.
  --ab_tol            <int>   Max. difference between AO+RO and DP, in reads,   [default: $ab_tol]
                                for a het to be judged on allele balance.
                                0 restricts the filter to exact AO+RO == DP.
  --drop_multiallelic <bool>  Remove multiallelic positions                     [default: $drop_multiallelic]
                                After 'norm --multiallelics -any --atomize',
                                a multiallelic site or an MNV with a 1/2 genotype
                                leaves >1 record at the same position. These
                                cannot be rendered correctly by bcftools consensus:
                                the records are applied independently, so e.g. a
                                site that is truly hom-ALT gets an IUPAC het code.
                                With 'true', such positions are dropped and become
                                'N' in the consensus.
  --more_opts         <str>   Quoted string with one or more additional options
                                for $TOOL_NAME
    
UTILITY OPTIONS:
  --env_type          <str>   Whether to use a Singularity/Apptainer container  [default: $env_type]
                                ('container') or a Conda environment ('conda') 
  --container_url     <str>   URL to download a container from                  [default (if any): $container_url]
  --container_dir     <str>   Dir to download a container to                    [default: $container_dir]
  --container_path    <file>  Local container image file ('.sif') to use        [default (if any): $container_path]
  --conda_path        <dir>   Full path to a Conda environment to use           [default (if any): $conda_path]
  -h/--help                   Print this help message
  -v/--version                Print script and $TOOL_NAME versions
    
TOOL DOCUMENTATION:
  $TOOL_DOCS
"
}

# Function to source the script with Bash functions
source_function_script() {
    # NOTE: the argument is optional - some call sites pass none and rely on
    #       the IS_SLURM global instead
    local is_slurm=${1:-${IS_SLURM:-false}} candidate

    # Determine the location of this script, and based on that, the function script
    if [[ "$is_slurm" == true ]]; then
        script_path=$(scontrol show job "$SLURM_JOB_ID" | awk '/Command=/ {print $1}' | sed 's/Command=//')
        script_dir=$(dirname "$script_path")
        SCRIPT_NAME=$(basename "$script_path")
    else
        script_dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
        SCRIPT_NAME=$(basename "$0")
    fi
    function_script_name="$(basename "$FUNCTION_SCRIPT_URL")"

    # Look for a local copy first, in order of preference, and only download as a
    # last resort: the download writes into the working dir, which many jobs share
    for candidate in "$script_dir"/../dev/"$function_script_name" \
                     "$script_dir"/../core-scripts/dev/"$function_script_name" \
                     "$function_script_name"; do
        if [[ -s "$candidate" ]]; then
            source "$candidate"
            check_functions_loaded
            return 0
        fi
    done

    # Download to a temp file, then move into place, so that concurrent jobs
    # can never source a half-written file
    echo "Can't find script with Bash functions ($function_script_name), downloading from GitHub..."
    tmp_script=$(mktemp "$function_script_name".XXXXXX)
    if ! wget -q "$FUNCTION_SCRIPT_URL" -O "$tmp_script"; then
        rm -f "$tmp_script"
        echo "ERROR: Failed to download $FUNCTION_SCRIPT_URL" >&2
        exit 1
    fi
    mv -f "$tmp_script" "$function_script_name"
    source "$function_script_name"
    check_functions_loaded
}

# Make sure the function script really provided the functions we rely on
check_functions_loaded() {
    if ! declare -F log_time check_val die load_env >/dev/null; then
        echo "ERROR: Sourced $function_script_name but its functions are missing" >&2
        echo "       (an outdated or truncated copy may be in the way - try deleting it)" >&2
        exit 1
    fi
}

# Check if this is a SLURM job, then load the Bash functions
if [[ -z "${SLURM_JOB_ID:-}" ]]; then IS_SLURM=false; else IS_SLURM=true; fi
source_function_script "$IS_SLURM"

# ==============================================================================
#                          PARSE COMMAND-LINE ARGS
# ==============================================================================
# Initiate variables
version_only=false  # When true, just print tool & script version info and exit
vcf=
ref_fa=
outdir=
bed=
sample_id=
more_opts=
threads=

# Parse command-line options
all_opts="$*"
all_opts_q=$(printf '%q ' "$@")   # Shell-quoted, so it can be re-run exactly
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i | --vcf )        check_val "$1" "${2:-}"; shift; vcf=$1 ;;
        --ref )             check_val "$1" "${2:-}"; shift; ref_fa=$1 ;;
        --bed )             check_val "$1" "${2:-}"; shift; bed=$1 ;;
        --sample_id )       check_val "$1" "${2:-}"; shift; sample_id=$1 ;;
        --snp_gap )         check_val "$1" "${2:-}"; shift; snp_gap=$1 ;;
        --min_qual )        check_val "$1" "${2:-}"; shift; min_qual=$1 ;;
        --min_ab )          check_val "$1" "${2:-}"; shift; min_ab=$1 ;;
        --max_ab )          check_val "$1" "${2:-}"; shift; max_ab=$1 ;;
        --ab_tol )          check_val "$1" "${2:-}"; shift; ab_tol=$1 ;;
        --drop_multiallelic ) check_val "$1" "${2:-}"; shift; drop_multiallelic=$1 ;;
        --filter_expr )     check_val "$1" "${2:-}" lax; shift; filter_expr=$1 ;;
        -o | --outdir )     check_val "$1" "${2:-}"; shift; outdir=$1 ;;
        --more_opts )       check_val "$1" "${2:-}" lax; shift; more_opts=$1 ;;
        --env_type )        check_val "$1" "${2:-}"; shift; env_type=$1 ;;
        --conda_path )      check_val "$1" "${2:-}"; shift; conda_path=$1 ;;
        --container_dir )   check_val "$1" "${2:-}"; shift; container_dir=$1 ;;
        --container_url )   check_val "$1" "${2:-}"; shift; container_url=$1 ;;
        --container_path )  check_val "$1" "${2:-}"; shift; container_path=$1 ;;
        -h | --help )       script_help; exit 0 ;;
        -v | --version)     version_only=true ;;
        * )                 die "Invalid option $1" "$all_opts" ;;
    esac
    shift
done

# ==============================================================================
#                          INFRASTRUCTURE SETUP
# ==============================================================================
# Load software
load_env "$env_type" "$conda_path" "$container_dir" "$container_path" "$container_url"
[[ "$version_only" == true ]] && print_version "$VERSION_COMMAND" && exit 0

# Check options provided to the script
[[ -z "$vcf" ]] && die "No input file specified, do so with -i/--vcf" "$all_opts"
[[ -z "$ref_fa" ]] && die "No reference FASTA file specified, do so with --ref" "$all_opts"
[[ -z "$bed" ]] && die "No BED file specified, do so with --bed" "$all_opts"
[[ -z "$outdir" ]] && die "No output dir specified, do so with -o/--outdir" "$all_opts"
[[ -z "$sample_id" ]] && die "No sample ID specified, do so with --sample_id" "$all_opts"
[[ ! -f "$vcf" ]] && die "Input file $vcf does not exist"
[[ ! -f "$ref_fa" ]] && die "Reference FASTA file $ref_fa does not exist"
[[ ! -f "$bed" ]] && die "BED file $bed does not exist"

# The reference index must already exist: this script is run as many concurrent
# per-sample jobs, which would otherwise race to create the same '.fai' file
[[ ! -f "$ref_fa".fai ]] &&
    die "Reference FASTA index $ref_fa.fai does not exist - create it first with 'samtools faidx $ref_fa'"

indir=$(realpath -m "$(dirname "$vcf")")
[[ "$indir" == "$(realpath -m "$outdir")" ]] && die "Input and output directories must be different"

# Define outputs based on script parameters
# NOTE: LOG_DIR is made absolute so log paths keep resolving if the working dir changes
LOG_DIR=$(realpath -m "$outdir")/logs
mkdir -p "$LOG_DIR" "$outdir"/fa "$outdir"/vcf "$outdir"/stats
vcf_out="$outdir"/vcf/"$sample_id".vcf.gz
fasta="$outdir"/fa/"$sample_id".fasta
region_counts_file="$outdir"/stats/"${sample_id}"_snp-counts-per-region.tsv
regions_file="$outdir"/stats/"$sample_id"_regions.txt
# NOTE: per-sample name - many jobs share one outdir, so a fixed name would be clobbered
vcf_stats_file="$outdir"/stats/"${sample_id}"_vcf-stats.txt

# Record how this script was called (and, under Slurm, which job ran it)
log_provenance "$LOG_DIR"

# ==============================================================================
#                         REPORT PARSED OPTIONS
# ==============================================================================
log_time "Starting script $SCRIPT_NAME, version $SCRIPT_VERSION"
echo "=========================================================================="
echo "All options passed to this script:        $all_opts"
echo "Working directory:                        $PWD"
echo
echo "Input VCF file:                           $vcf"
echo "Input reference FASTA file:               $ref_fa"
echo "Input BED file:                           $bed"
echo
echo "Filter expr. for low-quality genotypes:   $filter_expr"
echo "Minimum distance of SNPs to indels:       $snp_gap bp"
echo "Minimum QUAL at ALT-bearing sites:        $min_qual"
echo "Allele-balance range for het genotypes:   $min_ab - $max_ab"
echo "Allele-balance tolerance, |AO+RO-DP|:     $ab_tol"
echo "Drop multiallelic positions:              $drop_multiallelic"
echo "Output dir:                               $outdir"
echo "Output filtered VCF file:                 $vcf_out"
echo "Output consensus FASTA file:              $fasta"
[[ -n $more_opts ]] && echo "Additional options for $TOOL_NAME:        $more_opts"
log_time "Listing the input file(s):"
ls -lh "$vcf" "$ref_fa" "$bed"
echo "Number of regions in the BED file:        $(grep -c . "$bed" || true)"
set_threads "$IS_SLURM"
[[ "$IS_SLURM" == true ]] && slurm_resources

# ==============================================================================
#                               RUN
# ==============================================================================
# Getting basic stats on the input VCF file
log_time "Getting basic stats on the input VCF file..."
$TOOL_BINARY bcftools stats "$vcf" > "$vcf_stats_file"

# Build the QUAL filter expression.
# IMPORTANT: restricted to ALT-bearing records. The VCF is called with
# --report-monomorphic, and FreeBayes assigns invariant sites QUAL ~ 0
# (QUAL = P(site is polymorphic)). An unrestricted 'QUAL < x' filter therefore
# matches essentially every reference call and would blank out the consensus.
if [[ $(awk -v q="$min_qual" 'BEGIN{print (q>0)}') == 1 ]]; then
    qual_expr='ALT!="." && QUAL < '"$min_qual"
else
    qual_expr=
fi

# Build the SnpGap filter options.
# IMPORTANT: two passes are needed. '--SnpGap' hard-filters by default, which drops the
# records and makes 'bcftools consensus' emit the reference base for them -- a false
# hom-ref call rather than missing data, and one that is biased towards the reference
# because it is precisely the ALT-bearing calls that get erased. '--set-GTs' has no
# effect on '--SnpGap' failures, so the first pass only tags them (FILTER=SnpGap) and
# the second pass blanks the tagged genotypes.
# NOTE: '--SnpGap' only flags ALT-bearing SNP records; monomorphic records within the
# same window keep their reference base.
snpgap_opts=
snpgap_expr=
if [[ "$snp_gap" -gt 0 ]]; then
    snpgap_opts="--SnpGap $snp_gap -s SnpGap -m +"
    snpgap_expr='FILTER~"SnpGap"'
fi

# Build the allele-balance filter expression.
# IMPORTANT: restricted to records where AO+RO is within ab_tol reads of DP. After
# 'bcftools norm --multiallelics -any' a true 1/2 genotype becomes two records,
# each showing almost no support for the reference, so an unrestricted AB filter
# would throw away legitimate multiallelic het calls. Those records fall short of
# DP by the depth of the other allele, far more than a read or two, while exact
# equality also exempts the many records that are off by just one or two reads.
# The difference is bounded on both sides as sums, without subtracting: bcftools
# 1.23 evaluates a chain like 'DP - AO - RO' wrongly, so 'DP - AO - RO <= 2'
# matched no record at all and the filter silently blanked nothing.
# NOTE: AO has Number=A, so it must be indexed explicitly (FMT/AO[0:0]) for
# arithmetic in bcftools filter expressions to work.
if [[ $(awk -v lo="$min_ab" -v hi="$max_ab" 'BEGIN{print (lo>0 || hi<1)}') == 1 ]]; then
    ab='FMT/AO[0:0]/(FMT/AO[0:0]+FMT/RO)'
    ab_expr='GT="het" && (FMT/AO[0:0]+FMT/RO) <= (FMT/DP+'"$ab_tol"') && FMT/DP <= (FMT/AO[0:0]+FMT/RO+'"$ab_tol"') && ('"$ab"' < '"$min_ab"' || '"$ab"' > '"$max_ab"')'
else
    ab_expr=
fi

# Filter the VCF
log_time "Filtering the VCF file..."
# Step 0: Split multi-bp monomorphic block records (REF several bp, ALT '.') into
# one record per base. FreeBayes occasionally merges a run of invariant positions
# with identical stats into a single such record; left intact, one can span past
# the end of the small per-region reference slice 'bcftools consensus' is given
# further down and crash it with 'Requested allele outside valid range'. Splitting
# does not change what the consensus looks like: GT/DP/etc. are carried over
# unchanged to every base of the run, exactly as they applied to the block as a whole.
n_block=$($TOOL_BINARY bcftools view -H "$vcf" | awk -F'\t' '$5=="." && length($4)>1' | wc -l)
echo "Multi-bp monomorphic block records found: $n_block"
# Step 1: Soft-filter bad genotypes to missing (./.)
# Step 2: Set genotypes at low-QUAL variant sites to missing (ALT-bearing records only)
# Step 3: Tag SNPs within x bp of an indel with FILTER=SnpGap, then set those
#         genotypes to missing in a second pass, so they become 'N' in the consensus
# Step 4: Normalize the VCF to split multiallelics and MNVs, and left-align indels, for the next step
# Step 5: Remove indels, keeping only SNPs
# Step 6: Set allele-balance outlier het genotypes to missing (biallelic records only)
#? https://samtools.github.io/bcftools/bcftools.html#filter
$TOOL_BINARY bcftools view "$vcf" |
    awk -F'\t' -v OFS='\t' '
        /^#/ { print; next }
        $5 == "." && length($4) > 1 {
            ref = $4; pos = $2
            for (i = 0; i < length(ref); i++) { $2 = pos + i; $4 = substr(ref, i+1, 1); print }
            next
        }
        { print }' |
    $TOOL_BINARY bcftools view -Ou |
    runstats $TOOL_BINARY bcftools filter -e "$filter_expr" --set-GTs . -O u |
    runstats $TOOL_BINARY bcftools filter ${qual_expr:+-e "$qual_expr" --set-GTs .} -O u |
    runstats $TOOL_BINARY bcftools filter ${snpgap_opts} -O u |
    runstats $TOOL_BINARY bcftools filter ${snpgap_expr:+-e "$snpgap_expr" --set-GTs .} -O u |
    runstats $TOOL_BINARY bcftools norm --multiallelics -any --atomize -O u |
    runstats $TOOL_BINARY bcftools view -V indels -O u |
    runstats $TOOL_BINARY bcftools filter ${ab_expr:+-e "$ab_expr" --set-GTs .} \
        --threads "$threads" -O z -o "$vcf_out"

# Step 7: Drop positions that still carry more than one record.
# 'norm --multiallelics -any --atomize' splits a multiallelic site, and also
# decomposes an MNV carrying a 1/2 genotype, into several records at the same
# position. bcftools consensus applies those records independently and therefore
# renders them wrongly -- e.g. an MNV 'CTATG -> TTATA,CTATA' with GT 1/2 yields
# two 'G>A' records (1/0 and 0/1) at the last position, where the true genotype
# is A/A with zero reference reads, yet the consensus gets the IUPAC code 'R'.
# Genuine 1/2 SNPs are mis-rendered for the same reason. Dropping these positions
# makes them 'N' in the consensus (via 'bcftools consensus --absent N').
if [[ "$drop_multiallelic" == true ]]; then
    log_time "Dropping multiallelic positions..."
    dup_pos_file="$outdir"/stats/"$sample_id"_multiallelic_pos.txt
    $TOOL_BINARY bcftools query -f '%CHROM\t%POS\n' "$vcf_out" | uniq -d > "$dup_pos_file"
    n_dup=$(grep -c . "$dup_pos_file" || true)
    echo "Multiallelic positions found: $n_dup"

    if [[ "$n_dup" -gt 0 ]]; then
        vcf_tmp="$vcf_out".tmp.vcf.gz
        runstats $TOOL_BINARY bcftools filter -T ^"$dup_pos_file" \
            --threads "$threads" -O z -o "$vcf_tmp" "$vcf_out"
        mv "$vcf_tmp" "$vcf_out"
    fi
fi

echo "# Resulting file:"
ls -lh "$vcf_out"

log_time "Counting sites in the raw and filtered VCF file..."
raw_count=$($TOOL_BINARY bcftools view -H "$vcf" | wc -l)
filt_count=$($TOOL_BINARY bcftools view -H "$vcf_out" | wc -l)
echo "Raw / filtered VCF: $raw_count / $filt_count sites (removed: $((raw_count - filt_count)))"
if [[ -n "$snpgap_expr" ]]; then
    n_snpgap=$($TOOL_BINARY bcftools view -H -i "$snpgap_expr" "$vcf_out" | wc -l)
    echo "Genotypes blanked by the SnpGap filter: $n_snpgap"
fi

# Break the SNP records of the filtered VCF down by genotype.
# NOTE: a record whose genotype was set to missing by the filters above still
# carries an ALT allele, so it is counted by 'bcftools view -v snps' even though
# it becomes an 'N' in the consensus rather than a variant. Count genotypes, not
# records, or the number of variable sites is overstated by a wide margin.
log_time "Summarizing genotypes in the filtered VCF file..."
$TOOL_BINARY bcftools view -v snps "$vcf_out" |
    $TOOL_BINARY bcftools query -f '[%GT]\n' |
    awk '
        { n_rec++ }
        /\./                        { n_miss++; next }
        $1 == "0/0" || $1 == "0|0"  { n_ref++; next }
        { n_alt++; split($1, gt, /[\/|]/); if (gt[1] != gt[2]) n_het++ }
        END {
            printf "SNP records in the filtered VCF:   %d\n", n_rec
            printf "Variable sites (ALT genotypes):    %d\n", n_alt
            printf "  of which het (IUPAC-coded):      %d\n", n_het
            printf "Hom-ref genotypes:                 %d\n", n_ref
            printf "Missing genotypes (become N):      %d (%.1f%%)\n",
                   n_miss, n_rec ? 100 * n_miss / n_rec : 0
        }'

# Count number of SNPs per region in the BED file (Add a column with the sample ID with awk)
# NOTE: restricted to ALT-bearing genotypes, for the reason given just above
log_time "Counting variable sites per region..."
runstats $TOOL_BINARY bedtools intersect -c -a "$bed" \
    -b <($TOOL_BINARY bcftools view -v snps -i 'GT="alt"' "$vcf_out") |
    awk -v OFS='\t' -v sample="$sample_id" '{print $0, sample}' > "$region_counts_file"
echo "# Resulting file:"
ls -lh "$region_counts_file"
echo "# Regions with 0 variable sites: $(awk -F'\t' '$(NF-1) == 0' "$region_counts_file" | wc -l)"

# Set the VCF's sample name to --sample_id. Sarek takes it from the sample sheet,
# so without this the name inside the VCF can differ from the file name used from
# here on, and the merged VCF built by 04d would carry the Sarek names.
log_time "Setting the VCF sample name to '$sample_id'..."
samplename_file=$(mktemp)
echo "$sample_id" > "$samplename_file"
vcf_reheadered="$vcf_out".reheader.vcf.gz
$TOOL_BINARY bcftools reheader -s "$samplename_file" -o "$vcf_reheadered" "$vcf_out"
mv -f "$vcf_reheadered" "$vcf_out"
rm -f "$samplename_file"

# Index the new VCF ('-f' so that re-runs don't fail on an existing index)
log_time "Indexing the filtered VCF file..."
runstats $TOOL_BINARY bcftools index -f --threads "$threads" "$vcf_out"

# Create a regions-file for samtools (1-based coordinates)
log_time "Formatting BED coordinates..."
awk '{print $1":"$2+1"-"$3}' "$bed" > "$regions_file"
echo "# Resulting file:"
ls -lh "$regions_file"

# Generate consensus sequence
# 1. samtools faidx grabs the reference sequence just for the focal loci
# 2. bcftools consensus reads the >chr:start-end header, maps it to the VCF, and applies the variants/Ns
log_time "Generating consensus sequences per locus..."
runstats $TOOL_BINARY samtools faidx -r "$regions_file" "$ref_fa" |
    runstats $TOOL_BINARY bcftools consensus \
        --haplotype I --missing N --absent N $more_opts "$vcf_out" \
        > "$outdir"/fa/"${sample_id}"_init.fa

#! Note: when bcftools reports 'Applied X variants', this includes the counts of
#! missing genotypes (set to 'N'), so it will be higher than the variant count above.
#? --missing 'N' => set missing genotypes to 'N' in the consensus sequence
#? --haplotype I => for heterozygous sites, use IUPAC codes

# Rename the FASTA entries from region to locus ID
# NOTE: many per-sample jobs share one output dir. Every job derives identical
# content from the same BED file, and a same-filesystem 'mv' is an atomic rename,
# so a concurrent reader always sees either no file or a complete one.
log_time "Renaming the FASTA entries from region to locus ID..."
lookup_file="$outdir"/stats/lookup.tsv
if [[ ! -f "$lookup_file" ]]; then
    lookup_tmp="$lookup_file".tmp.$$
    awk -v OFS="\t" '{print $1":"$2+1"-"$3, $4}' "$bed" > "$lookup_tmp"
    mv -f "$lookup_tmp" "$lookup_file"
fi

$TOOL_BINARY seqkit replace \
    -p '^(\S+)(.*)$' \
    -r '{kv} $1' \
    -k "$lookup_file" \
    "$outdir"/fa/"${sample_id}"_init.fa > "$fasta"

rm "$outdir"/fa/"${sample_id}"_init.fa

echo -e "\n# Resulting FASTA file:"
ls -lh "$fasta"
echo -e "\n# First couple of headers of the FASTA file:"
grep -m5 "^>" "$fasta" || log_time "WARNING: no FASTA headers found in $fasta"

# Summarize missing data in the consensus, which is the key QC signal here:
# a locus with poor coverage comes out as (mostly) N rather than as a short sequence
log_time "Summarizing the consensus FASTA file..."
$TOOL_BINARY seqkit fx2tab -nl -B N "$fasta" |
    awk -F'\t' -v n_bed="$(grep -c . "$bed" || true)" '
        { n_loci++; bp += $(NF-1); pct_n_sum += $NF; if ($NF >= 100) n_all_n++ }
        END {
            printf "Loci in the consensus FASTA:       %d (BED regions: %d)\n", n_loci, n_bed
            printf "Total consensus length:            %d bp\n", bp
            printf "Mean %% N per locus:                %.1f%%\n",
                   n_loci ? pct_n_sum / n_loci : 0
            printf "Loci that are entirely N:          %d\n", n_all_n + 0
            if (n_loci != n_bed)
                print "WARNING: number of FASTA records differs from the number of BED regions"
        }'

# ==============================================================================
#                               WRAP-UP
# ==============================================================================
log_time "Listing files in the output dir:"
ls -lhd "$(realpath "$outdir")"/*

# Copy Slurm output file to the logs directory
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    slurm_file="slurm-${SLURM_JOB_NAME}-${SLURM_JOB_ID}.out"
    if [[ -f "$slurm_file" ]]; then
        cp "$slurm_file" "$LOG_DIR"/ 2>/dev/null || true
    fi
fi

final_reporting "$LOG_DIR"
