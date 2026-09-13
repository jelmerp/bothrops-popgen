# Prepare QC summary tables for the filtered consensus sequences and VCFs.
#
# Reads the per-sample/per-locus output of scripts 02c, 03c, 03d and 04a-04h,
# and writes four tables to results/stats/summaries. The per-sample stats of the
# flagged contaminant samples (run/1_main.md, 4I) fill in their own rows of the
# per-sample table and are used nowhere else.
#
#   sample-by-locus-stats.tsv   Long format, one row per sample x locus x metric:
#                               n_het / n_alt / n_diff as a percentage of called
#                               sites, n_miss as a percentage of locus length.
#                               Read by 05b_qc-loci.qmd and 05c_genotype-counts.qmd.
#   per-locus-stats.tsv         One row per locus, spanning all loci that entered
#                               step 3 so that the ones dropped along the way are
#                               still visible, with a `stage` column saying where
#                               each was dropped.
#   per-sample-stats.tsv        One row per sample, spanning all samples in the
#                               metadata, likewise with a `stage` column.
#   filtering-summary.tsv       Loci and samples remaining after each filter step.
#
# Loci dropped by 04b (sex-linked) and 04c (excess heterozygosity) are kept in the
# tables and flagged rather than removed, so that the downstream docs can show why
# they went. Use `retained_final` to restrict to the loci in the final data set.

# SETUP ------------------------------------------------------------------------
library(tidyverse)
library(here)

# --- Output files
outdir <- here("results/stats/summaries")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
all_file <- file.path(outdir, "sample-by-locus-stats.tsv")
locus_file <- file.path(outdir, "per-locus-stats.tsv")
sample_file <- file.path(outdir, "per-sample-stats.tsv")
summary_file <- file.path(outdir, "filtering-summary.tsv")

# --- Input files
# Sample metadata (the TSV is the canonical copy -- see run/1_main.md)
meta_file <- here("metadata/metadata-final.tsv")
# Locus filtering stats (script 03c)
filt_dir <- here("results/stats/locus_filtering")
locus_len_file <- file.path(filt_dir, "per_locus_length.tsv")
locus_miss_file <- file.path(filt_dir, "per_locus_missingness.tsv")
sample_miss_file <- file.path(filt_dir, "per_sample_mean_missingness.tsv")
short_loci_file <- file.path(filt_dir, "removed_short_loci.txt")
removed_samples_file <- file.path(filt_dir, "removed_samples.txt")
kept_bed_file <- file.path(filt_dir, "retained_loci.bed")
# Per-locus genotype counts and depth, one file per sample (script 04a)
vcf_dir <- here("results/stats/vcf_qc")
# Sex calls and sex-linked locus classes (script 04b)
sex_dir <- here("results/stats/sex_z")
sex_calls_file <- file.path(sex_dir, "sex_calls.tsv")
z_class_file <- file.path(sex_dir, "z_locus_classes.tsv")
z_exclude_file <- file.path(sex_dir, "loci_to_exclude.txt")
# Excess-heterozygosity stats and the final locus set (script 04c)
het_dir <- here("results/stats/excess_het")
het_stats_file <- file.path(het_dir, "per_locus_het_stats.tsv")
excess_het_file <- file.path(het_dir, "excess_het_loci.txt")
ab_skew_file <- file.path(het_dir, "ab_skew_loci.txt")
final_bed_file <- file.path(het_dir, "retained_loci_final.bed")
# Per-locus ClipKIT trimming (script 03d)
clipkit_file <- file.path(filt_dir, "per_locus_clipkit.tsv")
# Per-locus variable / parsimony-informative site counts (script 04e)
inform_file <- here("results/stats/informativeness/per_locus_informativeness.tsv")
# Per-sample mapping and coverage metrics from Sarek (script 02c)
sarek_file <- here("results/stats/sarek/sarek_stats.tsv")
# Per-sample sequencing run and flowcell (script 02d)
seq_batch_file <- here("results/stats/seq_batch/seq_batch.tsv")
# Per-sample Ts/Tv and singleton counts from the merged VCF (script 04f)
vcf_stats_file <- here("results/stats/vcf_tstv/per_sample_vcf_stats.tsv")
# Het allele balance per sample (script 04h)
ab_dir <- here("results/stats/allele_balance")
# Per-sample stats for the flagged contaminant samples, kept out of every locus
# decision (run/1_main.md, 4I): 04a and 04h run on them alone, 04b and 04f on them
# together with the retained samples, and only the flagged rows are used
vcf_flagged_dir <- here("results/stats/vcf_qc_flagged")
ab_flagged_dir <- here("results/stats/allele_balance_flagged")
sex_flagged_file <- here("results/stats/sex_z_flagged/sex_calls.tsv")
vcf_stats_flagged_file <- here("results/stats/vcf_tstv_flagged/per_sample_vcf_stats.tsv")

# --- Helpers
# Locus IDs carry the reference-genome prefix in the 03c output but not in the
# 04a-04c output; the short form is used throughout here
strip_prefix <- function(x) sub("^Bothrops_I9999_", "", x)

# Read a set of per-sample files and take the sample ID from the file name
read_per_sample <- function(files, suffix, ...) {
  read_tsv(files, id = "src_file", ...) |>
    mutate(sample = sub(suffix, "", basename(src_file)), .before = 1) |>
    select(!src_file)
}

# 02c, 03d, 04e and 04f are later additions; the script still produces its tables
# without them rather than failing, so that a partial rerun is possible
read_optional <- function(path, ...) {
  if (!file.exists(path)) {
    warning("Optional input not found, columns will be NA: ", path)
    return(NULL)
  }
  read_tsv(path, ...)
}

join_optional <- function(df, tbl, by) {
  if (is.null(tbl)) df else left_join(df, tbl, by = by)
}

read_lines_or_empty <- function(path) {
  if (!file.exists(path)) {
    warning("File not found, treating as empty: ", path)
    return(character(0))
  }
  x <- readLines(path)
  x[nzchar(x)]
}

# READ INPUT FILES -------------------------------------------------------------
meta <- read_tsv(meta_file, show_col_types = FALSE) |>
  mutate(contam = as.logical(contam))

# --- Per-locus genotype counts and depth (04a), one file per sample
read_counts <- function(files) {
  read_per_sample(files, suffix = "_counts-per-locus\\.tsv$", col_types = "ciiiidi") |>
    rename(dp_mean = mean_dp) |>
    mutate(
      locus = strip_prefix(locus),
      # 04a counts every position of the locus interval, so these sum to its length
      n_sites = n_het + n_alt + n_ref + n_miss,
      n_called = n_het + n_alt + n_ref,
      call_frac = n_called / n_sites
    )
}

# The retained samples: only these feed the per-locus and sample-by-locus tables,
# and they define `retained`
count_files <- list.files(vcf_dir, pattern = "_counts-per-locus\\.tsv$", full.names = TRUE)
stopifnot(length(count_files) > 0)
counts <- read_counts(count_files)

# The flagged contaminant samples (4I), used for their own per-sample rows only
flagged_samples <- meta$sample[meta$contam %in% TRUE]
count_files_flagged <- list.files(vcf_flagged_dir, pattern = "_counts-per-locus\\.tsv$",
                                  full.names = TRUE)
counts_flagged <- NULL
if (length(count_files_flagged) > 0) {
  counts_flagged <- read_counts(count_files_flagged) |>
    filter(sample %in% flagged_samples)
} else {
  warning("No 04a tables for the flagged samples in ", vcf_flagged_dir,
          "; their per-sample columns will be NA")
}

# --- Locus lengths and pre-filter missingness (03c)
locus_lens <- read_tsv(locus_len_file, col_types = "ci") |>
  mutate(locus = strip_prefix(locus)) |>
  rename(length = length_bp)

locus_miss <- read_tsv(locus_miss_file, col_types = "cd") |>
  mutate(locus = strip_prefix(locus)) |>
  rename(pct_N_prefilt = mean_pct_N)

short_loci <- strip_prefix(read_lines_or_empty(short_loci_file))

kept_bed <- read_tsv(
  kept_bed_file,
  col_names = c("chrom", "start", "end", "locus"),
  col_types = "ciic"
) |>
  mutate(locus = strip_prefix(locus))

# --- Sample missingness and the samples dropped by 03c
sample_miss <- read_tsv(sample_miss_file, col_types = "cd") |>
  rename(pct_N = mean_pct_N)
removed_samples <- read_lines_or_empty(removed_samples_file)

# --- Per-sample mean depth (04a)
dp_files <- c(list.files(vcf_dir, pattern = "_depth\\.tsv$", full.names = TRUE),
              list.files(vcf_flagged_dir, pattern = "_depth\\.tsv$", full.names = TRUE))
sample_dp <- read_tsv(dp_files, col_names = c("sample", "mean_dp"), col_types = "cd")

# --- Sex calls and sex-linked locus classes (04b)
sex_calls <- read_tsv(sex_calls_file, col_types = "cddddiccc")
# The flagged samples' calls come from the separate 04b run in 4I, which scored
# them together with the retained samples; only their own rows are taken from it
sex_calls_flagged <- read_optional(sex_flagged_file, col_types = "cddddiccc")
if (!is.null(sex_calls_flagged)) {
  sex_calls <- bind_rows(
    sex_calls,
    filter(sex_calls_flagged, sample %in% flagged_samples, !sample %in% sex_calls$sample)
  )
}
sex_calls <- sex_calls |>
  select(sample, sex, z_depth_ratio = depth_ratio,
         z_het_density = het_density_hemi, n_het_gametolog,
         sex_confident = confident, sex_flags = flags) |>
  mutate(sex_confident = sex_confident == "yes")

z_classes <- read_tsv(z_class_file, col_types = cols(locus = "c", class = "c")) |>
  select(locus, z_class = class) |>
  mutate(locus = strip_prefix(locus))

# Which of the classed Z loci are actually dropped depends on --z_policy
z_exclude_loci <- strip_prefix(read_lines_or_empty(z_exclude_file))

# --- Excess-heterozygosity stats and the final locus set (04c)
het_stats <- read_tsv(het_stats_file, show_col_types = FALSE) |>
  # The allele-balance columns are absent from stats written before 04c scored it
  select(locus, n_sites_het_scored = n_sites, n_excess_sites,
         excess_site_frac = site_frac, max_site_het_frac = max_het_frac,
         total_het, het_density,
         any_of(c("n_ab_samples", "n_ab_skewed", "ab_skew_frac"))) |>
  mutate(locus = strip_prefix(locus))

excess_het_loci <- strip_prefix(read_lines_or_empty(excess_het_file))
ab_skew_loci <- strip_prefix(read_lines_or_empty(ab_skew_file))

final_loci <- read_tsv(
  final_bed_file,
  col_names = c("chrom", "start", "end", "locus"),
  col_types = "ciic"
) |>
  pull(locus) |>
  strip_prefix()

# --- ClipKIT trimming per locus (03d). The pre_n_* site classes are left in the
# source file: ClipKIT scored them on the untrimmed alignment over all 191
# samples, so they do not describe the final data set (04e does).
clipkit <- read_optional(clipkit_file, show_col_types = FALSE)
if (!is.null(clipkit)) {
  clipkit <- clipkit |>
    select(locus, n_sites_pre, n_trimmed, pct_trimmed, trim_failed) |>
    mutate(locus = strip_prefix(locus), trim_failed = as.logical(trim_failed))
}

# --- Variable and parsimony-informative sites per locus (04e). Only the loci in
# the final set are scored, so the rest stay NA.
inform <- read_optional(inform_file, show_col_types = FALSE)
if (!is.null(inform)) {
  inform <- inform |>
    select(locus, n_sites_scored, n_variable, n_pi, n_singleton,
           prop_variable, prop_pi) |>
    mutate(locus = strip_prefix(locus))
}

# --- Sarek mapping and coverage metrics per sample (02c)
sarek <- read_optional(sarek_file, show_col_types = FALSE)

# --- Sequencing run and flowcell per sample (02d). The two runs differ about
# two-fold in capture efficiency, so `seq_batch` is a confounder for every
# depth- and missingness-based metric here and belongs in any such comparison.
# Only `seq_batch` is kept: instrument and run are constant or redundant with it.
seq_batch <- read_optional(seq_batch_file, col_types = "ccccc")
if (!is.null(seq_batch)) seq_batch <- select(seq_batch, sample, seq_batch)

# --- Ts/Tv and singleton counts per sample (04f). `mean_dp` is dropped: bcftools
# reports it over the merged records, while 04a's is over the target intervals.
vcf_stats <- read_optional(vcf_stats_file, show_col_types = FALSE)
# As for the sex calls, the flagged samples' rows come from the 4I run. Its merged
# VCF also holds the retained samples, so singletons are counted against both.
vcf_stats_flagged <- read_optional(vcf_stats_flagged_file, show_col_types = FALSE)
if (!is.null(vcf_stats) && !is.null(vcf_stats_flagged)) {
  vcf_stats <- bind_rows(
    vcf_stats,
    filter(vcf_stats_flagged, sample %in% flagged_samples, !sample %in% vcf_stats$sample)
  )
}
if (!is.null(vcf_stats)) {
  vcf_stats <- vcf_stats |>
    select(sample, n_transitions, n_transversions, ts_tv, n_indels,
           n_singletons, n_missing_vcf = n_missing)
}

# --- Allele balance of het genotypes per sample (04h), retained and flagged, over
# the final loci and the hets 02b judged ('exempt' ones have no balance). The `ab_`
# prefix keeps these apart from the columns 05b derives from the same tables.
ab_files <- c(list.files(ab_dir, pattern = "_het-ab\\.tsv\\.gz$", full.names = TRUE),
              list.files(ab_flagged_dir, pattern = "_het-ab\\.tsv\\.gz$", full.names = TRUE))
ab_sample <- NULL
if (length(ab_files) > 0) {
  ab_sample <- read_per_sample(ab_files, suffix = "_het-ab\\.tsv\\.gz$",
                               na = c("NA", "."), col_types = "cicdddnc") |>
    filter(locus %in% final_loci, status != "exempt") |>
    summarize(
      ab_n_het = n(),
      ab_pct_blanked = 100 * mean(status == "blanked"),
      ab_median = median(ab),
      ab_pct_below_0.35 = 100 * mean(ab < 0.35),
      .by = sample
    )
} else {
  warning("No 04h allele-balance tables found; the ab_* columns will be NA")
}

# PER-LOCUS TABLE --------------------------------------------------------------
# Aggregate the per-sample counts over samples
locus_from_counts <- counts |>
  summarize(
    n_samples = n(),
    n_samples_called = sum(n_called > 0),
    mean_dp = mean(dp_mean, na.rm = TRUE),
    min_dp = ifelse(all(is.na(dp_mean)), NA_real_, min(dp_mean, na.rm = TRUE)),
    mean_call_frac = mean(call_frac),
    sd_call_frac = sd(call_frac),
    min_call_frac = min(call_frac),
    mean_n_het = mean(n_het),
    mean_n_alt = mean(n_alt),
    .by = locus
  )

# Anchored on the pre-filter locus list so that loci dropped by 03c are included
per_locus <- locus_miss |>
  full_join(locus_lens, by = "locus") |>
  left_join(kept_bed, by = "locus") |>
  left_join(locus_from_counts, by = "locus") |>
  left_join(z_classes, by = "locus") |>
  left_join(het_stats, by = "locus") |>
  join_optional(clipkit, by = "locus") |>
  join_optional(inform, by = "locus") |>
  mutate(
    retained_03c = locus %in% kept_bed$locus,
    retained_final = locus %in% final_loci,
    # Where each locus left the pipeline
    stage = case_when(
      retained_final ~ "retained",
      locus %in% excess_het_loci ~ "dropped: excess het",
      locus %in% ab_skew_loci ~ "dropped: allele balance",
      locus %in% z_exclude_loci ~ paste0("dropped: Z (", z_class, ")"),
      locus %in% short_loci ~ "dropped: too short",
      !retained_03c ~ "dropped: missingness",
      TRUE ~ "dropped: other"
    )
  ) |>
  relocate(locus, chrom, start, end, length, stage, retained_03c, retained_final) |>
  arrange(locus)

# PER-SAMPLE TABLE -------------------------------------------------------------
# Anchored on the metadata so that the samples dropped by 03c are included
per_sample <- meta |>
  left_join(sample_miss, by = "sample") |>
  left_join(sample_dp, by = "sample") |>
  left_join(sex_calls, by = "sample") |>
  left_join(
    # The only place the flagged samples' 04a tables are used
    bind_rows(counts, counts_flagged) |>
      filter(locus %in% final_loci) |>
      summarize(
        n_loci_called = sum(n_called > 0),
        mean_pct_N_final = mean(100 * (1 - call_frac)),
        het_pct_final = 100 * sum(n_het) / sum(n_called),
        .by = sample
      ),
    by = "sample"
  ) |>
  join_optional(sarek, by = "sample") |>
  join_optional(ab_sample, by = "sample") |>
  join_optional(seq_batch, by = "sample") |>
  join_optional(vcf_stats, by = "sample") |>
  mutate(
    retained = sample %in% counts$sample,
    stage = case_when(
      retained ~ "retained",
      contam ~ "dropped: flagged contaminant",
      sample %in% removed_samples ~ "dropped: missingness",
      TRUE ~ "dropped: no data"
    )
  ) |>
  relocate(sample, stage, retained) |>
  relocate(any_of("seq_batch"), .after = "retained") |>
  arrange(sample)

# `enrichment`, `pct_on_target` and `target_mean_dp` come from 02c, which gets
# them by intersecting the mosdepth windows with the probe BED, so they exist for
# every sample rather than only the retained ones.
#
# `target_mean_dp` and `mean_dp` are not the same quantity and are expected to
# differ: the first averages over every position of every target including the
# uncovered ones, the second averages FORMAT/DP over the sites where 04a found a
# called genotype. Their ratio says how unevenly coverage sits within a target.
if (all(c("mean_dp", "target_mean_dp") %in% names(per_sample))) {
  per_sample <- per_sample |>
    mutate(dp_called_vs_target = mean_dp / target_mean_dp, .after = "target_mean_dp")
}

# SAMPLE-BY-LOCUS TABLE --------------------------------------------------------
# Long format, as consumed by 05b and 05c. Loci dropped by 04b/04c are kept and
# flagged; filter on `retained_final` to get the final data set.
all <- counts |>
  left_join(
    per_locus |> select(locus, z_class, retained_final),
    by = "locus"
  ) |>
  rename(length = n_sites) |>
  mutate(n_diff = n_alt + (n_het / 2)) |>
  select(sample, locus, length, z_class, retained_final,
         n_called, n_het, n_alt, n_diff, n_miss) |>
  pivot_longer(
    cols = c(n_het, n_alt, n_diff, n_miss),
    names_to = "genotype",
    values_to = "count"
  ) |>
  mutate(geno_pct = ifelse(
    genotype != "n_miss", (count / n_called) * 100, (count / length) * 100
  )) |>
  select(!c(n_called, count))

# FILTERING SUMMARY ------------------------------------------------------------
# The step is kept to a short token and the prose moved to its own column, so
# that the table still lines up under `column -t -s$'\t'`
n_retained <- length(unique(counts$sample))
filtering_summary <- tribble(
  ~step,   ~filter,                                    ~n_loci,                          ~n_samples,
  "input", "consensus sequences from step 2",          nrow(per_locus),                  nrow(meta),
  "03c",   "flagged contaminant samples removed",      nrow(per_locus),                  nrow(meta) - sum(meta$contam, na.rm = TRUE),
  "03c",   "samples with excess missingness removed",  nrow(per_locus),                  n_retained,
  "03c",   "loci too short or with excess missingness", sum(per_locus$retained_03c),     n_retained,
  "04b",   "sex-linked loci removed",                  sum(per_locus$retained_03c & !per_locus$locus %in% z_exclude_loci), n_retained,
  "04c",   "paralog loci removed (excess het or allele balance)", length(final_loci),               n_retained
)

# CHECKS AND REPORT ------------------------------------------------------------
# The locus lengths from 03c and the interval lengths counted by 04a must agree
len_check <- counts |>
  distinct(locus, n_sites) |>
  inner_join(locus_lens, by = "locus") |>
  filter(n_sites != length)
if (nrow(len_check) > 0) {
  warning("Locus length mismatch between 03c and 04a for ", nrow(len_check), " loci")
}

message("# Samples: ", length(unique(counts$sample)), " retained of ", nrow(meta))
message("# Flagged contaminant samples with per-sample stats (4I): ",
        length(unique(counts_flagged$sample)), " of ", length(flagged_samples))
message("# Loci:    ", length(final_loci), " retained of ", nrow(per_locus))
if ("n_pi" %in% names(per_locus)) {
  message("# Parsimony-informative sites: ",
          sum(locus_final_pi <- per_locus$n_pi[per_locus$retained_final], na.rm = TRUE),
          " over ", sum(!is.na(locus_final_pi)), " loci; ",
          sum(locus_final_pi == 0, na.rm = TRUE), " loci with none")
}
if ("ts_tv" %in% names(per_sample)) {
  ts_tv <- per_sample$ts_tv[per_sample$retained]
  message("# Ts/Tv: mean ", round(mean(ts_tv, na.rm = TRUE), 3),
          ", range ", paste(round(range(ts_tv, na.rm = TRUE), 3), collapse = "-"))
}
print(per_locus |> count(stage))
print(per_sample |> count(stage))
print(filtering_summary)

# WRITE OUTPUT FILES -----------------------------------------------------------
write_tsv(all, all_file)
write_tsv(per_locus, locus_file)
write_tsv(per_sample, sample_file)
write_tsv(filtering_summary, summary_file)
message("# Wrote: ", all_file, ", ", locus_file, ", ", sample_file, ", ", summary_file)
