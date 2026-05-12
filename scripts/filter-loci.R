# Load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# Define command-line options
option_list <- list(
  make_option(
    c("--blast-file"),
    type = "character",
    default = "results/blast_probes/blast_out_final.tsv",
    help = "Path to BLAST output file [default: %default]",
    metavar = "FILE"
  ),
  make_option(
    c("--outdir"),
    type = "character",
    default = "results/blast_probes",
    help = "Output directory for BED files [default: %default]",
    metavar = "DIR"
  ),
  make_option(
    c("--buffer"),
    type = "numeric",
    default = 500,
    help = "Buffer size in bp on each side of locus [default: %default]",
    metavar = "INT"
  )
)

# Parse arguments
parser <- OptionParser(option_list = option_list)
opts <- parse_args(parser)

# Extract options
blast_file <- opts$`blast-file`
output_dir <- opts$outdir
BUFFER_SIZE <- opts$buffer

# Validate buffer size
if (is.na(BUFFER_SIZE) || BUFFER_SIZE < 0) {
  stop("Buffer size must be a non-negative number")
}

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define output files with buffer size in filename
bedfile_nobuffer <- file.path(output_dir, "probes_no-buffer.bed")
bedfile_buffer <- file.path(output_dir, paste0("probes_", BUFFER_SIZE, "bp-buffer.bed"))

# Read input files
blast_init <- read_tsv(blast_file, show_col_types = FALSE) |>
  select(q = qseqid, pident, qcov = qcovhsp, chrom = sacc, start = sstart, end = send)

# Get multi-hits
multi_hits <- blast_init |>
  group_by(q) |>
  add_count(name = "n_hits") |>
  filter(n_hits > 1)

# Filter the loci
blast_filt <- blast_init |>
  filter(any(pident >= 95), .by = q) |>
  filter(any(qcov >= 90), .by = q) |>
  filter(!q %in% unique(multi_hits$q)) |>
  slice_head(n = 1, by = q)

# Change to BED format
bed <- blast_filt |>
  mutate(
    start_final = ifelse(start < end, start, end) - 1,
    end_final = ifelse(start < end, end, start)
  ) |>
  select(chrom, start = start_final, end = end_final, probe = q) |>
  arrange(chrom, start, end)
#any(duplicated(bed$probe)) # FALSE

# BED with 500-bp buffer
bed_buffer <- bed |>
  mutate(
    start = ifelse(start - BUFFER_SIZE >= 0, start - BUFFER_SIZE, 0),
    end = end + BUFFER_SIZE
  )

# Write the output files
write_tsv(bed, bedfile_nobuffer, col_names = FALSE)
write_tsv(bed_buffer, bedfile_buffer, col_names = FALSE)

# Report
message("\n# Input BLAST file: ", blast_file)
message("# Output directory: ", output_dir)
message("# Buffer size: ", BUFFER_SIZE, " bp")
message("# Number of loci before filtering: ", length(unique(blast_init$q)))
message("# Number of loci after filtering:  ", nrow(bed))
system(paste("ls -lh", bedfile_buffer, bedfile_nobuffer))
