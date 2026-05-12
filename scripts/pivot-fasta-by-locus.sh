#!/usr/bin/env bash
#SBATCH --account=PAS0471
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=pivot-fasta
#SBATCH --output=slurm-pivot-fasta-%j.out
set -euo pipefail

usage() {
  cat <<'EOF'
Pivot per-sample multi-locus FASTA files into per-locus FASTA files.

Usage:
  pivot-fasta-by-locus.sh <input_dir> <output_dir> [extension]

Arguments:
  input_dir   Directory containing per-sample FASTA files.
  output_dir  Directory where per-locus FASTA files will be written.
  extension   File extension to match (default: .fa)

Output:
  - One FASTA file per locus in output_dir
  - Entry headers are sample names (input filename without extension)
  - locus_name_map.tsv mapping original locus names to output filenames
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

input_dir="$1"
output_dir="$2"
ext="${3:-.fasta}"

echo "Starting script pivot-fasta-by-locus.sh"
date
echo "Input directory:              $input_dir"
echo "Output directory:             $output_dir"
echo "File extension to match:      $ext"
echo "----------------------------------------"
echo

if [[ ! -d "$input_dir" ]]; then
  echo "ERROR: Input directory not found: $input_dir" >&2
  exit 1
fi

mkdir -p "$output_dir"/logs

# Use a temp staging area to avoid partial output on failure.
tmp_out="$(mktemp -d "${output_dir%/}/.pivot_tmp.XXXXXX")"
cleanup() {
  rm -rf "$tmp_out"
}
trap cleanup EXIT

shopt -s nullglob
fasta_files=("$input_dir"/*"$ext")
shopt -u nullglob

if [[ ${#fasta_files[@]} -eq 0 ]]; then
  echo "ERROR: No input files matching *$ext in $input_dir" >&2
  exit 1
fi

map_file="$tmp_out/locus_name_map.tsv"
printf "locus\tfile\n" > "$map_file"

for fasta in "${fasta_files[@]}"; do
  sample="$(basename "$fasta")"
  sample="${sample%$ext}"

  awk \
    -v sample="$sample" \
    -v out_dir="$tmp_out" \
    -v map_file="$map_file" '
function sanitize(raw,   s) {
  s = raw
  gsub(/[^A-Za-z0-9._-]/, "_", s)
  if (s == "") s = "locus"
  return s
}

function flush_record(   safe_name, out_fa, n, candidate) {
  if (locus == "") return

  safe_name = sanitize(locus)
  candidate = safe_name
  n = 1
  while ((candidate in file_to_locus) && file_to_locus[candidate] != locus) {
    n++
    candidate = safe_name "__" n
  }

  file_to_locus[candidate] = locus
  if (!(candidate in seen_map)) {
    print locus "\t" candidate ".fa" >> map_file
    close(map_file)
    seen_map[candidate] = 1
  }

  out_fa = out_dir "/" candidate ".fa"
  print ">" sample >> out_fa
  print seq >> out_fa
  close(out_fa)
}

/^>/ {
  flush_record()
  locus = substr($0, 2)
  seq = ""
  next
}

{
  gsub(/[[:space:]]/, "", $0)
  seq = seq $0
}

END {
  flush_record()
}
' "$fasta"
done

# Move staged outputs into final destination.
shopt -s nullglob
for produced in "$tmp_out"/*.fa "$tmp_out"/locus_name_map.tsv; do
  mv "$produced" "$output_dir/"
done
shopt -u nullglob

# Report
echo "Done"
date
echo "Listing the files in the output directory: $output_dir"
ls -lh "$output_dir"
