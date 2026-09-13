
Lisle sent a tree from mtDNA sequences that Felipe made, saying:
> This data is only concatenated coding sequences from mtDNA – no tRNAs or non-coding regions (snakes have 2 d loops).

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
