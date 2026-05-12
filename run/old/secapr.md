
## SECAPR Clean Reads

```bash
sbatch -A PAS1533 -c 8 -t 400 --wrap="
    echo 'Running SECAPR clean_reads'
    secapr clean_reads --input data/fastq/P0021_FG_I12277 --output TEST --sample_annotation_file data/test2/sample_annotation.tsv
    echo 'Done'
"
```

## SECAPR Reference-assembly

<https://raw.githubusercontent.com/AntonelliLab/seqcap_processor/refs/heads/master/secapr/reference_assembly.py>

```bash
ref=data/ref/Binsularis_primary_chromosomes.fasta
outdir=results/secapr_ref-asm
sbatch scripts/secapr_ref-asm.sh data/test2 "$ref" "$outdir"

secapr reference_assembly --reads data/test2 --output results/secapr_ref-asm --reference data/ref/Binsularis_primary_chromosomes.fasta --reference_type user-ref-lib --cores 1
```

## Paralogous regions

```bash
secapr paralogs_to_ref -h
#file:///Users/poelstra.1/Library/CloudStorage/Dropbox/mcic/assist/2026-01_lisle/decapr/align_paralogs.html
```
