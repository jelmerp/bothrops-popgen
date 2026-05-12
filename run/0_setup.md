
## Create project dir structure (2026-02-11)

```bash
cd /fs/ess/PAS1533/users/jelmer
mkdir -p data/fastq notes run scripts software
git clone https://github.com/mcic-osu/mcic-scripts.git
```

## Copy the sequence data (2026-02-11)

```bash
ddir=/fs/ess/PAS1533/users/andrewmason/projects/Bothrops_jararaca_anchored
cp -v "$ddir"/00_raw_data.tar.bz2 data/
tar -xjvf data/00_raw_data.tar.bz2 -C data/fastq
```

## Download the reference genome (2026-02-11)

- Reference genome <https://academic.oup.com/gbe/article/18/1/evaf243/8379231>
- Downloaded from FigShare on 2025-02-16 (<https://figshare.com/projects/Bothrops_insularis_genome/237995>)
  Referred to as version 1

- Create a FASTA file with chromosomes and unplaced scaffolds:

```bash
primary=data/ref/Binsularis_primary_chromosomes.fasta
unplaced=data/ref/Binsu_unplaced.fasta
cat "$primary" "$unplaced" > data/ref/binsu_all.fna
```

## Remove sample with much too few reads (2026-04-01)

```bash
rm -rv data/fastq/P0021_FG_I1765_1
```
