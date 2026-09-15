## 2026-09-13 -- 05c merged into 05b; four more samples excluded

**I1495, I1507, I1517 and I1624 excluded as likely contaminants** (not flagged in
the metadata). All four have low allele balance and many hets blanked for their
batch; I1495, I1507 (females, raised Z depth ratio) and I1624 (male) also have a
sex call that is not confident. Added to 03c's a priori removal via
`metadata/contam_samples_qc.txt` and treated as flagged contaminants in 4I and 05a.
Leaves 121 samples (66 on 605, 55 on 705). Pipeline needs a rerun from 3C.
I1671 fits the pattern more weakly and is kept.

**05c genotype-count plots moved into 05b** (new last section, *Genotype counts
by group*); 05c removed. Problems found in 05c along the way: its x axis,
`contam`, is FALSE for all 125 samples in the sample-by-locus table (flagged
samples are not in it), so every plot had one box per panel; the y label said
"sites per 100 bp" but `geno_pct` is a percentage of called sites; only 2
*jararaca* are from islands (vs 103 mainland), and only *jararaca* has more than
2 males. The merged version plots all three metrics per figure, by species (split
by land type) and by *jararaca* state (north to south), with points colored by
batch and group sizes on the axis. The land-type-within-*jararaca* and by-sex
plots were dropped; sex is covered by *Heterozygosity per sample*.

**SC heterozygosity follows city, not batch.** Within *jararaca* in SC,
heterozygosity splits into a low group (0.25-0.32%: Florianópolis, Fraiburgo,
Capinzal; 6 of 7 samples from 605) and a high group (0.44-0.59%: mostly northern
SC cities, also mostly 605). Florianópolis has one high sample, I12294 (705). The
two island *jararaca* (I1767, I1768) are both Florianópolis and low. So the split
inside one batch supports locality over batch for the state differences. The
`locality` column is NA for all SC samples; `city` is filled in. Added to 05b as a
by-city plot. Also in 05b: species labels italic (ggtext), title now "Locus and
sample QC and genotype counts", and *Heterozygosity per sample* moved to directly
before *Genotype counts by group*.

## 2026-09-12 -- Sequencing batch effect; two 05b outliers explained; locus length vs. depth

**Rerun from 02b (ab_tol = 2, AB locus test in 04c): first checks.** 04c drops
L1306, L1813, L473, L663 (excess het) and L1813, L1826, L473 (allele balance), so
L1826 is the only new drop: 1697 final loci. 04b now flags 7 sex calls as not
confident (was 5): females I1671, I20937, I1507, I1795, I1495 on depth ratio,
I1624 (depth ratio + gametolog hets), I12279 (het density). Of the 14 samples 03c
dropped for missingness, 10 are from the 705 batch. 05b proofread and shortened;
out-of-date numbers are now inline R where possible, and a flagged-samples section was
added.

**Flagged contaminant samples (4I) vs. retained.** 05a had four stale input
paths (`results/stats/vcf`, `vcf_stats`, `vcf_flagged`, `vcf_stats_flagged`; now
`vcf_qc`, `vcf_tstv`, `vcf_qc_flagged`, `vcf_tstv_flagged`) and would have failed
on the first; fixed and run. Singletons and sex-call confidence depend on the
cohort (I12277: 165 singletons among the 125 retained, 138 among all 177; I12279
and I1671 are not confident only in the 125-sample 04b run), so 05b compares both
groups on the 4I runs. Medians, flagged (52) vs retained (125): allele balance
0.495 vs 0.494, hets below 0.35 9.3% vs 6.8%, hets blanked 1.6% vs 0.9%, het
0.43% vs 0.46%, depth 38x vs 47x; sex call not confident in 23% vs 4%. The
medians hide a clear subset: 13 flagged samples have median AB < 0.45 with
21-54% of hets below 0.35, lowest I1613, I1681, I1622, I1683, I1625, I1609
(0.33-0.36). A second group of not-confident flagged samples (I1554, I1557,
I1639, I18405; *germanoi* I1774, I1775) has median AB 0.5 but very low
heterozygosity (0.07-0.28%) and mostly low depth. 17 of the 52 flagged samples
have > 10% missing data over the final loci, so the missingness filter would
likely have removed them regardless of the flag. 05b rendered against the new
summaries; the hand-typed sample claims in the allele-balance, I20937, and
'Samples to watch' text were checked and updated (I20937 now 3,583 singletons,
9.9%; I1495/I1507 depth ratios 0.57/0.53 vs female median 0.43; I1671 added as a
weaker low-AB + sex-flag case).

**What 04a's `mean_dp` averages over.** 04a runs on `retained_loci.bed`, whose
coordinates 03c trims with the ClipKIT log, so the depth intervals are not probe
+ the full 600 bp of buffer. In the final loci, a median of 353 bp of buffer
remains on short loci (553-bp intervals) and 310 bp on long ones; buffer is ~34%
of all counted bases. 108 loci also lost some probe bases to trimming. Rough
estimate from the 04g region depths: ~51x over retained intervals against ~64x on
probes only. Per sample, buffer depth / probe depth (full 300-bp buffer, 04g
logic) is 0.36-0.48 (median 0.40; r = -0.16 with insert size), so buffer share
barely affects between-sample comparisons, but does affect between-locus ones.
Corrected the depth-section intro in 05b accordingly.

**Per-sample stats for the flagged contaminant samples (workflow change, not yet
run).** The 52 samples flagged a priori were removed in 03c, so the only stats for
them came from 02c. Checking the flags against the data (allele balance, sex-call
anomalies, het rate, singletons) needs the same per-sample stats as the retained
samples. They still have to be removed before any cohort-level decision, since
almost every step after 03c pools samples (locus missingness, Z classes, excess
het and AB skew, informativeness, depth profile, merged VCF). So 03c keeps
removing them, and a new step 4I computes their stats on the side, feeding nothing
back:

- 04a and 04h run on them alone, into separate directories (04b/04c/04e read every
  table in their input directory).
- 04b (sex calls) and 04d + 04f (Ts/Tv, singletons) depend on the other samples;
  both rerun on retained + flagged together, and 05a takes only the flagged rows.
  04f got a `--bed` option to restrict the stats to the final loci, since the
  flagged VCFs cover all loci.
- 05a adds these to the flagged samples' rows only, plus new per-sample columns
  for everyone: `het_pct_final` and `ab_*` (allele balance over the final loci).
- Caveats: the 04b/04f reference cohort for the flagged rows includes the flagged
  samples; ClipKIT (3B) still trims from all 191 samples, as before.
- Also fixed the 04h sample loop in 4C: the `*"$smp"` glob could match a longer ID
  (I1674 -> I11674).

**Paralog filter and AB tolerance changed (not yet rerun).** Two follow-ups to the
allele-balance findings below:

- *04c now also drops loci on allele balance.* For each sample with >= 5 evaluated
  hets at a locus, the locus is skewed there when their median AB is < 0.4; a
  locus is dropped when it is skewed in > 50% of >= 20 such samples. A locus
  failing either this or the excess-het test goes. Reads the 04h tables, so 04h
  now runs before 04c (run.md step 4C; later steps relettered). Tested on the
  current data: flags L1826 and L1813 (already excess-het), so the only new drop
  is L1826 (1698 -> 1697 loci). L10 (0.486) sits just under 0.5; a cutoff of
  0.45 would drop it too. 05a labels these 'dropped: allele balance'.
- *02b's AB filter now judges hets with |AO + RO - DP| <= 2* (`--ab_tol`, was
  exact equality). Simulated on the current 04h tables: unevaluated hets fall
  from 17.8% to 9.2%, and the blanked share barely moves (2.09% -> 2.12% of
  evaluated), since the newly judged hets are mostly balanced (median 0.48).
- *Bug caught while testing:* bcftools 1.23 evaluates a chained subtraction like
  'DP - AO - RO' wrongly, so the first version of the tolerance expression matched
  no record and blanked nothing. Rewritten as sums without subtraction. Also, a
  record can keep a spanning-deletion ALT ('*') after `-V indels`, giving AO like
  '6,.'; 04h now reads the first value, as 02b does.
- *Validated on I12277:* at tolerance 2, 02b and 04h agree position by position
  (6,138 kept, 466 exempt, 91 blanked, all missing in 02b's output); at tolerance
  0, 02b reproduces the existing filtered VCF with no genotype differences.

The pipeline has to be rerun from 02b for either change to reach the results.

**Allele balance (new script 04h).** 02b's allele-balance filter was applied but
never reported, because a blanked het looks like any other missing genotype. 04h
replays the 02b filter chain on the raw VCFs and records every het that reaches
the AB step. It reproduces the filtered VCFs exactly (checked position by position
for I12277). Figures below are for retained loci and samples.

- *Filter effect is small.* 2.1% of evaluated hets are blanked; median sample 0.9%
  (IQR 0.7-1.2%), unrelated to depth (r = 0.05).
- *The exemption is far broader than intended.* 17.8% of hets are never evaluated
  because AO + RO != DP. The case the exemption was written for, RO = 0 after
  splitting a 1/2 genotype, is only 27% of these; 48% are off by 1-2 reads. A
  tolerance (e.g. |AO + RO - DP| <= 2) would bring most of them under the filter.
  Not changed.
- *Pooled AB is reference-biased, not paralog-contaminated.* Single mode at 0.5;
  8.7% of hets below 0.35 vs 5.8% above 0.65, strongest at 11-40x. No mode near 0.25.
- *L1826 looks like a paralog that 04c missed.* Per-sample measure: share of samples
  (with >= 5 hets) whose median AB at the locus is < 0.4. L1826 is at 66% of 114
  samples, level with the dropped paralogs L1839 and L1813 (68%); the 99th
  percentile of retained loci is 30%. It also has 2x the long-class median depth
  (129x vs 64x) and 7 excess-het sites, under the 04c threshold. Balance of
  0.25-0.35 at the same sites across dozens of samples. Decision pending: drop
  it, or add an AB criterion to 04c. L10 is next (49% of 72) but weaker: normal
  depth, no excess-het sites.
- *Possible contamination in I1495, I1507, I1624 (and I1517).* Low median AB for
  their batch (0.37-0.44), and they are the samples that fail sex calling: 3 of the
  6 low-AB samples are 'sex not confident' vs 2 of 119 others (Fisher p = 6e-4).
  I1495 and I1507 are females with raised Z depth ratios (0.61, 0.58 vs ~0.48);
  I1624 is a male with gametolog hets. The three lowest-AB females (I1495, I1507,
  I1520) have the three highest female Z het densities. None is in the metadata
  contamination flags. Both AB flags added to 'Samples to watch'.
- *Unexplained: I1768, I1670, I1684.* ~20% of hets blanked, but balance centred on
  0.5, normal library metrics, among the fewest transitions and very few singletons.
  Not a one-sided contamination pattern.

**Locus length vs. depth: the earlier explanation in 05b was wrong.** It claimed
longer loci "carry more probes" and so capture better per base. Each locus is
exactly one probe (01d keeps one BLAST hit per query and drops multi-hit probes),
and the pre-trimming lengths are bimodal at ~800 and ~2,500 bp, not ~800/~1,950.
The bimodality is in the bait set itself: 1,474 baits of ~200 bp and ~500 of
~1,900-2,000 bp.

New script `04g_depth-profile.sh` bins every called position by its distance from
the nearest probe edge (125 samples). Capture efficiency ramps up over ~300 bp
inward from a bait edge and plateaus at ~73x; both bait classes trace the *same*
ramp. A ~200-bp bait is entirely inside it -- its most interior base is 100 bp
from an edge -- so it tops out at 47x and averages 45x on probe, against 70x for
the long class. Buffer positions sit at ~25x and are 97.6% called, so they do
enter `mean_dp`. Buffer share (75% vs 23% of the interval) plus the ramp together
give the whole ~2x gap in per-locus depth. Extra bait length past the ramp buys
nothing: within the long class r(probe_len, depth) = -0.07.

The trimming half of the story survives and is separable: within the short class
alone, final length still tracks depth at r = 0.82 (long class: 0.18).

**Consequence for the depth figures.** Every depth number in 05b -- per-sample
`mean_dp`, per-locus `mean_dp`, and 02c's `target_mean_dp`/`pct_on_target`/
`enrichment`, which also use the buffered BED -- averages over intervals that are
mostly buffer for the short class. They understate on-target depth and should not
be compared to a per-probe figure; between-sample comparisons are unaffected.
Caveat added to the per-sample depth section.


**Batch scan over the rest of 05b.** Tested every per-sample metric in the
report against batch (Wilcoxon, BH-adjusted), then asked of each significant
one whether `state` absorbs it. The split is clean and diagnostic:

- *Technical, batch is the cause.* `pct_on_target`/`enrichment` (+84%),
  `pct_duplication` (+86%), `pct_mq0` (-22%), `error_rate` (-9.5%),
  `insert_size_avg` (419 vs 390 bp), depth, read count, missingness. For
  `mean_dp` and `mean_pct_N_final`, adding batch on top of state is still
  highly significant (p = 5e-8, p = 7e-4), i.e. not geography.
- *Biological, geography is the cause.* Per-sample heterozygosity differs by
  batch within *jararaca* (0.492% vs 0.445%, p = 1.1e-4) and so do singletons
  (354 vs 254, p = 0.003) -- but `state` absorbs both completely (batch adds
  nothing: p = 0.28 and p = 0.20). Depth does not explain het either
  (coefficient n.s., and r = 0.12 within *jararaca*).
- *No batch effect at all.* Ts/Tv, `pct_mapped`, `z_depth_ratio` within sex
  (0.480/0.483 female, 0.851/0.832 male -- the ratio self-normalises, so sex
  calling is batch-immune), sex ratio, `sex_confident`.

**Batch is confounded with locality** (Cramer's V = 0.583, p = 1.1e-6): RS
(13/14), RJ (10/11) and PR (4/4) are essentially 705-only, SC (16/19) and ES
(11/13) mostly 605. Samples were evidently batched by collection region. So the
het/singleton batch signal is geographic structure, not an artifact -- but the
confounding also means the two cannot be fully separated, and 05c's *state*
comparison is the one to treat with care (species p = 0.007, land type p = 0.14,
sex p = 0.25).

**The I20937 capture argument in 05b is wrong as written.** It says 9.7% on
target and enrichment 68.6 are "near the bottom of the panel" -- true only
because I20937 is in the weak-capture 705 batch. Within 705 it ranks 14/57 on
both, just under the batch medians (10.2%, 72.0), well inside the batch range
(8.7-11.9%). And I1773 (*germanoi*, the n=1 control with "normal capture") is in
the 605 batch, so that comparison is cross-batch and inflated. The conclusion
still holds on the other evidence -- the singleton count itself, lowest-in-panel
het (0.19%), normal Ts/Tv (2.07), and an error rate that *is* genuinely highest
(rank 57/57 within its own batch) -- but the capture evidence is much weaker
than stated. I1767 (depth outlier) is unaffected: 122x against a 605 median of
57x is still a 2x outlier within its own batch.

**The quantile flags in "Samples to watch" mostly select batch**, since they cut
the tail of a bimodal distribution: low depth 7/7 in 705, high duplication 7/7
in 605, many singletons 6/7 and low Ts/Tv 6/7 in 605. The composite list
(>= 2 flags) happens to come out balanced 5/5, but the individual flags should be
computed within batch. Switched 05b to within-batch quantiles: flags now split
4/3 per batch, and the list goes from 10 samples to 11 (I12294 and I1723 in,
I1495 out). `sex not confident` is left panel-wide, since 04b already z-scores
within the assigned sex group.

**The excess-het filter is not batch-driven** -- all four dropped loci are
elevated in both batches (L1306 1.03/1.04, L1813 1.49/1.56, L663 1.40/1.38,
L473 1.32/0.90, against panel means of 0.38/0.42). Per-locus ranking is well
preserved across batches too (het r = 0.92, missingness r = 0.86), so locus
filtering decisions are robust; only the global offset differs (67% of all 1698
loci are higher-het in 605, p < 2e-16).

**Two sequencing batches, ~2x apart in capture efficiency.** The read-count vs.
depth plot in 05b splits into two lines with slopes 17.0 and 8.6 depth per
million reads. The cause is the sequencing run, read from the FASTQ read names:
run 605 (flowcell HCTGLBCXX, 68 retained samples) against run 705 (H2HLLBCXY,
57 retained).

The mechanism is the on-target rate and nothing else: depth-per-read correlates
with `pct_on_target` at r = 0.998 overall and 0.988 within each batch, and the
batches do not overlap on that axis (605: 13.3-25.4%, mean 18.8; 705: 8.7-11.9%,
mean 10.2). Duplication is *not* the driver, though it is confounded with batch
(6.2% vs 3.3%), which is why colouring that plot by duplicate rate looked
explanatory -- the gradient was re-encoding the batch. Within batch 605 the
duplication/depth-per-read correlation is only -0.10 (705: -0.56).

Partly self-cancelling: 705 was sequenced ~30% deeper (3.63M vs 2.81M reads) and
has the lower duplicate rate, so it sits down and to the right rather than
simply below. Net effect on the final data is modest -- mean missingness 4.13%
(705) vs 3.21% (605), and loci called per sample is 1696 of 1698 either way.

Batch is near-even for *jararaca* (55/50) and for mainland samples (53/51), but
not for the island endemics: *insularis* (6) and *sazimai* (5) are 605-only,
*alcatraz* is 6 of 7 on 705. Check any coverage-derived between-species
comparison against this.

One crossover sample: I1619 is on 605 but sits with the 705 cloud -- highest
duplicate rate in the dataset (14.4%) plus the lowest on-target rate of its own
batch (13.3%).

Added as script 02d (step 2D), joined into the per-sample table by 05a as
`seq_batch`, and 05b now colours by it and has a "Sequencing batch" subsection.

**I20937, the singleton outlier.** 3,587 singletons vs a median of 307 -- 10.1%
of all singleton sites in the dataset, from one of 125 samples. It is the only
*B. jabrensis* and the probes are *jararaca*-based.

Being alone in the panel is not enough on its own: *B. germanoi* (I1773) is also
n=1 but has 143 singletons and normal capture. I20937 captures badly (9.7% on
target, enrichment 68.6, both near the bottom of the panel), and the same
divergence turns every fixed difference into a singleton. Not contamination:
that would raise heterozygosity, and its heterozygosity is among the lowest in
the panel (0.19% vs 0.46% median) with Ts/Tv at 2.07. Its highest-in-panel error
rate (0.0197) and missingness (7.7%) follow from the capture failure. Usable,
but expect it to look artificially diverged and homozygous against the jararaca
bulk.

**I1767, the depth outlier.** 122x vs a median of 47x, ~30x above the next
sample. Just sequenced deeper: 2nd-highest read count (6.2M vs 3.1M median) at
ordinary capture efficiency, normal reads-to-depth conversion (19.6x per million
vs 17.8). Lowest missingness in the panel (2.1%), 31 singletons, normal Ts/Tv
and error rate. Duplication at 10.8% (2nd highest) is the expected cost of
sequencing one library this deep. No action.

**Locus length vs. depth** (new plot). r = 0.85, and the causation runs both
ways:

- Loci with a larger capture target carry more probes, so they are both longer
  pre-trimming and better covered. Pre-trim length correlates with depth at
  r = 0.79 and is itself bimodal at ~800 vs ~1950 bp.
- ClipKIT trims terminal mostly-`N` columns, so a thinly covered locus loses its
  ends and finishes short: `pct_trimmed` vs depth is r = -0.84.

Colouring by `pct_trimmed` separates the two. Median depth is 29.8 for loci
under 1,000 bp (n=1,265) vs 63.8 for the longer ones (n=433), so the short half
of the locus set is where genotyping is weakest.

**Heterozygosity per sample** is now faceted by sex. The two panels are
near-identical, which is the expected result: with `z_policy=all` every
sex-linked locus is gone, so female hemizygosity can no longer depress female
heterozygosity (medians 0.44% female vs 0.49% male overall; 0.47% vs 0.49%
within jararaca). The bimodality visible in both panels is species, not sex --
the low cluster is all 20 non-*jararaca* samples (medians 0.13-0.24%) plus 29
jararaca.

Note: the global x-axis break doubling collided on the locus-length histogram
(0-2500 bp), so that one plot keeps 250-bp breaks.

## 2026-09-12 -- Sex-call confidence rewritten; 05b plot/table fixes

The `confident` flag from `04b_sex-and-z-class.sh` was close to meaningless. It
was `auto_depth >= 15 AND het-signal agrees`, where:

- the depth check never fired (lowest autosomal depth in the run is 29.6), and
- the female branch was `het_density < 0.0005 OR n_het_gametolog > 0`, and every
  female has 132-190 heterozygous gametolog sites, so no female could ever be
  flagged at any depth ratio.

The net effect was one flagged sample (I12279) and three obviously intermediate
samples (ratios 0.579, 0.585, 0.613) passing as confident.

Replaced with within-group robust z-scores. Sex is assigned by a gap split, so
the only thing that can be judged is how well a sample sits inside the group it
was put in. Three statistics -- depth ratio, hemizygous het density, and
gametolog het count -- are z-scored within the assigned sex group using the
median and a MAD-based SD (`--max_group_z`, default 3); a statistic whose group
MAD is zero is skipped, which is the usual case for female het density. A new
`flags` column in `sex_calls.tsv` records which check fired.

Now flags 5 of 125, each for a specific reason:

| sample | sex | ratio | flag |
|---|---|---|---|
| I1507  | female | 0.579 | depth_ratio (z=5.4) |
| I1795  | female | 0.585 | depth_ratio (z=5.7) |
| I1495  | female | 0.613 | depth_ratio (z=7.2) |
| I1624  | male   | 0.745 | gametolog_hets (z=5.4) -- 24 sites vs 0-16 for all other males |
| I12279 | male   | 0.880 | het_density (z=3.1) -- almost no het over the hemizygous Z |

Sex assignments and locus classes are unchanged; only the confidence flag moved.

Group summaries (written to stderr by 04b): female n=89, ratio median 0.4807
(SD 0.0183); male n=36, ratio median 0.8385 (SD 0.0324).

Other 05b fixes:

- `fmt_markdown()` was parsing table cells like `"0. Block splitting"` as a
  markdown ordered list, so the first column showed the text twice. Escaped the
  period.
- Table font 13 -> 15 px.
- `theme_bw` -> `theme_classic`, and x-axis break density doubled everywhere by
  shadowing `scale_x_continuous` in the global env with an `n.breaks = 10`
  default -- ggplot resolves a plot's default scale by looking up
  `scale_<aes>_<type>` in the global environment before its own namespace.
- Excess-het plot: the log1p y axis had default (linear-spaced) breaks, now
  explicit; x switched from site count to percentage of scored sites, with the
  0.5% threshold drawn in. On the percentage scale the four dropped loci
  (0.79-1.59%) sit well clear of the retained distribution, which tops out at
  0.42% -- so the threshold is not cutting into a continuous tail.

## 2026-09-11 -- Pipeline/filtering overview moved into the 05b QC report

`doc/filtering.md` is now reproduced as a "Pipeline overview" section at the top
of `scripts/05b_qc-loci.qmd`, so the parameter tables sit next to the numbers
they produced. Two things in `doc/filtering.md` were out of date and were
corrected on the way in:

- Section 4 still described the old per-locus het-density percentile filter.
  The current 04c filter is the per-site one (`min_call_frac` 0.25,
  `max_het_frac` 0.6, `max_site_frac` 0.005, `min_het_sites` 3), and the sex-linked
  locus removal by 04b was missing entirely.
- The count table in section 5 predates the rerun (1,750 loci / 126 samples).
  The current numbers are 1,698 loci / 125 samples, and they now come from
  `filtering-summary.tsv` rather than being typed out.

`doc/filtering.md` itself was left as is, so it will drift again unless it is
cut down to a pointer at the report.

All 05b tables are now `gt` (helpers `qc_gt()`/`doc_gt()` in the setup chunk):
integer-valued columns without decimals, everything else at 2, natural width
instead of full width, and code hidden (`echo: false`).

## 2026-09-11 -- The 6.47% variable sites from 04e is a sample-size artefact, not high diversity

The 04e summary line (1.62 Mb scored, 6.47% variable, 4.28% parsimony-informative)
looks far too high for what is essentially a within-*jararaca* dataset. It isn't.

`prop_variable` is the fraction of sites *segregating in the sample*, which grows
with the number of chromosomes sampled (~the harmonic number), so it is not a
diversity estimate. Rerunning 04e's counting on random subsets of the 105
*jararaca*, same loci, same thresholds:

| n   | variable | PI    |
|-----|----------|-------|
| 5   | 1.60%    | 0.61% |
| 10  | 2.46%    | 1.21% |
| 25  | 3.77%    | 2.11% |
| 50  | 4.88%    | 3.01% |
| 105 | 6.14%    | 4.11% |

The island endemics contribute almost nothing to the headline number: *jararaca*
alone gives 6.14% of the 6.47%. The 4.28% is within-species polymorphism.

The n-independent measure is per-sample heterozygosity (from the 04a
`*_counts-per-locus.tsv`, retained loci only):

- *jararaca* 0.44%, *sazimai* 0.24%, *insularis* 0.19%, *jabrensis* 0.19%,
  *alcatraz* 0.17%, *germanoi* 0.12%

0.44% is an unremarkable pi for a widespread squamate, and the 2-4x lower island
values are the expected bottleneck signature. Nothing here says the calls are
inflated: mean Ts/Tv is 2.01 (04f) and the paralog loci are already gone (04c).

Worth noting: theta_W for *jararaca* is 0.0614/5.92 = 1.04%, ~2.4x pi, i.e. a
large excess of rare alleles. Expected given 105 animals sampled across many
localities (structure + expansion), but singletons are also where residual error
would sit. The PI count excludes singletons by construction, so it is the more
robust of the two columns.

Also: "parsimony-informative" is a weak bar at n=125 -- under neutrality ~80% of
segregating sites clear MAC>=2 on both sides. A high PI% here is not evidence of
tree-like signal.

## 2026-09-11 -- 02c now intersects the mosdepth windows with the probe BED

Previously 02c only took the genome-wide mean out of the mosdepth summary and the
on-target side was approximated in 05a as `04a mean_dp / genome_mean_dp`. It now
does the intersection properly and emits `pct_on_target`, `target_mean_dp` and
`enrichment` itself.

Done in awk rather than with bedtools: the targets are indexed into 100 kb
buckets so each 500 bp window is only tested against the few that could contain
it, the same trick 04a uses. 4.4 s per sample, ~18 min for all 191. Checked
against `bedtools intersect -wo` on I12278: 62.52x both ways. The window-summed
total also agrees with mosdepth's own exact total to 0.001%, so treating depth as
uniform within a 500 bp window costs nothing here.

Results: mean 14.2% of sequenced bases on target (range 6.2-25.4%), mean
on-target depth 34.2x, mean enrichment 100x. The targets are 2.2 Mb of a 1.57 Gb
genome, so 14% on 0.14% of the genome is the enrichment working.

Two things this fixed rather than just refined:

- The old enrichment could only be computed for the 125 retained samples, since
  04a never ran on the rest. The mosdepth-based one exists for all 191, which is
  what made the flagged-contaminant comparison in 05b possible at all.
- `genome_mean_dp` was being read from mosdepth's own mean column, which is
  rounded to two decimals (0.31). At that magnitude the rounding is ~1.5% of the
  value and it is the denominator of the enrichment ratio. Now computed as
  bases/length from the exact counts in the same row.

`target_mean_dp` (34.2x mean) and 04a's `mean_dp` are different quantities and
both are kept: the first averages over every target position including uncovered
ones, the second averages FORMAT/DP over called sites only. They correlate at
r = 0.9999 with a tight ratio of 1.34-1.40, so the called-site mean runs ~36%
high, which is what excluding uncovered positions should do. 05a carries the
ratio as `dp_called_vs_target`; a sample off that trend has patchy rather than
low capture.

On the contamination flag: the 52 flagged samples do sit lower on capture
(median 10.9% vs 13.3% on target, enrichment 77x vs 94x) while mapping rate and
duplicate rate are indistinguishable. A shift in medians, not a separation, so it
is weak corroboration at best -- the flag is still mostly carrying information
from outside the sequencing.

## 2026-09-11 -- Proofreading 05b turned up two silent data bugs

Both were invisible in the rendered report, and one of them had produced a wrong
conclusion.

**`z_class` was being read as all-NA.** In `sample-by-locus-stats.tsv` the column
is NA for every autosomal locus, which is 851,000 of 881,000 rows, so readr
guessed logical from the leading rows and turned each real class ("hemizygous",
"gametolog", ...) into NA -- 30,000 parsing problems behind a warning nobody
reads. Nothing in 05b or 05c used the column yet, so it broke nothing visible.
Fixed with an explicit `col_types` in both docs.

**`geno_pct` is NA where a locus is uncalled in a sample** (`n_called == 0`; 897
rows). 05c already passed `na.rm = TRUE`; 05b did not, so a single uncalled pair
made the whole locus or sample mean NA. This mattered: the per-sample
heterozygosity correlations were computed over the biased subset that survived,
giving r(depth, het) = -0.30 and r(missingness, het) = +0.27. I had written those
up as a possible contamination signature. With `na.rm = TRUE` the same
correlations are +0.10 and -0.12 -- both near zero, meaning heterozygosity is
*not* set by coverage and the 05c comparisons are not confounded. The report now
says that.

Two smaller fixes. The excess-het table was pulling 04e's `n_sites_scored`, which
is NA for the four excess-het loci because 04e only scores the final set, so the
table showed NA for exactly the loci it was about -- it wants 04c's
`n_sites_het_scored`. And `max_het_frac` was the name of both the threshold
parameter (0.6) and the observed per-locus maximum in 05a's output; the column is
now `max_site_het_frac`.

Added to 05b: per-sample heterozygosity with the depth check above; a comparison
of the 52 flagged contaminants against the rest on the Sarek metrics, which are
the only statistics that exist for them (everything else is computed after they
are removed, so `enrichment` and any depth metric is empty for that group); a
cumulative-informativeness curve (343 loci carry half the parsimony-informative
sites, 1171 carry 90%); and a consolidated "samples to watch" list, which turns
up 7 samples failing two or more checks, `I20937` failing three.

## 2026-09-11 -- Four QC metrics added: informativeness, trimming, mapping, Ts/Tv

New scripts `02c_sarek-stats.sh`, `03d_clipkit-stats.sh`,
`04e_locus-informativeness.sh` and `04f_vcf-stats.sh`, all joined into the 05a
tables and plotted in 05b. No earlier step had to be rerun -- every one of these
reads output that was already on disk.

**Locus informativeness (04e).** The real gap. Every statistic in the pipeline
measured distance from the *reference*; nothing measured whether a locus
separates the *samples*. A locus can carry hundreds of differences from the
reference and still be constant across all 125 samples. Reads 04a's
`*_sites.tsv.gz` rather than the VCFs, per the rule set in 04a's header. Per
site, over samples passing the call-rate threshold: `carriers_alt = n_het+n_alt`,
`carriers_ref = n_called-n_alt`; variable when both > 0, parsimony-informative
when both >= 2. Result over the final 1698 loci: 1.62M scored sites, 6.47%
variable, 4.28% parsimony-informative, mean 41 informative sites/locus. Only
`L379` has none. Caveat: 04a collapses multi-allelics to het/alt/miss, so this is
a biallelic approximation.

**Ts/Tv (04f).** `bcftools stats -s -` on the merged VCF. 1.90-2.11 across the
125 samples, mean 2.012 -- no sample is calling noise, which is the first real
evidence that the genotype filtering is not just producing plausible-looking
output. Singletons average 289/sample but `I20937` has 3,587, a 12x outlier
worth a look. Ts/Tv had to come from the VCF: 04a records only whether a genotype
is het/alt/missing, never which bases are involved.

**ClipKIT trimming (03d).** Parses the per-locus `.fa.log` files 03b already
wrote, replacing the log greps in step 3B. Mean 29.1% of each alignment trimmed;
38 loci trimmed by >75%, worst `L12` at 99.4% (800 bp -> 5 bp). Two things to
know: ClipKIT signals failure by writing an *empty* log (not a missing one) next
to an empty `.untrimmed` marker, so an awk pass keyed on `FNR == 1` silently
skipped all 34 failures until they were handled separately. And ClipKIT's own
site classes in those logs are computed on the untrimmed alignment over all 191
samples, so they do not describe the final data set -- they stay in the source
file under a `pre_` prefix and 04e's counts are the ones to use.

**Sarek mapping/coverage (02c).** From `samtools stats`, the Picard
markduplicates metrics, and the mosdepth summary. Mean 3.0M reads, 99.2% mapped,
4.6% duplicates (range 2.4-14.4%); lowest mapping rate 96.3% (`I1594`). The
MultiQC `general_stats` table was the obvious source and the wrong one: 956 rows,
one per fastq/lane, with the useful columns mostly empty. The per-sample report
files parse cleanly instead.

Two notes on this run's Sarek output: its sample IDs are still the long form
(`P0021_FG_I1603`), so 02c strips to the short form; and mosdepth was run over
the genome rather than the probe intervals, so `genome_mean_dp` is genome-wide
(~0.3x). That turns out to be useful as the denominator of a capture enrichment
ratio against 04a's on-target depth: 83-243x across samples. An exact on-target
read fraction would need the mosdepth windows intersected with the probe BED for
every sample, which is not worth it for a QC ratio.

## 2026-09-11 -- Step 5 QC tables rebuilt on the 04a-04c output

`05a_prep-qc-tables.R` had fallen behind the 04x rewrite and could not run: it
read `per_locus_len.tsv` (now `per_locus_length.tsv`), and its sample-ID regexes
still expected the long `P0021_FG_I1603` form, so on short IDs `sub()` found no
match and silently left the *full file path* in the `sample` column. Its sanity
checks (161 samples, 1751 loci) were from an older run; it is 125 and 1762.

It also joined three BED-layout files per sample (`_het_`, `_hom-alt_`,
`_hom-ref_per-locus.tsv`) when 04a writes one `_counts-per-locus.tsv` holding all
of it plus `mean_dp`, and it re-derived `n_miss` as `length - n_called` when 04a
counts it directly. Now one read per sample instead of three joins.

Two correctness points, not just tidying:

- Only the 4 excess-het loci were excluded; the 60 sex-linked loci that 04b drops
  (`z_policy=all`) stayed in the QC tables, so the summaries described 1758 loci
  while the final FASTA/VCF have 1698.
- 05b/05c read `results/geno/consensus/stats/summaries/` and `sample-depths.tsv`;
  05a writes `results/stats/summaries/` and wrote `sample-dps.tsv`. Neither the
  directory nor the file name matched, so neither doc could have rendered.

Output is now four tables in `results/stats/summaries`: `sample-by-locus-stats.tsv`
(long, as before), `per-locus-stats.tsv`, `per-sample-stats.tsv`, and
`filtering-summary.tsv`. Dropped loci and samples are retained in the tables with
a `stage` column instead of being filtered out, so the reports can show what each
threshold removed rather than only what survived.

Metrics that were on disk but unused, now folded in: per-locus mean depth (04a),
sample sex and Z:autosome depth ratio (04b), Z locus class (04b), per-locus
excess-het stats (04c), pre-filter missingness (03c), and the spread of the
per-sample call rate per locus. The last one matters because a mean %N cannot
distinguish a locus called in 60/125 samples from one uniformly 40% N.

Still not collected, and worth doing: segregating and parsimony-informative sites
per locus (nothing currently measures locus informativeness -- all counts are
per-sample differences from the reference), ClipKIT trim fraction per locus
(currently only in slurm logs), and the Sarek mosdepth/markduplicates stats
(on-target and duplicate rate, the usual explanation for a low-depth sample).

## 2026-09-11 -- One counting step for 04b and 04c

04b and 04c were scoring excess heterozygosity in two different ways: 04b took a
female:male ratio of per-locus het *density*, 04c a per-site het *fraction*. The
per-locus statistic is the weaker one, and on this data it demonstrably fails:
`L1839` has a female:male het ratio of 1.00 and was classified PAR-like, but it
carries 17 sites with `f > 0.6` in each sex group. It is a Z-linked collapsed
duplication affecting both sexes, which a sex-contrast of densities cannot see by
construction.

Counting now happens once, in `04a_qc-vcf.sh`, and 04b and 04c both read its
output instead of opening the VCFs themselves.

**04a** derives everything from a single pass (previously three `bcftools norm`
pipelines plus three `bedtools intersect` calls) and adds `<sample>_sites.tsv.gz`:
one row per site whose genotype is not hom-ref, as chrom, pos, locus, code
(`het` / `alt` / `miss`). Hom-ref is the implicit default and is left out. Per
sample that is ~6.5k het, ~5k alt and ~86k missing rows against 1.68M sites, i.e.
~5% of the full sample-by-site matrix, so the whole matrix is reconstructable from
a file small enough to read repeatedly: at a site `n_het` is the number of samples
listing it as het, and `n_called` is the sample count minus those listing it as
missing.

`miss` also covers positions inside a locus with no VCF record at all -- ~5% of
positions, since the VCFs do not cover every base. Without them the denominators
would be wrong, so 04a emits the complement of the positions it saw per locus. The
per-locus counts now sum exactly to the locus length for all 1762 loci. The three
BED-layout count files 05a reads are byte-for-byte what the old implementation
produced, checked against it directly; runtime is unchanged at ~19 s per sample.

**04b** takes `--qc_dir` instead of `--vcf_in` and `--het_dir`, and classifies Z
loci on the per-site statistic within each sex group. That splits the collapsed
copies into two classes:

| class      | n  | meaning                                                     |
| ---------- | -- | ----------------------------------------------------------- |
| gametolog  | 4  | excess-het sites in the heterogametic group only -- W on Z  |
| paralog    | 1  | excess-het sites in both groups -- `L1839`                  |
| hemizygous | 50 | female:male het density <= --hemi_ratio                      |
| PAR-like   | 5  | none of the above                                            |

The four gametologs are the same ones the ratio found (`L1970`, `L489`, `L918`,
`L879`), now with 116 / 29 / 23 / 12 excess-het sites in females and 0 in males.
`--z_policy` gained a `paralogs` value (gametolog + paralog).

**04c** takes `--qc_dir` and `--bed` and aggregates the same per-site records, so
its own counting step is gone. Because 04a's rows are sparse this is a sort over
~11M lines rather than the 281M genotype records the previous version streamed out
of the VCFs.

**Ordering.** 04a has to stay ahead of 04b and 04c, but now because it is the
counting step rather than because 04b needed a side-effect of it. The per-locus QC
tables step 5 wants can be derived from the same per-site records after 04c, so the
earlier idea of running 04a a second time on the final locus set is unnecessary.

## 2026-09-11 -- The 99th-percentile excess-het cutoff has no support in the data

`04c_filter-excess-het.sh` dropped the top 1% of the per-locus het-density
distribution. A percentile cutoff removes a fixed number of loci whether or not
any locus is anomalous, and here the distribution gives no reason to cut at 1%:
after the Z is removed it is smooth and unimodal, median 0.00349, p99 0.00937,
max 0.01827 -- the most extreme autosomal locus is only 5.2x the median and
there is no gap anywhere in the tail.

The 4 Z gametologs are known collapsed paralogs and give a calibration point:
their het density is 0.029-0.049, i.e. 8-14x the autosomal median. The p99
cutoff sits at 2.7x the median, far below any real paralog signature.

**Per-site het separates paralogs cleanly; per-locus het density does not.**
For each site, the fraction of called samples that are het (`f`). Under HWE the
maximum is 0.5, so `f > 0.6` is already an excess. Counted over all 1762 loci
(sites with >= 50 calls): 1682 loci have zero such sites, 65 have 1-2, and only
15 have >= 3. Median frac of sites with `f > 0.6` is 0, p95 is 0, p99 is 0.002.
This is signal vs no signal, unlike the density.

Of the 18 loci the p99 filter currently drops, 14 have **zero** sites with
`f > 0.6` and a max `f` of 0.38-0.65 -- indistinguishable from random control
loci. They are simply diverse. Only `L1813`, `L473`, `L1306` and `L663` carry the
paralog signature. Meanwhile loci at density ranks 85, 113, 184, 566, 880 and
1077 (`L1826`, `L1914`, `L476`, `L1727`, `L1720`, `L1652`) do carry it and are
currently retained.

**BLAST cross-check works and needs no new BLAST run.** `blast_out_raw.tsv` is
unfiltered (down to 70.6% identity, 2% qcov), while `01d_filter-loci.R` already
drops probes with >1 hit at >= 90% id and >= 80% qcov -- so any paralog still in
the set is a diverged one that only the raw file sees. Scoring each probe by the
bitscore of its best off-target hit relative to its best hit (`sec_ratio`,
off-target = different chrom or > 10 kb away): 3.5% of loci exceed 0.4 overall,
but 22% of the 18 flagged loci do, a 6-fold enrichment. `L473` (0.65), `L1813`
(0.52) and `L994` (0.41) are corroborated this way. BLAST alone is weak though --
the 59 loci with `sec_ratio > 0.4` have a mean het density of 0.00412 vs 0.00375
for the rest.

**Proposed replacement for the percentile cutoff:** flag a locus when the
fraction of its sites with `f > 0.6` exceeds ~0.005 (p99.5) with at least 3 such
sites, and use `sec_ratio` from the raw BLAST output as corroboration. That
flags ~6-8 autosomal loci instead of 18, and they are the right ones.

**Also:** `L1839`, classified PAR-like by `04b` (female:male het ratio 1.03), has
12 sites with `f > 0.6` and max `f` 0.845. It is a Z-linked collapsed duplication
affecting both sexes, which is why the sex-contrast misses it. Moot under
`--z_policy all`, relevant if the policy changes.

**Implemented.** `04c_filter-excess-het.sh` now takes `--vcf_in` and `--bed` and
computes the per-site statistic itself, so it no longer depends on `04a`'s
per-locus het counts. Defaults: `--min_call_frac 0.25`, `--max_het_frac 0.6`,
`--max_site_frac 0.005`, `--min_het_sites 3`. Het density is still reported per
locus in `per_locus_het_stats.tsv` but no longer filters. On the current data it
flags 5 autosomal loci (`L473`, `L1813`, `L1306`, `L476`, `L663`) instead of 18,
and 1617 of 1702 loci have no excess-het site at all.

`--min_call_frac` is deliberately low. The excess-het sites of `L473`, `L476` and
`L663` sit in windows called in only 56-86 of 177 samples, against a median of
174, and are tightly clustered -- `L473`'s ten sites span 41 bp. That is the
signature of a diverged paralog that maps poorly in most samples, so raising the
threshold to 0.5 discards exactly the loci worth catching. `L1813` is the other
pattern: normal call rates, sites spread over the whole locus.

**Bug found in `03c_filter-loci.sh` step 8.** `metadata/contam_samples.txt` holds
short IDs (`I1674`) but VCF filenames use full IDs (`P0021_FG_I1674`). Step 1
matches FASTA headers with `seqkit grep -r`, but step 8 uses `grep -qxF` on the
full sample ID, so the 52 contaminated samples are removed from the FASTA set and
kept in the VCF set: 125 samples in `locus_filtered/`, 177 in `vcf/`. `04a`, `04b`
and `04c` all ran on 177. It does not change the het ranking -- contaminated
samples have slightly *lower* mean het density (0.00353 vs 0.00385) and the
flagged set is identical when computed on the 125 clean samples only -- but the
merged VCF from `04d` will carry all 52.

**Fixed by moving to short sample IDs throughout.** The sample sheet now puts the
short ID in the `patient`/`sample` columns, which is what Sarek names its output
directories and output files after, so the short form propagates from there.
`03c` step 1 now matches sample IDs exactly instead of using an end-anchored regex
to bridge the two forms.

The sample name inside the VCF needed handling separately: Sarek sets it to
`<patient>_<sample>`, so the current VCFs carry `P0021_FG_I12277_P0021_FG_I12277`
and a short-ID sample sheet alone would only turn that into `I12277_I12277`. `02b`
now reheaders its output VCF to `--sample_id`, which is what gets plain short IDs
into the merged VCF from `04d`. Getting
this into the existing output needs a rerun from `02b` onward; the `smp` line in
the `02b` loop strips any prefix so that works without rerunning Sarek.

## 2026-09-10 -- The excess-het filter was mostly catching the Z chromosome

The top loci flagged by the per-locus het-density filter were all on `Binsu_Z`:
the Z holds 3.4% of loci (60/1762) but supplied 5 of the top 18, and the top 4
outright. These are not autosomal collapsed paralogs.

Snakes are ZW. Splitting samples by Z:autosome depth ratio and contrasting het
density between the two groups sorts the 60 Z loci cleanly into three classes:

| class      | n  | female:male het | interpretation                          |
| ---------- | -- | --------------- | --------------------------------------- |
| gametolog  | 4  | 7.4-22          | W copy collapsed onto Z; true paralogs  |
| PAR-like   | 6  | ~1.0            | undifferentiated, behaves autosomally   |
| hemizygous | 50 | 0.03-0.39       | differentiated Z; females haploid       |

The 4 gametologs are `L1970`, `L489`, `L918`, `L879`. Removing the Z drops the
maximum het density from 0.0329 to 0.0183, i.e. the gametologs were the entire
top of the distribution and were distorting the percentile cutoff.

**Sex calls.** Mean DP over the 50 hemizygous Z loci, divided by autosomal mean DP,
is sharply bimodal: 116 samples at ~0.43 (ZW female) and 61 at ~0.81 (ZZ male),
a clean 2:1 contrast with a gap at 0.57-0.66. Both modes sit below 1.0 because Z
capture efficiency is lower overall; the ratio between modes is the signal. Depth,
het over the hemizygous loci, and het at L1970 agree for 154/177 samples. Of the
23 discordant, the low-coverage ones (autosomal DP 8-14) falsely look hemizygous
by het, and contaminated samples are mildly enriched (9/52 vs 11/125).

Sex is called by splitting the sorted depth ratios at their largest internal gap,
so nothing is hardcoded. Since the ratio is cleanest over hemizygous loci, and
those can only be identified once sex is known, the two steps are iterated twice.
The first pass (all Z loci) mis-assigns 7 samples and recovers only 2 of the 4
gametologs; the second pass converges on 116/61 and all 4.

**Biggest downstream consequence:** female Z genotypes are currently called
diploid-homozygous rather than haploid, so any per-sample heterozygosity estimate
is deflated in proportion to the Z fraction, confounded with sex. Anything using
genotypes (He, ROH, inbreeding, PCA, admixture) needs the Z dropped or females
recalled as haploid. `04b_sex-and-z-class.sh` defaults to `--z_policy all`
(drop all 60), which costs 3.4% of loci.

## 2026-09-10 -- 32 consensus jobs died with "check_val: command not found"

Cause: `source_function_script()` looked for the Bash function script at
`scripts/../dev/bash_functions.sh`, which does not exist since the submodule moved to
`core-scripts/`. It therefore fell through to downloading `bash_functions.sh` into the
working dir. With 191 jobs launched at once, two of them started the `wget` while the
other 32 tested `[[ -f ... ]]`, found the half-written file and sourced it, so no
functions were defined. The other 159 jobs happened to source a complete file and ran.

Fix: the scripts now search `../dev/`, `../core-scripts/dev/` and the working dir for a
non-empty copy, download to a temp file and `mv` it into place only as a last resort, and
abort with a clear message if the expected functions are still undefined. Applied to
01a, 01b, 02a, 02b, 03b, to `core-scripts/dev/template.sh`, and to the 144 other scripts
in `core-scripts/` that carry this block (all now identical to the template). The stale
`bash_functions.sh` copies in the repo root and in `scripts/` were leftovers from the
failed downloads and were deleted.

## 2026-09-10 -- snp_gap: 10 bp is not supported by the FreeBayes data

Checked the distance from every callable SNP call to the nearest indel across 24
samples (raw FreeBayes VCFs, after the DP/GQ and QUAL filters), using
`n ALT SNP genotypes / n callable non-indel records` per distance so that coverage
differences are controlled for.

| dist to indel | SNPs/kb | rel. to bg | % het | het rel. |
| ------------- | ------- | ---------- | ----- | -------- |
| 1-3 bp        | 4.7     | 0.82       | 74%   | 1.24     |
| 4-5 bp        | 10.3    | 1.76       | 69%   | 1.16     |
| 6-10 bp       | 10.8    | 1.86       | 68%   | 1.15     |
| 11-20 bp      | 9.6     | 1.65       | 69%   | 1.15     |
| 21-50 bp      | 8.7     | 1.50       | 68%   | 1.14     |
| 51-100 bp     | 8.3     | 1.43       | 67%   | 1.12     |
| 101-200 bp    | 7.6     | 1.30       | 66%   | 1.11     |
| >200 bp (bg)  | 5.8     | 1.00       | 60%   | 1.00     |

Two things:

1. **There is no pile-up immediately next to indels** -- 1-3 bp is *below*
   background (0.74-0.91x). FreeBayes merges variants within its haplotype window
   into a single complex/MNV allele, so those SNPs are never emitted as separate
   records (note the sharp step between 3 and 4 bp -- that's the window edge, not
   biology). The BCFtools-call artefact from 2026-03-22 is genuinely gone, as
   predicted in that entry.
2. **The elevation at 4-30 bp is indistinguishable from the elevation at
   51-200 bp** and has not reached background by 200 bp. Same for the mild het
   excess. That is a regional effect -- indels and SNPs both accumulate in
   variable/repetitive regions -- not local misalignment. There is no shoulder in
   the curve, so no cutoff separates artefact from real variation.

Cost of the window, now that it produces N rather than reference bases:

| `snp_gap` | % of variable sites blanked | sites/sample |
| --------- | --------------------------- | ------------ |
| 3         | 0.6%                        | ~71          |
| 5         | 1.5%                        | ~172         |
| 10        | 3.8%                        | ~434         |
| 15        | 5.8%                        | ~669         |

10 bp discards ~3.8% of variable sites at the same artefact-to-signal ratio as
anywhere else in the locus, i.e. it is not buying anything specific. Set to 5,
which covers the edge of the FreeBayes haplotype window at ~1.5%. 0 would be
defensible too, relying on the DP/GQ, QUAL and allele-balance filters.

Caveat: the >200 bp background includes other loci, so part of the 1.3x at
101-200 bp may be a between-locus difference rather than a within-locus gradient.
Doesn't change the conclusion -- what matters is that 4-30 bp looks like 51-200 bp.

## 2026-09-10 -- SNPs near indels set to missing instead of reverted to reference

`02b_bcftools-consensus.sh` used `bcftools filter --SnpGap`, which hard-filters:
the records are dropped, so `bcftools consensus` falls back to the reference base.
That turns "we don't trust this call" into a positive hom-ref assertion, and a
biased one -- it is precisely the ALT-bearing calls that get erased, pushing
heterozygosity down and reference-allele frequency up. Because indels differ
between samples, it also makes sample-specific artefacts look like real differences.

`--SnpGap N -S .` does **not** work: tested, `--set-GTs` has no effect on `--SnpGap`
failures and the records are still dropped. Two passes are needed:

```bash
bcftools filter --SnpGap "$snp_gap" -s SnpGap -m + -O u |
bcftools filter -e 'FILTER~"SnpGap"' --set-GTs . -O u
```

Verified end-to-end on a toy 20-bp locus: a SNP 2 bp from a deletion now comes out
as `N` (was the reference base); a SNP 8 bp away is unaffected; `--snp_gap 0`
disables the filter cleanly.

Related findings:
- `--SnpGap` only flags ALT-bearing SNP records -- monomorphic records inside the
  same window are untouched and keep their reference base. Remaining (unaddressed)
  reference bias; the thorough fix would be a per-sample BED of indel spans passed
  to `bcftools consensus --mask ... --mask-with N`.
- Discarded indels do *not* produce false reference bases: FreeBayes emits no
  monomorphic records for positions spanned by a deletion (checked in
  P0021_FG_I12277 at Binsu_ma-1:18981111, a 20-bp deletion -- no records for
  18981112-18981131), so `--absent N` already makes those positions `N`.
- `snp_gap` is 10 bp in `run/1_main.md` vs. a script default of 5. The 10 bp came
  from the 2026-03-22 BCFtools-call check (excess SNPs near indels trailing off
  around 15 bp), but that was pre-FreeBayes. Now that the window costs data rather
  than silently reverting it, worth re-checking the `bcftools stats` IDD section on
  the current FreeBayes VCFs before the rerun.

## 2026-09-10 -- bcftools consensus crash on multi-bp monomorphic blocks

5 of 192 `02b_bcftools-consensus.sh` jobs failed with:
```
[E::bcf_get_variant_type] Requested allele outside valid range
```
All 5 used `--bed probes_300bp-buffer.bed` (the narrowest buffer); samples
P0021_FG_I1238, I1631, I1652, I1685, I1690.

Root cause: FreeBayes (with `--report-monomorphic`) occasionally merges a run of
invariant positions with identical stats into a single VCF record -- REF several
bp long, ALT `.` (e.g. `POS=8603771 REF=GAACAACT ALT=.`). Confirmed via bisection
that when such a record starts inside a BED region but its REF span extends past
the region's end, `bcftools consensus` crashes -- it's given only the small
per-region reference slice (via `samtools faidx -r regions_file | bcftools
consensus`), not the whole chromosome, so the record's span runs off the end of
what it was handed. Narrower BED buffers make this more likely; block sizes
observed across a few samples topped out at 45bp, well within a 300bp region, so
it's specifically about a block straddling the *edge* of a region, not fitting
within one.

Fix in `scripts/02b_bcftools-consensus.sh` (now step 0 of the filtering, before
everything else): split any REF>1bp/ALT=`.` block record into one record per
base via a short awk pass. GT/DP/etc. are carried over unchanged to every base
of the run -- verified this produces identical consensus output to the unsplit
record for an ordinary in-region block, and resolves the crash for the region
that triggered it. ~9 such records per sample, out of ~2M total.

To fix the previous run: just re-submit the 5 failed samples with `--bed
probes_300bp-buffer.bed` through the updated script.

## 2026-09-10 -- Variable-site counts in step 02b were overstated

`bcftools view -v snps` counts *records*, but a record whose genotype the filters
set to missing still carries its ALT allele. Those become `N` in the consensus,
not variants. For `P0021_FG_I12277`: 18,270 SNP records, but only 12,755 with an
ALT genotype -- 5,123 (28%) missing and 392 hom-ref. The logged "Variable sites"
number and the per-region counts in `*_snp-counts-per-region.tsv` were therefore
~43% too high. Both now count genotypes (`-i 'GT="alt"'`) rather than records.

Also in `scripts/02b_bcftools-consensus.sh`:

- `bcftools stats` wrote to a fixed `$outdir/vcf_stats.txt`. Since the runner
  submits one job per sample into a shared `-o`, every job clobbered that file.
  Now `stats/${sample_id}_vcf-stats.txt`.
- `bcftools index` lacked `-f`, so re-runs failed on the existing index.
- The reference `.fai` must now exist up front; concurrent jobs would otherwise
  race to build it.
- The lookup-file "atomic create" block was dead code: `mv` always overwrites, so
  the wait loops it guarded were unreachable. Simplified to a plain rename.
- `--more_opts` was parsed and documented but never passed to anything.
- 5 CPUs were requested but never used; `--threads` now goes to the compressing
  steps.

New per-sample QC lines in the log: genotype breakdown (ALT / het / hom-ref /
missing) and a consensus summary (loci, total bp, mean %N, all-N loci). The
latter is the useful one -- for `P0021_FG_I12277` at 500bp buffer: mean 49.4% N
and 24 loci entirely N.

## 2026-09-10 -- Multiallelic positions dropped in step 02

`bcftools consensus` cannot render a position that carries more than one record,
and `norm --multiallelics -any --atomize` produces those both for genuine
multiallelic sites and for MNVs with a `1/2` genotype. The records get applied
independently, so the consensus base is wrong. Traced example:

```
# Raw VCF - one MNV, genotype 1/2:
Binsu_ma-1  241914287  CTATG  TTATA,CTATA  1/2

# After atomization - two records at 241914291, both G>A:
Binsu_ma-1  241914291  G  A  1/0  AO=15 RO=0 DP=28
Binsu_ma-1  241914291  G  A  0/1  AO=13 RO=0 DP=28
```

Both alleles carry A and there are no reference reads, so the truth is `A/A`.
The consensus got `R`. Genuine `1/2` SNPs fail the same way (`A>C` + `A>G`
should give `S`).

Scale per sample: 626 affected positions (3.6% of ALT-bearing positions), of
which 74.4% are same-ALT MNV artifacts. Of positions where both records were
het, 79.5% had the same ALT and 94.3% of those had RO=0 - i.e. really hom-ALT
but counted as two het calls, which also inflated the input to the excess-het
locus filter.

Added `--drop_multiallelic` (default `true`) to `scripts/02_bcftools-consensus.sh`
as a step 7: positions still carrying >1 record are removed, so
`bcftools consensus --absent N` makes them `N`. The positions are logged per
sample to `stats/<sample>_multiallelic_pos.txt`.

Verified on one sample: 626 positions dropped, 0 duplicate positions remaining,
het calls 7,336 -> 6,623, monomorphic records retained (2,016,481), and the
example position now `N` instead of `R`. Also confirmed on the same run that all
four genotype filters are now active - no non-missing genotype survives with
DP<5, GQ<20, QUAL<30 at an ALT site, or biallelic het AB outside 0.25-0.75.

NOTE: step 02 had already been re-run for 185/191 samples before this change,
so it needs re-running from step 02 onwards.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-09-09 -- Audit of genotype/locus filtering; three problems found and fixed

Wrote up all filtering parameters in `doc/filtering.md`. Three problems found:

1. **The `FORMAT/GQ < 20` genotype filter never did anything.** FreeBayes only
   emits `FORMAT/GQ` with `--genotype-qualities`, which was not passed (neither
   Sarek's germline default nor `config/nfc-sarek.config` includes it). GQ was
   `.` for every genotype, and a bcftools comparison against a missing value is
   false. Only `DP >= 5` was ever applied.

   ```bash
   bcftools view -H -i 'FORMAT/GQ < 20' <raw>.vcf.gz | wc -l   # 0
   bcftools view -H -i 'FORMAT/DP < 5'  <raw>.vcf.gz | wc -l   # 698,776
   ```

2. **Sarek's `*.freebayes.filtered.vcf.gz` must NOT be used.** It comes from
   `vcffilter -f 'QUAL > 30'`. Because we call with `--report-monomorphic` and
   FreeBayes gives invariant sites QUAL ~ 0 (QUAL = P(site is polymorphic)), that
   filter keeps 16,645 of 2,409,758 records - 0/2,389,990 monomorphic sites survive.
   Feeding that to `bcftools consensus --absent N` would blank out the consensus.
   Using the unfiltered VCF was correct; keep doing it.

3. **ClipKit silently lost 34 of 1819 loci.** It crashes (`IndexError` in
   `warn_if_entry_contains_only_gaps`) when `--ends_only` trimming removes every
   column, which happens when all columns exceed the `--gaps` threshold. The 34
   loci were 83-99.99% N, so the `max_missing_locus` filter would have dropped
   them anyway - but the loss was invisible.

Changes made:

- `config/nfc-sarek.config`: added `--genotype-qualities` to the FreeBayes
  `ext.args`. Note this `ext.args` replaces Sarek's own FreeBayes defaults
  (`--min-alternate-fraction 0.1 --min-mapping-quality 1`); `--min-mapping-quality 1`
  is the FreeBayes default anyway, and 0.05 vs 0.1 min-alternate-fraction affects
  <1% of het calls, so the override is fine.
- `scripts/02_bcftools-consensus.sh`: added `--min_qual` (default 30), applied only
  to ALT-bearing records (`ALT!="." && QUAL < 30`) so monomorphic sites are exempt;
  and `--min_ab`/`--max_ab` (0.25/0.75) for heterozygote allele balance, applied
  only to clean biallelic records (`(FMT/AO[0:0]+FMT/RO)==FMT/DP`). The biallelic
  restriction matters: after `bcftools norm -m -any` a true 1/2 genotype is split
  across two records and looks like AB ~ 1.0. 17.2% of het calls are these, and
  98.6% of AB>=0.9 hets are split multiallelics - an unrestricted AB filter would
  discard them. `FMT/AO` has `Number=A` and must be indexed explicitly for
  arithmetic in bcftools expressions to work.
- `scripts/05b_filter-excess-het.sh` (new): drops loci whose per-locus het density
  (`total het / (n_samples * length)`) exceeds the 99th percentile - the collapsed-
  paralog signature. On the current data: median 0.0015/bp, p99 0.0054, max 0.0217;
  18 loci flagged, worst are L1970, L1813, L1839 at 8-14x the median.
  `scripts/05_qc-vcf.sh` renamed to `05a_qc-vcf.sh` to match.
- `scripts/07a_prep-qc-tables.R`: excludes the excess-het blacklist, since the
  per-locus count tables in `stats/vcf` are produced before that filter.
- `run/1_main.md`: new settings, a ClipKit fallback that copies the untrimmed
  alignment when ClipKit produces no output (so step 04 filters it explicitly
  instead of it vanishing), and the new 05b step before merging.

Effect of the new genotype filters, tested on one sample/chromosome: het calls
2821 -> 2281, monomorphic records preserved (685,877), record count unchanged.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
## 2026-03-23 -- Coverage drop-off around probes

```bash
ref=data/ref/binsu_all.fna && ls -lh "$ref"
bed=results/blast_probes/probes_sorted.bed && ls -lh "$bed"

# Convert CRAM to BAM
cram=results/geno/nfc-sarek/preprocessing/markduplicates/P0021_FG_I1485/P0021_FG_I1485.md.cram
samtools view -b -T "$ref" -o P0021_FG_I1485.bam "$cram"

# Subtract 1000 from column 2 and add 1000 to column 3, if column 2 is below 0 after subtracting 1000, change it to 0:
awk -F'\t' -v OFS='\t' '{print $1, ($2-1000<0?0:$2-1000), $3+1000, $4}' "$bed" > buffer.bed

# Quantify depth
bedtools coverage -d -a buffer.bed -b P0021_FG_I1485.bam > depth.tsv
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-23 -- Checking the VCF filtering steps manually

```bash
# Manually
bcftools filter -e 'FORMAT/DP < 5 || FORMAT/GQ < 20' --set-GTs . "$vcf" -O u |
    bcftools filter --SnpGap 5 -O u |
    bcftools view -v snps -O z -o TEST.vcf.gz

bcftools filter -e 'FORMAT/DP < 5 || FORMAT/GQ < 20' --set-GTs . -O z -o filt1.vcf.gz "$vcf"
bcftools filter --SnpGap 10 -O z -o filt2.vcf.gz filt1.vcf.gz
bcftools norm -m -any --atomize -O u filt2.vcf.gz | bcftools view -V indels -O z -o filt3.vcf.gz

#? -a, --atomize                   Decompose complex variants (e.g. MNVs become consecutive SNVs)

zgrep -vc "^#" "$vcf" filt1.vcf.gz filt2.vcf.gz filt3.vcf.gz

bcftools index TEST.vcf.gz
bcftools consensus --haplotype I --missing N -f "$ref" TEST.vcf.gz > consensus.fasta
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-22 -- SNPs near indels with BCFtools call

I was calling SNPs with BCFtools call, but when I checked the distribution
of SNPS around indels, there were losts of SNPs near indels.

```bash
# Run BCFtools call
ref=data/ref/binsu_all.fna && ls -lh "$ref"
vcf=results/geno/bcftools/"$smp".vcf.gz && ls -lh "$vcf"
sbatch mcic-scripts/popgenom/bcftools-call.sh \
    --ref_fasta "$ref" \
    --bam_dir results/geno/filt-bam \
    --vcf "$vcf" \
    --allsites

# Check if many SNPs occur near indels - section `IDD` has indel distance
bcftools stats "$vcf" > vcf_stats.txt

grep IDD vcf_stats.txt
```

- Many more SNPs near indels, only trailing off at around 15 bp distance.
  Therefore it's better to use FreeBayes (or GATK), which will perform local
  assembly around indels and thus should be more robust to this issue.

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-03-01

- Andrew's Python script was run like so:

```bash
python ~/python_scripts/Anchored_mapping_consensus.py \
    -r ../../BothropsjararacaGenomeAssemblyFinalFinal.normalized.fasta \
    -bt ../../probes_blast_tab \
    -b ../../02_mapping/P0021_FG_I1698/P0021_FG_I1698_mark_dups.bam \
    -o P0021_FG_I1698.fasta
```

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

## 2026-02-23 - SECAPR => Nextflow workflow

- Made a Nextflow workflow based on `SECAPR` `reference_assembly.py`.
  But later decided to not use this, because it doesn't really do any variant calling.

```bash
reads=data/test2
outdir=results/nf-refasm
workdir=/fs/scratch/PAS1533/jelmer/nf-refasm
sbatch workflows/ref_assembly/main.sh \
    --local_wf workflows/nf-probes \
    --reads "$reads" \
    --reference "$ref" \
    --outdir "$outdir" \
    --workdir "$workdir"
```
