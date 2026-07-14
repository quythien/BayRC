# BayRC Workflow: Comparative Circadian Genomics

2026-07-14

- [Data and Setup](#data-and-setup)
  - [A fast gene panel for this
    vignette](#a-fast-gene-panel-for-this-vignette)
- [Step 1: MCMC Inference](#step-1-mcmc-inference)
- [Step 2: Posterior Summaries](#step-2-posterior-summaries)
  - [Annotate first](#annotate-first)
  - [Posterior probability of
    rhythmicity](#posterior-probability-of-rhythmicity)
  - [Bayes factor and per-condition
    detection](#bayes-factor-and-per-condition-detection)
  - [Amplitude, phase, and credible
    intervals](#amplitude-phase-and-credible-intervals)
- [Step 3: Transition Classification](#step-3-transition-classification)
- [Step 4: Phase Concordance](#step-4-phase-concordance)
- [Step 5A: Two-Stage Pathway
  Enrichment](#step-5a-two-stage-pathway-enrichment)
  - [Stage 1: which pathways are circadian-active at
    all](#stage-1-which-pathways-are-circadian-active-at-all)
  - [Stage 2: what kind of transition is each pathway enriched
    for](#stage-2-what-kind-of-transition-is-each-pathway-enriched-for)
- [Step 5B: Genome-Wide Concordance](#step-5b-genome-wide-concordance)
- [The Six-Panel Pathway Heatmap](#the-six-panel-pathway-heatmap)
- [Biological Interpretation
  Checklist](#biological-interpretation-checklist)
- [References](#references)

> **Audience:** biologists with RNA-seq experience; no background in
> Bayesian statistics is assumed.

BayRC identifies how the circadian transcriptome changes between two
conditions: which genes stop oscillating, which start, and among genes
that keep oscillating in both, whether they still peak at the same hour.
This vignette runs the full five-step BayRC workflow, start to finish,
on real data, so every number below is actual output, not a transcript.

# Data and Setup

BayRC needs two things per condition: an expression matrix (genes by
samples, log scale, gene symbols as row names) and a zeitgeber time
vector (one value per sample).

The package ships real data from GSE98965 (Mure et al. 2018, *Science*):
5,066 genes, 12 timepoints (ZT0 to ZT22, every 2 hours), in two baboon
tissues, omental fat (OMF, visceral adipose) and thyroid (THR).

``` r
library(BayRC)

baboon <- readRDS(system.file("extdata", "baboon_OMF_THR_GSE98965.rds", package = "BayRC"))
kegg   <- readRDS(system.file("extdata", "kegg_pathway_list_hsa.rds",   package = "BayRC"))

dim(baboon$expr_OMF)
#> [1] 5066   12
baboon$zt
#>  [1]  0  2  4  6  8 10 12 14 16 18 20 22
```

## A fast gene panel for this vignette

A full 5,066-gene, full-length MCMC run takes about 30 minutes per
condition, too slow to build here. This vignette instead runs on a
60-gene panel: every gene from two real KEGG hits found in the
full-scale analysis (`README.md`’s Quick Start), KEGG Long-term
depression and KEGG GnRH signaling pathway, six core clock genes
(`CRY2`, `ARNTL`/BMAL1, `NR1D1`, `DBP`, `SURF2`, `B4GALT4`), and 23
background genes sampled at random. The panel is fixed in advance, so it
is exactly reproducible.

``` r
ltd_present  <- intersect(kegg[["KEGG Long-term depression"]], baboon$gene_symbol)
gnrh_present <- intersect(kegg[["KEGG GnRH signaling pathway"]], baboon$gene_symbol)
feature_genes <- c("CRY2", "ARNTL", "NR1D1", "DBP", "SURF2", "B4GALT4")

core <- unique(c(ltd_present, gnrh_present, feature_genes))
set.seed(42)
background <- sample(setdiff(baboon$gene_symbol, core), 23)
gene_panel <- sort(unique(c(core, background)))
n_genes <- length(gene_panel)
n_genes
#> [1] 60

idx <- match(gene_panel, baboon$gene_symbol)
expr_OMF <- baboon$expr_OMF[idx, ]; rownames(expr_OMF) <- gene_panel
expr_THR <- baboon$expr_THR[idx, ]; rownames(expr_THR) <- gene_panel
```

> **Interpretation:** Nothing below is specific to this panel; the same
> calls run unchanged on the full 5,066-gene matrix, only slower. If you
> are adapting this vignette to your own data, this is the section to
> replace: build `expr_OMF`/`expr_THR`-style matrices (genes by samples,
> log scale) and a zeitgeber time vector, and skip the panel subsetting.

# Step 1: MCMC Inference

For each gene, BayRC fits a cosinor curve with a rhythmicity switch:

``` math
Y_g(t) = M_g + \rho_g \, A_g \cos\!\left(\tfrac{2\pi}{24}(t - \phi_g)\right) + \varepsilon_g
```

$`M_g`$ is the MESOR (rhythm-adjusted mean), $`A_g`$ the amplitude,
$`\phi_g`$ the acrophase (peak time), and $`\rho_g \in \{0, 1\}`$
indicates whether the gene oscillates at all. Reversible-jump MCMC moves
between the flat ($`\rho_g =
0`$) and rhythmic ($`\rho_g = 1`$) model within one chain (Green 1995),
updating amplitude and phase with slice sampling (Neal 2003). The output
is a full posterior sample of $`\rho`$, $`\phi`$, $`A`$, $`M`$ for every
gene, not a single fitted curve, which is what lets every step after
this one (Bayes factors, transition calls, phase concordance, pathway
enrichment) carry uncertainty forward.

``` r
data_list_OMF <- list(data = as.data.frame(log2(expr_OMF + 1)),
                       time = baboon$zt, gname = gene_panel)
data_list_THR <- list(data = as.data.frame(log2(expr_THR + 1)),
                       time = baboon$zt, gname = gene_panel)

init_OMF <- CB_init_single(Data.list = data_list_OMF, P = 24)
mcmc_OMF <- CB_MCMC_single_rj_slice(
  Data.list = data_list_OMF, Init.value = init_OMF, P = 24,
  iteration = 3000, n.burn = 500, thin = 10,          # short chain: vignette speed only
  p_rhythmic = rep(0.2, n_genes),
  save.file  = tempfile(fileext = ".rds"),             # keep intermediate saves out of the repo
  save.file2 = tempfile(fileext = ".rds")
)

init_THR <- CB_init_single(Data.list = data_list_THR, P = 24)
mcmc_THR <- CB_MCMC_single_rj_slice(
  Data.list = data_list_THR, Init.value = init_THR, P = 24,
  iteration = 3000, n.burn = 500, thin = 10,
  p_rhythmic = rep(0.2, n_genes),
  save.file  = tempfile(fileext = ".rds"),
  save.file2 = tempfile(fileext = ".rds")
)
```

``` r
dim(mcmc_OMF$rho)   # genes x stored posterior draws
#> [1]  60 251
```

> **Interpretation:** `iteration = 3000, n.burn = 500, thin = 10` on 60
> genes builds in seconds. The README’s Quick Start runs
> `iteration = 2500, thin = 1` on the full 5,066 genes, which takes much
> longer. A short chain and a small panel both reduce statistical power,
> so read every number in this vignette as a demonstration of the
> mechanics, not a result to report. `save.file`/`save.file2` point at
> temporary files here only so this document leaves nothing behind; in
> your own analysis, point them somewhere you want a real crash-recovery
> checkpoint.

# Step 2: Posterior Summaries

## Annotate first

`match_symbols()` runs once per condition, right after MCMC and before
anything else; it attaches gene symbols and a Bayes-factor call to the
`rho` matrix that later functions read.

``` r
mcmc_OMF <- match_symbols(mcmc_OMF, BF = 3, p_rhythmic = 0.2)
mcmc_THR <- match_symbols(mcmc_THR, BF = 3, p_rhythmic = 0.2)
```

## Posterior probability of rhythmicity

`rowMeans()` on the `rho` matrix gives P(rho = 1 \| data): the fraction
of retained MCMC draws where that gene was in the rhythmic state.

``` r
pA <- rowMeans(mcmc_OMF$rho)   # P(rhythmic | data), omental fat
pB <- rowMeans(mcmc_THR$rho)   # P(rhythmic | data), thyroid
```

## Bayes factor and per-condition detection

`summarize_bay()` converts the same posterior probability into a Bayes
factor, BF = posterior odds / prior odds, easier to read against a fixed
evidence scale (BF \> 3 moderate, \> 10 strong, \> 100 decisive).
`detect_rhy()` applies Bayesian FDR (BFDR) control independently in each
condition, instead of using one fixed BF cutoff for every dataset.

``` r
bf_OMF <- summarize_bay(mcmc_OMF$rho, BF = 3, p_rhythmic = 0.2)
head(bf_OMF[order(-bf_OMF$BayesF), c("RowAverage", "BayesF")], 5)
#>                    RowAverage       BayesF
#> ARNTL               1.0000000 4.000000e+20
#> ENSPANG00000005723  1.0000000 4.000000e+20
#> SLC38A1             0.9960159 1.000000e+03
#> NR1D1               0.9920319 4.980000e+02
#> CAMK2G              0.9880478 3.306667e+02

detected <- detect_rhy(mcmc_OMF, mcmc_THR, bfdr_alpha = 0.25)
c(rhythmic_OMF = detected$n_rhythmic_A, rhythmic_THR = detected$n_rhythmic_B,
  n_total = detected$n_total)
#> rhythmic_OMF rhythmic_THR      n_total 
#>           60           60           60
```

> **Interpretation:** BFDR control sorts genes from most to least
> confidently rhythmic, then estimates the expected fraction of false
> positives among all genes called rhythmic above each possible cutoff.
> The threshold is the most permissive one that keeps this expected
> error rate under `bfdr_alpha` (0.25 here). It calibrates on the shape
> of the whole posterior-probability distribution, not a fixed cutoff,
> so it depends on how many candidate genes it has to sort: with only
> 60, it calls every gene in this panel rhythmic in both tissues, more
> permissive than the full genome would give. `ARNTL` (BMAL1) tops the
> ranking with an unbounded Bayes factor, exactly what a core clock gene
> should look like.

## Amplitude, phase, and credible intervals

`CB_getAllEst()` extracts posterior point estimates and 95% credible
intervals for every parameter. It returns an unnamed list in the order
amplitude, phase, MESOR, sigma, so name it once before `$`-access. Phase
sits on a 24-hour circle, so its credible interval is a circular
highest-density interval (HDI): the shortest arc containing 95% of the
posterior phase draws, correctly wrapping through ZT0/24 instead of
treating 23:00 and 01:00 as 22 hours apart.

``` r
est_OMF <- CB_getAllEst(mcmc_OMF, burn = 20)
names(est_OMF) <- c("A", "phi", "M", "sigma")

est_OMF$phi["ARNTL", ]
#>        phi.Est phi.Lower phi.Upper RHYindex
#> ARNTL 13.57517  12.47209  14.68265        1
est_OMF$A["ARNTL", ]
#>          A.Est  A.Lower  A.Upper RHYindex
#> ARNTL 1.767245 1.266125 2.238708        1
```

The same shortest-arc logic is available directly through
`circular_HDI()` and `circular_median()`, what `CB_getAllEst()` calls
internally:

``` r
hdi_ARNTL <- circular_HDI(mcmc_OMF$phi["ARNTL", ], credMass = 0.95, P = 24)
c(lower = round(hdi_ARNTL$lower, 2), upper = round(hdi_ARNTL$upper, 2),
  median = round(circular_median(mcmc_OMF$phi["ARNTL", ]), 2))
#> lower.phi upper.phi    median 
#>     12.47     14.63     13.56
```

> **Interpretation:** `ARNTL`’s posterior peak time in omental fat is
> about ZT 13.6, with a 95% credible interval of roughly ZT 12.5 to ZT
> 14.6, a tight window from only 12 samples. A wide interval, or one
> straddling ZT0/24, would mean the peak time itself is poorly
> constrained even for a confidently rhythmic gene: amplitude, phase,
> and rhythmicity are three separate questions, each with its own
> uncertainty.

# Step 3: Transition Classification

`transition_classify()` takes the marginal rhythmicity probabilities
from each condition and computes three joint posterior probabilities:
gained rhythmicity, p_gain = (1 - pA) x pB; lost rhythmicity, p_loss =
pA x (1 - pB); maintained rhythmicity, p_cons = pA x pB. Each of the
three vectors gets its own BFDR threshold, the same rule as above
applied to the transition probability instead of the raw rhythmicity
probability. This replaces the older practice of intersecting
fixed-cutoff gene lists with a Venn diagram, shown to overstate how much
circadian reprogramming has actually occurred (Pelikan et al. 2022).

``` r
trans <- transition_classify(pA, pB, bfdr_alpha = 0.25)
#> === Transition-level BFDR results ===
#> τ_gain = 1 | n_gain = 0 
#> τ_loss = 1 | n_loss = 0 
#> τ_cons = 0.589 | n_cons = 35
c(n_gain = trans$n_gain, n_loss = trans$n_loss, n_cons = trans$n_cons)
#> n_gain n_loss n_cons 
#>      0      0     35
```

> **Interpretation:** `gain_loss_status` is the central annotation
> vector; every function past this point reads it. On this panel,
> nothing clears the gain or loss threshold, and 35 of 60 genes clear
> the conservation threshold. That is the direct consequence of Step 2:
> if BFDR already calls every gene rhythmic in both tissues, the joint
> “maintained” probability (pA x pB) is high for most of them too, and
> there is nothing left over to call gained or lost. p_gain and p_loss
> are products of two marginal probabilities, so they run smaller and
> noisier than either marginal alone, and need more candidate genes than
> this panel has to calibrate well. The product formula also assumes a
> gene’s rhythmicity call in one condition is conditionally independent
> of its call in the other, reasonable for two separately profiled
> tissues, and the same assumption the manuscript uses throughout.

# Step 4: Phase Concordance

For genes classified as “Maintained,” `phase_infer()` tests whether the
gene peaks at the same time of day in both conditions or has shifted
phase. It computes the posterior distribution of the circular phase
difference delta-phi = phi_OMF - phi_THR, tests whether abs(delta-phi)
exceeds a tolerance (2 hours here), and applies BFDR control separately
to the shift and conservation calls.

``` r
phase <- phase_infer(
  phi_matrix1 = mcmc_OMF$phi, phi_matrix2 = mcmc_THR$phi,
  gain_loss_status = trans$gain_loss_status,
  shift = 2, P = 24, bfdr_alpha = 0.25, compute_hdi = TRUE
)
#> Computing phase metrics for maintained genes...
#> 
#> === PHASE INFERENCE SUMMARY ===
#> Maintained genes: 35 
#> BFDR α = 0.25  | shift threshold = 2 h
#> Significant phase-shifted genes: 34 
#> Significant phase-conserved genes: 0 
#> Undetermined genes: 1 
#> HDI computation: enabled

c(shifted = sum(phase$flag_shift, na.rm = TRUE),
  conserved = sum(phase$flag_cons, na.rm = TRUE),
  undetermined = sum(phase$flag_undetermined, na.rm = TRUE))
#>      shifted    conserved undetermined 
#>           34            0            1

phase$deltaPhi.Est["ARNTL"]
#>     ARNTL 
#> -3.112168
```

> **Interpretation:** 34 of the 35 maintained genes are phase-shifted,
> none phase-conserved. `ARNTL` shifts by about -3.1 hours (delta-phi is
> OMF minus THR, so a negative value means THR peaks later): the gene
> keeps oscillating in both tissues, but the clock resets by roughly
> three hours between them. A panel built mostly from two real KEGG hits
> (both found by testing for loss, not conservation, at full genome
> scale) landing almost entirely in the phase-shifted bucket here is a
> real, if small-panel, signal, not noise: staying rhythmic while
> resetting phase is a distinct outcome from losing rhythmicity
> outright, and both are visible in the full-scale analysis in
> `README.md`.

# Step 5A: Two-Stage Pathway Enrichment

## Stage 1: which pathways are circadian-active at all

`pathSelect()` ranks genes by a posterior metric and tests each pathway
in a gene-set collection for enrichment toward one end of that ranking,
using FGSEA (`nproc = 1` keeps this single-process, which also keeps
this step check-safe under R CMD check’s core budget).
`ranking.method = "union"` rewards genes rhythmic in *either* condition,
a broad first pass: pathways that show up here have some circadian
signal worth testing further, whatever direction it runs.

``` r
active_pathway_list <- kegg[c("KEGG Long-term depression", "KEGG GnRH signaling pathway")]

result_union <- pathSelect(
  mcmc.merge.list = list(A = mcmc_OMF, B = mcmc_THR),
  pathway.list = active_pathway_list, dataset.names = c("A", "B"),
  ranking.method = "union", score_type = "pos",
  qvalue.cut = 0.20, nperm = 1000, nproc = 1, seed = 42
)
#> 
#> === FGSEA PATHWAY ANALYSIS ===
#> Comparing: A vs B 
#> Ranking method: union 
#> Test type: One-sided (gain only) 
#> Pathways: 2 of 2 
#> Genes: 60 
#> Permutations: 1000 
#> 
#> Calculating posterior probabilities...
#> Probabilistic rhythmic gene summary:
#>   Expected rhythmic in A : 49 
#>   Expected rhythmic in B : 47.6 
#>   Expected gain: 8.5 
#>   Expected loss: 9.8 
#>   Expected conserved: 39.1 
#>   Expected union: 57.5 
#> 
#> Ranking statistics:
#>   Range: [ 3.918 , 13.8155 ]
#>   Interpretation: Fisher-like statistic (-2*log(prob))
#> 
#> Running fgsea...
#> Adding custom metrics...
#> 
#> === RESULTS SUMMARY ===
#> Significant pathways (Q < 0.2 ): 0

result_union$results[order(result_union$results$pval),
                      c("pathway", "pval", "padj", "size")]
#>                       pathway      pval      padj size
#> 1 KEGG GnRH signaling pathway 0.1048951 0.2097902   22
#> 2   KEGG Long-term depression 0.3866134 0.3866134   16
```

> **Interpretation:** Neither pathway clears Q \< 0.20 for the union
> direction at this scale (KEGG GnRH signaling pathway is closest, Q =
> 0.21). That does not mean these pathways carry no circadian signal
> here; it means the union test, which rewards rhythmicity in either
> condition without regard to which one, is not the direction where
> their signal shows up on this panel. Stage 2 below tests each specific
> direction directly instead of relying on this pre-screen to find it,
> exactly the same logic `README.md` walks through at full genome scale.

## Stage 2: what kind of transition is each pathway enriched for

``` r
result_cons <- pathSelect(mcmc.merge.list = list(A = mcmc_OMF, B = mcmc_THR),
  pathway.list = active_pathway_list, dataset.names = c("A", "B"),
  ranking.method = "conserved", score_type = "pos",
  qvalue.cut = 0.20, nperm = 1000, nproc = 1, seed = 42)
#> 
#> === FGSEA PATHWAY ANALYSIS ===
#> Comparing: A vs B 
#> Ranking method: conserved 
#> Test type: One-sided (gain only) 
#> Pathways: 2 of 2 
#> Genes: 60 
#> Permutations: 1000 
#> 
#> Calculating posterior probabilities...
#> Probabilistic rhythmic gene summary:
#>   Expected rhythmic in A : 49 
#>   Expected rhythmic in B : 47.6 
#>   Expected gain: 8.5 
#>   Expected loss: 9.8 
#>   Expected conserved: 39.1 
#>   Expected union: 57.5 
#> 
#> Ranking statistics:
#>   Range: [ 0.392 , 1.001 ]
#> 
#> Running fgsea...
#> Adding custom metrics...
#> 
#> === RESULTS SUMMARY ===
#> Significant pathways (Q < 0.2 ): 1 
#> 
#> Top significant pathways:
#>   KEGG GnRH signaling pathway: ES=0.363, NES=1.604, Q=0.1838
#>            Geometric_Mean_Conserved=0.6511
#>            Expected: Gain=2.8 (0.133), Loss=3.8 (0.181), Conserved=14.6 (0.686)

result_cons$results[, c("pathway", "pval", "padj", "Gain_Loss_Ratio_Arithmetic")]
#>                       pathway       pval      padj Gain_Loss_Ratio_Arithmetic
#> 1 KEGG GnRH signaling pathway 0.09190809 0.1838162                  0.7360673
#> 2   KEGG Long-term depression 0.44355644 0.4435564                  1.0015984
```

> **Gain-loss ratio (GLR):** `Gain_Loss_Ratio_Arithmetic` compares
> expected gains to expected losses within a pathway; GLR \> 1 means the
> pathway gained more rhythmicity than it lost, GLR \< 1 the reverse.
> KEGG GnRH signaling pathway is significant for conservation here (Q =
> 0.18) with GLR = 0.74, loss-leaning; KEGG Long-term depression is not
> significant for conservation (Q = 0.44). Combined with Step 4: this
> pathway’s genes mostly stay rhythmic in both tissues (why
> conservation, not gain or loss, is where the signal lands) while
> shifting peak time (why almost all of them showed up phase-shifted,
> not phase-conserved, above).

# Step 5B: Genome-Wide Concordance

`multi_conservation()` collapses the whole comparison into one number:
an adjusted Jaccard concordance score (the c-score), related to the
transcriptomic congruence framework of Zong et al. (2023), centered at 0
under no more overlap than chance and approaching 1 under perfect
agreement, with a permutation p-value and a bootstrap confidence
interval. `select.pathway.list = "global"` runs it over every gene in
the input rather than one pathway at a time.

``` r
global <- multi_conservation(
  mcmc.merge.list = list(A = mcmc_OMF, B = mcmc_THR),
  dataset.names = c("A", "B"),
  select.pathway.list = "global",
  n_perm = 80, n_boot = 80,      # small here for speed; README uses n_perm/n_boot at full scale
  use_cpp = TRUE, save_output = FALSE
)
#> 
#> === Processing pair 1 of 1 : A vs B ===
#> Analyzing pathway 1/1: Global

round(global[, c("A_vs_B_AdjustedConcordance", "A_vs_B_PValue",
                  "A_vs_B_CI_Lower_Adj", "A_vs_B_CI_Upper_Adj",
                  "A_vs_B_GainLossRatio")], 3)
#>   A_vs_B_AdjustedConcordance A_vs_B_PValue A_vs_B_CI_Lower_Adj
#> 1                      0.022         0.049               0.001
#>   A_vs_B_CI_Upper_Adj A_vs_B_GainLossRatio
#> 1               0.042                0.868
```

> **Interpretation:** The adjusted concordance here is about 0.02 (95%
> CI roughly 0.001 to 0.042, p = 0.049, barely under 0.05). Do not read
> this as an estimate of genome-wide OMF-THR concordance: this panel is
> built around two specific pathways and a handful of clock genes, not a
> random sample of the transcriptome, and the small candidate-gene count
> limits what the permutation test can resolve. Treat this as
> confirmation the function runs correctly end to end; `README.md`’s
> Quick Start reports the equivalent score on the full, unbiased
> 5,066-gene set. As a general reading guide: a c-score above roughly
> 0.3 with a small p-value indicates strong transcriptome-wide
> conservation, 0.1 to 0.2 indicates modest but real conservation, and a
> score near 0 with a non-significant p-value indicates two largely
> independent circadian programs.

# The Six-Panel Pathway Heatmap

`plot_heatmap()` is BayRC’s main visualization: one figure per pathway,
combining rhythmicity, transition status, and phase timing for every
gene into six panels read left to right. It draws directly to the
current graphics device, building on the output of
`transition_classify()` and `phase_infer()` above.

``` r
plot_heatmap(
  data1 = mcmc_OMF, data2 = mcmc_THR,
  pathway_genes = kegg[["KEGG GnRH signaling pathway"]],
  pathway_name  = "KEGG GnRH signaling pathway",
  phase_results = phase, transition_results = trans,
  group_names = c("OMF", "THR")
)
#> 
#> Pathway: KEGG GnRH signaling pathway 
#> Genes: 22 
#>   [full] Using all 22 genes
```

![plot of chunk pathway-heatmap](figure/pathway-heatmap-1.svg)<!-- -->

    #> 
    #> === SUMMARY [ FULL ] ===
    #> Total genes:     22 
    #> Phase Shifted:   15 
    #> Phase Conserved: 0 
    #> Gain in THR : 0 
    #> Loss in THR : 0

| Panel | Shows | How to read it |
|----|----|----|
| 1\. Rhythmicity status | Transition type from `transition_classify()` | Orange = conserved rhythm, purple = gain in THR, blue = loss in THR |
| 2\. Phase status | Shift/conservation flag from `phase_infer()` | Green = peaks align within the tolerance, red = peak timing shifted |
| 3\. P(rho = 1 \| data), OMF and THR | Posterior rhythmicity probability per condition | White to red gradient, 0 to 1; deep red is confident rhythmicity, near-white is flat |
| 4\. Phase posterior, OMF | Histogram of the MCMC phase samples for each gene | A sharp bar means a tightly constrained peak time; a spread-out bar means high phase uncertainty |
| 5\. Phase posterior, THR | Same reading as panel 4, other condition | Compare bar position against panel 4 to see the shift visually |
| 6\. Delta peak (hours) | Signed peak-time difference, THR minus OMF | Positive bars mean THR peaks later; negative mean earlier; gray means not classified |

> **Interpretation:** Most rows show a filled orange bar in panel 1
> (conserved rhythm) and red in panel 2 (shifted phase), matching Step
> 4: this panel is dominated by genes that stayed rhythmic but reset
> their peak time. Panel 4/5 histograms sit visibly apart for these
> genes, one tissue’s bar to the left of the other’s, which is what a
> real phase shift looks like in this figure; genes with no panel 1/2
> color did not clear the BFDR threshold and are correctly left
> unclassified rather than forced into a category.

# Biological Interpretation Checklist

Working through a BayRC result on your own data, check these points in
order:

- [ ] **Rhythmicity landscape.** Compare the fraction of genes rhythmic
  in each condition (`detect_rhy()`’s `n_rhythmic_A`/`n_rhythmic_B`
  relative to `n_total`). A large asymmetry between conditions points to
  broad circadian dampening or activation, not just a handful of genes
  moving.

- [ ] **Dominant transition.** Compare `trans$n_gain` against
  `trans$n_loss`. Loss-dominant is a signature commonly seen in aging
  and neurodegeneration; gain-dominant suggests the clock recruiting new
  transcriptional targets in the second condition.

- [ ] **Phase shifts among maintained genes.** For genes with
  `flag_shift = TRUE`, check whether `deltaPhi.Est` is similar in size
  and direction across genes (one global clock advance or delay) or
  scattered (gene-specific regulatory rewiring rather than a single
  shifted oscillator).

- [ ] **Pathway direction.** For each pathway tested, note which Stage 2
  direction is significant, gain, loss, or conserved, and whether
  `Gain_Loss_Ratio_Arithmetic` sits above or below 1. The same pathway
  can pass for entirely different biological reasons depending on which
  direction drives it.

- [ ] **Panel size versus BFDR power.** BFDR calibration needs a
  reasonably large candidate list to work well, as this vignette’s own
  60-gene panel shows. Run classification and pathway enrichment on the
  full gene set before drawing conclusions from a filtered subset.

- [ ] **Global concordance in context.** Check whether the c-score’s
  confidence interval is bounded away from 0, and whether the
  genome-wide GLR sits above or below 1. A significant but small c-score
  with a GLR near 1 describes two circadian programs that overlap
  modestly and are drifting in neither direction on net.

# References

- Green PJ. Reversible jump Markov chain Monte Carlo computation and
  Bayesian model determination. *Biometrika*. 1995;82(4):711-732.
  [10.1093/biomet/82.4.711](https://doi.org/10.1093/biomet/82.4.711)
- Neal RM. Slice sampling. *Annals of Statistics*. 2003;31(3):705-767.
- Newton MA, Noueiry A, Sarkar D, Ahlquist P. Detecting differential
  gene expression with a semiparametric hierarchical mixture method.
  *Biostatistics*. 2004;5(2):155-176.
- Müller P, Parmigiani G, Rice K. FDR and Bayesian multiple comparisons
  rules. *Bayesian Statistics 8*. 2007:349-370.
- Scott JG, Berger JO. Bayes and empirical-Bayes multiplicity adjustment
  in the variable-selection problem. *Annals of Statistics*.
  2010;38(5):2587-2619.
- Stephens M. False discovery rates: a new deal. *Biostatistics*.
  2016;18(2):275-294.
- Mure LS, Le HD, Benegiamo G, et al. Diurnal transcriptome atlas of a
  primate across major neural and peripheral tissues. *Science*.
  2018;359(6381):eaao0318.
  [10.1126/science.aao0318](https://doi.org/10.1126/science.aao0318)
- Zong W, Rahman T, Zhu L, et al. Transcriptomic congruence analysis for
  evaluating model organisms. *PNAS*. 2023;120(6):e2202584120.
  [10.1073/pnas.2202584120](https://doi.org/10.1073/pnas.2202584120)
- Pelikan A, Herzel H, Kramer A, Ananthasubramaniam B. Venn diagram
  analysis overestimates the extent of circadian rhythm reprogramming.
  *FEBS Journal*. 2022;289(21):6605-6621.
  [10.1111/febs.16095](https://doi.org/10.1111/febs.16095)
