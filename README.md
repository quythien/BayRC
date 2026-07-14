# BayRC

**Bayesian Rhythmicity Comparison** (BayRC) is a unified statistical framework for comparing and interpreting circadian rhythms across biological conditions: age, disease state, tissue, species, or sex.

BayRC jointly infers gene-level rhythmicity and phase, computes posterior probabilities of rhythmic and phase concordance, and classifies rhythmic gain, loss, conservation, and direction-specific phase shifts under Bayesian false discovery rate (BFDR) control. It further supports pathway-level enrichment and genome-wide concordance scoring, providing a unified, uncertainty-aware framework for comparative circadian analysis across tissues, species, and disease contexts.

---

## Biological Questions Answered

BayRC answers a multifaceted set of biological questions by moving through three
levels of resolution: the gene, the pathway, and the genome.

| Biological Question | Statistical Output |
|---|---|
| Which genes oscillate with a 24-hour rhythm? | Posterior P(rhythmic) + Bayes Factor per gene |
| How confident are we in the amplitude and peak timing? | Posterior estimates + 95% credible intervals (circular HDI for phase) |
| Which genes gained, lost, or conserved rhythmicity? | BFDR-controlled transition classification |
| Do conserved genes peak at the same time of day? | Circular phase concordance + 95% HDI of Δφ |
| Which biological pathways are remodeled? | Two-stage pathway enrichment + gain-loss ratio (GLR) |
| How similar are two transcriptomes globally? | Adjusted Jaccard c-score + permutation p-value + bootstrap CI |

---

## Installation

```r
install.packages("devtools")
devtools::install("path/to/BayRC")

# During development:
devtools::load_all("path/to/BayRC")
```

**Dependencies:** `Rcpp`, `circular`, `ggplot2`, `dplyr`  
**Suggested:** `ComplexHeatmap`, `KEGGREST`, `biomaRt`, `edgeR`, `DESeq2`, `parallel`

---

## Quick Start

BayRC ships a real dataset: baboon putamen (PUT) and substantia nigra (SUN)
expression, 5,066 genes, 12 zeitgeber timepoints each, from the diurnal
transcriptome atlas of Mure et al. (GEO accession
[GSE98965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98965)),
ortholog-matched to the same gene backbone used throughout the manuscript.
It's the same PUT-vs-SUN comparison as the paper's Case Study 2, run fresh
through the current package code, not simulated data.

MCMC on the full 5,066 genes takes about 30 minutes per condition, so the
block below is realistic to run but not something you'd casually re-run
while reading. Every number in the comments is real output from that run,
including the modest ones: at these settings (2,000 iterations, 500 burn-in)
most genes carry weak evidence either way, and only 8 genes clear the BFDR
threshold for conservation, with zero clearing gain or loss. That's an
honest result at this chain length and sample size (12 samples per gene),
not a bug; the paper's actual results use longer chains on the same data.

```r
library(BayRC)

# ── STEP 1: Load the bundled data and run MCMC ────────────────────────────
baboon <- readRDS(system.file("extdata", "baboon_PUT_SUN_GSE98965.rds", package = "BayRC"))
n_genes <- nrow(baboon$expr_PUT)

data_list_PUT <- list(data = as.data.frame(log2(baboon$expr_PUT + 1)),
                      time = baboon$zt, gname = baboon$gene_symbol)
data_list_SUN <- list(data = as.data.frame(log2(baboon$expr_SUN + 1)),
                      time = baboon$zt, gname = baboon$gene_symbol)

init_PUT <- CB_init_single(Data.list = data_list_PUT, P = 24)
mcmc_PUT <- CB_MCMC_single_rj_slice(Data.list = data_list_PUT, Init.value = init_PUT,
                                    P = 24, iteration = 2000, n.burn = 500,
                                    p_rhythmic = rep(0.2, n_genes))
init_SUN <- CB_init_single(Data.list = data_list_SUN, P = 24)
mcmc_SUN <- CB_MCMC_single_rj_slice(Data.list = data_list_SUN, Init.value = init_SUN,
                                    P = 24, iteration = 2000, n.burn = 500,
                                    p_rhythmic = rep(0.2, n_genes))
# mcmc_PUT$rho — posterior samples of rhythmicity (0/1), one row per gene
# mcmc_PUT$phi — posterior samples of phase (hours, ZT scale)

# ── STEP 2: Annotate and summarize ────────────────────────────────────────
# match_symbols() must run before anything downstream, or the row names and
# Bayes factors below won't line up.
mcmc_PUT <- match_symbols(mcmc_PUT, BF = 3, p_rhythmic = 0.2)
mcmc_SUN <- match_symbols(mcmc_SUN, BF = 3, p_rhythmic = 0.2)

pA <- rowMeans(mcmc_PUT$rho)   # Pr(rhythmic | data), PUT
pB <- rowMeans(mcmc_SUN$rho)   # Pr(rhythmic | data), SUN
# mean P(rhythmic): 0.154 in PUT, 0.147 in SUN

est_PUT <- CB_getAllEst(mcmc_PUT, burn = 50)
names(est_PUT) <- c("A", "phi", "M", "sigma")   # returned unnamed; name once for $-access

# ── STEP 3: Transition classification (gain / loss / conserved) ──────────
trans <- transition_classify(pA, pB, bfdr_alpha = 0.25)
# tau_gain = 1, n_gain = 0 | tau_loss = 1, n_loss = 0 | tau_cons = 0.68, n_cons = 8

# ── STEP 4: Phase concordance among conserved genes ───────────────────────
phase <- phase_infer(phi_matrix1 = mcmc_PUT$phi, phi_matrix2 = mcmc_SUN$phi,
                     gain_loss_status = trans$gain_loss_status,
                     shift = 2, P = 24, bfdr_alpha = 0.25, compute_hdi = TRUE)
# of the 8 conserved genes: 0 phase-conserved, 8 phase-shifted

# ── STEP 5A: Two-stage pathway enrichment ─────────────────────────────────
kegg <- readRDS("data/kegg_pathway_list_hsa.rds")   # 354 KEGG pathways, real gene sets

result_union <- pathSelect(mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
                           pathway.list = kegg, dataset.names = c("A", "B"),
                           ranking.method = "union", score_type = "pos",
                           qvalue.cut = 0.20, nperm = 500)
active <- dplyr::filter(result_union$results, pval < 0.05)$pathway
# 11 of 220 testable pathways pass Stage 1 (pval < 0.05); top hit is KEGG
# Oxidative phosphorylation (pval = 1.3e-4), which lines up with the paper's
# own SUN-PUT finding that oxidative phosphorylation is enriched here

# ── STEP 5B: Genome-wide concordance ──────────────────────────────────────
global <- multi_conservation(mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
                             dataset.names = c("A", "B"),
                             select.pathway.list = "global",
                             n_perm = 200, n_boot = 200, use_cpp = TRUE,
                             save_output = FALSE)
# AdjustedConcordance = 0.044 (95% CI 0.038-0.049), p = 0.005 -- significant,
# though genuinely small given the short chain length above
```

---

## The BayRC Pathway Heatmap

A key deliverable of BayRC is an integrated pathway heatmap (Figure 5 in the manuscript) that reads **across five panels from left to right** for each gene in a pathway of interest:

| Panel | Shows | How to read it |
|---|---|---|
| 1. Rhythmicity status | Transition type (what happened biologically) | Orange = conserved rhythm, blue = loss in B, purple = gain in B |
| 2. Phase status | Whether peak timing shifted | Green = conserved (peaks align within ±δ hours), red = shifted (the clock resets) |
| 3. P(ρ=1 \| data), A and B | Posterior probability of oscillation in each condition | White to red gradient, 0 to 1; deep red = confident rhythmic, near white = flat. Gain genes are red only in B, loss genes only in A |
| 4. Phase posterior, condition A | MCMC posterior distribution of peak time (ZT −6 to 18) | A sharp blue bar means confident peak timing; a spread bar means high uncertainty. Low-rhythmicity genes show low intensity |
| 5. Phase posterior, condition B | Same as panel 4, for condition B | Same reading as panel 4 |

This design lets you read the entire circadian landscape of a pathway (which genes oscillate, when they peak, and whether that timing is preserved) in a single glance.

`plot_pathway_integrated()` builds this figure from `transition_classify()` and
`phase_infer()` output. Below is a real one, the KEGG Circadian rhythm pathway
(23 genes) in baboon putamen vs. substantia nigra, from the same GSE98965 data
used elsewhere in this README. At this run's MCMC scale no gene crosses the
BFDR threshold for gain/loss/conserved (hence "Non-classified" throughout the
left two panels), but the rhythmicity and phase-posterior columns show real,
gene-by-gene structure. A full-scale run, like the ones behind the manuscript
figures, is what resolves the left two panels into calls.

![KEGG Circadian rhythm pathway heatmap, baboon PUT vs SUN](man/figures/pathway_heatmap_demo.png)

---

## Key Functions

BayRC exports 21 functions, grouped below the same way the paper's Methods
section is organized (§2.1 through §2.4).

### 1. MCMC Core (paper §2.1)
| Function | Purpose |
|---|---|
| `CB_init_single()` | Initialize MCMC chain from cosinor fit or random draws |
| `CBt_init_single()` | Initialization for the pipeline's time-error-aware variant |
| `CB_MCMC_single_rj_slice()` | Core Reversible Jump MCMC sampler — main engine |
| `CB_getAllEst()` | Posterior point estimates + 95% credible intervals; uses `circular_HDI()` for phase |
| `CBt_sim_data()` | Simulate circadian data for testing and tutorials |
| `circular_HDI()` | Shortest-arc 95% credible interval for a phase posterior |
| `circular_median()` | Circular median of a phase posterior |

`circular_HDI()` matters because phase is periodic: a gene peaking near ZT23
and one peaking near ZT01 are one hour apart, not 22. A naive linear
credible interval would miss that. The panel below shows real posterior
phase distributions for four core clock genes (baboon putamen, real
GSE98965 data), with the 95% circular HDI drawn as the red arc. `DBP` and
`PER1` show an HDI that sits entirely within one day; `CSNK1E` and `NR1D2`
show the arc correctly wrapping through ZT0/ZT24, which is exactly the case
a linear interval gets wrong.

![Posterior phase distribution and 95% circular HDI for four core clock genes](man/figures/circular_hdi_demo.png)

### 2. Gene-Level Biomarker Detection with BFDR Control (paper §2.2)

> **Workflow position:** `match_symbols()` runs once per condition immediately
> after MCMC and before any downstream analysis. Skipping it causes silent
> failures further down the pipeline.

| Function | Purpose |
|---|---|
| `match_symbols()` | Annotate MCMC output with gene symbols; required before classification |
| `bfdr_from_posterior()` | BFDR threshold τ from a vector of posterior probabilities (paper Eq. 2) |
| `detect_rhy()` | Condition-specific rhythmic gene sets with BFDR control |
| `summarize_bay()` | Per-gene Bayes Factor: `BF = posterior_odds / prior_odds` |
| `transition_classify()` | Joint posterior BFDR for gain / loss / conservation |
| `transition_classify_marginal()` | Marginal per-condition BFDR (alternative to joint) |
| `phase_infer()` | Phase-shift vs. conservation classification + 95% circular HDI on Δφ |

### 3. Pathway-Level Rhythmic Enrichment and Directionality (paper §2.3)
| Function | Purpose |
|---|---|
| `pathSelect()` | Stage 1: `ranking.method="union"` (active pathways); Stage 2: `"gain"`, `"loss"`, `"conserved"` — one function, all stages |
| `plot_pathway_integrated()` | The five-panel pathway heatmap described above (Figure 5 in the manuscript) |
| `multi_conservation_pathway()` | Pathway-level concordance score for a chosen gene set |
| `multi_conservation_pathway_bootstrap()` | Pathway-level concordance with bootstrap confidence intervals |

### 4. Genome-Wide Concordance Summary (paper §2.4)
| Function | Purpose |
|---|---|
| `multi_conservation()` | Full pipeline: c-score + GLR + permutation p-value + bootstrap CI |

### 5. Cross-Species Alignment

Needed only when comparing across species (e.g. the baboon-human lung
comparison in the manuscript); skip for same-species comparisons.

| Function | When to use | What it does |
|---|---|---|
| `match_homologs()` | Automated, needs internet | Uses biomaRt to find 1:1 orthologs; aligns all datasets to reference gene space |
| `merge_mcmc()` | Reproducible, needs a pre-built ortholog table | Alternative to `match_homologs()` using an explicit ortholog database; more reproducible than live biomaRt queries |

---

## Reproducing the Manuscript Figures

Each figure in the paper is generated by a specific script in `inst/analysis/`.
The table below summarizes the mapping; see `inst/analysis/README_figures.md` for
full details on inputs and parameters.

| Figure | Description | Generating script |
|---|---|---|
| 1 | Framework overview flowchart | none (created externally as an illustration) |
| 2 | Genome-wide concordance heatmaps across 26 baboon tissues | `circa_concordance.plots.R` |
| 3 | Within-species phase concordance scatter (SCN-HIP and SUN-PUT) | `Baboon_SCN_HIP.R`, `Baboon_SUN_PUT.R` |
| 4 | SUN-PUT pathway enrichment dotplots | `Baboon_SUN_PUT.R`, `plot_enrich_SUN_PUT.R` |
| 5 | SUN-PUT heatmaps | `Baboon_SUN_PUT.R` |
| 6 | Cross-species lung circadian analysis (baboon vs human) | `Baboon_Human_LUN.R` |

These scripts depend on pre-computed MCMC `.RData` outputs for the baboon and
human tissue data, which are controlled-access or too large to bundle with the
package. The scripts are included here for transparency so the figures can be
traced to their source, but they are not runnable out of the box without that
underlying data.

---

## Citation

Pham T, Kauffman K. *BayRC: A Bayesian framework for comparative circadian genomics with FDR control and concordance scoring.* (manuscript in preparation)
