# BayRC

**Bayesian Rhythmicity Comparison** — a unified statistical framework for comparative circadian genomics.

BayRC enables rigorous comparison of circadian transcriptomic programs across conditions (tissues, species, age groups, disease states). It combines Reversible Jump MCMC for gene-level Bayesian inference with a hierarchical pipeline for biomarker classification, pathway enrichment, and genome-wide concordance scoring — all with calibrated uncertainty quantification.

---

## Biological Questions Answered

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

```r
library(BayRC)

# ── STEP 1: MCMC Inference ───────────────────────────────────────────────────
# Run once per condition. Each call produces posterior samples for one dataset.
# Inputs: expression matrix (G genes × N samples), Zeitgeber time vector.
# Core functions require a named list: data (data.frame), time, gname.
data_list_A <- list(data  = as.data.frame(expr_A),
                    time  = time_A,
                    gname = rownames(expr_A))
data_list_B <- list(data  = as.data.frame(expr_B),
                    time  = time_B,
                    gname = rownames(expr_B))

init_A <- CB_init_single(Data.list = data_list_A, P = 24)
mcmc_A <- CB_MCMC_single_rj_slice(Data.list  = data_list_A,
                                   Init.value = init_A,
                                   P          = 24,
                                   iteration  = 3000,
                                   n.burn     = 1000)
init_B <- CB_init_single(Data.list = data_list_B, P = 24)
mcmc_B <- CB_MCMC_single_rj_slice(Data.list  = data_list_B,
                                   Init.value = init_B,
                                   P          = 24,
                                   iteration  = 3000,
                                   n.burn     = 1000)
# Outputs (per gene, G × K matrices where K = floor(iteration / thin)):
#   mcmc_A$rho  — MCMC samples of ρ ∈ {0,1} (rhythmicity indicator)
#   mcmc_A$phi  — MCMC samples of φ (phase/acrophase, in hours, ZT scale)
#   mcmc_A$A    — amplitude posterior samples
#   mcmc_A$M    — MESOR (midline) posterior samples
# All matrices have rownames = gname from data_list_A.

# ── STEP 1.5: Annotate MCMC Output (REQUIRED before any downstream step) ─────
# match_symbols() adds gene symbols as rownames on all matrices,
# sets attr(rho, "symbols") and attr(rho, "RHYindex"), resolves duplicate symbols.
# If rownames are Ensembl IDs, pass ensemble = useMart(...) for auto-conversion.
mcmc_A <- match_symbols(mcmc_A, BF = 3)   # BF threshold for binary rhythmicity call
mcmc_B <- match_symbols(mcmc_B, BF = 3)
# For cross-species (e.g. human vs baboon): call match_homologs() AFTER match_symbols().
# See vignette section 3.5 for the full cross-species workflow.

# ── STEP 2: Posterior Summaries — Rhythmicity, Amplitude, Phase ──────────────
# 2a. Posterior rhythmicity probability
pA <- rowMeans(mcmc_A$rho)   # Pr(rhythmic | data), one value per gene
pB <- rowMeans(mcmc_B$rho)

# 2b. Bayes Factor for rhythmicity: BF = posterior odds / prior odds
#     BF >> 1 → strong evidence for rhythmicity (BF > 3 is a common threshold)
bf_A <- summarize_bay(mcmc_A$rho, BF = 3, p_rhythmic = 0.5)
# bf_A$BayesF       — per-gene Bayes Factor
# bf_A$Rhythmicity  — binary call (1 = rhythmic at BF threshold)
# bf_A$RowAverage   — same as pA (posterior mean)

# 2c. Point estimates and 95% credible intervals for all parameters
est_A <- CB_getAllEst(mcmc_A, burn = 100)
est_B <- CB_getAllEst(mcmc_B, burn = 100)
# est_A$phi — data frame: phi.Est (peak hour), phi.Lower, phi.Upper
#             95% circular HDI computed via sliding-window arc algorithm
# est_A$A   — A.Est, A.Lower, A.Upper  (amplitude, linear scale)
# est_A$M   — M.Est, M.Lower, M.Upper  (MESOR)

# ── STEP 3: BFDR-Controlled Transition Classification ────────────────────────
# Identifies which genes gained, lost, or maintained rhythmicity.
# Uses joint transition posteriors: p_gain = (1−pA)·pB, etc.

trans <- transition_classify(pA, pB, bfdr_alpha = 0.25)
# trans$gain_genes         — indices of genes that gained rhythm in B
# trans$loss_genes         — indices of genes that lost rhythm in B
# trans$cons_genes         — indices of genes with conserved rhythm
# trans$gain_loss_status   — character vector ("Gain"/"Loss"/"Maintained"/"Non-rhythmic")
# trans$tau_gain/loss/cons — BFDR thresholds used

# ── STEP 4: Phase Concordance within Conserved Genes ─────────────────────────
# For each gene with conserved rhythmicity, estimates the phase difference Δφ
# and classifies whether peak timing is conserved or shifted.

phase <- phase_infer(
  phi_matrix1      = mcmc_A$phi,          # G × K posterior phase samples
  phi_matrix2      = mcmc_B$phi,
  gain_loss_status = trans$gain_loss_status,
  shift            = 2,                    # tolerance δ in hours (e.g., 2 or 4)
  P                = 24,
  bfdr_alpha       = 0.25,
  compute_hdi      = TRUE
)
# phase$flag_cons       — logical, TRUE if gene is phase-conserved
# phase$flag_shift      — logical, TRUE if gene is phase-shifted
# phase$deltaPhi.Est    — posterior median of Δφ (hours, circular)
# phase$HDI_lower/upper — 95% circular HDI of Δφ (shortest arc on circle)
# phase$peak1, $peak2   — peak hour estimates in each condition

# ── STEP 5A: Two-Stage Pathway Enrichment ────────────────────────────────────
# Stage 1 — Union test: which pathways contain rhythmically active genes?
result_union <- pathSelect(
  mcmc.merge.list = list(A = mcmc_A, B = mcmc_B),
  pathway.list    = kegg_pathway_list,        # named list of gene vectors
  dataset.names   = c("A", "B"),
  ranking.method  = "union",                  # S_g = −2·log[(1−pA)(1−pB)]
  score_type      = "pos",
  qvalue.cut      = 0.20,
  nperm           = 1000
)
active_pathways     <- dplyr::filter(result_union$results, pval < 0.05)$pathway
active_pathway_list <- kegg_pathway_list[active_pathways]

# Stage 2 — Transition enrichment within active pathways:
result_gain <- pathSelect(mcmc.merge.list = list(A = mcmc_A, B = mcmc_B),
                          pathway.list = active_pathway_list,
                          dataset.names = c("A", "B"),
                          ranking.method = "gain",      nperm = 1000)
result_loss <- pathSelect(mcmc.merge.list = list(A = mcmc_A, B = mcmc_B),
                          pathway.list = active_pathway_list,
                          dataset.names = c("A", "B"),
                          ranking.method = "loss",      nperm = 1000)
result_cons <- pathSelect(mcmc.merge.list = list(A = mcmc_A, B = mcmc_B),
                          pathway.list = active_pathway_list,
                          dataset.names = c("A", "B"),
                          ranking.method = "conserved", nperm = 1000)
# Each result$results: pathway name, pval, qval, GLR (gain-loss ratio),
#                      n_gain, n_loss, n_cons

# ── STEP 5B: Genome-Wide Concordance ─────────────────────────────────────────
global <- multi_conservation(
  mcmc.merge.list = list(A = mcmc_A, B = mcmc_B),
  dataset.names   = c("A", "B"),
  n_perm          = 1000,
  n_boot          = 1000,
  use_cpp         = TRUE
)
# global$c_score   — adjusted Jaccard concordance (centred at 0; range ≈ −1 to 1)
# global$GLR       — genome-wide gain-loss ratio (log GLR > 0 → activation dominant)
# global$pvalue    — permutation p-value
# global$CI_lower / $CI_upper — 95% bootstrap confidence interval on c-score
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

---

## Key Functions

### MCMC Core
| Function | Purpose |
|---|---|
| `CB_MCMC_single_rj_slice()` | Core Reversible Jump MCMC sampler — main engine |
| `CB_init_single()` | Initialize MCMC chain from cosinor fit or random draws |
| `CBt_MCMC_single()` | MCMC with spike-and-slab prior on time-of-measurement error |
| `CBt_init_single()` | Initialization for the time-error variant |
| `CB_getAllEst()` | Posterior point estimates + 95% credible intervals; uses `circular_HDI()` for phase |
| `Cosinor_fit()` | Fast OLS cosinor fit (initialization baseline / comparison) |

### Posterior Summaries
| Function | Purpose |
|---|---|
| `summarize_bay()` | Per-gene Bayes Factor: `BF = posterior_odds / prior_odds` |
| `bfdr_from_posterior()` | BFDR threshold τ from a vector of posterior probabilities |

### Biomarker Classification
| Function | Purpose |
|---|---|
| `detect_rhy()` | Condition-specific rhythmic gene sets with BFDR control |
| `transition_classify()` | Joint posterior BFDR for gain / loss / conservation |
| `transition_classify_marginal()` | Marginal per-condition BFDR (alternative to joint) |
| `phase_infer()` | Phase-shift vs. conservation classification + 95% circular HDI |

### Pathway Enrichment (Two-Stage)
| Function | Purpose |
|---|---|
| `pathSelect()` | Stage 1: `ranking.method="union"` (active pathways); Stage 2: `"gain"`, `"loss"`, `"conserved"` — one function, all stages |

### Genome-Wide Concordance
| Function | Purpose |
|---|---|
| `multi_conservation()` | Full pipeline: c-score + GLR + permutation p-value + bootstrap CI |
| `concordance()` | Pairwise adjusted Jaccard c-score from two MCMC output objects |
| `perm_conservation_global()` | Permutation null distribution for c-score |
| `p_conservation_global()` | Empirical p-values from permutation output |

### Annotation & Gene Alignment

> **Workflow position:** These functions run once per condition immediately after
> MCMC and before any downstream analysis. Skipping them causes silent failures.

| Function | When to use | What it does |
|---|---|---|
| `match_symbols()` | **Always** — after every MCMC run | Sets gene symbols as rownames + attributes on all matrices; computes binary RHYindex; resolves duplicate symbols by Bayes Factor; converts Ensembl IDs if needed |
| `match_homologs()` | Cross-species only — after `match_symbols()` | Uses biomaRt to find 1:1 orthologs; aligns all datasets to reference gene space |
| `merge_mcmc()` | Cross-species with pre-built ortholog table | Alternative to `match_homologs()` using an explicit ortholog database (e.g. `hw_orth.RData`); more reproducible than live biomaRt queries |

> **`match_homologs` vs `merge_mcmc`:** Both do cross-species alignment.
> Use `match_homologs()` for automated queries via biomaRt (requires internet).
> Use `merge_mcmc()` when you have a curated ortholog table and want reproducible
> results independent of database updates.

### Concordance Scoring

| Function | Use case | Mathematical basis |
|---|---|---|
| `compute_iteration_jaccard()` | **Canonical c-score** — propagates full MCMC uncertainty | Per-iteration Jaccard on binary ρ samples; mean over K iterations |
| `compute_concordance_minimal()` | Quick wrapper for `compute_iteration_jaccard()` | Returns scalar obs/adj/null Jaccard |
| `congruence()` | Analytical expected c-score (no permutation) | E[intersection] / E[union] using continuous posterior probabilities |
| `concordance2()` | Diagnostic / comparison | Conditional concordance + hard-threshold Jaccard on binarized posteriors |
| `perm_conservation_global()` | Full pipeline: c-score + p-value + bootstrap CI | Combines `compute_iteration_jaccard()` with permutation null |

> **Which to use:** `compute_iteration_jaccard()` is the paper-canonical statistic.
> `congruence()` is a fast analytical approximation. `concordance2()` is a legacy
> function for comparison with published point-estimate methods.

### Utilities
| Function | Purpose |
|---|---|
| `CBt_sim_data()` | Simulate circadian data with time-of-measurement error |
| `adjust.to.2pi()` | Map angles to [0, 2π] (used internally by MCMC) |
| `circular_HDI()` | Shortest-arc 95% credible interval for phase posteriors |
| `bfdr_from_posterior()` | Compute BFDR threshold from a vector of posterior probabilities |

---

## Repository Layout

```
BayRC/
├── R/                    # Package functions
│   ├── mcmc.R            # CB_MCMC_single_rj_slice (RJMCMC engine)
│   ├── mcmc_time.R       # CBt variant (time-of-measurement error)
│   ├── init.R / init_time.R
│   ├── estimates.R       # CB_getAllEst, circular_HDI
│   ├── internal.R        # summarize_bay, detect_rhy, transition_classify,
│   │                     # phase_infer, concordance, bfdr_from_posterior
│   ├── path_select.R     # pathSelect (two-stage enrichment)
│   ├── pathway.R         # multi_conservation_pathway
│   ├── global.R          # multi_conservation, permutation, bootstrap
│   ├── merge.R           # ortholog matching across species
│   ├── simulate.R        # CBt_sim_data
│   └── utils.R           # adjust.to.2pi, circular utilities
├── src/                  # C++ acceleration (congruence, permutations)
├── inst/analysis/        # Reproducibility scripts for paper figures
│   ├── Baboon_SUN_PUT.R  # Fig 2–5: within-species (SUN vs PUT)
│   ├── Baboon_Human_LUN.R # Fig 6: cross-species (lung)
│   ├── Human_Y_O.R       # aging comparison
│   └── pipeline/         # batch MCMC pipeline scripts
├── data/                 # Pathway annotation data (KEGG, GO)
├── vignettes/            # Workflow vignettes
└── paper/                # Manuscript (local only)
```

---

## Citation

Pham T, Kauffman K. *BayRC: A Bayesian framework for comparative circadian genomics with FDR control and concordance scoring.* (manuscript in preparation)
