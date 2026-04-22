# BayRC

**Bayesian Rhythmicity Comparison** — a unified framework for comparative circadian genomics.

BayRC detects conserved and altered circadian rhythms across conditions, tissues, or species using Bayesian FDR control, circular phase concordance analysis, and genome-wide concordance scoring.

## Overview

BayRC has two layers:

1. **Core MCMC** (`R/mcmc.R`, `R/init.R`, `R/estimates.R`, ...): Gene-level Bayesian cosinor model with Reversible Jump MCMC. Run once per condition to produce posterior samples of amplitude, phase, MESOR, and rhythmicity indicator.

2. **Comparative analysis** (`R/internal.R`, `R/pathway.R`, `R/global.R`, ...): Classifies rhythm transitions (gain/loss/conservation) with Bayesian FDR, computes circular phase concordance, pathway enrichment (GLR), and genome-wide adjusted Jaccard c-score with permutation p-values and bootstrap CIs.

## Installation (from source)

```r
# Install devtools if needed
install.packages("devtools")

# Load from local clone
devtools::load_all("/path/to/BayRC")

# Or install
devtools::install("/path/to/BayRC")
```

## Quick start

```r
library(BayRC)

# Step 1: Run MCMC for each condition (run once per dataset)
mcmc_A <- CB_MCMC_single_rj_slice(Y = expr_A, t = time_A, ...)
mcmc_B <- CB_MCMC_single_rj_slice(Y = expr_B, t = time_B, ...)

# Step 2: Extract posterior rhythmicity probabilities
pA <- colMeans(mcmc_A$l)   # posterior P(rhythmic) per gene, condition A
pB <- colMeans(mcmc_B$l)

# Step 3: Classify transitions with BFDR control
trans <- transition_classify(pA, pB, bfdr_alpha = 0.05)
# trans$gain_genes, trans$loss_genes, trans$cons_genes

# Step 4: Phase concordance for conserved genes
phi_A <- mcmc_A$phi   # G x K posterior phase samples
phi_B <- mcmc_B$phi
phase <- phase_infer(phi_A, phi_B, trans$gain_loss_status)

# Step 5: Genome-wide concordance
cscores <- congruence(mcmc_A$l, mcmc_B$l)
pval    <- p_conservation_global(cscores, n_perm = 10000)
```

## Key functions

| Function | Purpose |
|---|---|
| `CB_MCMC_single_rj_slice()` | Core RJMCMC sampler |
| `CB_getAllEst()` | Extract posterior estimates + credible intervals |
| `detect_rhy()` | Identify rhythmic genes with BFDR |
| `transition_classify()` | Gain/loss/conservation classification (joint posteriors) |
| `transition_classify_marginal()` | Same, using marginal per-condition BFDR |
| `phase_infer()` | Phase concordance classification (conserved/shifted) |
| `congruence()` | Adjusted Jaccard c-score between two MCMC outputs |
| `perm_conservation_global()` | Permutation null for c-score |
| `p_conservation_global()` | Empirical p-values from permutations |
| `multi_conservation_pathway()` | Pathway-level enrichment and GLR |
| `pathSelect()` | Stage 1: select rhythmically active pathways |

## Repository layout

```
BayRC/
├── R/                  # All package functions
├── src/                # C++ acceleration (congruence, permutations)
├── inst/analysis/      # Reproducibility scripts (paper figures)
├── data/               # Pathway annotation data
├── vignettes/          # Workflow vignettes
├── tests/              # Unit tests
└── paper/              # LaTeX manuscript
```

## Data

Large pathway data files (`hw_orth.RData`, `human.pathway.list.RData`) are stored at:
`/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/pathway_data/`

These are referenced by absolute path in `inst/analysis/` scripts.

## Citation

Pham T, Kauffman K. *BayRC: A Bayesian framework for comparative circadian genomics.* (manuscript in preparation)

## GitHub

`git@github.com:quythien/BayRC.git`
