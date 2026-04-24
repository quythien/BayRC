# BayRC Workflow: A Complete Guide to Comparative Circadian Genomics

**Package:** BayRC — Bayesian Rhythmicity Comparison  
**Audience:** Biologists with RNA-seq experience; no Bayesian statistics background required  
**Estimated reading time:** 30–45 minutes

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Data Requirements and Setup](#2-data-requirements-and-setup)
3. [Step 1: MCMC Inference](#3-step-1-mcmc-inference)
4. [Step 2: Posterior Summaries](#4-step-2-posterior-summaries)
5. [Step 3: Transition Classification](#5-step-3-transition-classification)
6. [Step 4: Phase Concordance](#6-step-4-phase-concordance)
7. [Step 5A: Pathway-Level Analysis](#7-step-5a-pathway-level-analysis-two-stage)
8. [Step 5B: Genome-Wide Concordance](#8-step-5b-genome-wide-concordance)
9. [Visualization — The BayRC Pathway Heatmap](#9-visualization--the-bayrc-pathway-heatmap)
10. [Interpreting Your Results — A Biological Checklist](#10-interpreting-your-results--a-biological-checklist)
11. [Function Quick Reference](#11-function-quick-reference)

---

## 1. Introduction

Every cell in the body runs on a roughly 24-hour molecular clock, coordinating
metabolism, stress responses, immune function, and hundreds of other processes to
the appropriate time of day. When you compare two biological states — two brain
regions, two tissues, two age groups, diseased versus healthy — a fundamental
question emerges: **has the circadian transcriptome been remodeled?** Which genes
switched from oscillating to flat? Which gained rhythmicity? Do genes that kept
oscillating still peak at the same hour?

Standard differential expression tools answer questions about *mean* expression
levels. BayRC is built to answer a different question: the **gain, loss,
conservation, and timing** of 24-hour oscillations between two conditions.

The statistical foundation is Bayesian inference. Where a classical method returns
a single yes/no decision for each gene, a Bayesian approach returns a **full
probability distribution** over all possible parameter values — amplitude, phase,
and whether the gene even oscillates. A gene with P(rhythmic) = 0.95 deserves a
different interpretation than one at P(rhythmic) = 0.51, even if both cross the
same threshold. BayRC propagates that uncertainty through every downstream step,
so the confidence you have in a single gene's rhythmicity directly feeds into
pathway enrichment scores and concordance calculations.

The deliverables of a full BayRC analysis are: (1) per-gene posterior probability
of rhythmicity and Bayes Factor; (2) amplitude and phase point estimates with
95% credible intervals; (3) BFDR-controlled transition calls (Gain / Loss /
Maintained / Non-rhythmic); (4) phase-shift versus conservation classification for
rhythmically maintained genes; (5) two-stage pathway enrichment with gain-loss
ratio; and (6) a genome-wide concordance c-score with permutation p-value and
bootstrap confidence interval.

BayRC operates in two layers. The **gene layer** runs per-condition Reversible Jump
MCMC to generate posterior samples of all oscillation parameters for every gene.
The **comparison layer** then uses those posteriors jointly — never collapsing them
to point estimates prematurely — to classify transitions, infer phase shifts, and
test pathway enrichment.

---

## 2. Data Requirements and Setup

### Input format

BayRC needs two objects per condition:

- **Expression matrix** (`Y`): a numeric matrix with **G rows (genes) × N columns
  (samples)**. Row names must be gene symbols (or any consistent identifiers).
  Values should be log-normalized counts: log-CPM from edgeR, variance-stabilized
  counts from DESeq2, or similar log-scale data. Do not use raw counts.
- **Zeitgeber time vector** (`t`): a numeric vector of length N giving the
  collection time of each sample in hours on a 24-hour scale (e.g., ZT0 = lights
  on, ZT12 = lights off). The order must match the column order of the expression
  matrix.

```r
# Check your input dimensions before proceeding
dim(expr_A)   # should be [G genes, N samples], e.g., [15000, 48]
length(time_A) # must equal ncol(expr_A), e.g., 48

head(time_A)
# [1]  0  2  4  6  8 10 ...
# Time points do not need to be evenly spaced or balanced.
# Biological replicates at the same ZT are entered as separate columns.
```

> **Note on preprocessing:** BayRC fits each gene independently using a cosinor
> model. Low-count genes with near-zero variance will have poor MCMC mixing. Filter
> out genes with very low expression (e.g., CPM < 1 in fewer than half the samples)
> before running MCMC.

### Installation

```r
# Install from source directory
install.packages("devtools")
devtools::install("/path/to/BayRC")

# During active development, load without a formal install:
devtools::load_all("/path/to/BayRC")

# Required dependencies
install.packages(c("Rcpp", "circular", "ggplot2", "dplyr", "HDInterval"))

# For pathway enrichment:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("fgsea")

# For concordance output (Excel):
install.packages("openxlsx")
```

### Loading your data

```r
library(BayRC)

# Load expression matrices and time vectors for two conditions.
# Example: Condition A = Putamen (PUT), Condition B = Substantia Nigra (SUN),
# mirroring the baboon brain region comparison from the BayRC paper.
load("expression_data.RData")  # loads expr_PUT, expr_SUN

# Zeitgeber time: typically collected across the full 24-hour cycle
time_PUT <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22,   # replicate 1
              0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)    # replicate 2
time_SUN <- time_PUT  # same sampling design in both conditions here
```

---

## 3. Step 1: MCMC Inference

### What is being estimated?

For each gene g, BayRC fits a cosinor (sinusoidal) model to the expression time
series:

```
Y_g(t) = M_g + A_g * cos( (2π/24) * (t − φ_g) ) + ε_g
```

- **M** (MESOR): the midline expression level — the 24-hour average
- **A** (amplitude): half the peak-to-trough distance; larger = more pronounced
  daily swing in expression
- **φ** (phi, acrophase): the time of peak expression in hours on the ZT scale
- **ρ** (rho): a binary switch — ρ = 1 means the cosine term is active (the gene
  oscillates); ρ = 0 means the gene is flat and A is effectively zero

The key deliverable of MCMC is not just one estimate for each parameter but a full
**posterior distribution**: thousands of plausible combinations of (M, A, φ, ρ)
that are consistent with the data. The proportion of MCMC samples where ρ = 1
directly estimates the posterior probability P(rhythmic | data).

### What is Reversible Jump MCMC doing?

Standard MCMC explores a fixed parameter space. But the flat model (ρ = 0) and the
rhythmic model (ρ = 1) have different numbers of parameters — the flat model has no
amplitude or phase, while the rhythmic model does. **Reversible Jump MCMC (RJMCMC)**
handles this by allowing the chain to jump between these two model spaces in a
single unified run, proposing moves that add or remove the oscillation parameters
in a way that respects detailed balance. The result is a single MCMC chain that
naturally estimates both the probability of rhythmicity *and* the parameter
distributions conditional on rhythmicity. Slice sampling is used within the
rhythmic model to ensure efficient exploration of the amplitude-phase space.

### Running MCMC for both conditions

```r
# ── Condition A (Putamen) ────────────────────────────────────────────────────
# Core functions require input as a named list: data (data.frame), time, gname.
data_list_PUT <- list(
  data  = as.data.frame(expr_PUT),  # G × N expression as data.frame
  time  = time_PUT,                  # length-N Zeitgeber time vector (hours)
  gname = rownames(expr_PUT)         # length-G gene identifiers
)

# Initialize the chain using a fast OLS cosinor fit for starting values.
# This prevents the chain from starting too far from the likely region.
init_PUT <- CB_init_single(Data.list = data_list_PUT, P = 24)

# Run the RJMCMC sampler.
# iteration: total post-thinning iterations stored per gene (K = iteration / thin)
# n.burn:    warm-up iterations discarded before storing (improves stationarity)
# thin:      thinning interval (default = 20; every 20th sample is kept)
mcmc_PUT <- CB_MCMC_single_rj_slice(
  Data.list  = data_list_PUT,  # data list from above
  Init.value = init_PUT,       # initialization from CB_init_single
  P          = 24,             # period in hours
  iteration  = 3000,           # 3000 thinned iterations stored (K = iteration)
  n.burn     = 1000,           # 1000 warm-up iterations discarded
  thin       = 20,             # keep 1 sample every 20 raw steps
  # p_rhythmic: prior Pr(rho=1) per gene — the prior probability that each gene
  # oscillates. Used in the RJMCMC log prior odds: log(p/(1-p)).
  # rep(0.2, G): uniform 20% prior rhythmicity (conservative — appropriate for most datasets).
  # Can be gene-specific (e.g. derived from a cosinor pre-screen p-value).
  p_rhythmic = rep(0.2, nrow(data_list_PUT$data))
)

# ── Condition B (Substantia Nigra) ───────────────────────────────────────────
data_list_SUN <- list(
  data  = as.data.frame(expr_SUN),
  time  = time_SUN,
  gname = rownames(expr_SUN)
)
init_SUN <- CB_init_single(Data.list = data_list_SUN, P = 24)

mcmc_SUN <- CB_MCMC_single_rj_slice(
  Data.list  = data_list_SUN,
  Init.value = init_SUN,
  P          = 24,
  iteration  = 3000,
  n.burn     = 1000,
  p_rhythmic = rep(0.2, nrow(data_list_SUN$data))
)

# Save immediately — MCMC is the most expensive computation in the pipeline
saveRDS(mcmc_PUT, "mcmc_PUT.rds")
saveRDS(mcmc_SUN, "mcmc_SUN.rds")
```

### Understanding the output

Each `mcmc_*` object is a list with four G × K matrices of posterior samples,
where G = number of genes and K = iteration (stored iterations after warm-up):

| Field | Dimensions | What each number represents |
|---|---|---|
| `mcmc_PUT$rho` | G × K | ρ indicator: 1 = rhythmic at iteration k, 0 = flat |
| `mcmc_PUT$phi` | G × K | Peak time (acrophase) of gene g at iteration k, in hours |
| `mcmc_PUT$A` | G × K | Amplitude of gene g at iteration k |
| `mcmc_PUT$M` | G × K | Midline (MESOR) expression of gene g at iteration k |

For a confidently rhythmic gene, `mcmc_PUT$rho[g, ]` will be a vector of nearly all
1s. For a clearly flat gene, nearly all 0s. The mixture across iterations is the
raw material for every downstream inference step.

> **On run time:** At 5,000 iterations for ~15,000 genes, expect 20–60 minutes per
> condition on a single core. Increasing to 10,000 iterations improves convergence
> for genes near the rhythmicity boundary, at the cost of proportionally more
> compute time. The `parallel` package can distribute genes across cores; see
> `?CB_MCMC_single_rj_slice` for options.

---

## 3.5. Step 1.5: Annotate the MCMC Output (Required)

> **This step is mandatory.** Raw MCMC output has gene identifiers only as
> `rownames`. All downstream BayRC functions expect the MCMC object to carry
> additional attributes — gene symbols, a binary rhythmicity index, and (for
> cross-species analyses) Ensembl IDs. Skipping this step will cause silent
> failures in `detect_rhy()`, `phase_infer()`, pathway enrichment, and concordance
> functions that rely on gene name matching.

### What annotation does

`match_symbols()` enriches every matrix in the MCMC list with three pieces of metadata:

| Attribute | Where stored | What it is |
|---|---|---|
| Gene symbols | `rownames(mcmc$rho)`, `rownames(mcmc$phi)`, etc. | Canonical gene symbol (e.g. `"BMAL1"`) used for matching across conditions |
| `attr(mcmc$rho, "symbols")` | Named vector of length G | Same symbols, used by `detect_rhy()` and visualization functions |
| `attr(mcmc$rho, "RHYindex")` | Integer vector (0/1) | Binary rhythmicity call at the chosen Bayes Factor threshold |
| `attr(mcmc$rho, "ensembl_gene_ids")` | Character vector (only when input was Ensembl IDs) | Required for cross-species matching via `match_homologs()` |

It also resolves **duplicate gene symbols**: when multiple Ensembl IDs map to the same gene symbol, the one with the highest Bayes Factor is kept.

### Same-species workflow (gene symbols already as rownames)

When your expression matrix row names are already HGNC/gene symbols (the common case for human/mouse data from processed databases):

```r
library(BayRC)

# After running MCMC: annotate both conditions.
# BF = 3 is the Bayes Factor threshold for the binary rhythmicity call stored
# in attr(mcmc$rho, "RHYindex"). This does NOT filter genes — all G genes are kept.
# p_rhythmic = 0.2 must match the prior used in CB_MCMC_single_rj_slice (default 0.2).
mcmc_PUT <- match_symbols(mcmc_PUT, BF = 3)
mcmc_SUN <- match_symbols(mcmc_SUN, BF = 3)

# What changed:
# - rownames(mcmc_PUT$rho)  now confirmed = gene symbols (e.g., "BMAL1", "CLOCK")
# - rownames(mcmc_PUT$phi)  now = gene symbols
# - rownames(mcmc_PUT$A)    now = gene symbols
# - attr(mcmc_PUT$rho, "symbols")   = character vector of gene symbols
# - attr(mcmc_PUT$rho, "RHYindex") = integer vector, 1 = BF > 3, 0 otherwise
```

### Same-species workflow (Ensembl IDs as rownames)

When your data uses Ensembl IDs (e.g., `"ENSG00000134057"`), `match_symbols()` 
calls `biomaRt` to convert them and also stores the original IDs for cross-species
matching:

```r
library(biomaRt)
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

mcmc_PUT <- match_symbols(mcmc_PUT, BF = 3, ensemble = ensembl_mart)
mcmc_SUN <- match_symbols(mcmc_SUN, BF = 3, ensemble = ensembl_mart)

# After this step:
# - rownames(mcmc_PUT$rho) = HGNC gene symbols
# - attr(mcmc_PUT$rho, "ensembl_gene_ids") = original Ensembl IDs (used by match_homologs)
# - Duplicate symbols resolved by keeping highest-BF Ensembl entry
```

> **Note:** `match_symbols()` detects whether rownames look like Ensembl IDs (start with
> `ENSG`/`ENSMUSG`, long strings) or gene symbols (short, mostly letters). If your
> data uses Ensembl IDs but you do not pass `ensemble`, the function will stop with
> an informative error.

### Cross-species workflow (human vs. baboon, human vs. mouse, etc.)

When comparing conditions from different species, gene identifiers must be mapped
to a common reference space before any joint analysis. BayRC requires:

1. Both MCMC objects have been through `match_symbols()` with Ensembl IDs input
   (so `attr(rho, "ensembl_gene_ids")` is populated)
2. `match_homologs()` maps each non-reference species to the reference space and
   computes the intersection of 1:1 orthologous genes

```r
# Example: comparing human Putamen (condition A) vs. baboon Putamen (condition B)
# Both MCMC objects must already be annotated with match_symbols().

# match_homologs() arguments:
#   input_dfs   — named list of annotated MCMC objects
#   species_from — character vector, same length, indicating each dataset's species
#   ref          — reference species for gene ID space (default: "human")
#
# Supported species: "human" (homo_sapiens), "mouse" (mus_musculus)
# Additional species can be added to the internal species_map in match_homologs().

aligned <- match_homologs(
  input_dfs    = list(human = mcmc_PUT, baboon = mcmc_SUN),
  species_from = c("human", "baboon"),
  ref          = "human"
)

# aligned is a list of two MCMC objects with matching rownames in human gene space.
# Only genes with 1:1 orthologs found in both species are retained.
mcmc_PUT_aligned <- aligned[[1]]   # human, subset to shared orthologs
mcmc_SUN_aligned <- aligned[[2]]   # baboon, re-indexed to human gene IDs

# Verify alignment:
identical(rownames(mcmc_PUT_aligned$rho), rownames(mcmc_SUN_aligned$rho))  # must be TRUE
cat("Shared orthologs:", nrow(mcmc_PUT_aligned$rho), "\n")
```

> **Warning:** `match_homologs()` requires the `"ensembl_gene_ids"` attribute to be
> set on both datasets. This attribute is only created by `match_symbols()` when the
> input row names are Ensembl IDs. If you call `match_homologs()` on data whose row
> names are already gene symbols, it will stop with an error about missing Ensembl
> IDs — go back and re-run your MCMC pipeline from expression data with Ensembl IDs
> as row names, or use `merge_mcmc()` with an ortholog table (see `?merge_mcmc`)
> as an alternative cross-species matching approach.

### What breaks if you skip annotation

| Downstream function | What fails | Error symptom |
|---|---|---|
| `detect_rhy(dat1, dat2)` | Gene names in output are `Gene_1, Gene_2, ...` instead of symbols | Output has no biological meaning |
| `phase_infer(phi1, phi2, ...)` | `rownames(phi_matrix)` is NULL → phase results have no gene IDs | Results cannot be linked back to genes |
| `pathSelect(...)` | Gene-pathway overlap uses rownames; NULL names = zero overlap for all pathways | All p-values → 1.0 |
| `match_homologs(...)` | `attr(rho, "ensembl_gene_ids")` is NULL | Hard stop with error |

---

## 4. Step 2: Posterior Summaries

Once MCMC is complete, extract three complementary summaries of rhythmicity
evidence per gene.

### 4a. Posterior rhythmicity probability

The most direct summary is `rowMeans(mcmc$rho)` — the fraction of MCMC iterations
where ρ = 1, which estimates P(rhythmic | data) for each gene.

```r
# pA[g]: posterior probability that gene g oscillates in condition A
pA <- rowMeans(mcmc_PUT$rho)   # one value per gene, range [0, 1]
pB <- rowMeans(mcmc_SUN$rho)

# Visualize the distribution to understand the landscape of rhythmicity
hist(pA, breaks = 50,
     main = "Posterior P(rhythmic | data) — Putamen",
     xlab = "Posterior probability",
     col  = "steelblue")
abline(v = 0.5, lty = 2, col = "red")
# Most genes will cluster near 0 (flat) or near 1 (confidently rhythmic).
# A bimodal distribution indicates good separation between the two classes.
```

> **Interpretation:** A value of pA = 0.92 for a gene means that in 92% of MCMC
> samples, the data were better explained by the rhythmic model than the flat model.
> This is a direct probability statement about the gene. A gene with pA = 0.15
> almost certainly does not oscillate in this condition. Unlike a p-value, this
> probability has a literal biological meaning: it is your updated belief about
> rhythmicity after seeing the expression data.

### 4b. Bayes Factor for rhythmicity

The **Bayes Factor (BF)** quantifies how much the data updated your belief about
rhythmicity relative to a neutral starting point. `summarize_bay()` computes this
from the posterior samples of ρ and provides binary rhythmicity calls.

```r
# BF = (posterior odds of rhythmicity) / (prior odds of rhythmicity)
# p_rhythmic = 0.5 is a neutral prior: equal weight to rhythmic and flat models
# p_rhythmic must match the prior used in CB_MCMC_single_rj_slice (default 0.2).
# Using a different value here will produce miscalibrated Bayes Factors.
bf_PUT <- summarize_bay(mcmc_PUT$rho, BF = 3, p_rhythmic = 0.2)
bf_SUN <- summarize_bay(mcmc_SUN$rho, BF = 3, p_rhythmic = 0.2)

# Output fields:
# bf_PUT$BayesF      — per-gene Bayes Factor (numeric vector)
# bf_PUT$Rhythmicity — binary call: 1 if BayesF > threshold (here, BF > 3)
# bf_PUT$RowAverage  — posterior mean of rho; identical to pA computed above

cat("PUT rhythmic genes (BF > 3): ", sum(bf_PUT$Rhythmicity), "\n")
cat("PUT rhythmic genes (BF > 10):", sum(bf_PUT$BayesF > 10), "\n")
```

**Reading Bayes Factors on the Jeffreys scale:**

| BF | Interpretation |
|---|---|
| < 1 | Evidence *against* rhythmicity |
| 1–3 | Anecdotal or weak |
| 3–10 | Moderate — gene is likely rhythmic |
| 10–100 | Strong — gene is confidently rhythmic |
| > 100 | Very strong / decisive |

> **Interpretation:** BF > 3 is the standard reporting threshold in BayRC, matching
> accepted Bayesian evidence standards. The BF naturally accounts for data quantity
> and quality: a gene measured at only two time points near its peak will not achieve
> a high BF even if expression appears to rise and fall, because sparse data do not
> strongly constrain the model. This is a feature — it prevents overconfident calls
> from low-information genes.

### 4c. Amplitude and phase estimates with credible intervals

`CB_getAllEst()` extracts posterior point estimates and 95% credible intervals for
all circadian parameters. Use this when you want to report specific peak times,
amplitudes, or MESOR values for genes of interest.

```r
# CB_getAllEst: posterior estimates and credible intervals from stored MCMC
# burn: additional burn-in to discard from the stored chain (fine-tuning)
est_PUT <- CB_getAllEst(mcmc_PUT, burn = 100)
est_SUN <- CB_getAllEst(mcmc_SUN, burn = 100)

# est_PUT is a list; access individual parameter data frames by index:
#   est_PUT[[1]] — phi (acrophase / peak time): phi.Est, phi.Lower, phi.Upper
#   est_PUT[[2]] — A  (amplitude):              A.Est,   A.Lower,   A.Upper
#   est_PUT[[3]] — M  (MESOR):                  M.Est,   M.Lower,   M.Upper
#   est_PUT[[4]] — sigma (noise variance)

# Phase estimates for the most confident rhythmic genes:
phase_df <- est_PUT[[1]]
head(phase_df[order(bf_PUT$BayesF, decreasing = TRUE), ])
#         phi.Est phi.Lower phi.Upper RHYindex
# ARNTL      6.1       4.2       8.0        1
# PER1       2.8       1.1       5.5        1
# NR1D1     14.3      11.2      17.0        1
```

**Why is the phase credible interval circular?** Phase wraps around the 24-hour
clock: a gene peaking at ZT23 and one peaking at ZT1 are only 2 hours apart, not
22. `CB_getAllEst()` uses a **circular Highest Density Interval (cHDI)** algorithm
that finds the shortest arc on the 24-hour circle containing 95% of posterior phase
samples. This correctly handles the ZT 0/24 boundary where a standard linear
interval would give nonsensical widths.

**What does HDI width tell you?** A narrow interval (e.g., ZT 6.1 with bounds
[4.2, 8.0]) means the data tightly constrain the peak time. A wide interval
spanning 8+ hours means peak timing is uncertain — possibly because the gene is
only weakly rhythmic, or because few samples were collected near its actual peak.

---

## 5. Step 3: Transition Classification

### The biological question

With rhythmicity probabilities in both conditions in hand, the central comparative
question is: **which genes changed their oscillatory status?** A gene rhythmic in
condition A but flat in condition B has *lost* circadian regulation. One that was
flat in A but rhythmic in B has *gained* it. A gene oscillating in both has
*maintained* rhythmicity and is the subject of phase concordance analysis (Step 4).

### The statistical approach: joint transition posteriors with BFDR

Rather than first applying a threshold to call each gene rhythmic or not — which
collapses all uncertainty into a binary decision — BayRC computes **joint posterior
transition probabilities** directly from pA and pB:

- **P(gain)**       = (1 − pA) × pB — probability flat in A, rhythmic in B
- **P(loss)**       = pA × (1 − pB) — probability rhythmic in A, flat in B
- **P(maintained)** = pA × pB — probability rhythmic in both

A gene with pA = 0.6 and pB = 0.8 has P(maintained) = 0.48 — moderate evidence
for conservation. This propagates the uncertainty in both conditions, rather than
making a hard call at each step.

**BFDR control** then determines a threshold τ for each transition type such that
the expected fraction of false discoveries among all classified genes does not
exceed α. This is the Bayesian analog of the classical FDR.

```r
# transition_classify: BFDR-controlled transition classification
# bfdr_alpha = 0.25 is appropriate for exploratory discovery;
# use 0.05 or 0.10 for confirmatory analyses
trans <- transition_classify(pA, pB, bfdr_alpha = 0.25)

# The function prints its thresholds and counts to the console:
# === Transition-level BFDR results ===
# τ_gain = 0.142 | n_gain = 312
# τ_loss = 0.118 | n_loss = 523
# τ_cons = 0.198 | n_cons = 677

# Summary of transitions
cat("Rhythmicity gains in B:   ", trans$n_gain, "\n")
cat("Rhythmicity losses in B:  ", trans$n_loss, "\n")
cat("Conserved in both:        ", trans$n_cons, "\n")

# The master classification vector — one label per gene
table(trans$gain_loss_status)
# Non-rhythmic     Gain     Loss  Maintained
#        11934      312      523         677

# Inspect gene-level details in the results data frame
head(trans$results)
# gene      pA     pB   p_gain  p_loss  p_cons  classification
# ARNTL   0.97   0.95     0.05    0.05    0.92  Maintained
# PER1    0.93   0.89     0.06    0.10    0.83  Maintained
# GENE5   0.03   0.76     0.74    0.02    0.02  Gain
```

> **Interpretation:** The `gain_loss_status` vector is the central annotation
> object for all remaining analysis steps. Every downstream step — phase
> concordance, pathway enrichment, heatmap construction — uses these labels. A gene
> receives only one label; genes that cross thresholds for multiple transition types
> are resolved by the BFDR procedure.

**Biological meaning of each class:**

| Label | What happened biologically |
|---|---|
| **Gain** | Arrhythmic in condition A, oscillates in condition B. Circadian regulation was *acquired* — a new clock-output pathway activated, or chromatin opened at a clock-response element. |
| **Loss** | Oscillated in condition A, flat in condition B. Circadian regulation was *lost* — common in disease states, aging, or chronic stress that dampens clock-driven transcription. |
| **Maintained** | Oscillates in both conditions. This is the substrate for phase concordance analysis. |
| **Non-rhythmic** | No evidence of rhythmicity in either condition at the chosen BFDR level. |

---

## 6. Step 4: Phase Concordance

### The biological question

Knowing that a gene oscillates in both conditions (the "Maintained" class) is only
half the story. **Does it peak at the same time of day in both?** A gene peaking at
ZT6 in condition A and ZT14 in condition B oscillates in both, but the timing has
shifted by 8 hours — a fundamentally different biological program than one peaking
at ZT6 in both. Phase concordance analysis asks: for genes that kept oscillating,
did they keep the same peak time?

### What is Δφ and why is it circular?

The phase difference Δφ = φ_A − φ_B is the signed difference in peak timing
between the two conditions, computed on the circle using each gene's full MCMC
posterior samples of phase from both conditions. The uncertainty in peak timing
from both chains propagates directly into the uncertainty about Δφ.

Because time is circular, Δφ is interpreted on the circle [−12, +12) hours: a
difference of +23 hours is equivalent to −1 hour (the gene peaks 1 hour earlier in
condition A). The posterior distribution of Δφ lives on this circle.

- **P(shift)** = fraction of posterior samples of |Δφ| ≥ δ (the tolerance in hours)
- **P(conserved)** = 1 − P(shift)

BFDR is applied to these probabilities to make formal, multiple-testing-controlled
classifications.

### Running phase inference

```r
# phase_infer: classifies Maintained genes as phase-conserved or phase-shifted.
# Only genes labeled "Maintained" in gain_loss_status are analyzed;
# all others receive NA in the output vectors.
#
# shift: the tolerance δ in hours
#   δ = 2 h  — tight definition (appropriate for densely sampled, low-noise data)
#   δ = 4 h  — permissive (appropriate for noisier or more sparsely sampled data)
# compute_hdi = TRUE: also computes 95% circular HDI of Δφ (slower but gives CIs)

phase <- phase_infer(
  phi_matrix1      = mcmc_PUT$phi,            # G × K phase matrix, condition A
  phi_matrix2      = mcmc_SUN$phi,            # G × K phase matrix, condition B
  gain_loss_status = trans$gain_loss_status,  # from transition_classify()
  shift            = 2,                       # δ = 2 hours (paper primary threshold)
  P                = 24,
  bfdr_alpha       = 0.25,
  compute_hdi      = TRUE
)

# Console output (example):
# === PHASE INFERENCE SUMMARY ===
# Maintained genes: 677
# BFDR α = 0.25  | shift threshold = 2 h
# Significant phase-shifted genes: 39
# Significant phase-conserved genes: 335
# Undetermined genes: 303
# HDI computation: enabled

# Per-gene output (all vectors of length G; NA for non-Maintained genes):
# phase$flag_cons         — TRUE if gene is phase-conserved (BFDR controlled)
# phase$flag_shift        — TRUE if gene is phase-shifted (BFDR controlled)
# phase$flag_undetermined — TRUE if insufficient evidence for either class
# phase$deltaPhi.Est      — posterior median of Δφ in hours (circular)
# phase$deltaPhi.Lower    — lower bound of 95% circular HDI of Δφ
# phase$deltaPhi.Upper    — upper bound of 95% circular HDI of Δφ
# phase$peak1             — posterior median peak hour in condition A
# phase$peak2             — posterior median peak hour in condition B

cat("Phase-conserved:", sum(phase$flag_cons,  na.rm = TRUE), "genes\n")
cat("Phase-shifted:  ", sum(phase$flag_shift, na.rm = TRUE), "genes\n")
```

> **Interpretation:**
> - `flag_cons = TRUE`: the posterior distribution of Δφ is strongly concentrated
>   near zero — the peak time agrees between conditions within ±δ hours. The gene
>   oscillates and peaks at the same time of day in both conditions.
> - `flag_shift = TRUE`: posterior mass is concentrated away from zero — there is a
>   genuine shift in peak timing. Even though the gene is rhythmic in both
>   conditions, the timing has been reset.
> - **Reading the HDI of Δφ:** An interval entirely above zero, e.g.,
>   [deltaPhi.Lower = 3.1, deltaPhi.Upper = 7.4], means the gene consistently peaks
>   3–7 hours earlier in condition A than in condition B. An interval straddling
>   zero, e.g., [−1.8, 2.3], is consistent with no timing change.

---

## 7. Step 5A: Pathway-Level Analysis (Two-Stage)

### Why two stages?

Testing all pathways for gain enrichment, all for loss enrichment, and all for
conservation would include pathways with zero circadian activity in either
condition. A pathway where every gene has pA ≈ pB ≈ 0.02 contributes only noise
to a transition enrichment test. The two-stage design prevents this:

- **Stage 1 (Union test):** Identifies pathways containing *any* rhythmically
  active genes in either condition. Only these "circadian-active" pathways proceed.
- **Stage 2 (Transition tests):** Within active pathways, tests whether genes are
  enriched for Gain, Loss, or Conservation using gene set enrichment analysis
  (FGSEA) with posterior-derived gene scores.

### Stage 1: identifying circadian-active pathways

The union score for a gene is −2 × log[(1 − pA)(1 − pB)] — a Fisher-like
combination statistic that is large whenever a gene has high rhythmicity
probability in *at least one* condition. Pathways with genes scoring highly on
this metric are enriched for circadian activity in the combined sense.

```r
# Load pathway annotations (KEGG gene sets as a named list of gene symbols)
kegg_pathway_list <- readRDS("kegg_pathway_list_hsa.rds")

# Stage 1: union test — which pathways contain circadian-active genes?
result_union <- pathSelect(
  mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
  pathway.list    = kegg_pathway_list,
  dataset.names   = c("PUT", "SUN"),
  ranking.method  = "union",    # gene score = -2*log[(1-pA)(1-pB)]
  score_type      = "pos",      # one-sided: enrichment for high scores only
  qvalue.cut      = 0.20,
  nperm           = 1000,
  seed            = 42
)

# Filter to pathways that pass Stage 1
active_names        <- result_union$results$pathway[result_union$results$padj < 0.20]
active_pathway_list <- kegg_pathway_list[active_names]
cat(length(active_names), "pathways pass Stage 1 (circadian-active)\n")
```

### Stage 2: transition enrichment within active pathways

```r
# Stage 2a: which active pathways are enriched for rhythm GAIN in condition B?
result_gain <- pathSelect(
  mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
  pathway.list    = active_pathway_list,   # only Stage-1 active pathways
  dataset.names   = c("PUT", "SUN"),
  ranking.method  = "gain",      # gene score = P(gain) = pB * (1 - pA)
  score_type      = "pos",
  qvalue.cut      = 0.20,
  nperm           = 1000
)

# Stage 2b: which are enriched for rhythm LOSS?
result_loss <- pathSelect(
  mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
  pathway.list    = active_pathway_list,
  dataset.names   = c("PUT", "SUN"),
  ranking.method  = "loss",      # gene score = P(loss) = pA * (1 - pB)
  score_type      = "pos",
  qvalue.cut      = 0.20,
  nperm           = 1000
)

# Stage 2c: which are enriched for rhythm CONSERVATION?
result_cons <- pathSelect(
  mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
  pathway.list    = active_pathway_list,
  dataset.names   = c("PUT", "SUN"),
  ranking.method  = "conserved", # gene score = P(conserved) = pA * pB
  score_type      = "pos",
  qvalue.cut      = 0.20,
  nperm           = 1000
)

# Each result$results data frame contains (selected key columns):
# pathway                     — pathway name
# pval                        — permutation p-value from fgsea
# padj                        — BH-adjusted q-value
# NES                         — normalized enrichment score
# Gain_Index                  — expected proportion of pathway genes that are gains
# Loss_Index                  — expected proportion that are losses
# Conserved_Index             — expected proportion that are maintained
# Expected_N_Gain             — expected count of gain genes in the pathway
# Expected_N_Loss             — expected count of loss genes
# Expected_N_Conserved        — expected count of conserved genes
# Gain_Loss_Ratio_Arithmetic  — expected gains / expected losses (the GLR)
# Top_Gain_Genes              — top 5 candidates by P(gain), comma-separated
# Top_Loss_Genes              — top 5 candidates by P(loss)
# Top_Conserved_Genes         — top 5 candidates by P(conserved)

sig_gain <- result_gain$results[result_gain$results$padj < 0.20, ]
sig_loss <- result_loss$results[result_loss$results$padj < 0.20, ]
sig_cons <- result_cons$results[result_cons$results$padj < 0.20, ]

cat(nrow(sig_gain), "pathways enriched for gain\n")
cat(nrow(sig_loss), "pathways enriched for loss\n")
cat(nrow(sig_cons), "pathways enriched for conservation\n")
```

### Interpreting the GLR (Gain-Loss Ratio)

The **Gain-Loss Ratio** (GLR = Expected_N_Gain / Expected_N_Loss) for a pathway
summarizes whether circadian remodeling in that pathway is predominantly activating
or suppressing. Because BayRC uses continuous posterior probabilities rather than
hard calls, these expected counts are sums of P(gain) and P(loss) across all
pathway genes.

- **GLR > 1:** More gains than losses in this pathway — the circadian clock has
  *recruited* this biological process in condition B.
- **GLR < 1:** More losses than gains — this pathway has been *uncoupled* from
  circadian control in condition B.
- **GLR ≈ 1:** Gains and losses are balanced — the pathway may be reorganized
  without net change in total rhythmicity.

> **Interpretation example:** A result showing the "Oxidative phosphorylation"
> pathway with NES = 2.4, padj = 0.003, and GLR = 3.9 means: this pathway is
> significantly enriched for rhythmicity in condition B, with roughly 4× more gain
> genes than loss genes. Mitochondrial metabolism has come under stronger circadian
> control in this condition — a finding with direct implications for energy
> homeostasis and time-of-day metabolic capacity.

---

## 8. Step 5B: Genome-Wide Concordance

### The biological question

Beyond individual genes and pathways, you may want a single-number summary: **how
similar are the two circadian transcriptomes overall?** The **c-score** (adjusted
Jaccard concordance, centered at zero) provides this.

### What is the c-score?

The raw Jaccard index measures overlap as: (shared rhythmic genes) / (union of
rhythmic genes). Under the null hypothesis of no biological conservation, some
overlap occurs by chance — if both conditions have 15% of genes as rhythmic, random
overlap gives a Jaccard of roughly 0.08. The c-score subtracts this chance
expectation (estimated by gene-label permutation) and rescales so that:

- **c = 0:** concordance no better than random
- **c > 0:** more rhythmic gene overlap than expected by chance
- **c < 0:** less overlap than expected (rare; suggests anti-concordant programs)

A permutation p-value tests whether the observed c-score exceeds the null
distribution. A bootstrap confidence interval quantifies uncertainty in the
c-score estimate.

### Running genome-wide concordance

```r
# multi_conservation with select.pathway.list = "global" runs across all genes.
# This is the genome-wide concordance analysis, not pathway-level.
global <- multi_conservation(
  mcmc.merge.list     = list(A = mcmc_PUT, B = mcmc_SUN),
  dataset.names       = c("PUT", "SUN"),
  select.pathway.list = "global",   # genome-wide (not pathway-specific)
  n_perm              = 1000,       # permutation iterations for p-value
  n_boot              = 500,        # bootstrap iterations for CI
  output.dir          = "results/PUT_SUN_concordance",
  use_cpp             = TRUE,       # C++ acceleration (strongly recommended)
  compute_pvalue      = TRUE,
  compute_ci          = TRUE
)

# Key output: results are stored in matrices indexed by pathway name and pair name.
# For global analysis there is one row ("Global") and one column ("PUT&SUN").
cscore   <- global$Adjusted["Global", "PUT&SUN"]
pvalue   <- global$PValue["Global", "PUT&SUN"]
ci_lower <- global$CI_lower_adj["Global", "PUT&SUN"]
ci_upper <- global$CI_upper_adj["Global", "PUT&SUN"]
glr      <- global$Gain_loss_ratio["Global", "PUT&SUN"]

cat("Adjusted c-score:", round(cscore, 3), "\n")
cat("Permutation p-value:", round(pvalue, 4), "\n")
cat("95% CI: [", round(ci_lower, 3), ",", round(ci_upper, 3), "]\n")
cat("Genome-wide GLR:", round(glr, 3), "\n")
```

> **Interpretation:** A c-score of 0.42 (p < 0.001) means the two circadian
> transcriptomes share 42% more rhythmic gene overlap than expected by chance. A
> c-score near 0 means the two conditions have largely non-overlapping circadian
> programs — strong evidence for circadian remodeling. A negative c-score can occur
> when the two conditions have very different rhythmic gene counts and the smaller
> set falls largely outside the larger one.

**When is the c-score most meaningful?** The c-score is most informative when both
conditions have substantial rhythmicity (> 300–500 rhythmic genes each). For very
sparse circadian programs, permutation variance is high and confidence intervals
will be wide.

---

## 9. Visualization — The BayRC Pathway Heatmap

A central deliverable of BayRC is a **five-panel pathway heatmap** (Figure 5 in
the manuscript) that integrates every layer of circadian evidence for every gene in
a pathway on a single page. Rather than a separate table of gains, a separate phase
plot, and a separate probability bar chart, this heatmap places all information in
a single aligned visualization. Each gene row tells a complete biological story
read from left to right.

### The five-panel design

```
┌─────────────────────┬──────────────┬────────────────────┬─────────────────┬─────────────────┐
│  Rhythmicity Status │ Phase Status │  P(ρ=1 | data)     │ Phase posterior │ Phase posterior │
│  (transition type)  │ (timing)     │  Condition A | B   │ Condition A     │ Condition B     │
└─────────────────────┴──────────────┴────────────────────┴─────────────────┴─────────────────┘
```

**Panel 1 — Rhythmicity Status (leftmost annotation bar):**  
A discrete color annotation keyed to `trans$gain_loss_status`.

- **Orange** = Maintained (rhythmic in both conditions)
- **Purple** = Gain (newly rhythmic in condition B)
- **Blue** = Loss (lost rhythmicity in condition B)

This panel answers: *Did this gene's oscillatory status change between conditions?*

**Panel 2 — Phase Status (second annotation bar):**  
Keyed to `phase$flag_cons` and `phase$flag_shift`. Only meaningful for Maintained
genes; Gain and Loss genes receive a neutral color here.

- **Green** = Phase-conserved (peak timing held within ±δ hours)
- **Red** = Phase-shifted (timing changed significantly)
- **Gray** = Undetermined (insufficient evidence to classify)

This panel answers: *For genes that kept oscillating, did they keep the same peak
time?*

**Panel 3 — Rhythmicity probability heat columns:**  
Two side-by-side columns showing pA (condition A) and pB (condition B) for each
gene, encoded as a white-to-deep-red gradient (white = P = 0; deep red = P = 1).

- A **Gain** gene row: white in pA column, deep red in pB — rhythmic only in B
- A **Loss** gene row: deep red in pA column, white in pB — rhythmic only in A
- A **Maintained** gene row: deep red in both — confidently rhythmic in both

This panel answers: *How confident are we in each gene's rhythmicity call, in each
condition independently?*

**Panel 4 — Phase posterior distribution, Condition A:**  
A bar spanning ZT −6 to ZT +18 showing the MCMC posterior density of peak time for
each gene in condition A. A sharp, tall bar at one ZT point means tightly
constrained peak timing. A flat, spread bar means high phase uncertainty (common
for genes with low pA). Genes with very low pA are shown at reduced intensity to
avoid overinterpreting uncertain phase estimates.

**Panel 5 — Phase posterior distribution, Condition B:**  
The same representation for condition B. Comparing Panels 4 and 5 side by side
within a single gene row reveals directly whether the peak distribution shifted —
two sharp bars at the same ZT means phase-conserved; two sharp bars at different
ZTs means phase-shifted.

### Reading a gene row: two worked examples

> **Example 1 — A robustly maintained clock gene (ARNTL / BMAL1):**  
> Reading left to right:
> - **Orange Rhythmicity Status** → maintained oscillator in both brain regions
> - **Green Phase Status** → peak timing conserved within ±2 h
> - **Deep red in both pA and pB columns** → pA ≈ 0.97, pB ≈ 0.95: highly
>   confident rhythmicity in both conditions
> - **Sharp blue bar at ZT6 in Panel 4** → ARNTL peaks around ZT6 in the putamen,
>   low phase uncertainty
> - **Sharp blue bar at ZT6 in Panel 5** → same peak in the substantia nigra
>
> Conclusion: ARNTL is a robustly maintained circadian oscillator with conserved
> ZT6 peak timing across both brain regions. Core clock genes should show exactly
> this pattern, making ARNTL a useful positive control for the entire analysis.

> **Example 2 — A maintained but phase-shifted gene (CRY1):**  
> Reading left to right:
> - **Orange Rhythmicity Status** → maintained oscillator in both conditions
> - **Red Phase Status** → peak timing has shifted significantly
> - **Deep red in both pA and pB columns** → CRY1 is confidently rhythmic in both
>   (pA ≈ 0.93, pB ≈ 0.91) — the rhythmicity itself is not in question
> - **Sharp blue bar at ZT14 in Panel 4** → CRY1 peaks mid-night in the putamen
> - **Sharp blue bar at ZT22 in Panel 5** → CRY1 peaks late-night in the
>   substantia nigra, ~8 hours later
>
> Conclusion: CRY1 oscillates confidently in both brain regions, but its peak has
> shifted by ~8 hours. The timing of this circadian repressor has been reset between
> conditions despite maintained oscillation — a biologically meaningful dissociation
> that a simple gain/loss analysis would completely miss. The five-panel heatmap
> makes this visible at a glance without separate tables or plots.

---

## 10. Interpreting Your Results — A Biological Checklist

After completing a BayRC analysis, work through these questions to connect numbers
to biological conclusions.

- [ ] **What fraction of genes are rhythmic in each condition?**  
  Compare the distributions of `pA` and `pB`. How many genes have pA > 0.5? Use
  `detect_rhy()` for BFDR-controlled rhythmic gene counts per condition separately.
  A large difference (e.g., 25% rhythmic in A vs. 8% in B) suggests broad
  circadian dampening in condition B — the molecular clock is still running, but it
  has uncoupled from most of its transcriptional targets.

- [ ] **What is the dominant transition type?**  
  Compare `trans$n_gain` vs. `trans$n_loss` vs. `trans$n_cons`.  
  Gain-dominant (n_gain >> n_loss): condition B reflects circadian *activation* —
  the clock has recruited new gene targets.  
  Loss-dominant (n_loss >> n_gain): condition B reflects circadian *dampening* —
  clock output has weakened. This is a common signature of aging, neurodegeneration,
  and metabolic stress.  
  Balanced gain/loss with large n_cons: the circadian core is preserved but the
  peripheral output has been broadly reorganized.

- [ ] **Are phase shifts systematic or gene-specific?**  
  Examine the distribution of `phase$deltaPhi.Est` across phase-shifted genes. If
  values cluster around the same Δφ (e.g., most between +3 and +5 h), the entire
  clock output has advanced or delayed uniformly — a systematic timing change
  possibly reflecting altered feeding schedules, light exposure, or clock gene
  levels. If Δφ values scatter across all magnitudes and signs, phase remodeling is
  gene-specific, possibly driven by altered transcription factor binding patterns.

- [ ] **Which pathways are remodeled? What do GLR values tell you?**  
  Cross-reference Stage 2 results: which pathways appear in the gain list but not
  the loss list? Inspect `Gain_Loss_Ratio_Arithmetic` for each significant pathway.  
  High GLR (> 3) in a metabolic pathway: circadian clock newly controls this
  metabolism in condition B.  
  Low GLR (< 0.5) in an immune pathway: immune rhythmicity has been lost — relevant
  to time-of-day infection susceptibility and immune function.  
  Check `Top_Gain_Genes` and `Top_Loss_Genes` to identify specific driver genes
  within each enriched pathway for follow-up validation.

- [ ] **What is the global concordance (c-score)? Is it significantly different from
  chance?**  
  Check `global$Adjusted` and `global$PValue`.  
  c-score > 0.3, p < 0.01: strong transcriptome-wide circadian conservation.  
  c-score ≈ 0.1–0.2, p < 0.05: modest but significant conservation — a meaningful
  subset of the rhythmic program is shared.  
  c-score ≈ 0, p > 0.1: the two circadian programs are largely independent.  
  Also examine the genome-wide GLR: does condition B have systematically more gains
  or losses across the entire transcriptome? This single number summarizes the net
  direction of the global circadian shift.

---

## 11. Function Quick Reference

| Step | Function | Key parameters | Output |
|---|---|---|---|
| Initialize chain | `CB_init_single(Data.list, P=24)` | `Data.list=list(data, time, gname)`; `data` is `as.data.frame(Y)` | Init object |
| Run MCMC | `CB_MCMC_single_rj_slice(Data.list, Init.value, P=24, iteration=3000, n.burn=1000)` | `iteration=3000`, `n.burn=1000`, `thin=20` typical | List of G×K posterior sample matrices |
| Posterior P(rhythmic) | `rowMeans(mcmc$rho)` | — | Numeric vector of length G |
| Bayes Factor | `summarize_bay(mcmc$rho, BF, p_rhythmic)` | `BF=3`, `p_rhythmic=0.5` | BayesF, Rhythmicity, RowAverage |
| All estimates + HDI | `CB_getAllEst(mcmc, burn)` | `burn=100`, `credMass=0.95` | List of data frames: phi, A, M, sigma |
| BFDR threshold | `bfdr_from_posterior(probs, alpha)` | `alpha=0.25` | threshold, rhythmic_genes, n_rhythmic |
| Per-condition detection | `detect_rhy(dat1, dat2, bfdr_alpha)` | `bfdr_alpha=0.25` | Rhythmic gene sets per condition |
| Transition classification | `transition_classify(pA, pB, bfdr_alpha)` | `bfdr_alpha=0.25` | gain_loss_status vector, counts, thresholds |
| Marginal classification | `transition_classify_marginal(pA, pB, bfdr_alpha)` | — | Alternative to joint classification |
| Phase inference | `phase_infer(phi1, phi2, gain_loss_status, shift, bfdr_alpha, compute_hdi)` | `shift=2` or `4`, `compute_hdi=TRUE` | flag_cons, flag_shift, deltaPhi.Est, HDI |
| Pathway enrichment | `pathSelect(mcmc.merge.list, pathway.list, ranking.method, nperm)` | `ranking.method`: "union", "gain", "loss", "conserved" | fgsea results + GLR + top genes per pathway |
| Genome-wide concordance | `multi_conservation(..., select.pathway.list="global")` | `n_perm=1000`, `n_boot=500` | c-score, p-value, CI, GLR |
| Multi-dataset concordance | `multi_conservation_global(mcmc.merge.list, dataset.names)` | `B=1000` | Pairwise c-score matrix |
| Ortholog matching | `match_homologs(input_dfs, species_from, ref)` | — | Matched MCMC objects for cross-species analysis |

---

## Citation

Pham T, Kauffman K. *BayRC: A Bayesian framework for comparative circadian
genomics with FDR control and concordance scoring.* Manuscript in preparation, 2025.
