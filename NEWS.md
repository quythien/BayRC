# BayRC 0.3.0 (2026-07-14)

Public-release preparation: a curated public API, real bundled example data,
and documentation rewritten around it.

## Exported API curated from 61 to 22 functions

Every export was checked against two questions: does it answer one of the six
biological questions in the README, and does the paper's own methodology
describe it (not just "does something call it")? Kept: the MCMC core,
single- and two-group biomarker detection with BFDR control, pathway
enrichment, genome-wide concordance, and cross-species alignment. Un-exported
(kept in `R/` as internal code, not deleted): `FoxWright`/`1`/`2`/`3`, an
11-function variance-validation cluster in `simulate_parallel.R`,
`CBt_MCMC_single`, `congruence()` and `concordance2()` (legacy diagnostics
predating the package's BFDR-based approach; `concordance2()` in particular
is the hard-threshold Venn-diagram method the paper's own introduction
argues against), `compute_adjusted_jaccard_analytical_pvalue()` (analytical
shortcut, not in the paper), `transition_classify_marginal()` (duplicates
`detect_rhy()`'s per-condition logic, not in the paper), and about a dozen
more per-iteration/permutation helpers with no real callers outside the
functions that already wrap them. `pairwise_concordance()`, initially cut
for apparent lack of callers, was restored after confirming real usage in an
untracked working script. `Cosinor_fit()`, initially cut, was restored: it
is the classical (non-Bayesian) cosinor method the paper's own comparisons
contrast BayRC's RJMCMC approach with, and its un-exported `@examples` was
calling itself directly, which is a real `R CMD check` ERROR for anyone
installing the package.

## Real bugs found and fixed

- `plot_heatmap()` (renamed from `plot_pathway_integrated()`) checked
  `requireNamespace("ComplexHeatmap")` but never attached it, so
  `rowAnnotation()`/`Heatmap()`/`colorRamp2()` failed for anyone who hadn't
  separately loaded `ComplexHeatmap`/`circlize`. Now attaches both; added to
  `Suggests` (were undeclared).
- `phase_infer()` returned `flag_cons`, `flag_shift`, `peak1`, `peak2`, and
  related per-gene vectors unnamed, unlike `transition_classify()`'s
  `gain_loss_status` which is named by gene symbol. Indexing the unnamed
  vectors by gene symbol silently returned `NA` instead of erroring. All
  per-gene outputs are now named.
- An unescaped `%*%` in `update_A_phi_slice()`'s roxygen `@param` line was
  corrupting that Rd file's rendering.
- `R/BayRC-package.R` added with `@importFrom`/`@useDynLib` roxygen tags.
  Without them, `devtools::document()` silently stripped NAMESPACE's
  `useDynLib`/`importFrom` lines on every regeneration.
- `pairwise_concordance()` called `concordance()`, a function that does not
  exist anywhere in the package or its pre-refactor snapshot (only a
  commented-out stub remained); it would have errored for any caller.
  Inlined the Jaccard computation the stub described.
- `transition_classify()` and `transition_classify_marginal()` silently
  accepted a raw posterior sample matrix (`mcmc$rho`) in place of the
  documented per-gene probability vector (`rowMeans(mcmc$rho)`). The
  arithmetic still ran, but building the returned `results` data.frame
  recycled the still-matrix `pA`/`pB` against a length-`G*K` vector,
  observed to balloon memory past 100 GB on a real 5,066-gene,
  2,001-iteration posterior before it was caught and killed. Both now
  validate input shape and fail fast with an actionable error instead.
- `match_homologs()` called `getHomologs()`, a function that has never
  existed anywhere in this codebase's history (confirmed against the
  pre-refactor `Thien/` snapshot); it would have errored for any caller
  doing cross-species alignment. Replaced with a real `biomaRt::getLDS()`
  ortholog query, reusing the existing `try_any_mirror()` mirror-fallback
  helper.
- `match_symbols()`, `try_any_mirror()`, and `match_homologs()` called
  `getBM()`/`useEnsembl()`/`getLDS()` unqualified; since `biomaRt` is a
  `Suggests` dependency loaded only via `requireNamespace()` (never
  attached), these would error unless a caller happened to have `biomaRt`
  separately attached. Now namespace-qualified (`biomaRt::getBM()`, etc.).
- Non-ASCII characters (curly quotes, em dashes, Greek letters, section
  signs) removed from `R/internal.R` and `R/plot_heatmap.R` source: ASCII
  substitutes in comments, `\uXXXX` escapes in the few `cat()`/`warning()`
  strings that print Greek letters to the user, so printed output is
  unchanged.

## Real bundled example data

- `inst/extdata/baboon_PUT_SUN_GSE98965.rds`: real baboon putamen and
  substantia nigra expression (Mure et al. 2018, GEO accession GSE98965),
  5,066 genes ortholog-matched to the same backbone used in the manuscript,
  12 real zeitgeber timepoints. README's Quick Start now runs this real
  PUT-vs-SUN comparison end to end through the current package code
  (mirroring the paper's Case Study 2), replacing an earlier simulated-data
  walkthrough.

## Documentation

- README rewritten: real Quick Start on bundled data; a "How BayRC
  Categorizes Every Gene" table with real counts at every BFDR-controlled
  stage; a real pathway heatmap figure (KEGG Parkinson disease, the pathway
  the manuscript's own SUN-PUT case study centers on); a circular-HDI
  demonstration figure; a plain-language glossary with citations; a
  References section.
- `BayRC-manual.pdf` rebuilt with a proper title page and a categorized,
  hyperlinked function index grouped by the paper's own section numbers,
  restricted to the 22 curated exports (not all ~74 documented topics).
  The index is further split into single-group vs two-group biomarker
  detection (paper Sec. 2.2), matching the two distinct function
  signatures. Its `\describe{\item{}}` structure rendered each category's
  first bullet inline with the heading instead of on its own line; now
  uses `\strong{}` + `\itemize{}`.
- `LICENSE.md` and `CITATION.cff` added. `Authors@R` updated to match the
  manuscript's current author list. Maintainer email updated to
  `qtp1@pitt.edu`.
- `.gitignore` fixed: `LICENSE.md`, `CITATION.cff`, and `BayRC-manual.pdf`
  were previously swallowed by the blanket ignore rule and would not have
  appeared on the public GitHub repo. `.Rbuildignore` added (was missing
  entirely), excluding dev-only files (`paper/`, `Thien/`, `REVIEW_*.md`)
  from the built package.
- `go.pathway.list_hsa.RData`, `kegg.pathway.list_cel_GeneNames.RData`, and
  `kegg_pathway_list_hsa.rds` moved from `data/` to `inst/extdata/`: none
  are referenced anywhere in `R/`, only by `inst/analysis/` scripts via
  hardcoded external paths, so they were undocumented lazy-load "public"
  data with no real caller, and the `.rds` file is not an allowed `data/`
  directory type under `R CMD check`.
- `vignettes/BayRC_workflow.Rmd` rewritten to be genuinely runnable (was
  entirely `eval = FALSE` with placeholder objects that don't exist,
  which is a real `R CMD check` ERROR since `knitr::purl()` extracts
  chunk code regardless of the `eval` option).
- Unused `Rfast`/`truncdist` removed from `DESCRIPTION` `Imports`.

# BayRC 0.2.0 (2026-04-23)

## Critical bug fixes

- **MCMC Markov property restored**: `set.seed()` was called inside the sampling
  loop on every iteration, resetting the RNG state and breaking the Markov
  property. Moved to a single call before burn-in. (#1)

- **RJMCMC acceptance ratio corrected**: `log.phi.prior` was assigned `1/P`
  (the density) instead of `log(1/P)` (the log-density). This biased every
  birth move by ~3.2 log-units, inflating posterior rhythmicity probabilities
  globally. (#2)

- **`try_save` re-throws on error**: Previously caught errors silently and
  returned `NULL`, allowing corrupted values to propagate through the chain.
  Now saves checkpoint then re-throws. (#3)

- **`thin` parameter now enforced**: Storage `cbind` calls were running every
  iteration regardless of `thin`, causing 20× memory overruns and storing
  unintended autocorrelated samples. (#4)

## Statistical API fixes

- `transition_classify()` now returns `tau_rhythmic_A`, `tau_rhythmic_B`
  (paper §2.2 condition-specific BFDR thresholds) and exposes `p_gain`,
  `p_loss`, `p_cons` at the top level of the return list. (#9, #10)

- `summarize_bay()` default `p_rhythmic` corrected from 0.5 to 0.2 to match
  the MCMC prior, preventing 4× miscalibration of reported Bayes Factors. (#18)

- `phase_infer()` default `shift` corrected from 4 to 2 hours, matching the
  paper's primary "±2 h" phase concordance claim. (#24)

- `congruence()` now warns when `delta` or `units` are passed (these parameters
  are unused by the threshold-free probabilistic c-score). (#7)

- GLR zero-denominator now returns `Inf` (mathematically correct per paper
  Eq. 4) instead of `1e10` throughout `congruence()` and `path_select.R`. (#17)

## New features

- `CB_MCMC_single_rj_slice()` gains a `diagnostics = FALSE` parameter.
  When `TRUE`, appends per-gene acceptance rates, ESS estimates, and posterior
  rhythmicity probabilities to the return list, and prints a summary. (#19)

- `multi_conservation_pathway_bootstrap()` added as a backward-compatible
  wrapper around `multi_conservation()` for analysis scripts that call it by
  the old name. (#14)

- `multi_conservation()` gains `save_output = TRUE` to suppress the
  unconditional Excel write in CI/test environments. (#27)

## Bug fixes

- `adjust.to.2pi(2*pi)` returned `NULL` (boundary case fell through all
  branches). Replaced with `((x %% (2*pi)) + 2*pi) %% (2*pi)`. (#5)

- `circular_median()` returned class `c("circular","numeric")` instead of
  plain `numeric`, causing unexpected dispatch downstream. (#11)

- `CBt_sim_data()` was completely broken: default argument referenced
  undefined `my.TOJR.grid`; `matrix(replicate(...))` failed on data.frame
  rows; `Params["A"]` returned a list instead of a scalar. All three fixed. (#8)

- `rand()` in `src/permutation_functions.cpp` replaced with `Rcpp::sample`
  to use R's seeded RNG. Eliminates platform-dependent non-reproducibility
  and modulo bias. (#16)

- `BHM_PRC.R`: `whPRC` typo fixed to `while`; `quantPRC()` replaced with
  `quantile()`; missing `}` closing brace fixed (trapped `source(heatmap.R)`
  inside a conditional block). (#15)

## Package infrastructure

- 17 missing packages added to `Imports:` in DESCRIPTION; `edgeR` and
  `DESeq2` removed (not used in `R/`). Passes `R CMD check` dependency scan. (#12)

- All `library()` / `require()` calls removed from exported function bodies
  (CRAN Policy §1.4 violation). NAMESPACE `importFrom` entries added for
  ggplot2, dplyr, tibble, purrr, RColorBrewer. `biomaRt` uses
  `requireNamespace()` guard. (#13)

- `R/permutation_test.R` (misplaced analysis script) removed from package
  source. (#28)

## Paper fixes

- Supplementary Section S7 added: MCMC convergence diagnostics table covering
  all datasets and tissues used in the three case studies. All chains confirmed
  converged (0% of genes with ESS < 100 in any tissue).

- Stage 1 pathway filter description corrected in paper §2.3: "FDR < 0.05"
  updated to "unadjusted p < 0.05 as a liberal pre-screen". (#23)

- Pathway count "12" on line 366 corrected to "14" (consistent with lines
  402 and 537). (#25)

## Reproducibility

- `inst/analysis/config.R` created: centralises all hardcoded paths behind
  environment variables (`BAYRC_DATA_DIR`, `BAYRC_WD_DIR`, etc.). (#21)

- `inst/analysis/README_figures.md` created: maps each paper figure (1–6) to
  its generating script, required data, and key parameters. (#22)

## Tests

- `tests/testthat/test-phase-infer.R` added: 8 known-answer tests for
  `phase_infer()` covering zero phase difference (all conserved), P/2
  difference (all shifted), empty maintained set, and threshold sensitivity. (#26)

# BayRC 0.1.0 (2025)

- Initial package release with RJMCMC engine, BFDR classification, circular
  phase inference, pathway enrichment, and genome-wide concordance scoring.
