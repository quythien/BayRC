# BayRC Paper Figures — Generating Scripts

This file maps each figure to the script that produces it and the data that
script needs. Everything lives under `inst/analysis/`.

**Paths.** Every script now begins by sourcing `config.R`, which reads its
directories from environment variables and falls back to the original
`/home/qtp1/...` locations. Nothing else in the scripts hardcodes a path. To
run elsewhere, set the variables before launching R:

```
export BAYRC_DATA_DIR=/path/to/Collaborative
export BAYRC_WD_DIR=/path/to/Circadian
export BAYRC_RESULT_DIR=/path/to/GTEXdata/result      # per-tissue MCMC output
export BAYRC_SUMMARY_DIR=$BAYRC_RESULT_DIR/summary/hb # rho/phi summaries
export BAYRC_OUTPUT_DIR=/path/to/analysis/output
export BAYRC_FIGURE_DIR=/path/to/figure/output
```

Scripts in `inst/analysis/` find `config.R` next to themselves; scripts in
`plots/` and `pipeline/` find it one directory up. Both work under `Rscript`
and when sourced from the directory the script lives in.

---

## Order of operations

The MCMC output is not in the repository, and neither are the two `.RData`
summaries every figure script loads. Build them in this order:

```
1. CAMO_h_b.R                          per-tissue MCMC, 26 human + 26 baboon
2. pipeline/summarize_rho_phi.R        -> mcmc_rho_BF3.RData
                                       -> phi/mcmc_phi_BF3.RData
3. the figure scripts below
4. assemble_figures.R                  panels -> Figure_2 ... Figure_6
```

Step 2 writes the file the figure scripts load. They all begin with
`load(.../mcmc_rho_BF3.RData)`, which `pipeline/summarize_rho_phi.R` produces
from the per-tissue chains; `--validate` compares its output against an
existing copy instead of overwriting it.

---

## Figure 1 — BayRC framework overview

Hand-drawn flowchart (`flowchart_resized.jpeg`). No generating script.

---

## Figure 2 — Genome-wide and circadian concordance heatmaps (26 baboon tissues)

**Scripts:**
- Genome-wide panel: `plots/heatmap_baboon.R` — writes
  `Baboon_Concordance_Heatmap_0.25_<method>.pdf` and
  `Baboon_Concordance_Heatmap_Dissim_<method>.pdf` to `BAYRC_FIGURE_DIR`,
  plus `Baboon_Concordance_Matrix.csv`.
- KEGG Circadian panel: `plots/heatmap_circadian_pairs.R within_baboon` —
  writes `*_Heatmap_<method>.pdf` under `BAYRC_OUTPUT_DIR/figure/`.

`heatmap_circadian_pairs.R` also takes `within_baboon_with_scn`,
`within_human` and `cross_species`.

**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData`.
**Pairs:** `choose(26, 2) = 325` within-species; cross-species is `26 x 26 = 676`.

> Earlier versions of this file attributed Figure 2 to
> `circa_concordance.plots.R`. That script does not draw the concordance
> heatmaps — it draws the phase-concordance scatters that make up Figures 3
> and 6, and it expects the objects from the case-study scripts to already be
> in the environment. It is a plotting companion, not a standalone entry point.

---

## Figure 3 — Within-species phase concordance scatter

**Scripts:**
- SCN-HIP panel: `Baboon_SCN_HIP.R` → `Baboon_SCN_HIP_Peak_Concordance_0.25_2h_new.pdf`
- SUN-PUT panel: `Baboon_SUN_PUT.R` → `Baboon_SUN_PUT_Peak_Concordance_0.25_2h_new.pdf`

**Key parameter:** `shift = 2` (the ±2h phase window, paper §2.2).
**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData`.

---

## Figure 4 — SUN-PUT pathway enrichment dotplots

**Scripts:** `Baboon_SUN_PUT.R` (stage-2 enrichment), then
`plot_enrich_SUN_PUT.R` for the final dotplots
(`SUN_PUT_shifted_GO_BP_dotplot.pdf`, `SUN_PUT_shifted_KEGG_dotplot.pdf`).
**Required data:** the stage-2 enrichment results written by `Baboon_SUN_PUT.R`.

---

## Figure 5 — SUN-PUT heatmaps

**Script:** `Baboon_SUN_PUT.R`, heatmap section.
**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData`.

---

## Figure 6 — Cross-species lung

**Script:** `Baboon_Human_LUN.R` → `Baboon_Human_LUN_Peak_Concordance_0.25_2h.pdf`.
**Key parameter:** `shift = 2`.
**Required data:** baboon and human lung chains from `mcmc_rho_BF3.RData` and
`mcmc_phi_BF3.RData`.

---

## Assembly

`assemble_figures.R` collects every panel above into `paper/subfigures/` under
a stable name and writes `paper/figures/Figure_N.pdf`. Figures 3 and 4 put two
panels side by side; the script merges them into a two-page PDF and says so,
since placing them side by side is the one step still done by hand.

---

## Notes

- `multi_conservation` writes Excel sheets named `"Results"` and
  `"Column_Definitions"`, so `read.xlsx` calls need `sheet = "Results"`.
- The phase concordance window is `shift = 2` hours throughout.
- Two copies of these scripts exist: `inst/analysis/` (this one, the one with
  `config.R`) and `R/v1/R/Thien/analysis/`. They have diverged;
  `circa_concordance.plots.R` is identical between them, the case-study
  scripts are not. Use this copy.

---

## Parameter search

`pipeline/` holds the scripts that survey the tunable parameters behind
Figures 3-5 before they are fixed. They read the same rho/phi summaries as the
figure scripts and write under `BAYRC_OUTPUT_DIR/param_search/`.

```
sunput_cache.R            SUN and PUT rho/phi slices -> cache_PUT_SUN.rds
sunput_phase_grid.R <a>   bfdr_alpha x phase window, one alpha per call,
                          with the OxPhos, proteasome and Parkinson gene sets
sunput_pathselect_grid.R  union/gain/loss/conserved enrichment per KEGG
                          release, size filter opened up
sunput_enrich_sweep.R     pathway list x size rule x stage-1 cut x stage-2 q
                          applied to the stored pathSelect tables
sunput_shifted_ora_grid.R over-representation of the shifted set per grid cell,
                          with and without the measured-gene background
case_study_counts.R a s   every gene-level count the Applications section
                          reports, for all three case studies
```

`case_study_counts.R` takes the alpha and window on the command line, so
pointing `BAYRC_RESULT_DIR` at one run or the other gives the same table for
both and makes the comparison direct.
