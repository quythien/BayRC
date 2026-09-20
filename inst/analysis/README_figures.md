# BayRC Paper Figures

Every paper figure maps to one script here, with the inputs it reads and the
files it writes. Start from the table; the sections below give the detail.

| Figure | Script | Writes under |
|---|---|---|
| 1 | none, hand-drawn flowchart | — |
| 2A | `plots/heatmap_baboon.R` | `BAYRC_FIGURE_DIR` |
| 2B | `plots/heatmap_circadian_pairs.R within_baboon` | `BAYRC_OUTPUT_DIR/heatmap_circadian_pairs/` |
| 3A | `applications/Baboon_SCN_HIP.R` | `BAYRC_FIGURE_DIR/baboon_SCN_HIP/` |
| 3B, 4, 5 | `applications/Baboon_SUN_PUT.R` | `BAYRC_FIGURE_DIR/baboon_SUN_PUT/` |
| 6 | `applications/Baboon_Human_LUN.R` | `BAYRC_FIGURE_DIR/baboon_human_LUN/` |
| S4 | `plots/Cosinor_residual_diagnostics_LUN*.R` | `BAYRC_FIGURE_DIR` |
| S5 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `BAYRC_FIGURE_DIR` |

## Layout

```
config.R             paths, all overridable by environment variable
applications/        the three paper case studies, each end to end
pipeline/            shared upstream steps and the parameter searches
plots/               shared drawing code: theme, palettes, heatmaps, scatters
exploratory/         other tissue pairs, not used by any paper figure
assemble_figures.R   collects the panels into Figure_N.pdf
```

Each application script runs standalone from a cold session and writes a
`run_record.txt` beside its figures naming the summary directory, the date its
`mcmc_rho_BF3.RData` was written, the parameters and the run date. That file,
not this one, answers which run produced a given figure.

## Paths

`config.R` reads every directory from an environment variable and falls back to
the original locations. It prints the summary directory and the date its rho
file was written, and stops if that file is missing, so the run behind a figure
is visible in the log.

```
export BAYRC_DATA_DIR=/path/to/Collaborative
export BAYRC_WD_DIR=/path/to/Circadian
export BAYRC_RESULT_DIR=/path/to/GTEXdata/result_fixed
export BAYRC_SUMMARY_DIR=$BAYRC_RESULT_DIR/summary/hb
export BAYRC_OUTPUT_DIR=/path/to/analysis/output
export BAYRC_FIGURE_DIR=/path/to/figure/output
export BAYRC_PATHWAY_DIR=/path/to/R/pathway_data
```

The default `BAYRC_RESULT_DIR` is the 2025 `result/` tree; point it at
`result_fixed/` for the corrected sampler.

## Order of operations

The MCMC output and the two `.RData` summaries are not in the repository.

```
1. CAMO_h_b.R                     per-tissue MCMC, 26 human + 26 baboon
2. pipeline/summarize_rho_phi.R   -> mcmc_rho_BF3.RData, phi/mcmc_phi_BF3.RData
3. the application and plot scripts
4. assemble_figures.R             panels -> Figure_2 ... Figure_6
```

## Analysis parameters

The three application scripts state these as plain assignments near the top.

```
bfdr_alpha     0.25     BFDR level, in transition_classify and phase_infer
shift          2        phase window in hours
pathway list   kegg_pathway_list_hsa.rds, gene sets cut to the measured genes
min_measured   15       pathways kept at 15 or more measured genes
stage 1        union enrichment, BH q < 0.10, a screen rather than a test
stage 2        gain, loss and conservation on the active set, q < 0.20
nperm          10000    fgsea permutations
```

Stage 1 is deliberately the looser of the two: it screens for pathways with
rhythmic signal in either tissue, and stage 2 carries the inference.

## Figures

**Figure 2** — genome-wide and KEGG-circadian concordance across 26 baboon
tissues. `plots/heatmap_baboon.R` writes
`Baboon_Concordance_Heatmap_0.5_<method>.pdf` and `Baboon_Concordance_Matrix.csv`;
the cap and ramp come from `plots/palette_concordance.R`.
`plots/heatmap_circadian_pairs.R` also takes `within_baboon_with_scn`,
`within_human` and `cross_species`. 325 within-species pairs, 676 cross-species.

**Figure 3** — within-species peak-phase concordance scatters, panel A from the
SCN-HIP script and panel B from SUN-PUT. Both draw
`<pair>_Peak_Concordance.pdf` through `plots/peak_concordance.R`.

**Figure 4** — SUN-PUT pathway transition enrichment, from the stage-2 output:
`SUN_PUT_transition_enrichment.pdf`, with `stage1_union.csv` and
`stage2_significant.csv` beside it.

**Figure 5** — SUN-PUT pathway heatmaps for KEGG Parkinson disease and KEGG
Oxidative phosphorylation. `plot_heatmap()` writes one file per pathway, plus a
`_rhythmic_only` version that drops the genes rhythmic in neither tissue.

**Figure 6** — cross-species lung: the concordance scatter and the KEGG
Circadian rhythm heatmap.

## Notes

- `multi_conservation` writes Excel sheets named `Results` and
  `Column_Definitions`, so `read.xlsx` calls need `sheet = "Results"`.
- The gene sets rename `ARNTL` to `BMAL1`, the symbol the atlas uses.
- `exploratory/` holds the other tissue pairs. Several of them read objects
  from the workspace a previous script left behind and do not run on their own.
