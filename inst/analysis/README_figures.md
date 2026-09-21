# BayRC Paper Figures

Every paper figure maps to one script here, with the inputs it reads and the
files it writes. Start from the table; the sections below give the detail.

| Figure | Script | Writes under |
|---|---|---|
| 1 | none, hand-drawn flowchart | — |
| 2A | `plots/heatmap_baboon.R`, drawn by `plots/replot_figure2.R` | `BAYRC_FIGURE_DIR/figure2/` |
| 2B | `plots/heatmap_circadian_pairs.R within_baboon`, drawn by `plots/replot_figure2.R` | `BAYRC_FIGURE_DIR/figure2/` |
| 2C | `plots/figure2_panelC.R` | `BAYRC_FIGURE_DIR/figure2/` |
| 3A, 3D clock panel | `applications/Baboon_SCN_HIP.R` | `BAYRC_FIGURE_DIR/baboon_SCN_HIP/` |
| 3B, 5A | `applications/Baboon_PUT_SUN.R` | `BAYRC_FIGURE_DIR/baboon_PUT_SUN/` |
| 3C, 5B | `applications/Baboon_PUT_VIC.R` | `BAYRC_FIGURE_DIR/baboon_PUT_VIC/` |
| 3D, 3E | `plots/summary_panels.R` | `<paper>/demos/figure3/` |
| 4 | `plots/pathway_transition_panels.R` | `<paper>/figures/` |
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
plots/layout_figures.py  reflows the assembled figures, see below
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
1. CAMO_h_b.R                        per-tissue MCMC, 26 human + 26 baboon
2. pipeline/summarize_rho_phi.R      -> mcmc_rho_BF3.RData, phi/mcmc_phi_BF3.RData
3. plots/heatmap_baboon.R            concordance matrices for Figure 2
   plots/heatmap_circadian_pairs.R within_baboon
4. applications/Baboon_SCN_HIP.R     the four case studies, each end to end
   applications/Baboon_PUT_SUN.R     add --replot to redraw from their tables
   applications/Baboon_PUT_VIC.R
   applications/Baboon_Human_LUN.R
5. plots/replot_figure2.R            Figure 2 panels A and B and the colour bar
   plots/figure2_panelC.R            Figure 2 panel C
   plots/summary_panels.R            Figure 3 panels D and E
6. assemble_figures.R                panels -> Figure_2 ... Figure_6
7. plots/pathway_transition_panels.R Figure 4, written straight to the figures
8. plots/layout_figures.py <paper>   reflows Figures 3, 5 and 6
```

Steps 3 to 8 run from `inst/analysis`. Step 7 comes after step 6 because it
writes `Figure_4.pdf` itself rather than handing panels to the assembler, and
step 8 comes last because it works on the assembled files. Steps 5 and 7 read
the tables and caches the earlier steps wrote, so neither needs the MCMC output.

Before re-running step 8 on a fresh build, empty
`<paper>/archive/before_legend_layout/`. That directory holds the assembler's
output as the reflow found it, which is what makes the reflow repeatable; the
script reads the copy there in preference to the live file.

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

The putamen scripts also print the phase offset over the whole maintained set
rather than over the shifted class alone. The shifted class is selected for
exceeding the phase window, so its mean is biased upward by the threshold that
defined it.

## Figures

**Figure 2** — genome-wide concordance across 26 baboon tissues and
KEGG-circadian concordance across the 25 excluding SCN. `plots/heatmap_baboon.R`
and `plots/heatmap_circadian_pairs.R` compute the concordance matrices;
`plots/replot_figure2.R` draws both panels from those matrices as
`Fig2A_genomewide.pdf` and `Fig2B_circadian.pdf`, with the cap and ramp from
`plots/palette_concordance.R`, and writes the run record.
`plots/heatmap_circadian_pairs.R` also takes `within_baboon_with_scn`,
`within_human` and `cross_species`. 325 within-species pairs, 676 cross-species.

**Figure 3** — within-species peak-phase concordance scatters: panel A from the
SCN-HIP script, panels B and C from the two putamen circuits. All three draw
`<pair>_Peak_Concordance.pdf` through `plots/peak_concordance.R`, with
condition A on the x axis and the same ±2 h band.

**Figure 4** — pathway transition enrichment in the two putamen circuits.
`plots/pathway_transition_panels.R` reads the `stage2_significant.csv` each
putamen script wrote and draws the whole figure in one go: panel A puts both
comparisons on one set of pathway rows, coloured by `-log10(q)` and sized by the
posterior expected gene count; panel B gives each enriched pathway's transition
composition, with a star on the transition it was enriched for. Because both
panels come out of one plot, the assembler has no Figure 4 panels to merge and
the script writes `Figure_4.pdf` directly.

**Figure 5** — KEGG Parkinson disease drawn for both circuits and stacked.
`plot_heatmap()` writes one file per pathway, plus a `_rhythmic_only` version
that drops the genes rhythmic in neither tissue. Panel A is drawn with
`show_legend = FALSE`; panel B carries the legend for the pair along its
bottom edge.

**Figure 6** — cross-species lung: the concordance scatter and the KEGG
Circadian rhythm heatmap. The phase-class key sits inside the scatter, so the
strip along the bottom carries only the heatmap's own keys.

## Sizing figures for print

The manuscript includes every figure at `\textwidth`, 488.5 pt. A panel is
first scaled into its slot by the assembler and then the whole figure is scaled
into the text block, so type set on a wide canvas is reduced twice and can reach
3 or 4 pt in print while looking correct at native size. Two controls exist for
this: `canvas_width` in `plot_heatmap()` narrows the page a heatmap is drawn on,
which raises its type without changing any font size, and `font_scale`
multiplies the type inside a heatmap when the canvas cannot narrow further. The
block titles and the heatmap title take a capped share of `font_scale`, since
they sit at fixed offsets and run into their neighbours otherwise.

Check a figure at the size it will print rather than on screen.

## Notes

- `multi_conservation` writes Excel sheets named `Results` and
  `Column_Definitions`, so `read.xlsx` calls need `sheet = "Results"`.
- The gene sets rename `ARNTL` to `BMAL1`, the symbol the atlas uses.
- `exploratory/` holds the other tissue pairs. Several of them read objects
  from the workspace a previous script left behind and do not run on their own.
