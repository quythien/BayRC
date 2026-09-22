# BayRC Paper Figures

Every paper figure maps to one script here, with the inputs it reads and the
files it writes. Start from the table; the sections below give the detail and
the order to run things in.

`README_numbers.md` maps every number the manuscript quotes to the script and
output it comes from.

Figure 1 is a hand-drawn flowchart. No script produces it.

| Figure | Script | Writes |
|---|---|---|
| 1 | none, hand-drawn flowchart | — |
| 2A | `plots/heatmap_baboon.R`, drawn by `plots/replot_figure2.R` | `figure2/Fig2A_genomewide_nature.pdf` |
| 2B | `plots/heatmap_circadian_pairs.R within_baboon`, drawn by `plots/replot_figure2.R` | `figure2/Fig2B_circadian_nature.pdf`, `figure2/Fig2_concordance_legend.pdf` |
| 2C | `plots/figure2_panelC.R` | `figure2/Fig2C_circadian_membership_nature.pdf`, `figure2/Fig2C_circadian_membership.csv` |
| 2D | `plots/figure2/export_pathway.R`, then `plots/figure2/explore_phase_groups.R` | `figure2/option8_balanced_global_nature.pdf` |
| 2 | `plots/figure2/assemble_figure2_nature.py` | `figure2/Figure_2.pdf` |
| 3A, 3D clock panel | `applications/Baboon_SCN_HIP.R` | `baboon_SCN_HIP/` |
| 3B, 5A | `applications/Baboon_PUT_SUN.R` | `baboon_PUT_SUN/` |
| 3C, 5B | `applications/Baboon_PUT_VIC.R` | `baboon_PUT_VIC/` |
| 3D, 3E | `plots/summary_panels.R` | `<paper>/demos/figure3/summary_panels.pdf` |
| 3, 5 | `assemble_figures.R`, then `plots/layout_figures.py` | `<paper>/figures/Figure_3.pdf`, `Figure_5.pdf`, `Figure_5_row.pdf` |
| 4 | `plots/pathway_transition_panels.R` | `<paper>/figures/Figure_4.pdf` |
| 6A, 6B, 6C | `plots/figure6/make_figure6.R` | `figure6/Fig6A_human_baboon.pdf`, `Fig6B_human_mouse.pdf`, `Fig6C_heatmap.pdf` |
| 6 | `plots/figure6/assemble_figure6.py` | `figure6/Figure_6_column.pdf` |
| S4 | `plots/Cosinor_residual_diagnostics_LUN*.R` | `BAYRC_FIGURE_DIR` |
| S5 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `BAYRC_FIGURE_DIR` |
| Table S8 | `pipeline/bfdr_calibration.R <seed>` for seeds 1 to 10, then `--summary` | `calibration/bfdr_calibration_summary.csv` |
| S8 | `pipeline/bfdr_calibration.R <seed> 2000 - <n> <A> 1` over n in 12, 24, 48, A in 0.5, 1, 1.5, 2, 3 and seeds 1 to 10, then `plots/figure_S_operating.R` | `Figure_S_operating.pdf`, `calibration/bfdr_operating_characteristics.csv`, `calibration/bayrc_cosinor_auc.csv` |
| Table 1 | `applications/Baboon_PUT_SUN.R`, `applications/Baboon_PUT_VIC.R`, then `gain_loss_ratios.R` | `<pair>/stage2_significant.csv`, `<pair>/pathway_metrics.csv`, `gain_loss_ratios.csv` |

Relative paths in the third column are under `BAYRC_FIGURE_DIR`. `<paper>` is
`$BAYRC_RESULT_DIR/paper` with the settings below.

## Layout

```
config.R                 paths, all overridable by environment variable
applications/            the four paper case studies, each end to end
pipeline/                shared upstream steps and the parameter searches
plots/                   shared drawing code: theme, palettes, heatmaps, scatters
plots/figure2/           Figure 2 panel D and the Figure 2 layout
plots/figure6/           the three-species Figure 6 panels and layout
plots/layout_figures.py  reflows the assembled Figures 3 and 5
mouse/                   the mouse chains and the human/mouse summaries
exploratory/             other tissue pairs, not used by any paper figure
assemble_figures.R       collects the panels and merges Figures 3 and 5
```

Each application script runs standalone from a cold session and writes a
`run_record.txt` beside its figures naming the summary directory, the date its
`mcmc_rho_BF3.RData` was written, the parameters and the run date. That file,
not this one, answers which run produced a given figure.

## Requirements

- R with the package installed (`R CMD INSTALL .` from the checkout), plus
  `ComplexHeatmap`, `circlize`, `cluster` and `ggplot2`.
- `pdflatex` and `pdfinfo` on the PATH, for `assemble_figures.R`.
- Python 3 with PyMuPDF, for every `.py` script: `pip install pymupdf`. A
  script stops with that instruction when the module is missing.

## Paths

`config.R` reads every directory from an environment variable and derives the
fallback from where the checkout sits. It prints the summary directory and the
date its rho file was written, and stops if that file is missing, so the run
behind a figure is visible in the log. Set these before anything else:

```
export BAYRC_RESULT_DIR=/path/to/GTEXdata/result_fixed   # the MCMC output
export BAYRC_OUTPUT_DIR=$BAYRC_RESULT_DIR/analysis        # tables and caches
export BAYRC_FIGURE_DIR=$BAYRC_RESULT_DIR/analysis/figures
export BAYRC_GTEX_DIR=/path/to/GTEXdata                   # holds data/CAMO.bab.hum.RData
```

Optional:

```
export BAYRC_SUMMARY_DIR=$BAYRC_RESULT_DIR/summary/hb   # the default
export BAYRC_PATHWAY_DIR=/path/to/pathway_data          # default: R/pathway_data beside the checkout, else inst/extdata
export BAYRC_PIPELINE_DIR=/path/to/Pipeline             # holds one_cosinor_OLS_new.R, cosinor scripts only
export BAYRC_DATA_DIR, BAYRC_WD_DIR, BAYRC_AGING_DIR    # the aging pipeline scripts only
```

The Python scripts read `BAYRC_FIGURE_DIR` for their panel directory and take it
as their first argument otherwise.

## Summaries

The MCMC output and the summary objects are not in the repository. Two summary
directories sit under `BAYRC_RESULT_DIR/summary`:

- `hb` holds the human and baboon summaries, `mcmc_rho_BF3.RData` and
  `phi/mcmc_phi_BF3.RData`, built by `pipeline/summarize_rho_phi.R`. Every
  within-species figure and the baboon side of Figure 6 read from here; this is
  `BAYRC_SUMMARY_DIR`.
- `hm` holds the human and mouse summaries under the same two names, built by
  `mouse/build_mouse_summaries.R` from the six mouse tissue chains. Only
  `plots/figure6/make_figure6.R` and `mouse/mouse_diagnostics.R` read from here,
  and `make_figure6.R` stops when it is absent.

The per-tissue chains come from `run_fixed.R` for human and baboon and
`mouse/run_mouse.R` for mouse. Both take one tissue per process, so a driver can
fan the runs out across cores, and both read their atlas from
`BAYRC_GTEX_DIR/data`.

## Order of operations

Everything below runs from `inst/analysis` with the variables above exported.
Each step is a fresh R process. Steps 3 to 7 read the tables and caches the
earlier steps wrote, so none of them needs the chains after step 1.

```
# 1. summaries
Rscript pipeline/summarize_rho_phi.R                       # -> summary/hb
Rscript mouse/build_mouse_summaries.R $BAYRC_RESULT_DIR $BAYRC_RESULT_DIR/summary/hb \
        $BAYRC_RESULT_DIR/summary/hm ../extdata/ensembl_symbol_map.csv       # -> summary/hm

# 2. concordance matrices behind Figure 2
Rscript plots/heatmap_baboon.R
Rscript plots/heatmap_circadian_pairs.R within_baboon

# 3. the four case studies; add --replot to redraw from their tables
Rscript applications/Baboon_SCN_HIP.R
Rscript applications/Baboon_PUT_SUN.R
Rscript applications/Baboon_PUT_VIC.R
Rscript applications/Baboon_Human_LUN.R    # cut to the genes mouse also measures

# 4. Figure 2
Rscript plots/replot_figure2.R                             # panels A, B and the colour bar
Rscript plots/figure2_panelC.R                             # panel C and its membership table
Rscript plots/figure2/export_pathway.R                     # phase table for panel D
FIG2_CAIRO=1 Rscript plots/figure2/explore_phase_groups.R  # panel D
python3 plots/figure2/assemble_figure2_nature.py           # -> figure2/Figure_2.pdf

# 5. Figures 3 and 5
Rscript plots/summary_panels.R                             # Figure 3 panels D and E
Rscript assemble_figures.R                                 # panels -> <paper>/figures/Figure_3, Figure_5
python3 plots/layout_figures.py $BAYRC_RESULT_DIR/paper    # reflows Figures 3 and 5

# 6. Figure 4
Rscript plots/pathway_transition_panels.R                  # -> <paper>/figures/Figure_4.pdf

# 7. Figure 6
Rscript plots/figure6/make_figure6.R                       # panels A, B, C
python3 plots/figure6/assemble_figure6.py                  # -> figure6/Figure_6_column.pdf
python3 plots/figure6/check_heatmap.py $BAYRC_FIGURE_DIR/figure6/Fig6C_heatmap.pdf
```

The manuscript's `Figure_2.pdf` is `figure2/Figure_2.pdf` from step 4 and its
`Figure_6.pdf` is `figure6/Figure_6_column.pdf` from step 7; copy both into
`<paper>/figures/`. `assemble_figures.R` collects the Figure 2 and 6 panels for
the record but never writes either figure, and `layout_figures.py` leaves
`Figure_6.pdf` alone, so a rerun of step 5 cannot replace them.

`exploratory/run_figure2_chain.R` runs step 2 in one process per script, and
`assemble_figures.R --dry-run` lists which panel each figure is still missing.

Before re-running `layout_figures.py` on a fresh build, empty
`<paper>/archive/before_legend_layout/`. That directory holds the assembler's
output as the reflow found it, which is what makes the reflow repeatable; the
script reads the copy there in preference to the live file.

## Analysis parameters

The application scripts state these as plain assignments near the top.

```
bfdr_alpha     0.25     BFDR level, in transition_classify and phase_infer
shift          2        phase window in hours
pathway list   kegg_pathway_list_hsa.rds, gene sets cut to the measured genes
min_measured   15       pathways kept at 15 or more measured genes, which is
                        229 pathways on the 5,066 baboon-human genes and 225 on
                        the 4,893 the three species share
stage 1        union enrichment, BH q < 0.20, a screen rather than a test
stage 2        gain, loss and conservation on the active set, q < 0.05
nperm          10000    fgsea permutations
```

Stage 1 is deliberately the looser of the two: it screens for pathways with
rhythmic signal in either tissue, and stage 2 carries the inference.

The putamen scripts also print the phase offset over the whole maintained set
rather than over the shifted class alone. The shifted class is selected for
exceeding the phase window, so its mean is biased upward by the threshold that
defined it.

## Figures

**Figure 2** — genome-wide concordance across 26 baboon tissues, KEGG-circadian
concordance across the 25 excluding SCN, the membership of the circadian
pathway in those tissues, and the peak timing of its genes.
`plots/heatmap_baboon.R` and `plots/heatmap_circadian_pairs.R` compute the
concordance matrices; `plots/replot_figure2.R` draws panels A and B from those
matrices, with the cap and ramp from `plots/palette_concordance.R`, and writes
the run record. `plots/figure2_panelC.R` draws panel C and writes the
tissue-by-gene membership table. `plots/figure2/export_pathway.R` summarises
the posterior phase of every called cell of that table from
`phi/mcmc_phi_BF3.RData`, and `plots/figure2/explore_phase_groups.R` draws
panel D from it, grouping the genes by phase on both tissue groups weighted
equally; `FIG2_CAIRO=1` draws it with `cairo_pdf`, which embeds the fonts.
Panels that have to carry type at one printed size are built on one native
width, the `_nature` files, and `assemble_figure2_nature.py` only places them.
`assemble_figure2_column.py` is an alternative layout from the same panels.
`plots/heatmap_circadian_pairs.R` also takes `within_baboon_with_scn`,
`within_human` and `cross_species`: 325 within-species pairs, and 625
cross-species, which excludes SCN from both species.
`pipeline/pairwise_concordance_all.R` keeps SCN and so runs 676.

**Figure 3** — within-species peak-phase concordance scatters: panel A from the
SCN-HIP script, panels B and C from the two putamen circuits. All three draw
`<pair>_Peak_Concordance.pdf` through `plots/peak_concordance.R`, with
condition A on the x axis and the same ±2 h band. `plots/summary_panels.R`
draws the transition counts and posterior phase-class proportions from each
pair's `plot_data.rds`, and `layout_figures.py` seats them under the scatters.

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
bottom edge. `Figure_5_row.pdf` holds the same two panels along a row.

**Figure 6** — the circadian pathway in lung across three species, with human
as the reference. `plots/figure6/make_figure6.R` reads the `hb` and `hm`
summaries, classifies human against baboon and human against mouse on the
frozen parameters, draws the two concordance scatters with human on the x axis
in both, and draws the three-condition heatmap through `plot_heatmap()` with
`data3`. `assemble_figure6.py` writes the stacked and the column arrangement
and reports the smallest type each reaches at text width; the paper carries the
column. `check_heatmap.py` reads the drawn heatmap back and fails if the two
offset blocks were built from different rules or any labels collide.
`demo_three_regions.R` runs the same three-condition heatmap on SCN,
hippocampus and putamen as a check on that code path.

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
  from the workspace a previous script left behind and do not run on their own,
  and some carry object and label names from the pair they were copied from, so
  read the names against the tissues at the top of the file.
- `pipeline/CAMO_Aging_*.R`, `pipeline/Pipleline.R` and `plots/figures.R`
  belong to the aging analysis and read from `BAYRC_AGING_DIR`; they are not
  behind any figure in this paper.

## Table 1

Table 1 is assembled by hand from two files each case study writes, and it is
the one item here without a script of its own:

- `stage2_significant.csv` supplies the pathway size, the expected gain, loss
  and conserved counts, and the Stage 2 q-values. The expected union is their
  sum.
- `pathway_metrics.csv` supplies the adjusted concordance.
- `gain_loss_ratios.R` supplies the gain-loss ratio, E[Gain] / E[Loss] as the
  Methods define it: per pathway from the expected counts in
  `stage2_significant.csv`, and genome-wide from the marginal posterior
  probabilities. It also prints the Spearman correlation over the 325 baboon
  tissue pairs that the Results quote.
