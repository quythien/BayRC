# Manuscript figure order and reproduction

The manuscript uses the following order. Analysis directories and output names
follow these numbers; cached output from an older checkout must be regenerated
or migrated before running the assemblers.

| Figure | Content | Main generating scripts |
|---|---|---|
| 1 | Bayesian framework and biomarker classes | Author-supplied artwork |
| 2 | Brain contrasts: rhythmic transitions, timing, then gene profiles | `plots/summary_panels.R`, the three `applications/Baboon_*` brain scripts, `assemble_figures.R`, `plots/layout_figures.py` |
| 3 | Pathway enrichment and transition composition | `plots/pathway_transition_panels.R` |
| 4 | Parkinson disease pathway in PUT–SUN and PUT–VIC | Brain application scripts, `assemble_figures.R`, `plots/layout_figures.py` |
| 5 | Human–baboon and human–mouse lung comparisons | `plots/figure5/make_figure5.R`, `plots/figure5/assemble_figure5.py` |
| 6 | Genome-wide and circadian pathway concordance; rhythmicity and phase profiles | `plots/replot_figure6.R`, `plots/figure6_panelC.R`, `plots/figure6/export_pathway.R`, `plots/figure6/explore_phase_groups.R`, `plots/figure6/assemble_figure6_nature.py` |

Figure 2 panels A/B contain transition counts and phase classifications, with
contrasts ordered PUT–SUN, PUT–VIC, SCN–HIP. Panels C/D/E show the corresponding
gene profiles in that order. The SCN clock-reference diagnostic is not a main
Figure 2 panel.

Figure 5 panels A/B compare human with baboon and mouse, respectively; panel C
is the three-species circadian pathway heatmap. Figure 6 panels A/B show
concordance, C posterior rhythmicity, and D phase profiles.

## Build sequence

Run commands from `inst/analysis`, with paths configured through `config.R`.
After the application scripts have produced their plot caches and panels:

```sh
Rscript plots/summary_panels.R
Rscript assemble_figures.R /path/to/paper
python3 plots/layout_figures.py /path/to/paper
Rscript plots/pathway_transition_panels.R
Rscript plots/figure5/make_figure5.R
python3 plots/figure5/assemble_figure5.py
Rscript plots/replot_figure6.R
Rscript plots/figure6_panelC.R
Rscript plots/figure6/export_pathway.R
FIG6_CAIRO=1 Rscript plots/figure6/explore_phase_groups.R
python3 plots/figure6/assemble_figure6_nature.py
```

Copy `figure5/Figure_5_column.pdf` to `paper/figures/Figure_5.pdf` and
`figure6/Figure_6.pdf` to `paper/figures/Figure_6.pdf`.
`layout_figures.py` writes the final Figure 2 and compact Figure 4 directly.
Always run `assemble_figures.R` immediately before `layout_figures.py` so the
layout step receives fresh panels. `assemble_figures.R --dry-run` checks input
locations without writing figures; incomplete sets are not assembled.

## Changes from the previous order

Previous Figures 2, 3, 4, 5 and 6 are now Figures 6, 2, 3, 4 and 5,
respectively. The previous Figure 3 summary panels D/E are now Figure 2 A/B;
its scatter panels B/C/A are now Figure 2 C/D/E. The compact row layout is the
main Figure 4. The previously numbered SCN reference diagnostic stays separate.
The old `FIG2_*` concordance settings are now `FIG6_*`; the old `FIG6_*`
cross-species settings are now `FIG5_*`.

## Supplementary plots

`plots/figure_S_operating.R` regenerates the BFDR simulation figure. Its
`--summaries` mode uses saved summary CSVs; panels D/E omit amplitude-to-noise
ratio 0.5. Replicate counts come from the saved summaries.

`pathway_phase_hierarchy()`, `plot_pathway_phase_hierarchy()` and
`plot_pathway_profiles()` implement the pathway clustering and paired profile
plots. Their package help pages describe coverage, tree cuts and the optional
sparse-data dissimilarity shrinkage. Paper-specific settings and thyroid
posterior contrasts are documented in Supplementary Sections S9 and S10.

The manuscript, submission figures and Overleaf files remain in the local
`paper/` directory, which is excluded from this repository.
