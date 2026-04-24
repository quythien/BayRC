# BayRC Paper Figures — Generating Scripts

This file maps each figure in `paper/figures/` to the analysis script that produced it.
All scripts are in `inst/analysis/`. Run each script from the BayRC root directory
after sourcing `config.R` to set data paths.

**Prerequisites:** All scripts require pre-computed MCMC `.RData` objects
(not in the repository — deposit location to be added). Set paths via `config.R`.

---

## Figure 1 — BayRC Framework Overview (flowchart)

**File:** `paper/figures/Figure_1.pdf`  
**Source:** Manually created diagram (`paper/figures/flowchart_resized.jpeg`).  
No generating R script. Created externally (illustration tool).

---

## Figure 2 — Genome-Wide Concordance Heatmaps (26 baboon tissues)

**File:** `paper/figures/Figure_2_concordance_heatmaps_combined_ward.D2.pdf`  
**Script:** `inst/analysis/circa_concordance.plots.R`  
**Inputs:** Multi-tissue MCMC outputs (26 baboon tissues from GTEx baboon data)  
**Key output calls:** `ggsave()` producing concordance heatmaps (genome-wide + KEGG Circadian pathway)  
**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData` from `BAYRC_GTEX_DIR`

---

## Figure 3 — Within-Species Phase Concordance Scatter (SCN-HIP + SUN-PUT)

**File:** `paper/figures/Figure_3_within_species_scatter.pdf`  
**Scripts:**
- SCN-HIP panel: `inst/analysis/Baboon_SCN_HIP.R` — saves `Baboon_SCN_HIP_Peak_Concordance_0.25_2h_new.pdf`
- SUN-PUT panel: `inst/analysis/Baboon_SUN_PUT.R` — saves `Baboon_SUN_PUT_Peak_Concordance_0.25_2h_new.pdf`

Both panels are combined into Figure 3.  
**Key parameter:** `shift = 2` (±2h phase concordance window)  
**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData`

---

## Figure 4 — SUN-PUT Pathway Enrichment Dotplots

**File:** `paper/figures/Figure_4_enrichment_combined.pdf`  
**Script:** `inst/analysis/Baboon_SUN_PUT.R` and `inst/analysis/plot_enrich_SUN_PUT.R`  
**Outputs:** `_barplot.pdf` and `_dotplot.pdf` files  
Also produces: `paper/figures/SUN_PUT_shifted_GO_BP_dotplot.pdf` and  
`paper/figures/SUN_PUT_shifted_KEGG_dotplot.pdf`  
**Required data:** Stage 2 enrichment results from `Baboon_SUN_PUT.R`

---

## Figure 5 — SUN-PUT Heatmaps

**File:** `paper/figures/Figure_5_sunput_heatmaps.pdf`  
**Script:** `inst/analysis/Baboon_SUN_PUT.R` (heatmap section)  
**Required data:** `mcmc_rho_BF3.RData`, `mcmc_phi_BF3.RData`

---

## Figure 6 — Cross-Species Lung Circadian Analysis

**File:** `paper/figures/Figure_6_lung_circadian.pdf`  
**Script:** `inst/analysis/Baboon_Human_LUN.R`  
**Key output:** `Baboon_Human_LUN_Peak_Concordance_0.25_2h.pdf`  
**Key parameter:** `shift = 2` (±2h phase concordance window)  
**Required data:** Baboon lung MCMC output + Human lung MCMC output

---

## Recommended Execution Order

```
1. Baboon_SCN_HIP.R     → Figure 3 (panel A)
2. Baboon_SUN_PUT.R     → Figures 3 (panel B), 4, 5
3. plot_enrich_SUN_PUT.R → Figure 4 (finalized dotplots)
4. Baboon_Human_LUN.R   → Figure 6
5. circa_concordance.plots.R → Figure 2
```

Figure 1 requires no R script.

---

## Notes

- All scripts currently contain hardcoded paths (`/home/qtp1/Projects/...`).
  Replace with `source("config.R")` at the top of each script and use
  `BAYRC_DATA_DIR`, `BAYRC_WD_DIR` etc. from `config.R`.
- Excel output sheet names in `multi_conservation` results are `"Results"` and
  `"Column_Definitions"` — downstream `read.xlsx` calls should use `sheet = "Results"`.
- Phase concordance threshold is `shift = 2` hours throughout (paper §2.2, "±2h").
