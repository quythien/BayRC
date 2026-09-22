# BayRC Paper Numbers

Every number the manuscript and its supplement quote maps to one row here: where
it is printed, the script that produces it, and the file, column or printed line
it comes from. `README_figures.md` covers the figures themselves and the order
the scripts run in; this file covers the numbers inside the text, the captions
and the tables.

Scripts are named relative to `inst/analysis`. Output paths are relative to
`BAYRC_OUTPUT_DIR` unless they start with `FIG/`, which is `BAYRC_FIGURE_DIR`.
`figures/<pair>/` is the directory each application script wrote its tables,
`plot_data.rds` cache and `run_record.txt` into. A value in quotes in the last
column is the label of a line the script prints.

`paper_numbers.R` prints every number that no other script prints, grouped
under the paper section that quotes it. It reads the stored summaries, caches
and tables, runs in about five minutes and writes nothing:

```
Rscript paper_numbers.R > paper_numbers.txt
```

The other scripts named below print or write their own numbers as part of the
pipeline in `README_figures.md`. `pipeline/case_study_counts.R` prints several
of the same transition counts and reads the older KEGG release
`kegg.pathway.list_hsa.RData` from `BAYRC_PATHWAY_DIR`, so it needs that
directory to hold both lists.

---

## Abstract

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Abstract | Spearman r = 0.974 on gene ranks | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "Spearman r (posterior vs p-value)" |
| Abstract | Pearson r = 0.9998 on phase | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "Pearson r (phase, n = 1917)" |
| Abstract | 5 of 5,066 genes after Benjamini-Hochberg | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv`, `cosinor_total` on the `BH q < 0.05` row |
| Abstract | 26 baboon tissues; 25 outside the SCN | `paper_numbers.R` | Case Study 1: "tissues, pairs"; "circadian matrix tissues, SCN included" |
| Abstract | nine peripheral tissues, mean 0.413; sixteen remaining, 0.157 | `pipeline/cluster_rhythmic_content.R` | "tight cluster, mean within-cluster c-score" and "remaining, mean within-group c-score" |
| Abstract | genome-wide concordance 0.133 and 0.132 | `genome_concordance_pvalues.R` | `genome_concordance_pvalues.csv`, `adjusted_concordance` |
| Abstract | gain-loss ratio 0.55 against 0.99 | `gain_loss_ratios.R` | `gain_loss_ratios.csv`, genome-wide rows: 0.554 and 0.990 |
| Abstract | about 2.6 h later | `paper_numbers.R` | "mean comparator minus reference (Figure 3E)", PUT_SUN and PUT_VIC |
| Abstract | 74% and 77% phase-shifted | `paper_numbers.R` | "Phase-shifted", PUT_SUN and PUT_VIC |
| Abstract | 4,893 genes in three species | `paper_numbers.R` | Applications: "three-species genes" |
| Abstract | 19 rhythm-conserved with baboon, 12 with mouse | `paper_numbers.R` | Case Study 3: "tau_cons, rhythm-conserved" |
| Abstract | human lung 64 rhythmic genes | `paper_numbers.R` | Case Study 3: "rhythmic human / baboon / mouse" |
| Abstract | 12 of 19 phase-conserved; none of 12 | `paper_numbers.R` | Case Study 3: "phase conserved / shifted / undetermined" |
| Abstract | circular mean 10.0 h in lung | `mouse/mouse_diagnostics.R` | "mouse minus human peak offset", LUN |
| Abstract | intervals agree to 1.47 h | `mouse/mouse_diagnostics.R` | "peak measured from BMAL1 within each species", LUN |
| Abstract | dropped from the shortened abstract; screen numbers below | `applications/Baboon_Human_LUN.R` | `figures/baboon_human_LUN/stage1_union.csv`; `paper_numbers.R` "baboon LUN: genes, tested, retained, rank, q" |

## Methods

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| 2.2, transition labels | one gene, ENDOD1, clears gain and conservation | `paper_numbers.R` | Case Study 3, mouse: "genes clearing two thresholds"; the three baboon pairs print 0 and baboon lung prints none |
| 2.3, expected counts | 5 genes of 5,066 after Benjamini-Hochberg | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv`, `BH q < 0.05` |
| 2.3, Stage 1 and 2 | FDR < 0.20; q < 0.05 | `applications/*.R` | `figures/<pair>/run_record.txt`, `stage1_q`, `stage2_q` |
| 2.4, permutations | B = 1,000 pathway, 500 genome-wide | `applications/Baboon_PUT_*.R`, `genome_concordance_pvalues.R` | `pathway_metrics.csv` `PValue` floor 1/1,001; `genome_concordance_pvalues.csv` `permutations` |

## Applications

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Applications | 64 tissues profiled by Mure et al. | none, from the cited paper | — |
| Applications | 12 timepoints 2 h apart | `paper_numbers.R` | "baboon lung samples, spacing" |
| Applications | 26 tissues with a GTEx counterpart | `paper_numbers.R` | "baboon tissues matched to GTEx" |
| Applications | mouse 6 h spacing, 8 timepoints, 4 distinct phases | `paper_numbers.R` | "mouse samples, spacing, distinct phases" |
| Applications | six mouse tissues AOR CER HEA KIC LIV LUN | `paper_numbers.R` | "mouse tissues" |
| Applications | 5,066 baboon-human genes; 8,295 mouse-human | `paper_numbers.R` | "baboon-human matched genes"; "mouse-human matched genes" |
| Applications | 4,893 in the intersection; 173 with no mouse counterpart, no core clock gene | `paper_numbers.R` | "three-species genes"; "baboon-human genes with no mouse counterpart"; "of them clock genes / KEGG circadian genes" |
| Applications | alpha = 0.25; delta = 2 h | `applications/*.R` | `figures/<pair>/run_record.txt`, `bfdr_alpha`, `shift` |
| Applications | 229 pathways on 5,066 genes, 225 on 4,893 | `paper_numbers.R` | "KEGG pathways with >= 15 of 5,066 genes"; the 4,893-gene line gives 225 |
| Applications | Stage 1 q < 0.20, Stage 2 q < 0.05, 10,000 permutations | `applications/*.R` | `run_record.txt`, `stage1_q`, `stage2_q`, `nperm` |

## Posterior Against Cosinor in Baboon Lung

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Results | 5,066 genes | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "genes compared" |
| Results | p_rhythmic = 0.2 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "prior p_rhythmic" |
| Results | BF > 1 detects 3,440 and contains all 1,969 at p < 0.05 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R`; `paper_numbers.R` | `FIG/S5_threshold_grid.csv`, `BF1_cosinor_only` = 0; "BF > 1: posterior cut, genes" |
| Results | BF > 3 detects 2,280; BF > 5 detects 1,744 | `paper_numbers.R` | "BF > 3 ..." and "BF > 5: posterior cut, genes" |
| Results | Spearman r = 0.974; 202 of the top 253 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "top 5% cutoff k = 253 of 5066; overlap = 202 (80%)" |
| Results | 1,917 detected by both; Pearson r = 0.9998 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_detection_counts.csv`, "Both methods"; "Pearson r (phase, n = 1917)" |
| Results | Benjamini-Hochberg leaves 5; 12 samples | `plots/S5_Bayes_Cosinor_Agreement_LUN.R`; `paper_numbers.R` | `S5_threshold_grid.csv`; "baboon lung samples, spacing" |

## Case Study 1

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| CS1 | 26 tissues over 5,066 genes | `paper_numbers.R` | "tissues, pairs"; "baboon-human matched genes" |
| CS1 | median 0.11, no pair above 0.33 | `paper_numbers.R` | "median / max genome-wide score" (0.1118 / 0.3259) |
| CS1 | putamen with skeletal muscle, stomach and omental fat | `paper_numbers.R` | "Figure 2A clade holding PUT" (MUA OMF PUT STF) |
| CS1 | cerebellum with aorta, amygdala and pituitary | `paper_numbers.R` | "Figure 2A clade holding CER" |
| CS1 | SCN-HIP 0.243, p = 0.002 | `genome_concordance_pvalues.R` | `genome_concordance_pvalues.csv`, `SCN-HIP` |
| CS1 | SCN-HIP in the top tenth | `paper_numbers.R` | "SCN-HIP score, rank among pairs" (7 of 325) |
| CS1 | 23 measured genes; 25 tissues other than SCN | `paper_numbers.R` | "KEGG pathways with >= 15 of 5,066 genes"; "circadian matrix tissues, SCN included" |
| CS1 | nine tissues named; k = 2 cut | `pipeline/cluster_rhythmic_content.R` | the tissue list under "tight cluster" |
| CS1 | 0.413 against 0.157 | `pipeline/cluster_rhythmic_content.R` | as in the abstract |
| CS1 | nine core-loop genes called in almost every peripheral tissue | `plots/figure2_panelC.R` | `FIG/figure2/Fig2C_circadian_membership.csv`, `called` (8 or 9 of 9 each) |
| CS1 | 13.1 of 23 against 8.5 | `plots/figure2_panelC.R` | "genes called rhythmic, mean per tissue" |
| CS1 | 78.7% against 59.7% | `plots/figure2_panelC.R` | "pairs in the same block agreeing on a gene's call" |
| CS1 | r = 0.19 | `pipeline/cluster_rhythmic_content.R` | "correlation of mean c-score with pathway rhythmic percent" |
| CS1 | three timing groups and their genes | `plots/figure2/explore_phase_groups.R` | `phase_group_assignments.csv`, mode `balanced_global` |
| CS1 | 43 gained, 81 conserved, 94 lost | `paper_numbers.R` | SCN_HIP "gain / conserved / loss" |
| CS1 | 81.5% within 2 h | `paper_numbers.R` | SCN_HIP "conserved within 2 h" |
| CS1 | 74% conserved, 14% shifted, 12% undetermined | `paper_numbers.R` | SCN_HIP "Phase-conserved", "Phase-shifted", "Undetermined" |
| CS1 | hippocampus later for 7, earlier for 4 | `paper_numbers.R` | SCN_HIP "shifted, comparator later / earlier" |
| CS1 | mean difference +0.22 h | `paper_numbers.R` | SCN_HIP "mean comparator minus reference (Figure 3E)" |

## Case Study 2

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| CS2 | 0.133 and 0.132, p = 0.002 in both | `genome_concordance_pvalues.R` | `genome_concordance_pvalues.csv` |
| CS2 | 110 / 227 / 465 against SUN; 277 / 679 / 291 against VIC | `paper_numbers.R` | "gain / conserved / loss", PUT_SUN and PUT_VIC |
| CS2 | gain-loss ratio 0.55 and 0.99 | `gain_loss_ratios.R` | `gain_loss_ratios.csv`, genome-wide rows |
| CS2 | 325 tissue pairs, Spearman 0.99 | `gain_loss_ratios.R` | "325 tissue pairs: Spearman correlation ..." (0.992) |
| CS2 | +2.57 h and +2.65 h | `paper_numbers.R` | "mean comparator minus reference (Figure 3E)" |
| CS2 | conserved sets differ threefold | `paper_numbers.R` | "conserved PUT-VIC over PUT-SUN" (2.99) |
| CS2 | 24.7% and 26.5% within 2 h | `paper_numbers.R` | "conserved within 2 h" |
| CS2 | 74% and 77% shifted, 10% and 13% conserved, 15% and 10% undetermined | `paper_numbers.R` | "Phase-shifted", "Phase-conserved", "Undetermined" |
| CS2 | 169 shifted, none earlier; 524, 2 earlier | `paper_numbers.R` | "shifted, comparator later / earlier" |
| CS2 | +0.22 h for SCN against hippocampus | `paper_numbers.R` | SCN_HIP "mean comparator minus reference (Figure 3E)" |
| CS2 | 229 pathways tested; Stage 1 retained 8 in each | `paper_numbers.R` | "pathways tested / Stage 1 retained (q < 0.20)" |
| CS2 | Stage 2 found 7 in each | `paper_numbers.R` | "Stage 2 pathways / pathway-transition pairs" |
| CS2 | five pathways for conserved and three for lost genes against SUN | `paper_numbers.R` | "enriched pairs by transition" and the Table 1 print |
| CS2 | every VIC set a conservation set; ribosome q = 7e-4; none for loss or gain | `paper_numbers.R` | "enriched pairs by transition"; Table 1 print, Ribosome |
| CS2 | loss the largest slice against SUN, 44% to 56% | `paper_numbers.R` | "largest slice in each enriched pathway"; "loss share across enriched pathways" |
| CS2 | conservation the largest slice in every VIC row | `paper_numbers.R` | "largest slice in each enriched pathway", VIC |
| CS2 | 227 of 5,066 conserved | `paper_numbers.R` | "conserved genes PUT-SUN, share of all" |
| CS2 | Parkinson q = 0.027, 0.004, 0.039 | `applications/Baboon_PUT_*.R` | `figures/baboon_PUT_{SUN,VIC}/stage2_significant.csv`, `q` |
| CS2 | 117 measured genes | `paper_numbers.R` | "measured / lost / conserved / gained" |
| CS2 | pathway c-score 0.116 and 0.141, p ≤ 0.001 | `applications/Baboon_PUT_*.R` | `pathway_metrics.csv`, `AdjustedConcordance`, `PValue`; `paper_numbers.R` "Parkinson disease pathway score, p" |
| CS2 | 25 lost, 6 conserved, none gained against SUN | `paper_numbers.R` | PUT_SUN "measured / lost / conserved / gained" |
| CS2 | NDUFS6 +2.1, DDIT3 +2.8, ADRM1 +3.1, NDUFA13 +3.5 | `paper_numbers.R` | PUT_SUN "Phase-shifted (4)" |
| CS2 | HSPA5 -0.4, NDUFB2 +0.1 | `paper_numbers.R` | PUT_SUN "Phase-conserved (2)" |
| CS2 | twelve lost proteasome subunits, as listed | `paper_numbers.R` | PUT_SUN "lost proteasome subunits" |
| CS2 | eight lost respiratory chain subunits | `paper_numbers.R` | PUT_SUN "lost oxidative phosphorylation genes" |
| CS2 | 31 conserved, 14 lost, 4 gained against VIC | `paper_numbers.R` | PUT_VIC "measured / lost / conserved / gained" |
| CS2 | 28 shifted, all positive, +1.7 to +3.8 h | `paper_numbers.R` | PUT_VIC "Phase-shifted (28)"; "shifted range, all positive" (+1.75 to +3.84) |
| CS2 | COX5B +1.3; 2 undetermined | `paper_numbers.R` | PUT_VIC "Phase-conserved (1)"; "Undetermined (2)" |
| CS2 | five of the 14 lost are proteasome subunits | `paper_numbers.R` | PUT_VIC "lost proteasome subunits" |
| CS2 | ADRM1, DDIT3, NDUFA13 shifted in both at nearly the same magnitude | `paper_numbers.R` | "shifted in both, SUN / VIC" |

## Case Study 3

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| CS3 | 4,893 genes | `paper_numbers.R` | Case Study 3: "genes" |
| CS3 | lung the largest conserved set of six against baboon, 19, of which 12 phase-conserved | `mouse/cross_species_mouse.R` | `cross_species_mouse/cross_species_six_tissues.csv`, `maintained`, `phase_conserved`; printed by `paper_numbers.R` |
| CS3 | 1,976 baboon, 731 mouse, 64 (1.3%) human | `paper_numbers.R` | "rhythmic human / baboon / mouse" |
| CS3 | baboon tau_cons = 0.554, 19 conserved, 63.2% within 2 h, 12 / 6 / 1 | `paper_numbers.R` | baboon: "tau_cons, rhythm-conserved"; "within 2 h of human"; "phase conserved / shifted / undetermined" |
| CS3 | mouse tau_cons = 0.455, 12 conserved, 25% within 2 h, 0 / 10 / 2 | `paper_numbers.R` | mouse: the same three lines |
| CS3 | the same 64 human genes | `mouse/mouse_diagnostics.R`; `paper_numbers.R` | "rhythmic genes of 4893", LUN human; "rhythmic human / baboon / mouse" |
| CS3 | Stage 1 retained two of 225: circadian rhythm 23 genes q = 7.3e-6, circadian entrainment 20 genes q = 0.093 | `applications/Baboon_Human_LUN.R` | `figures/baboon_human_LUN/stage1_union.csv`; `paper_numbers.R` "baboon LUN: genes, tested, retained, rank, q" |
| CS3 | mouse q = 5.0e-5 | `applications/Baboon_Human_LUN.R` run for mouse | `figures/mouse_human_LUN_genome/stage1_union.csv`; `paper_numbers.R` "mouse LUN: ..." |
| CS3 | of 12 comparisons, top-ranked in 10, present in 11 | `paper_numbers.R` | the twelve "genes, tested, retained, rank, q" lines |
| CS3 | circadian rhythm shared by both comparators in five of six tissues; none in kidney cortex | `applications/Baboon_Human_LUN.R` per tissue | `figures/{baboon,mouse}_human_<t>_genome/stage2_significant.csv` |
| CS3 | CRY2, DBP, NR1D2, PER2 phase-conserved in baboon; BMAL1 +2.03; NR1D1 +2.04 | `paper_numbers.R` | baboon: "pathway genes phase-conserved", "... phase-shifted", "... phase-undetermined" |
| CS3 | mouse: those five and NR1D1 shifted, none conserved | `paper_numbers.R` | mouse: "pathway genes phase-shifted" |
| CS3 | the other 17 genes | `paper_numbers.R` | "pathway genes conserved in neither comparator" |
| CS3 | human 0.008, 0.010, 0.024; baboon 0.999, 0.881, 0.978 | `paper_numbers.R` | "human, baboon posterior PER1 CUL1 CSNK1D" |
| CS3 | 4 against 12; 23 against 4,893 | `paper_numbers.R` | baboon "pathway genes phase-conserved"; "KEGG Circadian rhythm genes"; "genes" |
| CS3 | offsets positive, 5.7 to 11.7 h | `mouse/mouse_diagnostics.R` | "across all six" |
| CS3 | circular means 9.6, 8.2, 8.5, 8.2 (one gene), 7.0, 10.0 | `mouse/mouse_diagnostics.R` | "mouse minus human peak offset, rhythm-conserved clock genes" |
| CS3 | mouse offsets +10.83, +10.33, +11.69, +9.36, +10.23, +7.52 h | `paper_numbers.R` | mouse: "pathway genes phase-shifted" |
| CS3 | baboon offsets +2.03, +2.04, -1.55, +1.99, +0.37, +0.37 h | `paper_numbers.R` | baboon: the three "pathway genes phase-..." lines |
| CS3 | 1.47 h in lung (5 of the 7 intervals from BMAL1), 1.76 h in liver (4 of 7) | `mouse/mouse_diagnostics.R` | "peak measured from BMAL1 within each species"; the 7 are the other clock genes measured in both species |
| CS3, Discussion | 64 human lung genes clear the baboon conserved threshold, 95 the mouse one | `paper_numbers.R` | "human genes reaching tau_cons", over the 4,893 genes, at `tau_cons` 0.5544 and 0.4548 |

## Discussion

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Discussion | 0.133 and 0.132; 0.55 and 0.99 | as in Case Study 2 | — |
| Discussion | 227 at +2.57 h, 679 at +2.65 h; +0.22 h | `paper_numbers.R` | "mean comparator minus reference (Figure 3E)" |
| Discussion | 708 gene slots, 286 distinct genes | `pipeline/sunput_pathway_overlap.R` | "sum of individual sizes"; "genes in the union of all seven" |
| Discussion | the same 40 OxPhos genes in four sets; 24 of 41 endocannabinoid | `pipeline/sunput_pathway_overlap.R` | "shared measured genes", OxPhos row |
| Discussion | thyroid shares no gene with OxPhos, 15 measured | `pipeline/sunput_pathway_overlap.R` | "shared measured genes"; "measured genes per pathway" |
| Discussion | five of the seven VIC sets on the same block | `paper_numbers.R` | "VIC-enriched sets sharing >= 40 OxPhos genes" |
| Discussion | twelve proteasome subunits against five | `paper_numbers.R` | "lost proteasome subunits", both pairs |
| Discussion | HSPA5 and NDUFB2 lose rhythmicity against VIC | `paper_numbers.R` | "HSPA5, NDUFB2 against VIC" |
| Discussion | 0 in testis to 110 in oesophagus; lung 64, fifth of 26 | `paper_numbers.R` | Discussion: "fewest, most"; "lung calls, rank of 26" |
| Discussion | 19 and 12 conserved | `paper_numbers.R` | Case Study 3 |
| Discussion | strongest of two of 225; 11 of 12; top-ranked in 10 | as in Case Study 3 | — |
| Discussion | 10.0 h in lung; 7.0 to 10.0 h across six | `mouse/mouse_diagnostics.R` | "mouse minus human peak offset ..." |
| Discussion | 1.47 h; 0 phase-conserved of 12; rotation of 10 h | as in Case Study 3 | — |
| Limitations | realised FDR 0.33, power 0.78 at alpha = 0.25, 12 samples | `pipeline/bfdr_calibration.R --summary` | `calibration/bfdr_calibration_summary.csv`, row 0.25 |
| Limitations | 4 distinct phases, twice; 12 in baboon; 1 residual df against 3 parameters | `paper_numbers.R` | "mouse samples, spacing, distinct phases"; 4 phases less 3 parameters |
| Limitations | 64 in lung, at most 110 | `paper_numbers.R` | Discussion lines |

## Figure Captions

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Figure 2 | 26 tissues; 5,066 genes; 23 genes; 25 tissues | `paper_numbers.R` | as in Case Study 1 |
| Figure 2B | k = 2; nine tissues listed; 0.413; sixteen; 0.157 | `pipeline/cluster_rhythmic_content.R` | as in Case Study 1 |
| Figure 2C | call at BFDR 0.25 | `plots/figure2_panelC.R` | `Fig2C_circadian_membership.csv`, `called` |
| Figure 2D | three timing groups G1 to G3 | `plots/figure2/explore_phase_groups.R` | `phase_group_assignments.csv` |
| Figure 3A | R_c = 81; 81.5% within 2 h | `paper_numbers.R` | SCN_HIP lines |
| Figure 3B | R_c = 227; 24.7%; +2.57 h | `paper_numbers.R` | PUT_SUN lines |
| Figure 3C | R_c = 679; 26.5%; +2.65 h; three times as many | `paper_numbers.R` | PUT_VIC lines; "conserved PUT-VIC over PUT-SUN" |
| Figure 3D | counts at BFDR 0.25 | `paper_numbers.R` | "gain / conserved / loss", three pairs |
| Figure 3E | +0.22, +2.57, +2.65 h | `paper_numbers.R` | "mean comparator minus reference (Figure 3E)" |
| Figure 4A | q < 0.05 cut; eight pairs across seven pathways; 15 thyroid genes; seven VIC pathways, all conserved | `paper_numbers.R` | "Stage 2 pathways / pathway-transition pairs"; "enriched pairs by transition"; Table 1 print |
| Figure 5 | 117 genes | `paper_numbers.R` | "measured / lost / conserved / gained" |
| Figure 5A | 25 lost, 6 conserved; 4 shifted positive; HSPA5 -0.4, NDUFB2 +0.1 | `paper_numbers.R` | PUT_SUN Parkinson lines |
| Figure 5B | 31 conserved, 14 lost; 28 of 31 positive; COX5B +1.3 | `paper_numbers.R` | PUT_VIC Parkinson lines |
| Figure 6 | 4,893 genes | `paper_numbers.R` | Case Study 3: "genes" |
| Figure 6A | R_c = 19; 63.2%; BMAL1 +2.03; NR1D1 +2.04 | `paper_numbers.R` | baboon lines |
| Figure 6B | R_c = 12; 25%; about 10 h; circular mean +10.0 h | `paper_numbers.R`; `mouse/mouse_diagnostics.R` | mouse lines; "mouse minus human peak offset", LUN |
| Figure 6C | stronger of the two of 225 retained by Stage 1 | `applications/Baboon_Human_LUN.R` | `figures/baboon_human_LUN/stage1_union.csv` |
| Figure 6C | other 17 genes; CRY1 at 0.474 | `paper_numbers.R` | "pathway genes conserved in neither comparator"; "highest human posterior among those" |

## Table 1

`paper_numbers.R` prints both blocks of Table 1 under "Case Study 2: two-stage
enrichment and Table 1", in the table's row order and at its rounding.

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| Table 1 caption | 8 of 229 pass Stage 1 at FDR < 0.20; 15 measured genes | `paper_numbers.R` | "pathways tested / Stage 1 retained (q < 0.20)" |
| Table 1 caption | five of the seven pathways in both blocks | `applications/Baboon_PUT_*.R` | `stage2_significant.csv`, `pathway`, both pairs |
| Table 1 caption | thyroid the highest adjusted concordance against SUN | `applications/Baboon_PUT_SUN.R` | `pathway_metrics.csv`, `AdjustedConcordance` |
| Table 1 headers | genome-wide adj c 0.133, 0.132 | `genome_concordance_pvalues.R` | `genome_concordance_pvalues.csv` |
| Table 1 headers | genome-wide GLR 0.55, 0.99 | `gain_loss_ratios.R` | `gain_loss_ratios.csv` |
| Table 1, Size | 188 175 117 110 62 41 15; 188 175 117 115 110 69 62 | `applications/Baboon_PUT_*.R` | `stage2_significant.csv`, `size` |
| Table 1, Adj c | 14 cells | `applications/Baboon_PUT_*.R` | `pathway_metrics.csv`, `<pair>_AdjustedConcordance` |
| Table 1, GLR | 14 cells | `gain_loss_ratios.R` | `gain_loss_ratios.csv`, per-pathway rows |
| Table 1, E[Gain], E[Loss], E[Cons] | 42 cells | `applications/Baboon_PUT_*.R` | `stage2_significant.csv`, `Expected_N_Gain`, `Expected_N_Loss`, `Expected_N_Conserved` |
| Table 1, E[Union] | 14 cells | `paper_numbers.R` | Table 1 print, `union`, the sum of the three |
| Table 1, Enrichment | 15 q-values and stars | `applications/Baboon_PUT_*.R` | `stage2_significant.csv`, `direction`, `q` |

## Supplement

| Paper location | Value as printed | Script | Output file / column or printed line |
|---|---|---|---|
| S4 | six most rhythmic genes in each of LUN, PUT, SUN | `plots/Cosinor_residual_diagnostics_LUN_v2_CI95.R` | the figure's panels |
| S4 | 5,066 genes; Shapiro-Wilk 94.2%, 94.2%, 94.1% | `plots/Cosinor_residual_diagnostics_LUN_v2_CI95.R` | `FIG/Cosinor_residual_diagnostics_{LUN,PUT,SUN}.csv`, `shapiro.pval`; `paper_numbers.R` S4 lines |
| S4 | constant variance 94.8%, 94.1%, 96.2% | `plots/Cosinor_residual_diagnostics_LUN_v2_CI95.R` | the same files, `hetero.lm.pval`; `paper_numbers.R` S4 lines |
| S5 | n = 5,066 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "genes compared" |
| S5 | p_rhythmic = 0.2; posterior cuts 0.20, 0.43, 0.56 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "posterior equivalent of each cut" |
| S5 | 3,440 containing all 1,969; 1,744 against 1,969; 56 and 281 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv` |
| S5 | 1,097 and 2 at p < 0.01 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv`, `p < 0.01` row |
| S5 | Benjamini-Hochberg leaves 5 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv`, `BH q < 0.05` |
| Table S5 | 18 cells and 6 totals | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_threshold_grid.csv`; totals in "Bayesian totals" |
| S5 | k = 1 to 1,000; k = 253; 202 of 253 (80%) | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "top 5% cutoff k = 253 of 5066; overlap = 202 (80%)" |
| S5 | Spearman r = 0.974 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | "Spearman r (posterior vs p-value)" |
| S5 | 1,917 genes; Pearson r = 0.9998 | `plots/S5_Bayes_Cosinor_Agreement_LUN.R` | `FIG/S5_detection_counts.csv`, "Both methods"; "Pearson r (phase, n = 1917)" |
| S7 | 708 slots, 286 distinct; 40 OxPhos genes; 24 of 41; thyroid 0, at most six | `pipeline/sunput_pathway_overlap.R` | "sum of individual sizes"; "genes in the union of all seven"; "shared measured genes" |
| Table S7 | 49 cells | `pipeline/sunput_pathway_overlap.R` | "shared measured genes" |
| Table S7 caption | Jaccard 0.64 (89 genes), 0.60 (136 genes) | `pipeline/sunput_pathway_overlap.R` | "Jaccard"; "shared measured genes" |
| S8 | 12 time points; ten replicates of 2,000; 20% rhythmic, 0.199 realised over 20,000 | `pipeline/bfdr_calibration.R` | `calibration/run_record.txt`; `paper_numbers.R` "rhythmic share realised, kept draws" |
| S8 | amplitudes on [0.3, 1.2]; residual SD 0.5; period 24 h; 2,001 retained draws | `pipeline/bfdr_calibration.R` | `calibration/run_record.txt`, `amplitude`, `noise_sd`, `period`, `kept_draws` |
| Table S8 | 8 rows of realised FDR, SE and power | `pipeline/bfdr_calibration.R --summary` | `calibration/bfdr_calibration_summary.csv`, `realised`, `se`, `power` |
| S8 | factor of 1.89 falling to 1.22 | `pipeline/bfdr_calibration.R --summary` | `bfdr_calibration_summary.csv`, `ratio`; `paper_numbers.R` "realised over nominal, strictest to loosest" |
| S8 | 0.334 with power 0.775 at alpha = 0.25 | `pipeline/bfdr_calibration.R --summary` | `bfdr_calibration_summary.csv`, row 0.25 |
| S8 | n = 12, 24, 48; A/sigma 0.5, 1, 1.5, 2, 3; sigma = 1; mesor 5; 20% rhythmic | `pipeline/bfdr_calibration.R` | the sweep file names and the generator in the script |
| S8 | AUCs differ by at most 0.008 | `plots/figure_S_operating.R` | `calibration/bayrc_cosinor_auc.csv`; `paper_numbers.R` "largest \|AUC BayRC - AUC cosinor\|" |
| S8 | power above 0.95 from A/sigma = 1 at n = 48 | `plots/figure_S_operating.R` | `calibration/bfdr_operating_characteristics.csv`, `power`; `paper_numbers.R` "smallest power, n = 48 and A/sigma >= 1" |
| S8 | type I error below nominal in every cell | `plots/figure_S_operating.R` | `bfdr_operating_characteristics.csv`, `type1`; `paper_numbers.R` "largest type I error minus alpha" |
| S8 | FDR within 0.03 of nominal at n = 48 for A/sigma ≤ 1 | `plots/figure_S_operating.R` | `bfdr_operating_characteristics.csv`, `fdr`; `paper_numbers.R` "largest FDR minus alpha, n = 48, A/sigma <= 1" |
| S8 | settles at 0.10, 0.17, 0.33 | `plots/figure_S_operating.R` | `bfdr_operating_characteristics.csv`; `paper_numbers.R` "FDR at n = 48, A/sigma >= 1.5, alpha ..." |
| S9 | 8 timepoints 6 h apart; 4 distinct phases; 8 pooled samples | `paper_numbers.R` | "mouse samples, spacing, distinct phases" |
| S9 | 4 samples, 3 parameters, 1 residual df | arithmetic on the design | — |
| S9 | pooled 15-39%, cycle 1 30-76%, cycle 2 48-69% | `mouse/mouse_diagnostics.R` | "range across tissues" |
| S9 | arms within 12% in each tissue; aorta clock amplitude 1.658, 1.659, 1.657 | `mouse/mouse_diagnostics.R` | "inclusion and core-clock amplitude, by arm", AOR `clock A` |
| S9 | pooled clock posterior highest in 5 of 6 | `paper_numbers.R` | "tissues where pooled clock posterior is highest" |
| S9 | median absolute difference 0.52 to 1.09 h | `paper_numbers.R` | "median \|cycle 1 - cycle 2\| clock peak, range" |
