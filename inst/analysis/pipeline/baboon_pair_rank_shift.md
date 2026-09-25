# Stability of the Within-Baboon Pair Ranking

Figure 6 ranks all `choose(26, 2) = 325` within-baboon tissue pairs by adjusted
concordance. Across two runs of the same pipeline that ranking is not stable,
and the instability tracks a property of the metric rather than a change in the
underlying data. This note records the numbers, the mechanism, and what it
means for reading the figure.

Produced by `baboon_pair_rank_shift.R`; the per-pair table is
`baboon_pair_rank_shift.csv`.

---

## What Moves

Over all 325 pairs the two runs agree closely on the value of the adjusted
concordance and much less closely on the order:

| | |
|---|---|
| Pearson, adjusted concordance | 0.836 |
| Spearman, adjusted concordance | 0.877 |
| median absolute rank shift | 30 |
| maximum absolute rank shift | 145 |
| pairs moving fewer than 10 ranks | 79 of 325 (24%) |
| pairs moving more than 50 ranks | 78 of 325 |

The top of the ranking turns over almost completely. Of the ten most
concordant pairs, **two are shared between runs**. Eight leave and eight
enter:

**Leaving the top 10:** Heart-Thyroid (rank 1 to 41), Substantia
nigra-Cortex (2 to 26), Colon-Spleen (4 to 13), Spleen-Substantia nigra
(6 to 18), Heart-Ileum (7 to 47), Omental fat-Stomach (8 to 22),
Ileum-Stomach (9 to 28), Frontal cortex-Substantia nigra (10 to 33).

**Entering the top 10:** Amygdala-Pituitary (41 to 1), Amygdala-Aorta
(110 to 2), Amygdala-Cerebellum (84 to 3), Amygdala-Oesophagus (107 to 4),
Aorta-Pituitary (37 to 6), Hippocampus-Pituitary (28 to 8),
Aorta-Oesophagus (20 to 9), Amygdala-Hippocampus (101 to 10).

Every entering pair contains Amygdala or Pituitary.

---

## The Mechanism

The obvious reading is that pairs move because one of their two tissues
changed more than the others. That is not what the data show. Taking each
tissue in turn and averaging the rank shift over the 25 pairs it belongs to:

| Predictor of a tissue's mean rank shift | Pearson | Spearman |
|---|---|---|
| that tissue's rhythmic-gene count | **+0.845** | +0.832 |
| how much that tissue's count changed between runs | +0.173 | +0.228 |

The count itself predicts the movement; the change in the count does not.
Tissues with few rhythmic genes rise and tissues with many fall. Amygdala
(53 rhythmic genes) and Aorta (504) carry the largest negative mean shifts,
-44.9 and -50.0; Heart (3,372) and Thyroid (3,122) the largest positive,
+46.4 and +40.5.

The reason is that the adjusted concordance retains a dependence on how many
rhythmic genes the two tissues have, and the sign of that dependence is not
constant. Correlating each pair's adjusted concordance against the geometric
mean of its two rhythmic-gene counts:

| | previous run | current run |
|---|---|---|
| within baboon | **+0.231** | **-0.183** |
| within human | -0.693 | -0.638 |

In baboon the dependence **changes sign** between runs. In human it does not,
and the human ranking is correspondingly stable: the ten most concordant
human pairs are the same ten in both runs, in nearly the same order
(Spearman 0.975).

The two species sit in different regimes. Baboon tissues are rhythmic across
8% to 67% of the 5,066 genes tested, and their counts move by about 9%
between runs. Human tissues are rhythmic across 0.02% to 6%, and their counts
move by about 47%. Far from the sparse limit the residual size-dependence of
the adjusted concordance is weak, so a modest shift in counts is enough to
reverse it; in the sparse human regime it is strong and negative and stays
that way under a much larger shift.

---

## How Figure 6 Should Be Read

The adjusted concordance values are reproducible: they correlate at 0.836
across runs, and the median within-baboon value is essentially unchanged
(0.1042 against 0.1118). What is not reproducible is their order near the
top.

Three consequences follow.

The genome-wide structure of the heatmap is a real result and can be read as
one. The overall level of concordance, the spread across pairs, and the
contrast between within-species and cross-species concordance are all stable.

Individual top-ranked pairs are not a stable quantity and should not carry
interpretation on their own. A sentence naming a specific pair as the most
concordant in the baboon panel describes that run, not the biology, and any
text that names baboon pairs needs re-reading against the current table.

The ordering near the top is most fragile for the sparsest tissues, because
those sit where the metric's size-dependence is weakest and most easily
reversed. Amygdala, with 53 rhythmic genes, is the clearest case: it appears
in half the pairs entering the top 10.

A concordance statistic whose dependence on gene-set size changes sign between
two runs of the same pipeline is measuring something other than shared
rhythmicity alone. Pair rankings from this metric carry that dependence with
them, and a size-matched null, or a statistic with the dependence removed,
would be needed before a ranking is interpreted directly.
