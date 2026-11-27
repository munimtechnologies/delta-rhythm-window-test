# READING-RHYTHM: the 'conserved optimal rhythm' survives the measurement window (graded 2026-11-27)

MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human
expert has checked this reading. Text CC BY 4.0, code Apache-2.0.

Graded from job `delta-rhythm-v1` (run on 2026-11-27 16:32-17:36 UTC, python:3.12, seed
20261125), exactly as registered in [PREREG.md](PREREG.md) on 2026-11-25 before any data file
was opened. No amendment was needed: the loader matched on the first run. Raw outputs are in
`results/rhythm/` (`results.json`, `grade.json`, `RESULTS-RHYTHM.md`, `load_report.json`,
`download_manifest.json`, `log.txt`). The authors' data are not redistributed; the manifest
gives the SHA-256 of every file used (GitHub chundrac/phylo-acoustic-rhythm at commit 656fc88,
Zenodo 10.5281/zenodo.19816728).

## Headline: OPTIMUM_SURVIVES_WINDOW

Our registered hypothesis was that the 2.7 Hz Ornstein-Uhlenbeck (OU) optimum reported by
Piette et al. 2026 (PLoS Biology 24(6): e3003798) is what Brownian evolution looks like through
the rhythm method's own window (calls less than 2 s apart, envelope low-passed at 20 Hz). **It
was wrong.** A Brownian model with reflecting bounds at that window (Boucher and Demery 2016)
explains the data no better than plain Brownian motion, and the OU model beats it by far more
than the window alone ever produces: in none of 400 window-only data sets did the OU advantage
reach the real one (p 0.0025, the smallest value 400 simulations can give), and the test had
power 0.995. All three sensitivity windows and the all-taxa analysis give the same verdict.

## Numbers (birds and mammals, primary window 0.5-20 Hz)

| quantity | value |
|---|---|
| rows matched (of 718 in the authors' bird and mammal orders) | 694 |
| species matched to the MCC tree | 82 |
| within-species variance of log10 rate (pooled) | 0.02348 |
| species geometric mean rate | 3.115 Hz (2.5-97.5%: 1.026-9.839 Hz) |
| SD of species log10 rate / uniform-window SD | 0.2637 / 0.4625 |
| AIC: BM / OU / window-bounded BM | 28.4107 / 18.7897 / 28.6538 |
| OU optimum | 3.123 Hz |
| OU half-life | 0.0910 tree heights |
| logL(OU) - logL(BM); p under the fitted BM (400 sims) | 5.8105; 0.0025 |
| logL(OU) - logL(window BM); 95th percentile under the window null | 6.9320; 2.8186 |
| p under the fitted window BM (400 sims) | 0.0025 |
| power against the fitted OU (200 sims) | 0.995 |
| window-only data sets in which OU beats BM by AIC | 146 of 400 (0.365) |
| rows outside 0.5-20 Hz | 0.86% |

Sensitivity windows (150 null and 75 power simulations each):

| window | logL(OU) - logL(window BM) | null 95th pct | p | power | verdict |
|---|---|---|---|---|---|
| 0.4-20 Hz | 7.0033 | 2.2836 | 0.0066 | 1.000 | OPTIMUM_SURVIVES_WINDOW |
| 0.5-10 Hz (widened to 11.43 Hz: at least one species mean lies above 10 Hz) | 5.4243 | 1.5530 | 0.0066 | 0.987 | OPTIMUM_SURVIVES_WINDOW |
| 0.33-30 Hz | 7.1914 | 3.5565 | 0.0066 | 0.973 | OPTIMUM_SURVIVES_WINDOW |

All taxa (92 species, adding insects, frogs, a fish and a reptile): OU optimum 3.174 Hz,
logL(OU) - logL(window BM) 9.7667 against a null 95th percentile of 2.5402, p 0.0050 (200
sims), power 1.0; verdict OPTIMUM_SURVIVES_WINDOW. Here OU beat BM by AIC in 0.5 of
window-only data sets.

Arm A (cycles per sequence; rate against sequence length) was UNTESTABLE_NO_LENGTH: neither
the main file nor the Zenodo archive has a per-row sequence-length column our registered
loader recognises, so P10 and P11 are not graded.

## What this means

- The concentration of call rhythms near 3 Hz is not produced by the method's window. Rhythms
  are much more tightly clustered (SD 0.26 log10 units) than Brownian evolution inside the
  window allows, and the window model is not even preferred over unbounded Brownian motion.
- The published model choice (OU over BM) is reached by the window alone fairly often, in 36.5%
  of window-only worlds by AIC (our stand-in for the authors' LOO stacking), so on its own it
  would be weak evidence. The real OU advantage, however, is about 2.5 times the window null's 95th
  percentile, so the conclusion does not rest on that weak step.
- We reproduce an optimum at 3.12 Hz (birds and mammals) with maximum likelihood on species
  means, against the paper's 2.7 Hz (95% interval 0.45-4.99 Hz) from Bayesian fits to every
  recording on 50 trees. Both are inside the 0.5-4 Hz delta band.
- Together with the companion test of Amichay, Balasubramanian and Abrams 2026
  ([animal-tempo-window-test](https://github.com/munimtechnologies/animal-tempo-window-test/blob/main/READING-TEMPO.md),
  PREFERRED_TEMPO_BEYOND_WINDOW, p 0.0118), both 2026 PLoS Biology claims of a slow, conserved
  animal tempo pass a window test that we built expecting them to fail.

## What this does not show

- That the cause is stabilising selection tied to delta-band brain rhythms. An OU fit with a
  half-life of 0.09 tree heights means rhythm carries little phylogenetic signal; selection
  toward one optimum is one reading, convergent pressures or measurement choices shared across
  species (for example which sound unit is taken as the rhythm) are others. We tested only the
  window.
- Anything about species selection: which species and recordings were chosen is not modelled.
- That our ML-on-species-means design equals the authors' Bayesian design; the deviation is
  stated in PREREG.md.

## Forecasts (registered priors) and Brier score

| id | forecast | prior | outcome | Brier |
|---|---|---|---|---|
| P1 | at least 60 bird and mammal species match the tree | 0.80 | true | 0.0400 |
| P2 | OU beats BM by AIC | 0.80 | true | 0.0400 |
| P3 | OU over BM survives the BM bootstrap | 0.65 | true | 0.1225 |
| P4 | window BBM beats BM by AIC | 0.55 | false | 0.3025 |
| P5 | OU advantage over window BBM not significant | 0.45 | false | 0.2025 |
| P6 | power at least 0.5 | 0.50 | true | 0.2500 |
| P7 | OU optimum in 2.0-3.6 Hz | 0.75 | true | 0.0625 |
| P8 | OU half-life below 1 tree height | 0.75 | true | 0.0625 |
| P9 | fewer than 2% of rows outside 0.5-20 Hz | 0.70 | true | 0.0900 |
| P10 | median cycles per sequence below 10 | 0.50 | not graded | - |
| P11 | within-species slope negative | 0.55 | not graded | - |
| P12 | all-taxa verdict equals primary | 0.70 | true | 0.0900 |
| P13 | all three sensitivity windows agree | 0.45 | true | 0.3025 |
| P14 | species geometric mean in 2.025-3.375 Hz | 0.75 | true | 0.0625 |
| P15 | OU beats BM by AIC in at least half of window-only sims | 0.50 | false | 0.2500 |
| P16 | species SD below half the uniform-window SD | 0.45 | false | 0.2025 |

Mean Brier over the 14 graded forecasts: **0.1486** (a coin-flip forecaster scores 0.25). Our
two central forecasts about the window (P4, P5) were both wrong.

## Prior work (checked 2026-11-27)

Searches on 2026-11-27 for critiques or reanalyses of Piette et al. 2026 (bounded Brownian
motion, measurement window, OU bias) found the paper, its correction (e3003980, author-name
format only), the supportive PLoS Biology commentary (e3003814), the 2024 bioRxiv preprint
(2.9 Hz), and general methods work: Boucher and Demery 2016 (Systematic Biology 65:651) on
bounded evolution, a 2026 bioRxiv preprint on discretised bounded Brownian models
(10.64898/2026.04.20.719455), and a 2026 bioRxiv anuran study (10.64898/2026.07.20.739646)
noting that omitting bounded BM can inflate support for adaptive models. None tests the rhythm
optimum against its measurement window. This is, as far as we found, the first such test.
