# PREREG: is the "conserved optimal rhythm" of animal calls more than the window the method can see? (registered 2026-11-25)

MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human
expert has reviewed this plan. Text CC BY 4.0, code Apache-2.0. The data and code analysed
belong to their authors and are downloaded at run time from their public sources (GitHub
chundrac/phylo-acoustic-rhythm at commit 656fc883f231e0a00f6ae5e15fbeb22bf2960f8e, and
Zenodo 10.5281/zenodo.19816728, CC BY 4.0). Nothing of theirs is redistributed here.

Registered 2026-11-25 (about 06:00 UTC), before any data file of the study was opened. What was
read before registering: the paper's text and open peer review, and the authors' analysis
scripts (`run_model.R`, the Stan model, the wavelet, length and control scripts), to learn
file names, column names and the order filter. No data row, tree or result file was opened.

## The claim

Piette, Cathcart, Barbieri, Ming, Grandjean, Bickel, Deaux & Giraud 2026, "Animal acoustic
communication has a conserved optimal rhythm within the neural delta range", PLoS Biology
24(6): e3003798 (doi 10.1371/journal.pbio.3003798). From 2,000+ recordings of 98 species they
report that call rhythm evolved under an Ornstein-Uhlenbeck (OU) process toward an optimum of
2.7 Hz (95% interval 0.45-4.99 Hz), better than Brownian motion (BM), and read it as
stabilising selection tied to delta-band brain oscillations. A companion PLoS Biology primer
(Xu & Kalra 2026, e3003814) and a second paper (Amichay, Balasubramanian & Abrams 2026,
e3003735) frame a 0.5-4 Hz "universal tempo" across animals.

## Why it might be the window

The rhythm measurement can only return values inside a window, by construction:

1. **Lower edge.** A sequence is "more than two calls separated by less than two seconds of
   silence" and lasts more than one second. A species whose calls come more than ~2 s apart
   cannot form a sequence, so rates below about 1/(2 s + call length), roughly 0.33-0.5 Hz,
   are excluded. The paper's lower bound is 0.45 Hz.
2. **Upper edge.** The amplitude envelope is low-passed at 20 Hz (4th-order Butterworth)
   before the wavelet, so faster modulations (insect pulses, bird trills) are suppressed and
   the slower phrase rhythm is picked instead.
3. **Wavelet.** `upperPeriod = floor(nrow/3)` means the slowest analysable rhythm is ~3/T Hz
   for a sequence of T seconds.

The geometric middle of 0.5-20 Hz is 3.2 Hz. Bounded Brownian motion is known to be
mistaken for OU (Boucher & Demery 2016, Syst Biol 65:651), and OU is known to be favoured
spuriously with small trees and measurement error (Cooper et al. 2016, Biol J Linn Soc
118:64). Neither paper's analysis, nor its peer review, tests either issue.

**Question.** Once the method's own window is built into a no-optimum model, is an OU
optimum still needed to explain the rhythms of birds and mammals?

## Data and loader (fixed now)

- `Data_Rhythm_Final_English.csv` and `Coding_mtDNA_HKY_Yule_simple.tree` (the MCC tree made
  by the authors' `make_MCC_tree.R`) from the pinned commit; `run_model.R` saved for the
  record; Zenodo `Data.zip` (used only to find a sequence-length column if the main file has
  none). SHA-256 of every file is written to `download_manifest.json`.
- Separator detected from the header (`;` expected); comma decimals accepted. If the header
  is one field shorter than the rows (R row names), the first field is treated as row names.
- Rows exactly as the authors' `run_model.R`: Order in their list of 12 mammal and 34 bird
  orders (their spellings, exact match); rate = `Log_Frequency` (log10 Hz; if absent,
  log10 of `Frequency`); rows with missing rate dropped; taxon = cleaned (lower case, spaces
  and underscores removed) R column 8 if it is a tree tip, else column 7; their five manual
  name fixes; rows kept only if the taxon is a tip.
- Species value y = mean of log10 rate over its rows; within-species variance tau^2 pooled
  over species with 2+ rows; measurement-error variance tau^2/n per species, used in every
  model.
- Tree: MCC tree, depths measured from the kept species' common ancestor, scaled to
  root-to-tip height 1.

## Models (maximum likelihood, same data, same error)

- **BM**: y ~ N(mu, sig2 C + E). 2 parameters.
- **OU**: the authors' kernel, y ~ N(theta, s2 exp(-a D/maxD) + E). 3 parameters.
- **BBM (window)**: Brownian motion with reflecting bounds at the method's window, fixed in
  advance, uniform (stationary) root, likelihood on a 100-bin grid (Boucher & Demery 2016).
  1 parameter (sig2). Primary window 0.5-20 Hz. Sensitivity windows: 0.4-20, 0.5-10 and
  0.33-30 Hz. If any species mean lies outside a window, that window is widened to the data
  range plus 0.02 log10 units and flagged.

## Statistics

- AIC of each model.
- p_BM: parametric bootstrap (400 data sets from the fitted BM) of logL(OU) - logL(BM).
- **Primary p_BBM**: parametric bootstrap (400 from the fitted primary-window BBM) of
  logL(OU) - logL(BBM). Small p means the data need more than the window.
- Power: share of 200 data sets simulated from the fitted OU whose statistic exceeds the
  95th percentile of the BBM null.
- Window-only check: share of the BBM-null data sets in which OU beats BM by AIC (how often
  the window alone reproduces the published model choice).
- Sensitivity windows: 150 null and 75 power simulations each. All-taxa analysis (every
  matched species), primary window, 200/200/100 simulations.
- Arm A (descriptive): if a sequence-length column exists, cycles per sequence (rate x
  length), and the within-species slope of log rate on log length (species-cluster
  bootstrap CI).
- Seed 20261125.

## Verdict (birds and mammals, primary window), in this order

1. `LOADER_STOP`: data or tree cannot be read, or fewer than 40 species match.
2. `OU_NOT_PREFERRED`: AIC(OU) >= AIC(BM). The published model choice is not reproduced here.
3. `OPTIMUM_SURVIVES_WINDOW`: p_BBM < 0.05.
4. `WINDOW_EXPLAINS_OPTIMUM`: p_BBM >= 0.05 and power >= 0.5.
5. `UNDECIDED_LOW_POWER`: otherwise.

If the loader stops because columns or names differ from what the scripts imply, a dated
AMENDMENT-1 will fix the loader only, written after reading `load_report.json` and before
reading any model output, and the job will be re-run under a new version.

## Forecasts (probability that each comes true)

| id | forecast | p |
|---|---|---|
| P1 | at least 60 bird and mammal species match the tree | 0.80 |
| P2 | OU beats BM by AIC | 0.80 |
| P3 | OU over BM survives the BM bootstrap (p < 0.05) | 0.65 |
| P4 | window BBM (primary) beats BM by AIC | 0.55 |
| P5 | OU advantage over window BBM is not significant (p >= 0.05) | 0.45 |
| P6 | power of that test is at least 0.5 | 0.50 |
| P7 | fitted OU optimum in 2.0-3.6 Hz | 0.75 |
| P8 | OU half-life below 1 tree height | 0.75 |
| P9 | fewer than 2% of rows outside 0.5-20 Hz | 0.70 |
| P10 | median cycles per sequence below 10 | 0.50 |
| P11 | within-species slope of log rate on log length negative, CI excludes 0 | 0.55 |
| P12 | all-taxa verdict equals the birds-and-mammals verdict | 0.70 |
| P13 | all three sensitivity windows give the primary verdict | 0.45 |
| P14 | species geometric mean rate in 2.025-3.375 Hz | 0.75 |
| P15 | in at least half of window-only simulations OU beats BM by AIC | 0.50 |
| P16 | species SD of log10 rate below half the uniform-window SD | 0.45 |

## Deviations from the paper, stated now

The authors fit Bayesian models to every recording on 50 posterior trees and compare them
by LOO stacking; we use maximum likelihood on species means with the same within-species
error, on the MCC tree, and calibrate by simulation. We test only the BM/OU contrast and the
window, not the regressions on body size or habitat. The window edges are our reading of the
paper's methods; the sensitivity windows exist because they are approximate.

## What each outcome would mean

- `WINDOW_EXPLAINS_OPTIMUM`: the data are what Brownian evolution looks like through this
  method's window; the 2.7 Hz "optimum" and the OU preference need not reflect selection.
- `OPTIMUM_SURVIVES_WINDOW`: rhythms are more concentrated than the window alone produces;
  the conserved-rhythm reading passes a test it had not faced.
- `UNDECIDED_LOW_POWER`: 98 species cannot tell the two apart; the claim is not established
  by this design either way.
