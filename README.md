# delta-rhythm-window-test

A pre-registered, calibrated test of the claim that animal acoustic communication evolved
toward a conserved optimal rhythm of about 2.7 Hz in the neural delta range (Piette et al.
2026, PLoS Biology 24(6): e3003798).

The question: the paper's rhythm measurement can only return values inside a window (calls
must be less than 2 s apart; the envelope is low-passed at 20 Hz). Brownian motion between
hard bounds is known to look like Ornstein-Uhlenbeck stabilising selection. Once the method's
own window is built into a no-optimum model, is an optimum still needed?

## Result (graded 2026-11-27): OPTIMUM_SURVIVES_WINDOW

Yes, it is. On 82 bird and mammal species, a Brownian model bounded at the method's window
fits no better than plain Brownian motion (AIC 28.65 against 28.41), and the OU model's
advantage over it (6.93 log-likelihood units) is far beyond anything the window alone produced
(95th percentile 2.82; p 0.0025 over 400 simulations; power 0.995). Every sensitivity window
and the all-taxa analysis agree; the OU optimum is 3.12 Hz. Our registered hypothesis, that the
window manufactures the optimum, was not supported. See [READING-RHYTHM.md](READING-RHYTHM.md).

## Files

- `PREREG.md`: the plan, verdict rules and 16 forecasts, fixed before any data was opened.
- `rhythm_test.py`: downloads the authors' public data and tree at a pinned commit, fits BM,
  OU (the authors' kernel) and window-bounded BM (Boucher & Demery 2016) by maximum
  likelihood with within-species error, and calibrates by parametric bootstrap.
- `test_rhythm_test.py`: unit tests (including an analytic check of the bounded-BM
  likelihood) and an end-to-end run on synthetic data.
- `READING-RHYTHM.md`: the graded reading; `test_reading_rhythm.py` checks it against
  `results/rhythm/results.json` and `results/rhythm/grade.json` (stdlib only).
- `results/rhythm/`: the job's outputs (no data of the authors).

Run: `python -m unittest -v test_rhythm_test` then
`python rhythm_test.py --data-dir data --out-dir out` (needs numpy and network access).

The authors' data and code are not redistributed here; they are fetched from
github.com/chundrac/phylo-acoustic-rhythm (commit 656fc88) and Zenodo 10.5281/zenodo.19816728.

Companion test: [animal-tempo-window-test](https://github.com/munimtechnologies/animal-tempo-window-test)
(Amichay et al. 2026). Ledger of all of MV2's calibrated re-tests:
[calibrated-animal-claims](https://github.com/munimtechnologies/calibrated-animal-claims).

MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human
expert has checked this work yet. Code Apache-2.0; text CC BY 4.0.
