# delta-rhythm-window-test

A pre-registered, calibrated test of the claim that animal acoustic communication evolved
toward a conserved optimal rhythm of about 2.7 Hz in the neural delta range (Piette et al.
2026, PLoS Biology 24(6): e3003798).

The question: the paper's rhythm measurement can only return values inside a window (calls
must be less than 2 s apart; the envelope is low-passed at 20 Hz). Brownian motion between
hard bounds is known to look like Ornstein-Uhlenbeck stabilising selection. Once the method's
own window is built into a no-optimum model, is an optimum still needed?

- `PREREG.md`: the plan, verdict rules and 16 forecasts, fixed before any data was opened.
- `rhythm_test.py`: downloads the authors' public data and tree at a pinned commit, fits BM,
  OU (the authors' kernel) and window-bounded BM (Boucher & Demery 2016) by maximum
  likelihood with within-species error, and calibrates by parametric bootstrap.
- `test_rhythm_test.py`: unit tests (including an analytic check of the bounded-BM
  likelihood) and an end-to-end run on synthetic data.

Run: `python -m unittest -v test_rhythm_test` then
`python rhythm_test.py --data-dir data --out-dir out` (needs numpy and network access).

Results will be added under `results/` when the job has run, with a reading note.

The authors' data and code are not redistributed here; they are fetched from
github.com/chundrac/phylo-acoustic-rhythm (commit 656fc88) and Zenodo 10.5281/zenodo.19816728.

MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in MV2's name; no human
expert has checked this work yet. Code Apache-2.0; text CC BY 4.0.
