# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Checks that READING-RHYTHM.md matches results/rhythm/*.json (stdlib only).

Run: python -m unittest -v test_reading_rhythm
"""
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results", "rhythm")


def _load(name):
    with open(os.path.join(RES, name), encoding="utf-8") as f:
        return json.load(f)


def _reading():
    with open(os.path.join(HERE, "READING-RHYTHM.md"), encoding="utf-8") as f:
        return f.read()


class TestReadingRhythm(unittest.TestCase):
    def setUp(self):
        self.r = _load("results.json")
        self.g = _load("grade.json")
        self.text = _reading()
        self.bm = self.r["birds_mammals"]
        self.win = {w["name"]: w for w in self.bm["windows"]}

    def test_headline(self):
        self.assertEqual(self.r["headline"], "OPTIMUM_SURVIVES_WINDOW")
        self.assertIn("## Headline: OPTIMUM_SURVIVES_WINDOW", self.text)
        for w in self.bm["windows"]:
            self.assertEqual(w["verdict"], "OPTIMUM_SURVIVES_WINDOW")
        self.assertEqual(self.r["all_taxa"]["windows"][0]["verdict"], "OPTIMUM_SURVIVES_WINDOW")

    def test_primary_numbers(self):
        p = self.win["primary_0.5_20"]
        self.assertEqual(self.bm["n_species"], 82)
        self.assertEqual(self.r["load"]["birds_mammals"]["rows_matched"], 694)
        self.assertEqual(self.r["load"]["birds_mammals"]["rows_in_orders"], 718)
        checks = {
            "28.4107": self.bm["aic_BM"], "18.7897": self.bm["aic_OU"],
            "28.6538": p["aic_BBM"], "5.8105": self.bm["stat_ll_OU_minus_BM"],
            "6.9320": p["stat_ll_OU_minus_BBM"], "2.8186": p["null_BBM_q95"],
            "0.2637": self.bm["species_sd_log10"], "0.4625": p["uniform_window_sd_log10"],
            "0.0910": self.bm["OU_half_life_tree_heights"],
        }
        for s, v in checks.items():
            self.assertEqual(s, "%.4f" % v, s)
            self.assertIn(s, self.text)
        self.assertEqual("%.3f" % self.bm["OU_theta_hz"], "3.123")
        self.assertEqual("%.3f" % self.bm["species_geomean_hz"], "3.115")
        self.assertEqual("%.4f" % p["p_BBM_null"], "0.0025")
        self.assertEqual("%.4f" % self.bm["p_BM_null"], "0.0025")
        self.assertEqual(p["n_null"], 400)
        self.assertAlmostEqual(p["power_vs_OU"], 0.995)
        self.assertIn("146 of 400 (0.365)", self.text)
        self.assertEqual(round(p["share_BBM_null_sims_OU_beats_BM_by_AIC"] * p["n_null"]), 146)
        self.assertEqual("%.2f%%" % (100 * self.r["rows_share_outside_primary_window"]), "0.86%")

    def test_sensitivity_windows(self):
        rows = {"low_0.4_20": ("7.0033", "2.2836"), "high_0.5_10": ("5.4243", "1.5530"),
                "wide_0.33_30": ("7.1914", "3.5565")}
        for name, (stat, q95) in rows.items():
            w = self.win[name]
            self.assertEqual("%.4f" % w["stat_ll_OU_minus_BBM"], stat)
            self.assertEqual("%.4f" % w["null_BBM_q95"], q95)
            self.assertEqual("%.4f" % w["p_BBM_null"], "0.0066")
            self.assertIn(stat, self.text)
            self.assertIn(q95, self.text)
        self.assertTrue(self.win["high_0.5_10"]["widened"])
        self.assertEqual("%.2f" % self.win["high_0.5_10"]["upper_hz"], "11.43")

    def test_all_taxa(self):
        a = self.r["all_taxa"]
        w = a["windows"][0]
        self.assertEqual(a["n_species"], 92)
        self.assertEqual("%.3f" % a["OU_theta_hz"], "3.174")
        self.assertEqual("%.4f" % w["stat_ll_OU_minus_BBM"], "9.7667")
        self.assertEqual("%.4f" % w["null_BBM_q95"], "2.5402")
        self.assertEqual("%.4f" % w["p_BBM_null"], "0.0050")
        for s in ("9.7667", "2.5402", "3.174 Hz"):
            self.assertIn(s, self.text)

    def test_brier(self):
        graded = [x for x in self.g if x["brier"] is not None]
        self.assertEqual(len(graded), 14)
        mean = sum(x["brier"] for x in graded) / len(graded)
        self.assertEqual("%.4f" % mean, "0.1486")
        self.assertIn("**0.1486**", self.text)
        outcome = {x["id"]: x["outcome"] for x in self.g}
        self.assertFalse(outcome["P4"])
        self.assertFalse(outcome["P5"])
        self.assertIsNone(outcome["P10"])
        self.assertIsNone(outcome["P11"])
        self.assertEqual(self.r["arm_a"]["status"], "UNTESTABLE_NO_LENGTH")


if __name__ == "__main__":
    unittest.main()
