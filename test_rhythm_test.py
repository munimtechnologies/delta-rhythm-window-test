# Copyright 2026 MV2, Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Tests for rhythm_test.py (numpy + stdlib only). Made with AI agents."""
import json
import math
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rhythm_test as rt  # noqa: E402

TREE3 = "((A:1,B:1):1,C:2);"


def tips_of(tr):
    return [k for k in range(len(tr["names"])) if not tr["kids"][k]]


class TreeTests(unittest.TestCase):
    def test_newick(self):
        tr = rt.parse_newick(TREE3)
        tips = tips_of(tr)
        self.assertEqual(sorted(tr["names"][k] for k in tips), ["A", "B", "C"])
        d = rt.tree_depths(tr)
        for k in tips:
            self.assertAlmostEqual(d[k], 2.0)

    def test_nexus_translate(self):
        text = ("#NEXUS\n[R-package APE]\nBEGIN TAXA;\n\tTAXLABELS\n\t\t'Tree swallow'\n\t;\nEND;\n"
                "BEGIN TREES;\n\tTRANSLATE\n\t\t1\tPan_troglodytes,\n"
                "\t\t2\t'Homo sapiens',\n\t\t3\tCanis_lupus\n\t;\n"
                "\tTREE * UNTITLED = [&R] ((1:1,2:1):1,3:2);\nEND;\n")
        tr = rt.read_tree(text)
        names = sorted(tr["names"][k] for k in tips_of(tr))
        self.assertEqual(names, ["Canis_lupus", "Homo sapiens", "Pan_troglodytes"])

    def test_matrices(self):
        tr = rt.parse_newick(TREE3)
        tips = tips_of(tr)          # A, B, C in creation order
        C, Dn, H, dmax_h = rt.tip_matrices(tr, tips)
        self.assertAlmostEqual(H, 2.0)
        self.assertAlmostEqual(C[0, 0], 1.0)
        self.assertAlmostEqual(C[0, 1], 0.5)
        self.assertAlmostEqual(C[0, 2], 0.0)
        self.assertAlmostEqual(Dn[0, 1], 0.5)
        self.assertAlmostEqual(Dn[0, 2], 1.0)
        self.assertAlmostEqual(dmax_h, 2.0)


class ModelTests(unittest.TestCase):
    def test_gauss_ll(self):
        tr = rt.parse_newick(TREE3)
        C = rt.tip_matrices(tr, tips_of(tr))[0]
        y = np.array([0.1, -0.2, 0.3])
        V = 0.5 * C + 0.1 * np.eye(3)
        ll, mu = rt.gauss_ll(y, V)
        Vi = np.linalg.inv(V)
        one = np.ones(3)
        mu2 = (one @ Vi @ y) / (one @ Vi @ one)
        r = y - mu2
        ll2 = -0.5 * (3 * math.log(2 * math.pi) + np.linalg.slogdet(V)[1] + r @ Vi @ r)
        self.assertAlmostEqual(mu, mu2, places=10)
        self.assertAlmostEqual(ll, ll2, places=10)

    def test_reflect(self):
        self.assertAlmostEqual(rt.reflect(-0.1, 0.0, 1.0), 0.1)
        self.assertAlmostEqual(rt.reflect(1.1, 0.0, 1.0), 0.9)
        self.assertAlmostEqual(rt.reflect(2.1, 0.0, 1.0), 0.1)
        self.assertAlmostEqual(rt.reflect(0.5, 0.0, 1.0), 0.5)

    def test_bbm_uniform_limit(self):
        tr = rt.parse_newick(TREE3)
        bbm = rt.BBM(tr, tips_of(tr), 2.0, -0.3, 1.3, nb=100)
        ll = bbm.loglik(np.array([0.5, 0.4, 0.6]), [1e-4] * 3, 1e6)
        self.assertAlmostEqual(ll, 3 * math.log(1 / 1.6), delta=0.02)

    def test_bbm_matches_bm_far_from_bounds(self):
        tr = rt.parse_newick("((A:0.5,B:0.5):0.5,C:1.0);")
        tips = tips_of(tr)
        bbm = rt.BBM(tr, tips, 1.0, -3.0, 3.0, nb=300)
        y = np.array([0.1, -0.2, 0.3])
        me = np.array([0.01] * 3)
        sig2 = 0.05
        ll = bbm.loglik(y, me, sig2)
        C = rt.tip_matrices(tr, tips)[0]
        V = sig2 * C + np.diag(me + bbm.h ** 2 / 12.0)
        llbm, _ = rt.gauss_ll(y, V)
        one = np.ones(3)
        prec = one @ np.linalg.inv(V) @ one
        expect = llbm + 0.5 * math.log(2 * math.pi / prec) - math.log(6.0)
        self.assertAlmostEqual(ll, expect, delta=0.02)

    def test_bbm_simulate_inside(self):
        tr = rt.parse_newick(TREE3)
        bbm = rt.BBM(tr, tips_of(tr), 2.0, 0.0, 1.0, nb=50)
        sims = bbm.simulate(np.random.default_rng(1), 50.0, [0.0] * 3, 20)
        self.assertEqual(sims.shape, (20, 3))
        self.assertTrue(np.all(sims >= 0.0) and np.all(sims <= 1.0))

    def test_verdict(self):
        self.assertEqual(rt.verdict(10, 5, 0.5, 0.9), "OU_NOT_PREFERRED")
        self.assertEqual(rt.verdict(1, 5, 0.01, 0.1), "OPTIMUM_SURVIVES_WINDOW")
        self.assertEqual(rt.verdict(1, 5, 0.2, 0.6), "WINDOW_EXPLAINS_OPTIMUM")
        self.assertEqual(rt.verdict(1, 5, 0.2, 0.3), "UNDECIDED_LOW_POWER")


def balanced(names):
    if len(names) == 1:
        return names[0]
    h = len(names) // 2
    return "(%s:1,%s:1)" % (balanced(names[:h]), balanced(names[h:]))


class EndToEnd(unittest.TestCase):
    def test_pipeline_runs(self):
        rng = np.random.default_rng(7)
        names = ["Genus_sp%d" % i for i in range(16)]
        with tempfile.TemporaryDirectory() as d:
            tree_path = os.path.join(d, "tree.nwk")
            with open(tree_path, "w") as fh:
                fh.write(balanced(names) + ";")
            data_path = os.path.join(d, "data.csv")
            lines = ["Common_Name;Frequency;Log_Frequency;Length;Order;Environement;"
                     "Genus;Species_Name"]
            for i, nm in enumerate(names):
                base = math.log10(2.7) + rng.normal(0, 0.2)
                for _ in range(3):
                    lf = base + rng.normal(0, 0.05)
                    lines.append("%s;%s;%s;%s;%s;Forest;Genus;%s" % (
                        "Common %d" % i, ("%.4f" % 10 ** lf).replace(".", ","),
                        ("%.5f" % lf).replace(".", ","), "3,5",
                        "Passeriformes" if i < 8 else "Primates", nm.replace("_", " ")))
            with open(data_path, "w") as fh:
                fh.write("\n".join(lines) + "\n")
            out = os.path.join(d, "out")
            rc = rt.main(["--data-file", data_path, "--tree-file", tree_path,
                          "--out-dir", out, "--skip-gate", "--quick"])
            self.assertEqual(rc, 0)
            for f in ("results.json", "grade.json", "RESULTS-RHYTHM.md",
                      "load_report.json", "download_manifest.json", "log.txt"):
                self.assertTrue(os.path.exists(os.path.join(out, f)), f)
            res = json.load(open(os.path.join(out, "results.json")))
            self.assertIn(res["headline"], rt.VERDICTS, res.get("trace"))
            self.assertEqual(res["load"]["birds_mammals"]["n_species"], 16)
            self.assertEqual(len(json.load(open(os.path.join(out, "grade.json")))), 16)
            self.assertEqual(res["arm_a"]["status"], "OK")


if __name__ == "__main__":
    unittest.main()
