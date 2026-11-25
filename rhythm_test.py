#!/usr/bin/env python3
# Copyright 2026 MV2, Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""delta-rhythm-window-test.

Can the measurement window alone explain the 'conserved optimal rhythm' of animal
acoustic communication (Piette et al. 2026, PLoS Biol 24(6): e3003798)?
Pre-registered in PREREG.md before any data file was opened.

Models fitted to species mean log10 rhythm on the authors' MCC tree, all with the
same within-species measurement error:
  BM  - Brownian motion (mu, sig2)
  OU  - the authors' stationary OU kernel s2 * exp(-a * D / max D) (theta, a, s2)
  BBM - Brownian motion with reflecting bounds fixed IN ADVANCE at the method's
        admissible window (0.5-20 Hz primary; Boucher & Demery 2016) (sig2 only)
Calibration by parametric bootstrap. Made with AI agents working in MV2's name.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
import time
import urllib.request
import zipfile

import numpy as np

COMMIT = "656fc883f231e0a00f6ae5e15fbeb22bf2960f8e"
RAW = "https://raw.githubusercontent.com/chundrac/phylo-acoustic-rhythm/" + COMMIT + "/"
DATA_FILE = "Data_Rhythm_Final_English.csv"
TREE_FILE = "Coding_mtDNA_HKY_Yule_simple.tree"
CODE_FILE = "OU_models/birds_and_mammals/run_model.R"
ZENODO_ZIP = "https://zenodo.org/records/19816728/files/Data.zip?download=1"
PREREG_URL = ("https://raw.githubusercontent.com/munimtechnologies/"
              "delta-rhythm-window-test/main/PREREG.md")
PREREG_MARK = "registered 2026-11-25"
UA = {"User-Agent": "research-job/1.0 (python urllib)"}

# Order filter copied from the authors' OU_models/birds_and_mammals/run_model.R
# (spellings kept exactly as in their script).
MAMMAL_ORDERS = ["Artiodactyla", "Carnivora", "Chiroptera", "Cingulata", "Eulipotyphla",
                 "Hyracoidea", "Lagomorpha", "Perissodactyla", "Primates", "Rodentia",
                 "Scadentia", "Sciuromorpha"]
BIRD_ORDERS = ["Accipitriformes", "Anseriformes", "Apodiformes", "Bucerotiformes",
               "Cariamiformes", "Casuariiformes", "Charadriiformes", "Colliformes",
               "Columbiformes", "Coraciiformes", "Cuculiformes", "Eurypygiformes",
               "Falconiformes", "Galiiformes", "Gaviiformes", "Gruiformes",
               "Leptosomiformes", "Mesitornithiformes", "Musophagiformes",
               "Nyctibiiformes", "Opisthocomiformes", "Otidiformes", "Passeriformes",
               "Phaethontiformes", "Phoenicopteriformes", "Piciformes",
               "Podecipediformes", "Procellariiformes", "Psittaciformes",
               "Pterocliformes", "Sphenisciformes", "Strigiformes", "Tinamiformes",
               "Trogoniformes"]
NAME_FIX = {"burhinusoedicnemus": "burhinusbistriatus",
            "zaedyuspichyi": "zaedyuspichiy",
            "pantroglodytes": "pantroglodytestroglodytes",
            "chlorocebuspygerythus": "chlorocebuspygerythrus",
            "opisthocomushoazin": "ophisthocomushoazin"}
# (name, lower Hz, upper Hz); the first is the primary window.
WINDOWS = [("primary_0.5_20", 0.5, 20.0), ("low_0.4_20", 0.4, 20.0),
           ("high_0.5_10", 0.5, 10.0), ("wide_0.33_30", 1.0 / 3.0, 30.0)]
VERDICTS = ["OU_NOT_PREFERRED", "OPTIMUM_SURVIVES_WINDOW", "WINDOW_EXPLAINS_OPTIMUM",
            "UNDECIDED_LOW_POWER"]


# ----------------------------------------------------------------- trees
def parse_newick(s):
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    names, lens, parent, kids = [], [], [], []

    def new(p):
        names.append("")
        lens.append(0.0)
        parent.append(p)
        kids.append([])
        k = len(names) - 1
        if p >= 0:
            kids[p].append(k)
        return k

    cur = new(-1)
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            cur = new(cur)
            i += 1
        elif c == ",":
            cur = new(parent[cur])
            i += 1
        elif c == ")":
            cur = parent[cur]
            i += 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in ",();":
                j += 1
            lens[cur] = float(s[i + 1:j].strip())
            i = j
        elif c in "'\"":
            j = s.index(c, i + 1)
            names[cur] = s[i + 1:j]
            i = j + 1
        elif c.isspace():
            i += 1
        else:
            j = i
            while j < n and s[j] not in ",():;":
                j += 1
            names[cur] = s[i:j].strip()
            i = j
    return {"names": names, "lens": lens, "parent": parent, "kids": kids}


def read_tree(text):
    t = re.sub(r"\[[^\]]*\]", "", text)
    trans = {}
    if t.lstrip().upper().startswith("#NEXUS"):
        m = re.search(r"\btranslate\b(.*?);", t, re.I | re.S)
        if m:
            for part in m.group(1).split(","):
                bits = part.strip().split(None, 1)
                if len(bits) == 2:
                    trans[bits[0]] = bits[1].strip().strip("'\"")
        m = re.search(r"\btree\b[^=;]*=(.*?;)", t, re.I | re.S)
        if not m:
            raise ValueError("no tree statement in NEXUS file")
        nwk = m.group(1)
    else:
        nwk = t.strip()
    tr = parse_newick(nwk)
    for k, nm in enumerate(tr["names"]):
        if not tr["kids"][k] and nm in trans:
            tr["names"][k] = trans[nm]
    return tr


def tree_depths(tr):
    d = [0.0] * len(tr["names"])
    for i, p in enumerate(tr["parent"]):
        if p >= 0:
            d[i] = d[p] + tr["lens"][i]
    return d


def tip_matrices(tr, tips):
    """C: shared path below the tips' MRCA, scaled so root-to-tip height is 1.
    Dn: patristic distance / its maximum. H: height in original units.
    Dmax_h: maximum patristic distance in height units."""
    d = tree_depths(tr)
    paths = []
    for t in tips:
        p, k = [], t
        while k >= 0:
            p.append(k)
            k = tr["parent"][k]
        paths.append(p)
    sets = [set(p) for p in paths]
    m = len(tips)
    S = np.zeros((m, m))
    for i in range(m):
        S[i, i] = d[tips[i]]
        for j in range(i + 1, m):
            v = max(d[k] for k in paths[i] if k in sets[j])
            S[i, j] = S[j, i] = v
    root_depth = float(S.min())
    dt = np.array([d[t] for t in tips])
    H = float(dt.max() - root_depth)
    C = (S - root_depth) / H
    D = dt[:, None] + dt[None, :] - 2.0 * S
    Dmax = float(D.max())
    return C, D / Dmax, H, Dmax / H


# ------------------------------------------------------- gaussian models
def gauss_ll(y, V):
    n = len(y)
    try:
        Lc = np.linalg.cholesky(V)
    except np.linalg.LinAlgError:
        return -np.inf, float("nan")
    a = np.linalg.solve(Lc, np.ones(n))
    b = np.linalg.solve(Lc, y)
    mu = float(a @ b / (a @ a))
    r = b - mu * a
    ll = -0.5 * (n * math.log(2 * math.pi) + 2.0 * float(np.sum(np.log(np.diag(Lc))))
                 + float(r @ r))
    return ll, mu


def golden_max(f, lo, hi, iters=30):
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - g * (b - a)
    d = a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc >= fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (c, fc) if fc >= fd else (d, fd)


def grid_golden(f, lo, hi, ngrid=25, iters=30):
    xs = np.linspace(lo, hi, ngrid)
    fs = [f(float(x)) for x in xs]
    i = int(np.argmax(fs))
    x, fx = golden_max(f, float(xs[max(i - 1, 0)]), float(xs[min(i + 1, ngrid - 1)]), iters)
    if fs[i] >= fx:
        return float(xs[i]), fs[i]
    return x, fx


def aic(fit):
    return 2.0 * fit["k"] - 2.0 * fit["ll"]


def fit_bm(y, C, me):
    M = np.diag(me)
    f = lambda ls: gauss_ll(y, math.exp(ls) * C + M)[0]
    ls, ll = grid_golden(f, math.log(1e-5), math.log(1e2), 25, 30)
    mu = gauss_ll(y, math.exp(ls) * C + M)[1]
    return {"model": "BM", "ll": ll, "k": 2, "mu": mu, "sig2": math.exp(ls)}


def fit_ou(y, Dn, me, na=31):
    M = np.diag(me)

    def inner(la):
        E = np.exp(-math.exp(la) * Dn)
        f = lambda ls: gauss_ll(y, math.exp(ls) * E + M)[0]
        ls, ll = grid_golden(f, math.log(1e-5), math.log(1e2), 12, 25)
        return ll, ls

    las = np.linspace(math.log(0.01), math.log(300.0), na)
    res = [inner(float(x)) for x in las]
    i = int(np.argmax([r[0] for r in res]))
    best_la, (best_ll, best_ls) = float(las[i]), res[i]
    la, ll = golden_max(lambda x: inner(x)[0], float(las[max(i - 1, 0)]),
                        float(las[min(i + 1, na - 1)]), 15)
    if ll > best_ll:
        best_la = la
        best_ll, best_ls = inner(la)
    E = np.exp(-math.exp(best_la) * Dn)
    theta = gauss_ll(y, math.exp(best_ls) * E + M)[1]
    return {"model": "OU", "ll": best_ll, "k": 3, "theta": theta,
            "alpha": math.exp(best_la), "s2": math.exp(best_ls)}


def sim_mvn(rng, mean, V, B):
    Lc = np.linalg.cholesky(V + 1e-10 * np.eye(len(V)))
    return mean + rng.standard_normal((B, len(V))) @ Lc.T


# ------------------------------------------------ bounded Brownian motion
def reflect(x, L, U):
    W = U - L
    r = (x - L) % (2.0 * W)
    if r > W:
        r = 2.0 * W - r
    return L + r


class BBM:
    """Reflecting-bounds Brownian motion on a grid (Boucher & Demery 2016),
    uniform (stationary) root. Branch lengths divided by `scale`."""

    def __init__(self, tr, tips, scale, L, U, nb=100):
        self.L, self.U, self.nb = float(L), float(U), nb
        self.W = self.U - self.L
        self.h = self.W / nb
        self.x = self.L + (np.arange(nb) + 0.5) * self.h
        A = -2.0 * np.eye(nb) + np.eye(nb, k=1) + np.eye(nb, k=-1)
        A[0, 0] = -1.0
        A[-1, -1] = -1.0
        lam, V = np.linalg.eigh(A)
        self.lam = np.minimum(lam, 0.0)
        self.V = V
        self.parent = tr["parent"]
        self.lens = [l / scale for l in tr["lens"]]
        n = len(self.parent)
        mark = [False] * n
        for t in tips:
            k = t
            while k >= 0 and not mark[k]:
                mark[k] = True
                k = self.parent[k]
        self.post = [i for i in range(n - 1, -1, -1) if mark[i]]
        self.root = self.post[-1]
        self.tipj = {t: j for j, t in enumerate(tips)}

    def loglik(self, y, me, sig2):
        sd = np.sqrt(np.asarray(me, float) + self.h ** 2 / 12.0)
        c = sig2 / (2.0 * self.h ** 2)
        acc = {}
        logs = 0.0
        for i in self.post:
            if i in self.tipj:
                j = self.tipj[i]
                z = (y[j] - self.x) / sd[j]
                lv = -0.5 * z * z - math.log(sd[j] * math.sqrt(2.0 * math.pi))
                m = float(lv.max())
                v = np.exp(lv - m)
                logs += m
            else:
                v = acc.pop(i, None)
                if v is None:
                    return -np.inf
                m = float(v.max())
                if not m > 0:
                    return -np.inf
                v = v / m
                logs += math.log(m)
            if i == self.root:
                s = float(v.mean())
                if not s > 0:
                    return -np.inf
                return logs + math.log(s)
            w = self.V @ (np.exp(c * self.lens[i] * self.lam) * (self.V.T @ v))
            w = np.maximum(w, 0.0)
            s = float(w.max())
            if not s > 0:
                return -np.inf
            logs += math.log(s)
            w = w / s
            p = self.parent[i]
            acc[p] = w if p not in acc else acc[p] * w
        return -np.inf

    def simulate(self, rng, sig2, me, B):
        out = np.zeros((B, len(self.tipj)))
        pre = self.post[::-1]
        sd_me = np.sqrt(np.asarray(me, float))
        for b in range(B):
            val = {}
            for i in pre:
                if i == self.root:
                    val[i] = rng.uniform(self.L, self.U)
                else:
                    step = rng.normal(0.0, math.sqrt(max(sig2 * self.lens[i], 0.0)))
                    val[i] = reflect(val[self.parent[i]] + step, self.L, self.U)
            for t, j in self.tipj.items():
                out[b, j] = val[t] + rng.normal(0.0, sd_me[j])
        return out


def fit_bbm(bbm, y, me):
    f = lambda ls: bbm.loglik(y, me, math.exp(ls))
    ls, ll = grid_golden(f, math.log(1e-4), math.log(1e3), 22, 30)
    return {"model": "BBM", "ll": ll, "k": 1, "sig2": math.exp(ls),
            "L_log10": bbm.L, "U_log10": bbm.U}


# ------------------------------------------------------------------ data
def clean(s):
    return re.sub(r"[ _]", "", str(s)).lower()


def load_table(path):
    raw = open(path, "rb").read()
    text = None
    for enc in ("utf-8-sig", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    first = text.splitlines()[0] if text.strip() else ""
    sep = ";" if first.count(";") > first.count(",") else ","
    rows = list(csv.reader(io.StringIO(text), delimiter=sep))
    header = [h.strip() for h in rows[0]]
    body = [r for r in rows[1:] if any(c.strip() for c in r)]
    offset = 0
    if body and len(header) == len(body[0]) - 1:
        header = ["_rowname"] + header   # R write.table row names
        offset = 1
    return header, body, sep, offset


def num(s):
    s = (s or "").strip().strip('"')
    if s in ("", "NA", "NaN", "nan", "NULL", "-", "#N/A"):
        return None
    try:
        v = float(s.replace(",", "."))
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def col_index(header, name):
    if name in header:
        return header.index(name)
    low = [h.lower() for h in header]
    return low.index(name.lower()) if name.lower() in low else None


def length_index(header):
    for k, h in enumerate(header):
        hl = h.lower()
        if hl.startswith("length") or hl in ("duration", "sequence_length", "seq_length"):
            return k
    return None


def build_dataset(header, rows, offset, tr, orders):
    i_order = col_index(header, "Order")
    i_logf = col_index(header, "Log_Frequency")
    i_f = col_index(header, "Frequency")
    i7, i8 = 6 + offset, 7 + offset
    rep = {"order_col": i_order, "log_frequency_col": i_logf, "frequency_col": i_f,
           "name_cols_r_7_8": [header[i7] if i7 < len(header) else None,
                               header[i8] if i8 < len(header) else None]}
    tipname = {}
    for k, nm in enumerate(tr["names"]):
        if not tr["kids"][k]:
            tipname.setdefault(clean(nm), k)
    per, unmatched = {}, {}
    n_in, n_na, n_strip_only = 0, 0, 0
    all_lf = []
    for r in rows:
        if len(r) < len(header):
            r = r + [""] * (len(header) - len(r))
        if orders is not None:
            if i_order is None:
                rep["error"] = "no Order column"
                break
            if r[i_order] not in orders:
                if r[i_order].strip() in orders:
                    n_strip_only += 1
                continue
        n_in += 1
        if i_logf is not None:
            lf = num(r[i_logf])
        elif i_f is not None:
            f = num(r[i_f])
            lf = math.log10(f) if f is not None and f > 0 else None
        else:
            rep["error"] = "no Log_Frequency or Frequency column"
            break
        if lf is None:
            n_na += 1
            continue
        n8, n7 = clean(r[i8]), clean(r[i7])
        nm = n8 if n8 in tipname else n7
        nm = NAME_FIX.get(nm, nm)
        if nm not in tipname:
            unmatched[nm] = unmatched.get(nm, 0) + 1
            continue
        all_lf.append(lf)
        per.setdefault(tipname[nm], []).append(lf)
    tips = sorted(per)
    y = np.array([float(np.mean(per[t])) for t in tips])
    n = np.array([len(per[t]) for t in tips])
    ss, dfree = 0.0, 0
    for t in tips:
        v = np.array(per[t])
        if len(v) >= 2:
            ss += float(((v - v.mean()) ** 2).sum())
            dfree += len(v) - 1
    tau2 = ss / dfree if dfree > 0 and ss > 0 else 0.01
    rep.update({"rows_in_orders": n_in, "rows_missing_rate": n_na,
                "rows_matched": len(all_lf), "rows_excluded_only_by_whitespace": n_strip_only,
                "n_species": len(tips), "tau2_within_species_log10": tau2,
                "tau2_fallback_used": not (dfree > 0 and ss > 0),
                "unmatched_names": dict(sorted(unmatched.items(),
                                               key=lambda kv: -kv[1])[:40]),
                "tip_names": [tr["names"][t] for t in tips],
                "rows_per_species": [int(v) for v in n]})
    return {"tips": tips, "y": y, "n": n, "me": tau2 / n, "tau2": tau2,
            "row_lf": np.array(all_lf), "report": rep}


def within_slope(groups, rng, B=1000):
    def slope(gs):
        sxy = sxx = 0.0
        for x, yv in gs:
            xc, yc = x - x.mean(), yv - yv.mean()
            sxy += float(xc @ yc)
            sxx += float(xc @ xc)
        return sxy / sxx if sxx > 0 else float("nan")
    est = slope(groups)
    bs = []
    for _ in range(B):
        idx = rng.integers(0, len(groups), len(groups))
        bs.append(slope([groups[i] for i in idx]))
    bs = np.array([b for b in bs if math.isfinite(b)])
    if len(bs) == 0:
        return est, None, None
    return est, float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def arm_a(tables, orders, rng):
    """Descriptive: cycles per sequence and the within-species slope of log rate
    on log sequence length, from the first table that has a length column."""
    for label, (header, rows, offset) in tables:
        i_len = length_index(header)
        i_order = col_index(header, "Order")
        i_logf = col_index(header, "Log_Frequency")
        i_f = col_index(header, "Frequency")
        if i_len is None or i_order is None or (i_logf is None and i_f is None):
            continue
        i_name = col_index(header, "Common_Name")
        if i_name is None:
            i_name = 7 + offset
        recs = []
        for r in rows:
            if len(r) < len(header):
                r = r + [""] * (len(header) - len(r))
            if r[i_order] not in orders:
                continue
            if i_logf is not None:
                lf = num(r[i_logf])
            else:
                f = num(r[i_f])
                lf = math.log10(f) if f is not None and f > 0 else None
            T = num(r[i_len])
            if lf is None or T is None or T <= 0:
                continue
            recs.append((r[i_name], lf, T))
        if not recs:
            continue
        Ts = np.array([t for _, _, t in recs])
        out = {"source": label, "length_column": header[i_len], "n_rows": len(recs),
               "median_length": float(np.median(Ts))}
        if out["median_length"] > 100:
            out["status"] = "LENGTH_UNITS_UNCLEAR"
            return out
        cyc = np.array([10 ** lf * t for _, lf, t in recs])
        out.update({"status": "OK", "median_cycles": float(np.median(cyc)),
                    "share_cycles_lt3": float(np.mean(cyc < 3)),
                    "share_cycles_lt4": float(np.mean(cyc < 4))})
        by = {}
        for nm, lf, t in recs:
            by.setdefault(nm, []).append((math.log10(t), lf))
        groups = [(np.array([a for a, _ in v]), np.array([b for _, b in v]))
                  for v in by.values() if len(v) >= 3]
        groups = [g for g in groups if float(np.ptp(g[0])) > 0]
        out["n_species_slope"] = len(groups)
        if groups:
            est, lo, hi = within_slope(groups, rng)
            out.update({"within_slope": est, "within_slope_ci": [lo, hi]})
        return out
    return {"status": "UNTESTABLE_NO_LENGTH"}


# --------------------------------------------------------------- study
def verdict(aic_ou, aic_bm, p_bbm, power):
    if not aic_ou < aic_bm:
        return "OU_NOT_PREFERRED"
    if p_bbm < 0.05:
        return "OPTIMUM_SURVIVES_WINDOW"
    if power >= 0.5:
        return "WINDOW_EXPLAINS_OPTIMUM"
    return "UNDECIDED_LOW_POWER"


def study(ds, tr, rng, windows, B, log):
    y, me, tips = ds["y"], ds["me"], ds["tips"]
    C, Dn, H, Dmax_h = tip_matrices(tr, tips)
    fbm, fou = fit_bm(y, C, me), fit_ou(y, Dn, me)
    stat1 = fou["ll"] - fbm["ll"]
    out = {"n_species": len(y), "tree_height_original_units": H,
           "species_geomean_hz": 10 ** float(y.mean()),
           "species_hz_2.5_97.5": [10 ** float(np.quantile(y, 0.025)),
                                   10 ** float(np.quantile(y, 0.975))],
           "species_sd_log10": float(y.std(ddof=1)),
           "BM": fbm, "OU": fou, "aic_BM": aic(fbm), "aic_OU": aic(fou),
           "OU_theta_hz": 10 ** fou["theta"],
           "OU_half_life_tree_heights": math.log(2.0) * Dmax_h / fou["alpha"],
           "stat_ll_OU_minus_BM": stat1}
    log("obs fits done: dll(OU-BM)=%.3f" % stat1)
    sims = sim_mvn(rng, fbm["mu"], fbm["sig2"] * C + np.diag(me), B["bm"])
    s1 = []
    for s in sims:
        s1.append(fit_ou(s, Dn, me)["ll"] - fit_bm(s, C, me)["ll"])
    out["p_BM_null"] = (1 + sum(v >= stat1 for v in s1)) / (len(s1) + 1)
    out["null_BM_stat_quantiles"] = [float(q) for q in np.quantile(s1, [0.5, 0.95])]
    log("BM null done: p=%.4f" % out["p_BM_null"])
    out["windows"] = []
    for idx, (name, lo_hz, hi_hz) in enumerate(windows):
        L, U = math.log10(lo_hz), math.log10(hi_hz)
        widened = False
        if float(y.min()) < L:
            L, widened = float(y.min()) - 0.02, True
        if float(y.max()) > U:
            U, widened = float(y.max()) + 0.02, True
        bbm = BBM(tr, tips, H, L, U)
        fb = fit_bbm(bbm, y, me)
        stat2 = fou["ll"] - fb["ll"]
        b0 = B["bbm"] if idx == 0 else B["sens"]
        bp = B["pow"] if idx == 0 else max(B["sens"] // 2, 2)
        s2, ou_beats_bm = [], 0
        for s in bbm.simulate(rng, fb["sig2"], me, b0):
            o = fit_ou(s, Dn, me)
            s2.append(o["ll"] - fit_bbm(bbm, s, me)["ll"])
            if idx == 0 and aic(o) < aic(fit_bm(s, C, me)):
                ou_beats_bm += 1
        p_bbm = (1 + sum(v >= stat2 for v in s2)) / (len(s2) + 1)
        q95 = float(np.quantile(s2, 0.95))
        pw = sim_mvn(rng, fou["theta"], fou["s2"] * np.exp(-fou["alpha"] * Dn) + np.diag(me), bp)
        hits = 0
        for s in pw:
            if fit_ou(s, Dn, me)["ll"] - fit_bbm(bbm, s, me)["ll"] > q95:
                hits += 1
        power = hits / bp
        w = {"name": name, "lower_hz": 10 ** L, "upper_hz": 10 ** U, "widened": widened,
             "BBM": fb, "aic_BBM": aic(fb), "stat_ll_OU_minus_BBM": stat2,
             "p_BBM_null": p_bbm, "null_BBM_q95": q95, "n_null": b0,
             "power_vs_OU": power, "n_power": bp,
             "uniform_window_sd_log10": (U - L) / math.sqrt(12.0)}
        if idx == 0:
            w["share_BBM_null_sims_OU_beats_BM_by_AIC"] = ou_beats_bm / b0
        w["verdict"] = verdict(out["aic_OU"], out["aic_BM"], p_bbm, power)
        out["windows"].append(w)
        log("window %s: p=%.4f power=%.2f verdict=%s" % (name, p_bbm, power, w["verdict"]))
    return out


def analyse(data_path, tree_path, zip_path, cfg, log=print):
    rng = np.random.default_rng(cfg["seed"])
    res = {"headline": None, "config": cfg}
    with open(tree_path, encoding="utf-8", errors="replace") as fh:
        tr = read_tree(fh.read())
    header, rows, sep, offset = load_table(data_path)
    res["load"] = {"data_header": header, "sep": sep, "rowname_offset": offset,
                   "n_rows": len(rows),
                   "n_tree_tips": sum(1 for k in tr["kids"] if not k)}
    bm_orders = set(MAMMAL_ORDERS + BIRD_ORDERS)
    sets = {}
    for label, orders in (("birds_mammals", bm_orders), ("all_taxa", None)):
        sets[label] = build_dataset(header, rows, offset, tr, orders)
        res["load"][label] = sets[label]["report"]
    main_ds = sets["birds_mammals"]
    if main_ds["report"].get("error") or len(main_ds["tips"]) < cfg["min_species"]:
        res["headline"] = "LOADER_STOP"
        return res
    lf = main_ds["row_lf"]
    res["rows_share_outside_primary_window"] = float(
        np.mean((lf < math.log10(0.5)) | (lf > math.log10(20.0))))
    res["birds_mammals"] = study(main_ds, tr, rng, WINDOWS, cfg["B"], log)
    all_ds = sets["all_taxa"]
    if len(all_ds["tips"]) >= cfg["min_species"]:
        res["all_taxa"] = study(all_ds, tr, rng, WINDOWS[:1], cfg["B_all"], log)
    tables = [("main", (header, rows, offset))]
    if zip_path and os.path.exists(zip_path):
        try:
            with zipfile.ZipFile(zip_path) as z:
                for nm in z.namelist():
                    if nm.lower().endswith(".csv") and "__macosx" not in nm.lower():
                        tmp = os.path.join(cfg["tmp_dir"], re.sub(r"[^A-Za-z0-9_.-]", "_", nm))
                        with open(tmp, "wb") as fh:
                            fh.write(z.read(nm))
                        h2, r2, _, o2 = load_table(tmp)
                        tables.append(("zip:" + nm, (h2, r2, o2)))
        except (zipfile.BadZipFile, OSError) as e:
            res["zip_error"] = str(e)
    res["arm_a"] = arm_a(tables, bm_orders, rng)
    res["headline"] = res["birds_mammals"]["windows"][0]["verdict"]
    return res


# --------------------------------------------------------------- grading
PRIORS = [
    ("P1", "loader matches at least 60 bird and mammal species to the tree", 0.80),
    ("P2", "OU beats BM by AIC on the observed data", 0.80),
    ("P3", "OU over BM survives the BM parametric bootstrap (p < 0.05)", 0.65),
    ("P4", "window-bounded BM (primary) beats BM by AIC", 0.55),
    ("P5", "OU advantage over window-bounded BM is NOT significant (p >= 0.05)", 0.45),
    ("P6", "power of the OU-vs-window test is at least 0.5", 0.50),
    ("P7", "fitted OU optimum lies in 2.0-3.6 Hz", 0.75),
    ("P8", "OU half-life is below 1 tree height", 0.75),
    ("P9", "fewer than 2% of rows lie outside 0.5-20 Hz", 0.70),
    ("P10", "median cycles per sequence is below 10", 0.50),
    ("P11", "within-species slope of log rate on log length is negative, CI excludes 0", 0.55),
    ("P12", "all-taxa verdict equals the birds-and-mammals verdict", 0.70),
    ("P13", "all three sensitivity windows give the primary verdict", 0.45),
    ("P14", "species geometric mean rate lies in 2.025-3.375 Hz", 0.75),
    ("P15", "in at least half of window-only simulations OU beats BM by AIC", 0.50),
    ("P16", "species SD of log10 rate is below half the uniform-window SD", 0.45),
]


def grade(res):
    bmres = res.get("birds_mammals")
    out = {}
    nsp = res.get("load", {}).get("birds_mammals", {}).get("n_species")
    out["P1"] = None if nsp is None else nsp >= 60
    if bmres:
        w0 = bmres["windows"][0]
        out["P2"] = bmres["aic_OU"] < bmres["aic_BM"]
        out["P3"] = bmres["p_BM_null"] < 0.05
        out["P4"] = w0["aic_BBM"] < bmres["aic_BM"]
        out["P5"] = w0["p_BBM_null"] >= 0.05
        out["P6"] = w0["power_vs_OU"] >= 0.5
        out["P7"] = 2.0 <= bmres["OU_theta_hz"] <= 3.6
        out["P8"] = bmres["OU_half_life_tree_heights"] < 1.0
        out["P9"] = res["rows_share_outside_primary_window"] < 0.02
        a = res.get("arm_a", {})
        out["P10"] = a["median_cycles"] < 10 if a.get("status") == "OK" else None
        ci = a.get("within_slope_ci") if a.get("status") == "OK" else None
        out["P11"] = (ci[1] < 0) if ci and ci[1] is not None else None
        out["P12"] = (res["all_taxa"]["windows"][0]["verdict"] == w0["verdict"]
                      if "all_taxa" in res else None)
        out["P13"] = all(w["verdict"] == w0["verdict"] for w in bmres["windows"][1:])
        out["P14"] = 2.025 <= bmres["species_geomean_hz"] <= 3.375
        out["P15"] = w0["share_BBM_null_sims_OU_beats_BM_by_AIC"] >= 0.5
        out["P16"] = bmres["species_sd_log10"] < 0.5 * w0["uniform_window_sd_log10"]
    rows = []
    for pid, text, prior in PRIORS:
        o = out.get(pid)
        rows.append({"id": pid, "prediction": text, "prior": prior, "outcome": o,
                     "brier": None if o is None else (prior - (1.0 if o else 0.0)) ** 2})
    return rows


def fmt(v, nd=4):
    if isinstance(v, float):
        return ("%." + str(nd) + "g") % v
    return str(v)


def write_md(res, grades):
    L = ["# RESULTS: delta-rhythm-window-test", "",
         "Generated by rhythm_test.py (MV2, Munim, Inc.; made with AI agents; no human "
         "expert has checked this). Pre-registration: PREREG.md.", "",
         "**Headline: %s**" % res.get("headline"), ""]
    bmres = res.get("birds_mammals")
    if bmres:
        L += ["## Birds and mammals (primary)", "",
              "| quantity | value |", "|---|---|",
              "| species | %d |" % bmres["n_species"],
              "| species geometric mean (Hz) | %s |" % fmt(bmres["species_geomean_hz"]),
              "| species 2.5-97.5%% (Hz) | %s - %s |" % tuple(fmt(v) for v in bmres["species_hz_2.5_97.5"]),
              "| species SD of log10 rate | %s |" % fmt(bmres["species_sd_log10"]),
              "| AIC BM / OU | %s / %s |" % (fmt(bmres["aic_BM"], 6), fmt(bmres["aic_OU"], 6)),
              "| OU optimum (Hz) | %s |" % fmt(bmres["OU_theta_hz"]),
              "| OU half-life (tree heights) | %s |" % fmt(bmres["OU_half_life_tree_heights"]),
              "| p, OU vs BM under BM null | %s |" % fmt(bmres["p_BM_null"]),
              "| rows outside 0.5-20 Hz | %s |" % fmt(res["rows_share_outside_primary_window"]),
              ""]
        L += ["| window | lower-upper Hz | AIC BBM | p (OU vs window) | power | verdict |",
              "|---|---|---|---|---|---|"]
        for w in bmres["windows"]:
            L.append("| %s | %s-%s%s | %s | %s | %s | %s |" % (
                w["name"], fmt(w["lower_hz"]), fmt(w["upper_hz"]),
                " (widened)" if w["widened"] else "", fmt(w["aic_BBM"], 6),
                fmt(w["p_BBM_null"]), fmt(w["power_vs_OU"]), w["verdict"]))
        w0 = bmres["windows"][0]
        L += ["", "Window-only simulations in which OU beats BM by AIC: %s" %
              fmt(w0["share_BBM_null_sims_OU_beats_BM_by_AIC"]), ""]
    if res.get("all_taxa"):
        w = res["all_taxa"]["windows"][0]
        L += ["## All taxa", "", "species %d; p (OU vs window) %s; power %s; verdict %s" % (
            res["all_taxa"]["n_species"], fmt(w["p_BBM_null"]), fmt(w["power_vs_OU"]),
            w["verdict"]), ""]
    if res.get("arm_a"):
        L += ["## Arm A (descriptive)", "", "```", json.dumps(res["arm_a"], indent=1), "```", ""]
    L += ["## Predictions", "", "| id | prediction | prior | outcome | Brier |",
          "|---|---|---|---|---|"]
    for g in grades:
        L.append("| %s | %s | %.2f | %s | %s |" % (g["id"], g["prediction"], g["prior"],
                                                  g["outcome"], fmt(g["brier"], 3)
                                                  if g["brier"] is not None else "-"))
    return "\n".join(L) + "\n"


def dump(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1,
                  default=lambda o: o.item() if hasattr(o, "item") else str(o))


# ------------------------------------------------------------------ main
def fetch(url, path):
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=300) as r, open(path, "wb") as fh:
            fh.write(r.read())
    data = open(path, "rb").read()
    return {"url": url, "path": path, "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def gate():
    try:
        req = urllib.request.Request(PREREG_URL, headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            return PREREG_MARK in r.read().decode("utf-8", "replace")
    except Exception:
        return False


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--data-file", default=None)
    ap.add_argument("--tree-file", default=None)
    ap.add_argument("--skip-gate", action="store_true")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    logs = []

    def log(msg):
        line = time.strftime("%H:%M:%S ") + msg
        logs.append(line)
        print(line, flush=True)

    if not a.skip_gate and not gate():
        dump(os.path.join(a.out_dir, "results.json"), {"headline": "PREREG_NOT_LIVE"})
        log("PREREG not live; no data downloaded")
        return 0
    cfg = {"seed": 20261125, "min_species": 40, "tmp_dir": a.out_dir,
           "B": {"bm": 400, "bbm": 400, "pow": 200, "sens": 150},
           "B_all": {"bm": 200, "bbm": 200, "pow": 100, "sens": 2}}
    if a.quick:
        cfg.update({"min_species": 8, "B": {"bm": 3, "bbm": 3, "pow": 2, "sens": 2},
                    "B_all": {"bm": 2, "bbm": 2, "pow": 2, "sens": 2}})
    zip_path = None
    if a.data_file and a.tree_file:
        data_path, tree_path = a.data_file, a.tree_file
        manifest = {"local": [a.data_file, a.tree_file]}
    else:
        os.makedirs(a.data_dir, exist_ok=True)
        manifest = {"commit": COMMIT, "files": []}
        data_path = os.path.join(a.data_dir, DATA_FILE)
        tree_path = os.path.join(a.data_dir, TREE_FILE)
        zip_path = os.path.join(a.data_dir, "zenodo-19816728-Data.zip")
        for url, p in ((RAW + DATA_FILE, data_path), (RAW + TREE_FILE, tree_path),
                       (RAW + CODE_FILE, os.path.join(a.data_dir, "run_model.R")),
                       (ZENODO_ZIP, zip_path)):
            try:
                manifest["files"].append(fetch(url, p))
            except Exception as e:
                manifest["files"].append({"url": url, "error": str(e)})
        try:
            with zipfile.ZipFile(zip_path) as z:
                manifest["zip_members"] = [[i.filename, i.file_size] for i in z.infolist()]
        except Exception as e:
            manifest["zip_error"] = str(e)
    dump(os.path.join(a.out_dir, "download_manifest.json"), manifest)
    tmp = os.path.join(a.out_dir, "tmp")
    os.makedirs(tmp, exist_ok=True)
    cfg["tmp_dir"] = tmp
    try:
        res = analyse(data_path, tree_path, zip_path, cfg, log)
    except Exception as e:
        import traceback
        res = {"headline": "LOADER_STOP", "error": repr(e), "trace": traceback.format_exc()}
    dump(os.path.join(a.out_dir, "load_report.json"), res.get("load", {}))
    grades = grade(res) if res.get("headline") in VERDICTS else []
    dump(os.path.join(a.out_dir, "results.json"), res)
    dump(os.path.join(a.out_dir, "grade.json"), grades)
    with open(os.path.join(a.out_dir, "RESULTS-RHYTHM.md"), "w", encoding="utf-8") as fh:
        fh.write(write_md(res, grades))
    with open(os.path.join(a.out_dir, "log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(logs) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
