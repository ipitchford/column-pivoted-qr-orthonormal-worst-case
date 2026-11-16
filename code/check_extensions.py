#!/usr/bin/env python3
"""Checks for the coherence-free law and the finite-k guidance in the paper.

Run from the package root:
    python code/check_extensions.py --output results/extensions.json

Sections of the output (all deterministic, no random numbers):
  1. C(H): the closed-form constant (2.2) against a direct 50-digit minimisation
     of (tau+2)(1+H e^{-tau}) over tau >= 0 (Section 4), and the elementary
     bound C(H)^2 <= log H + log(1+log H) + 5 used in Corollary 2.6.
  2. Finite-k comparison of the new bound sqrt(e) exp(C(H) sqrt k) with the
     Drmac-Gugercin factor sqrt(4^k+6k-1)/3 (Table 2 of the paper).
  3. Explicit generalised-Kahan witnesses (Lemma 7.3, Theorem 2.7(b)) at small k, built and
     run through exact-order QRCP in 50-digit arithmetic: orthonormality, pivot
     order and margins, the inverse-norm lower bound, and the upper bounds.
  4. A grid check of the elementary inequality (7.5) behind Theorem 2.7(b).
  5. The two-sided law in (N, k) along N = k^alpha (Corollary 2.8).
  6. The constant implicit in the construction of Chen, Liu, He and Dong.

The script never uses `assert`; it records every check and exits with status
1 if any check fails (also under `python -O`).
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
FAILURES: list[str] = []


def check(flag: bool, label: str) -> bool:
    if not flag:
        FAILURES.append(label)
    return bool(flag)


def rho(H):
    """rho >= 1 solving rho - log rho = 1 + log H."""
    H = mp.mpf(H)
    if H == 1:
        return mp.mpf(1)
    target = 1 + mp.log(H)
    lo, hi = mp.mpf(1), 2 * target + 2
    for _ in range(400):
        mid = (lo + hi) / 2
        if mid - mp.log(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def C(H):
    r = rho(H)
    return mp.sqrt(r) + 1 / mp.sqrt(r)


def direct_min(H):
    """sqrt of min over tau >= 0 of (tau+2)(1+H e^{-tau}) by ternary search."""
    H = mp.mpf(H)
    f = lambda t: (t + 2) * (1 + H * mp.e ** (-t))
    lo, hi = mp.mpf(0), mp.mpf(10) + 2 * mp.log(H + 1)
    for _ in range(600):
        a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if f(a) <= f(b):
            hi = b
        else:
            lo = a
    return mp.sqrt(f((lo + hi) / 2))


def section_constant():
    rows = []
    for H in ("1", "1.5", "2", "4", "10", "100"):
        c_closed, c_direct = C(H), direct_min(H)
        ok = check(abs(c_closed - c_direct) < mp.mpf("1e-20"), f"C({H}) closed vs direct")
        rows.append({"H": H, "rho": mp.nstr(rho(H), 12), "C_closed": mp.nstr(c_closed, 15),
                     "C_direct_min": mp.nstr(c_direct, 15), "agree": ok})
    worst = mp.mpf(-10)
    H = mp.mpf(1)
    while H < mp.mpf(10) ** 12:
        gap = C(H) ** 2 - (mp.log(H) + mp.log(1 + mp.log(H)) + 5)
        worst = max(worst, gap)
        H *= mp.mpf("1.05")
    check(worst <= 0, "C(H)^2 <= log H + log(1+log H) + 5 on grid")
    asym = [{"H": "1e%d" % e, "C_squared_over_logH": mp.nstr(C(mp.mpf(10) ** e) ** 2 / mp.log(mp.mpf(10) ** e), 6)}
            for e in (2, 4, 8, 16, 32)]
    return {"table": rows, "max_of_C2_minus_bound_on_grid_1_to_1e12": mp.nstr(worst, 8),
            "C_squared_over_log_H": asym}


def dg_log(k):
    return 0.5 * mp.log(mp.mpf(4) ** k + 6 * k - 1) - mp.log(3)


def section_dg():
    out = {"crossover": [], "log10_factors": []}
    for H in ("1", "2", "10", "100"):
        c = C(H)
        k = 2
        while not (mp.mpf("0.5") + c * mp.sqrt(k) < dg_log(k)):
            k += 1
        out["crossover"].append({"H": H, "first_k_new_bound_smaller": k})
    for k in (10, 20, 50, 100, 200):
        row = {"k": k, "DG_log10": mp.nstr(dg_log(k) / mp.log(10), 4)}
        for H in ("1", "2", "10"):
            row["new_log10_H" + H] = mp.nstr((mp.mpf("0.5") + C(H) * mp.sqrt(k)) / mp.log(10), 4)
        out["log10_factors"].append(row)
    return out


def qrcp_order(cols, k):
    """Exact QRCP on the given column list (each a list of mpf). Returns the
    selected indices, per-step margins (best minus runner-up squared residual,
    relative to best) and the residual norms of the winners."""
    n = len(cols)
    basis = []
    chosen, margins, winners = [], [], []
    for _ in range(k):
        res = []
        for j in range(n):
            if j in chosen:
                res.append(mp.mpf(-1))
                continue
            v = list(cols[j])
            for q in basis:
                dot = mp.fsum(a * b for a, b in zip(q, v))
                v = [a - dot * b for a, b in zip(v, q)]
            res.append(mp.fsum(a * a for a in v))
        order = sorted(range(n), key=lambda j: res[j], reverse=True)
        best, second = order[0], order[1]
        chosen.append(best)
        margins.append((res[best] - res[second]) / res[best])
        winners.append(mp.sqrt(res[best]))
        v = list(cols[best])
        for q in basis:
            dot = mp.fsum(a * b for a, b in zip(q, v))
            v = [a - dot * b for a, b in zip(v, q)]
        nv = mp.sqrt(mp.fsum(a * a for a in v))
        basis.append([a / nv for a in v])
    return chosen, margins, winners


def kahan_witness(k, L):
    """Lemma 7.3 family: A = (2 sqrt k)^{-1} diag(s^{i-1}) T diag(f), c^2 = 1-e^{-L/(k-1)}."""
    c2 = 1 - mp.e ** (-mp.mpf(L) / (k - 1))
    c, s = mp.sqrt(c2), mp.sqrt(1 - c2)
    f = [1 - mp.mpf(j) / (2 * k) for j in range(k)]
    A = mp.matrix(k, k)
    for i in range(k):
        for j in range(k):
            if j == i:
                A[i, j] = s ** i * f[j]
            elif j > i:
                A[i, j] = -c * s ** i * f[j]
        # rows scaled below
    A = A / (2 * mp.sqrt(k))
    m = int(mp.ceil(64 * k * s ** (-2 * (k - 1))))
    G = mp.eye(k) - A * A.T
    E, Qv = mp.eigsy(G)
    if min(E[i] for i in range(k)) <= 0:
        return None
    S = Qv * mp.diag([mp.sqrt(E[i]) for i in range(k)]) * Qv.T
    gram_err = mp.mnorm(A * A.T + S * S - mp.eye(k), 1)
    cols = [[A[i, j] for i in range(k)] for j in range(k)]
    cols += [[S[i, j] / mp.sqrt(m) for i in range(k)] for j in range(k)]
    chosen, margins, winners = qrcp_order(cols, k)
    sv = mp.svd_r(A, compute_uv=False)
    gamma = 1 / min(sv[i] for i in range(k))
    N = k * (m + 1)
    lower = 2 * mp.sqrt(k) * s ** (-(k - 1)) * c * (1 + c) ** (k - 2)
    upper_coh = mp.sqrt(mp.e * (N - k + 1)) * mp.e ** (C(mp.mpf(N) / k) * mp.sqrt(k))
    upper_range = mp.sqrt(mp.e * (N - k + 1)) * mp.e ** mp.sqrt((k - 1) * mp.log(N - k + 1))
    rec = {
        "k": k, "L": L, "c": mp.nstr(c, 10), "copies_m": m, "N": N,
        "gram_error_1norm": mp.nstr(gram_err, 5),
        "selected": chosen, "min_relative_margin": mp.nstr(min(margins), 6),
        "gamma": mp.nstr(gamma, 12), "gamma_lower_formula": mp.nstr(lower, 12),
        "log_gamma_over_sqrtN_minus_k_plus_1": mp.nstr(mp.log(gamma / mp.sqrt(N - k + 1)), 8),
        "upper_coherence_free": mp.nstr(upper_coh, 8), "upper_diagonal_range": mp.nstr(upper_range, 8),
    }
    tag = f"witness k={k} L={L}"
    rec["ok"] = all([
        check(gram_err < mp.mpf("1e-40"), tag + " orthonormal rows"),
        check(chosen == list(range(k)), tag + " pivots are the core columns in order"),
        check(min(margins) > 0, tag + " strict pivot margins"),
        check(gamma >= lower * (1 - mp.mpf("1e-30")), tag + " gamma >= lower formula"),
        check(gamma <= upper_coh and gamma <= upper_range, tag + " gamma below both upper bounds"),
    ])
    return rec


def lower_bound_terms(k, L):
    x = mp.mpf(L) / (k - 1)
    c = mp.sqrt(1 - mp.e ** (-x))
    exact = mp.log(2 / mp.sqrt(66)) - mp.log(k) / 2 + mp.log(c) + (k - 2) * mp.log(1 + c)
    claimed = mp.sqrt(k * mp.mpf(L)) - L - mp.log(k) - 3
    return exact, claimed


def section_grid():
    worst = None
    count = 0
    ks = sorted(set([3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 80, 120, 200, 400, 1000, 3000, 10000, 100000, 10 ** 6]))
    for k in ks:
        Lmax = mp.mpf(k - 1) / 2
        if Lmax < 1:
            continue
        for t in range(41):
            L = 1 + (Lmax - 1) * mp.mpf(t) / 40
            exact, claimed = lower_bound_terms(k, L)
            slack = exact - claimed
            count += 1
            if worst is None or slack < worst[0]:
                worst = (slack, k, L)
    check(worst[0] >= 0, "inequality (7.5) on grid")
    return {"points": count, "min_slack": mp.nstr(worst[0], 8), "at_k": worst[1], "at_L": mp.nstr(worst[2], 8)}


def section_law():
    rows = []
    for alpha in (3, 4, 6):
        for e in (3, 4, 6, 8):
            k = 10 ** e
            logN = alpha * mp.log(k)
            upper = mp.mpf("0.5") + C(mp.e ** (logN - mp.log(k))) * mp.sqrt(k)
            L = logN - mp.log(66 * mp.mpf(k) ** 2)
            lower = mp.sqrt(k * L) - L - mp.log(k) - 3
            norm = mp.sqrt(k * mp.log(k))
            rows.append({"alpha": alpha, "k": "1e%d" % e, "upper_over_sqrt_klogk": mp.nstr(upper / norm, 6),
                         "lower_over_sqrt_klogk": mp.nstr(lower / norm, 6),
                         "limits": [mp.nstr(mp.sqrt(alpha - 2), 6), mp.nstr(mp.sqrt(alpha - 1), 6)]})
            check(lower <= upper, f"law alpha={alpha} k=1e{e} lower <= upper")
    return rows


def section_chen():
    rows = []
    for H in ("1.5", "2", "4", "10"):
        rows.append({"H": H, "Chen_et_al_coefficient_sup": mp.nstr(mp.sqrt(mp.log(mp.mpf(H))) / 2, 6),
                     "sharp_C": mp.nstr(C(H), 6)})
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", default="results/extensions.json")
    args = ap.parse_args()
    witnesses = []
    for k, L in ((6, 1), (8, 1), (8, 2), (12, 1), (12, 2), (12, 3), (16, 2)):
        rec = kahan_witness(k, L)
        if rec is None:
            check(False, f"witness k={k} L={L}: completion not positive definite")
        else:
            witnesses.append(rec)
    out = {
        "constant": section_constant(),
        "drmac_gugercin_comparison": section_dg(),
        "kahan_witnesses_lemma_7_3": witnesses,
        "inequality_7_5_grid": section_grid(),
        "two_sided_law_N_equals_k_power_alpha": section_law(),
        "chen_et_al_implicit_coefficient": section_chen(),
        "arithmetic": "mpmath %s, %d decimal digits" % (mp.__version__, mp.mp.dps),
        "python": platform.python_version(),
        "failures": FAILURES,
        "passed": not FAILURES,
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("drmac_gugercin_comparison", "inequality_7_5_grid", "failures", "passed")}, indent=1))
    print("ALL CHECKS PASSED" if not FAILURES else "FAILED: %d checks" % len(FAILURES))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
