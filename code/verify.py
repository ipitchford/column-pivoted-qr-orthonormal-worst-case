#!/usr/bin/env python3
"""Independent finite checks for sharp coherent-frame QRCP asymptotics.

This is numerical and symbolic supporting evidence, not a proof of the
unrestricted theorem.  The JSON records every failed check and the precision
actually used.  No external data, network access, or generated code is needed.

Run from the package root: python code/verify.py --output results/verification.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import itertools
import math
import platform
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.optimize as opt
import sympy as sp


SEED = 20261115


def parameter_r(h: float) -> float:
    """Unique r >= 1 with r-log(r)=1+log(h), using a stable equation."""
    if h < 1:
        raise ValueError("The coherence parameter must be at least one")
    if h == 1:
        return 1.0
    target = math.log(h)
    # In x=r-1 coordinates, x-log1p(x) avoids cancellation near h=1.
    hi = max(2.0, 2.0 * (1.0 + target))
    while hi - math.log1p(hi) < target:
        hi *= 2
    x = opt.brentq(lambda x: x - math.log1p(x) - target,
                  0.0, hi, xtol=1e-14, rtol=4*np.finfo(float).eps)
    return 1.0 + x


def sharp_constant(h: float) -> float:
    r = parameter_r(h)
    return math.sqrt(r) + 1.0 / math.sqrt(r)


def exact_rational_checks() -> dict:
    """Exact rational tail, taper, weighted-energy and inverse checks."""
    records = []
    for k in (2, 3, 4, 5):
        c, s = sp.Rational(3, 5), sp.Rational(4, 5)
        d = [s**i for i in range(k)]
        v = [1 - sp.Rational(j, 10*k) for j in range(k)]
        core = sp.Matrix(k, k, lambda i,j:
                         d[i] if i == j else (-c*d[i] if i < j else 0))
        tapered = core * sp.diag(*v)
        equality = all(sum(core[l,j]**2 for l in range(i,j+1)) == d[i]**2
                       for i in range(k) for j in range(i,k))
        strict = all(sum(tapered[l,j]**2 for l in range(i,j+1))
                     < tapered[i,i]**2
                     for i in range(k) for j in range(i+1,k))
        expected = c*(1+c)**(k-2)/d[-1]
        inverse_entry = core.inv()[0,k-1]
        inverse_ok = inverse_entry == expected
        energy_records = []
        dd = [tapered[i,i] for i in range(k)]
        for a in (1, 2, 3):
            weighted = sp.Matrix(k,k,lambda i,j:
                                 tapered[i,j]*dd[j]**a/dd[i]**(a+1))
            energies = [sum(weighted[i,j]**2 for i in range(j))
                        for j in range(k)]
            good = all(value <= sp.Rational(1,a) for value in energies)
            energy_records.append({
                "power": a,
                "column_offdiagonal_energy": [str(x) for x in energies],
                "upper_bound": str(sp.Rational(1,a)),
                "determinant": str(weighted.det()),
                "passed": bool(good and weighted.det() == 1),
            })
        records.append({
            "rank": k,
            "core_matrix": [[str(core[i,j]) for j in range(k)]
                            for i in range(k)],
            "all_tail_equalities": bool(equality),
            "all_tapered_core_pivots_strict": bool(strict),
            "inverse_entry": str(inverse_entry),
            "inverse_entry_formula": str(expected),
            "inverse_entry_identity": bool(inverse_ok),
            "weighted_energy_checks": energy_records,
            "passed": bool(equality and strict and inverse_ok
                           and all(x["passed"] for x in energy_records)),
        })
    return {"arithmetic": "exact SymPy rationals", "records": records,
            "passed": all(x["passed"] for x in records)}


def inverse_norm_and_bounds(q: np.ndarray) -> dict:
    n, k = q.shape
    _, all_r, piv = la.qr(q.T, mode="economic", pivoting=True)
    r = all_r[:, :k]
    d = np.abs(np.diag(r))
    selected = q[piv[:k], :]
    # An independent SVD of selected rows, rather than the QR factor itself.
    gamma = 1.0/la.svdvals(selected)[-1]
    coherence = n/k*np.max(np.sum(q*q, axis=1))
    upper_log = .5 + .5*math.log(n-k+1) + sharp_constant(coherence)*math.sqrt(k)
    spread = math.log(d[0]**2/d[-1]**2)
    range_upper_log = .5 - math.log(d[-1]) + math.sqrt(max(0.,(k-1)*spread))
    worst_tail_error = 0.0
    worst_floor_error = 0.0
    for i in range(k):
        # Includes every unselected column, not just columns of the core.
        tails = np.sum(all_r[i:,i:]**2, axis=0)
        worst_tail_error = max(worst_tail_error, float(np.max(tails)-d[i]**2))
        floor = (k-i)/(n-i)
        worst_floor_error = max(worst_floor_error, floor-d[i]**2)
    # Check the weighted-energy lemma at the independently chosen power 2.
    weighted = (r / d[:,None]**3) * d[None,:]**2
    column_energies = np.sum(np.triu(weighted,1)**2, axis=0)
    orth_error = la.norm(q.T@q-np.eye(k), 2)
    checks = {
        "orthonormal_columns": orth_error < 2e-12,
        "all_greedy_tail_inequalities": worst_tail_error < 2e-12,
        "all_parseval_pivot_floors": worst_floor_error < 2e-12,
        "weighted_energy_lemma": np.max(column_energies) <= .5+2e-12,
        "main_sharp_upper_bound": math.log(gamma) <= upper_log+2e-12,
        "dynamic_range_upper_bound": math.log(gamma) <= range_upper_log+2e-12,
    }
    return {
        "rows": n, "rank": k,
        "coherence": float(coherence), "inverse_norm": float(gamma),
        "log_inverse_norm": math.log(gamma),
        "main_upper_log": upper_log,
        "dynamic_range_upper_log": range_upper_log,
        "main_log_slack": upper_log-math.log(gamma),
        "dynamic_range_log_slack": range_upper_log-math.log(gamma),
        "orthonormality_spectral_error": float(orth_error),
        "maximum_tail_violation": worst_tail_error,
        "maximum_pivot_floor_violation": float(worst_floor_error),
        "maximum_weighted_offdiagonal_column_energy": float(np.max(column_energies)),
        "checks": checks, "passed": all(checks.values()),
    }


def random_parseval_checks() -> dict:
    rng = np.random.default_rng(SEED)
    records = []
    for n,k in ((7,2),(11,4),(32,8),(64,16),(128,32),(256,64)):
        for repeat in range(4):
            q,_ = la.qr(rng.standard_normal((n,k)), mode="economic")
            record = inverse_norm_and_bounds(q)
            record["replicate"] = repeat
            records.append(record)
    # Exact constant-coherence frame in floating arithmetic, with pivot ties.
    k = 16
    q = np.tile(np.eye(k),(4,1))/2
    record = inverse_norm_and_bounds(q)
    record["case"] = "four_repeated_scaled_identity_bases"
    records.append(record)
    return {"seed": SEED, "records": records,
            "passed": all(x["passed"] for x in records)}


def make_core(k: int, h: float, delta: float):
    hbar = h/(1+delta)
    if hbar <= 1:
        raise ValueError("Need h/(1+delta)>1")
    r = parameter_r(hbar)
    u = np.arange(k,dtype=float)/k
    split = 1-1/r
    t = np.where(u <= split, r*u, np.log(hbar/(1-u)))
    # Smallest power of two at least k^3; exact integer arithmetic.
    m = 1 << (k**3-1).bit_length()
    d = np.sqrt(h*k/m)*np.exp(-t/2)
    c = np.sqrt(-np.expm1(-np.diff(t)))
    taper = np.exp(-np.arange(k,dtype=float)/(k*k))
    core = np.diag(d)
    for i in range(k-1):
        core[i,i+1:] = -d[i]*c[i]
    core *= taper[None,:]
    return core, d, c, taper, m, r


def completed_family_check(k: int, h: float, delta: float, cap: float) -> dict:
    core,d,c,taper,m,r = make_core(k,h,delta)
    gram_missing = np.eye(k)-core@core.T
    vals, vecs = la.eigh(gram_missing)
    if vals[0] <= 0:
        return {"rank": k, "h": h, "delta": delta, "coherence_cap":cap,
                "smallest_completion_eigenvalue":float(vals[0]),
                "passed":False, "failure":"Identity completion is not positive definite"}
    b = (vecs*np.sqrt(vals)[None,:])@vecs.T
    # The first k rows of a Sylvester Hadamard matrix have only L distinct
    # columns, each repeated m/L times.  This representation is exact.
    ell = 1 << (k-1).bit_length()
    patterns = la.hadamard(ell, dtype=float)[:k,:]
    filler_unique = b@patterns/math.sqrt(m)
    repetitions = m//ell
    full_gram = core@core.T + repetitions*(filler_unique@filler_unique.T)
    n = k+m
    gram_error = la.norm(full_gram-np.eye(k),2)
    row_norm_sq_max = max(float(np.max(np.sum(core*core,axis=0))),
                          float(np.max(np.sum(filler_unique*filler_unique,axis=0))))
    coherence = n/k*row_norm_sq_max
    margins = []
    filler_margins = []
    replayed = []
    for i in range(k):
        # Independent direct Euclidean residual scores. The prescribed prefix
        # spans precisely the first i coordinates because the core is triangular.
        scores = np.sum(core[i:,i:]**2,axis=0)
        filler_scores = np.sum(filler_unique[i:,:]**2,axis=0)
        intended = float(scores[0])
        second_core = float(np.max(scores[1:])) if len(scores)>1 else -math.inf
        filler_max = float(np.max(filler_scores))
        competitor = max(second_core,filler_max)
        margins.append((intended-competitor)/intended)
        filler_margins.append((intended-filler_max)/intended)
        replayed.append(int(np.argmax(np.concatenate([scores,filler_scores]))))
    # Duplicates have identical residual scores. Thus ordinary QRCP on this
    # compressed array is a separate replay of all first k pivots of the full
    # matrix, provided each pivot is from the prescribed core.
    _,_,piv = la.qr(np.concatenate([core,filler_unique],axis=1),
                   mode="economic",pivoting=True)
    replay_ok = bool(np.array_equal(piv[:k],np.arange(k)))
    inv = la.solve_triangular(core,np.eye(k),lower=False)
    gamma = float(la.norm(inv,2))
    formula_log = math.log(c[0]) + float(np.sum(np.log1p(c[1:]))) - math.log(d[-1])
    entry = float(inv[0,-1])
    entry_rel_error = abs(math.expm1(math.log(entry)-formula_log))
    # The dimension here includes the repeated filler columns.
    upper_log = .5 + .5*math.log(n-k+1)+sharp_constant(coherence)*math.sqrt(k)
    actual_diag = np.diag(core)
    spread = math.log(actual_diag[0]**2/actual_diag[-1]**2)
    range_upper_log = .5-math.log(actual_diag[-1])+math.sqrt((k-1)*spread)
    normalized_entry_log = formula_log-.5*math.log(n-k+1)
    checks = {
        "completion_positive_definite": bool(vals[0] > 0),
        "compressed_parseval_identity": bool(gram_error < 2e-11),
        "coherence_below_cap": bool(coherence < cap),
        "all_direct_residual_pivot_margins_strict": bool(min(margins)>0),
        "all_prescribed_residual_maxima": bool(all(x==0 for x in replayed)),
        "independent_scipy_qrcp_replay": replay_ok,
        "inverse_entry_product_vs_direct_solve": bool(entry_rel_error < 2e-10),
        "sharp_upper_bound": bool(math.log(gamma) <= upper_log+2e-10),
        "dynamic_range_upper_bound": bool(math.log(gamma) <= range_upper_log+2e-10),
    }
    return {
        "rank": k,"h": h,"delta":delta,"coherence_cap":cap,
        "effective_profile_coherence":h/(1+delta),"profile_parameter_r":r,
        "profile_asymptotic_constant":sharp_constant(h/(1+delta)),
        "total_rows":n,"completion_columns":m,
        "distinct_completion_columns":ell,"completion_repetitions":repetitions,
        "coherence":coherence,"smallest_completion_eigenvalue":float(vals[0]),
        "parseval_spectral_error":float(gram_error),
        "minimum_relative_pivot_margin":min(margins),
        "minimum_relative_margin_over_fillers":min(filler_margins),
        "all_relative_pivot_margins":margins,
        "inverse_norm_via_triangular_solve":gamma,
        "inverse_entry_direct":entry,"inverse_entry_product_log":formula_log,
        "inverse_entry_relative_error":entry_rel_error,
        "log_normalized_inverse_entry_div_sqrt_rank":normalized_entry_log/math.sqrt(k),
        "sharp_upper_log":upper_log,"dynamic_range_upper_log":range_upper_log,
        "checks":checks,"passed":all(checks.values()),
    }


def completed_family_checks() -> dict:
    records = []
    negative_controls = []
    # Generous completion slack makes these finite matrices useful actual
    # algorithm checks rather than merely eventually valid asymptotic inputs.
    for cap,h,delta in ((2.,1.8,.2),(4.,3.6,.5),(1.5,1.4,.1)):
        for k in (16,32,64,128):
            record = completed_family_check(k,h,delta,cap)
            if cap == 1.5 and k == 16:
                # This actual failure was found in the initial run, which is
                # retained in results/initial-failed-run.json.  It is
                # outside the theorem's sufficiently-large-rank guarantee.
                # Keep every failed assertion visible, and test that two
                # independent pivot checks agree in rejecting this finite case.
                record["classification"] = "Observed finite negative control; not an admissible lower-bound witness"
                record["explanation"] = (
                    "At rank 16 the taper can erase the delta=0.1 margin: "
                    "exp(-1/k)*sqrt(1+delta)<1.  A completion column overtakes "
                    "the proposed final pivot.  Ranks 32, 64 and 128 are checked separately."
                )
                record["negative_control_passed"] = bool(
                    record["minimum_relative_pivot_margin"] < 0
                    and not record["checks"]["independent_scipy_qrcp_replay"]
                    and record["checks"]["compressed_parseval_identity"])
                negative_controls.append(record)
            else:
                records.append(record)
    return {"records":records,"negative_controls":negative_controls,
            "passed":all(x["passed"] for x in records)
                     and all(x["negative_control_passed"] for x in negative_controls)}


def high_precision_asymptotic_checks() -> dict:
    mp.mp.dps = 80
    ld = np.longdouble
    ranks = (128,512,2048,8192,32768,131072,524288,1_000_000)
    records = []
    for hbar_text in ("1.01","1.5","8.0"):
        hbar_mp = mp.mpf(hbar_text)
        target = 1+mp.log(hbar_mp)
        r_mp = mp.findroot(lambda z:z-mp.log(z)-target,
                          (mp.mpf("1.000001"),max(mp.mpf(2),2*target)))
        r = ld(str(r_mp)); hbar = ld(hbar_text)
        constant_mp = mp.sqrt(r_mp)+1/mp.sqrt(r_mp)
        for k in ranks:
            u = np.arange(k,dtype=ld)/ld(k)
            split = 1-1/r
            t = np.where(u<=split,r*u,np.log(hbar/(1-u)))
            c = np.sqrt(-np.expm1(-np.diff(t)))
            sum_c = np.sum(c,dtype=ld)
            log_product = np.log(c[0])+np.sum(np.log1p(c[1:]),dtype=ld)
            # Independent 80-digit evaluation of sum(c_i), using the constant
            # exponential segment and an exact finite inverse-square-root sum.
            # This uses Hurwitz zeta only as a stable summation formula.
            split_mp = 1-1/r_mp
            last_front_index = int(mp.floor(k*split_mp))
            last_front_index = max(0,min(k-1,last_front_index))
            front_terms = last_front_index
            front_c = mp.sqrt(-mp.expm1(-r_mp/k))
            exact_sum = front_terms*front_c
            boundary_i = last_front_index
            if boundary_i < k-1:
                t_left = r_mp*boundary_i/k
                t_right = mp.log(hbar_mp/(1-mp.mpf(boundary_i+1)/k))
                exact_sum += mp.sqrt(-mp.expm1(-(t_right-t_left)))
                tail_top = k-boundary_i-1
                if tail_top >= 2:
                    exact_sum += mp.zeta(mp.mpf("0.5"),2)-mp.zeta(mp.mpf("0.5"),tail_top+1)
            relative = abs(mp.mpf(str(sum_c))-exact_sum)/max(1,abs(exact_sum))
            records.append({
                "effective_profile_coherence":hbar_text,"rank":k,
                "parameter_r_80_digits":str(r_mp),
                "limit_constant_80_digits":str(constant_mp),
                "sum_c_longdouble":str(sum_c),
                "sum_c_80_digit_independent":str(exact_sum),
                "relative_disagreement":str(relative),
                "sum_c_over_sqrt_rank":str(sum_c/np.sqrt(ld(k))),
                "log_dk_inverse_entry_over_sqrt_rank":str(log_product/np.sqrt(ld(k))),
                "passed":bool(relative < mp.mpf("2e-14")),
            })
    return {
        "scope":"Analytic selected-core growth only; these large full frames are not materialised or QR-factorised.",
        "vector_arithmetic":str(np.finfo(np.longdouble)),
        "independent_sum_arithmetic":"mpmath, 80 decimal digits",
        "records":records,"passed":all(x["passed"] for x in records),
    }


def sign_vectors(k: int) -> np.ndarray:
    return np.array(list(itertools.product((-1.,1.),repeat=k)),dtype=float).T


def radial_sign_phi(t: np.ndarray, signs: np.ndarray) -> np.ndarray:
    k = t.shape[0]
    x = la.expm(t/2)@signs
    x /= np.sqrt(np.sum(x*x,axis=0))[None,:]
    return (k/signs.shape[1])*(x@x.T)


def endpoint_checks() -> dict:
    """Finite support for the analytic implicit-function endpoint argument."""
    rng = np.random.default_rng(SEED+1)
    derivatives = []
    for k in (3,4,5):
        signs = sign_vectors(k)
        raw = rng.standard_normal((k,k))
        s = (raw+raw.T)/2
        s -= np.trace(s)/k*np.eye(k)
        s /= la.norm(s,"fro")
        diagonal = np.diag(np.diag(s))
        expected = diagonal+(1-2/k)*(s-diagonal)
        errors = []
        for step in (1e-4,1e-5,1e-6):
            observed = (radial_sign_phi(step*s,signs)
                        -radial_sign_phi(-step*s,signs))/(2*step)
            errors.append({"step":step,
                           "frobenius_error":float(la.norm(observed-expected,"fro")),
                           "trace_derivative":float(np.trace(observed))})
        derivatives.append({"rank":k,"trace_zero_direction":s.tolist(),
                            "predicted_offdiagonal_multiplier":1-2/k,
                            "finite_differences":errors,
                            "passed":all(x["frobenius_error"]<2e-8 for x in errors)})
    # A true singular direction at rank two: radial normalisation preserves
    # all four sign directions for exp(t*S/2), so Phi is constant along it.
    s2 = np.array([[0.,1.],[1.,0.]])/math.sqrt(2)
    signs2 = sign_vectors(2)
    negative_derivative = (radial_sign_phi(1e-4*s2,signs2)
                           -radial_sign_phi(-1e-4*s2,signs2))/(2e-4)
    negative_phi_change = la.norm(radial_sign_phi(.2*s2,signs2)-np.eye(2),"fro")
    negative = {"rank":2,"direction":s2.tolist(),
                "finite_difference_derivative_norm":float(la.norm(negative_derivative,"fro")),
                "phi_change_at_parameter_point_two":float(negative_phi_change),
                "expected":"Zero offdiagonal derivative; local inverse theorem unavailable in this direction",
                "passed":bool(la.norm(negative_derivative,"fro")<1e-10
                               and negative_phi_change<1e-12)}

    k, a, copies = 4, .5, 128
    signs = sign_vectors(k)
    n = k+copies*signs.shape[1]
    d = np.sqrt((k/n)*(np.arange(k,0,-1)/k)**a)
    c = np.sqrt(1-(d[1:]/d[:-1])**2)
    core = np.diag(d)
    for i in range(k-1):
        core[i,i+1:] = -d[i]*c[i]
    target = (n/(n-k))*(np.eye(k)-core@core.T)
    pairs = [(i,j) for i in range(k) for j in range(i+1,k)]
    dimension = k-1+len(pairs)

    def unpack(x):
        t = np.zeros((k,k))
        t[np.arange(k-1),np.arange(k-1)] = x[:k-1]
        t[-1,-1] = -np.sum(x[:k-1])
        for value,(i,j) in zip(x[k-1:],pairs):
            t[i,j]=t[j,i]=value
        return t

    def independent_coordinates(matrix):
        return np.array([matrix[i,i] for i in range(k-1)]
                        +[matrix[i,j] for i,j in pairs])

    def equation(x):
        return independent_coordinates(radial_sign_phi(unpack(x),signs)-target)

    solution = opt.root(equation,np.zeros(dimension),method="hybr",
                        options={"xtol":1e-11,"maxfev":2000})
    t = unpack(solution.x)
    phi_residual = la.norm(radial_sign_phi(t,signs)-target,2)
    directions = la.expm(t/2)@signs
    directions /= np.sqrt(np.sum(directions*directions,axis=0))[None,:]
    fillers = math.sqrt(k/n)*directions
    gram_error = la.norm(core@core.T+copies*(fillers@fillers.T)-np.eye(k),2)
    all_column_norms = np.concatenate([np.sum(core*core,axis=0),
                                      np.sum(fillers*fillers,axis=0)])
    equal_norm_error = float(np.max(np.abs(all_column_norms-k/n)))
    stages = []
    for i in range(k):
        scores = np.sum(core[i:,i:]**2,axis=0)
        filler_scores = np.sum(fillers[i:,:]**2,axis=0)
        stages.append({"stage":i+1,
                       "normalized_core_residual":float(n*scores[0]),
                       "normalized_max_filler_residual":float(n*np.max(filler_scores)),
                       "normalized_filler_margin":float(n*(scores[0]-np.max(filler_scores))),
                       "relative_core_tie_spread":float((np.max(scores)-np.min(scores))/scores[0])})
    endpoint_witness = {
        "rank":k,"profile_power":a,"copies_per_sign":copies,"total_rows":n,
        "root_solver_success":bool(solution.success),"root_solver_message":str(solution.message),
        "root_solver_evaluations":int(solution.nfev),"solution_trace_zero_matrix":t.tolist(),
        "local_inverse_residual_spectral_norm":float(phi_residual),
        "compressed_parseval_spectral_error":float(gram_error),
        "maximum_equal_column_norm_error":equal_norm_error,
        "coherence":float(n/k*np.max(all_column_norms)),"stages":stages,
        "passed":bool(solution.success and phi_residual<2e-11 and gram_error<2e-11
                      and equal_norm_error<2e-15
                      and all(x["normalized_filler_margin"]>0 for x in stages[1:])
                      and all(x["relative_core_tie_spread"]<2e-13 for x in stages)),
    }
    return {
        "scope":"The unrestricted existence argument uses the analytic inverse function theorem. These are finite derivative checks and one floating-point local solve, not a validated numerical proof or a uniform solver guarantee.",
        "tie_convention":"Smallest original column label at mathematical ties; first-stage ties are forced by coherence one.",
        "derivative_checks":derivatives,"rank_two_negative_control":negative,
        "finite_endpoint_witness":endpoint_witness,
        "passed":all(x["passed"] for x in derivatives) and negative["passed"]
                 and endpoint_witness["passed"],
    }


def finite_field_sign_design(q: int, modulus: int) -> np.ndarray:
    """Exact polynomial-basis field arithmetic for the specified small field."""
    degree = q.bit_length()-1

    def multiply(left: int, right: int) -> int:
        value = 0
        while right:
            if right & 1:
                value ^= left
            left <<= 1
            if left & q:
                left ^= modulus
            right >>= 1
        return value

    def trace(value: int) -> int:
        answer = 0
        power = value
        for _ in range(degree):
            answer ^= power
            power = multiply(power,power)
        if answer not in (0,1):
            raise ArithmeticError("The specified modulus did not give a binary field trace")
        return answer

    cube = [multiply(multiply(x,x),x) for x in range(q)]
    base = np.array([
        [1-2*trace(multiply(a,x)^multiply(b,cube[x])) for a in range(q) for b in range(q)]
        for x in range(q)],dtype=np.int64)
    return np.concatenate((base,-base),axis=1)


def quantitative_endpoint_checks() -> dict:
    """Targeted new-design checks; the full huge completion is never formed."""
    rng = np.random.default_rng(SEED+1701)
    records = []
    for q,modulus in ((4,0b111),(8,0b1011)):
        k = q
        signs = finite_field_sign_design(q,modulus)
        m = signs.shape[1]
        moment_records = []
        for order in range(1,5):
            sums = [int(np.sum(np.prod(signs[list(subset)],axis=0,dtype=np.int64),
                                   dtype=np.int64))
                    for subset in itertools.combinations(range(k),order)]
            moment_records.append({"order":order,"subset_count":len(sums),
                                   "max_absolute_unnormalized_moment":max(map(abs,sums),default=0),
                                   "all_zero":all(value==0 for value in sums)})

        direction = rng.normal(size=(k,k))
        direction = (direction+direction.T)/2
        direction -= np.trace(direction)/k*np.eye(k)
        direction /= la.norm(direction,"fro")
        diagonal = np.diag(np.diag(direction))
        predicted = diagonal+(1-2/k)*(direction-diagonal)
        step = 1e-5
        observed = (radial_sign_phi(step*direction,signs)
                    -radial_sign_phi(-step*direction,signs))/(2*step)
        derivative_error = float(la.norm(observed-predicted,"fro"))

        a = 0.5
        copies = math.ceil(2304*k**3/((1-a)*m))
        n = k+copies*m
        d = np.sqrt((k/n)*(np.arange(k,0,-1)/k)**a)
        c = np.sqrt(1-(d[1:]/d[:-1])**2)
        core = np.diag(d)
        for i in range(k-1):
            core[i,i+1:] = -d[i]*c[i]
        target = n/(n-k)*(np.eye(k)-core@core.T)
        perturbation = target-np.eye(k)

        def inverse_derivative(value):
            diag = np.diag(np.diag(value))
            result = diag+(value-diag)/(1-2/k)
            # Analytically the map preserves trace zero; remove roundoff only.
            result -= np.trace(result)/k*np.eye(k)
            return result

        t = np.zeros((k,k))
        steps = []
        for iteration in range(12):
            phi = radial_sign_phi(t,signs)
            linear = np.diag(np.diag(t))+(1-2/k)*(t-np.diag(np.diag(t)))
            updated = inverse_derivative(perturbation-(phi-np.eye(k)-linear))
            steps.append({"iteration":iteration+1,"increment_frobenius_norm":float(la.norm(updated-t,"fro"))})
            t = updated
            if steps[-1]["increment_frobenius_norm"]<1e-14:
                break

        phi_residual = float(la.norm(radial_sign_phi(t,signs)-target,2))
        directions = la.expm(t/2)@signs
        directions /= np.sqrt(np.sum(directions*directions,axis=0))[None,:]
        fillers = math.sqrt(k/n)*directions
        gram_error = float(la.norm(core@core.T+copies*(fillers@fillers.T)-np.eye(k),2))
        norms = np.concatenate((np.sum(core*core,axis=0),np.sum(fillers*fillers,axis=0)))
        equal_norm_error = float(np.max(np.abs(norms-k/n)))
        later_margins = [float(n*(np.sum(core[i:,i]**2)-np.max(np.sum(fillers[i:,:]**2,axis=0))))
                         for i in range(1,k)]
        perturbation_norm = float(la.norm(perturbation,"fro"))
        radius_ratio = float(la.norm(t,"fro")/(6*perturbation_norm))
        margin_bound = 31*(1-a)/64
        record = {
            "rank":k,"field_order":q,"modulus_binary":bin(modulus),
            "design_columns_with_multiplicity":m,"moment_arithmetic":"exact int64 sums",
            "mixed_moments":moment_records,
            "derivative_step":step,"derivative_frobenius_error":derivative_error,
            "profile_power":a,"copies_per_design_column":copies,"total_rows":n,
            "perturbation_frobenius_norm":perturbation_norm,
            "contraction_threshold":1/(1152*k),
            "iterations":steps,"solution_over_certified_radius":radius_ratio,
            "local_inverse_residual_spectral_norm":phi_residual,
            "compressed_parseval_spectral_error":gram_error,
            "maximum_equal_column_norm_error":equal_norm_error,
            "later_scaled_squared_pivot_margins":later_margins,
            "proved_scaled_margin_lower_bound":margin_bound,
            "passed":bool(all(x["all_zero"] for x in moment_records)
                          and derivative_error<2e-8
                          and perturbation_norm<=1/(1152*k)
                          and radius_ratio<=1+1e-10
                          and phi_residual<2e-12 and gram_error<2e-12
                          and equal_norm_error<2e-15
                          and min(later_margins)>margin_bound),
        }
        records.append(record)
    return {
        "scope":"Exact mixed moments through order four for the F4 and F8 designs, followed by finite-difference and floating-point contraction checks. Compressed multiplicities preserve the full covariance exactly in the formula; only k by (k+m) columns are materialised. This does not validate the analytic constants for all dimensions.",
        "records":records,"passed":all(x["passed"] for x in records),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",default="results/verification.json")
    args = parser.parse_args()
    start = time.perf_counter()
    results = {
        "status":"Finite evidence only; no proof-assistant verification has been performed.",
        "run_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "seed":SEED,
        "versions":{"python":sys.version,"platform":platform.platform(),
                    "numpy":np.__version__,"scipy":scipy.__version__,
                    "sympy":sp.__version__,"mpmath":mp.__version__},
    }
    suites = (("exact_rational",exact_rational_checks),
              ("random_parseval",random_parseval_checks),
              ("completed_lower_family",completed_family_checks),
              ("large_rank_core_asymptotics",high_precision_asymptotic_checks),
              ("coherence_one_endpoint",endpoint_checks),
              ("quantitative_endpoint_design",quantitative_endpoint_checks))
    for name,fn in suites:
        suite_start = time.perf_counter()
        try:
            record = fn()
        except Exception as exc:
            record = {"passed":False,"error_type":type(exc).__name__,"error":str(exc)}
        record["elapsed_seconds"] = time.perf_counter()-suite_start
        results[name] = record
        print(f"{name}: {'PASS' if record['passed'] else 'FAIL'} ({record['elapsed_seconds']:.3f} s)",flush=True)
    results["passed"] = all(results[name]["passed"] for name,_ in suites)
    results["elapsed_seconds"] = time.perf_counter()-start
    output = Path(args.output)
    output.write_text(json.dumps(results,indent=2,allow_nan=False,
                                 default=lambda x:x.item() if isinstance(x,np.generic)
                                 else str(x))+"\n",encoding="utf-8")
    print(f"Recorded all outcomes in {output.resolve()}",flush=True)
    return 0 if results["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
