#!/usr/bin/env python3
"""Exact identification fingerprints for a known rational Kahan ingredient.

These finite fingerprints identify the classical triangular core used in the
construction. They are not evidence that the sharp asymptotic theorem is new.

Run from any directory:
    python path/to/code/fingerprints.py
Or specify a destination:
    python path/to/code/fingerprints.py --output fingerprints.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


C = sp.Rational(3, 5)
S = sp.Rational(4, 5)
X = sp.Symbol("x")


def kahan_core(k: int) -> sp.Matrix:
    return sp.Matrix(k, k, lambda i, j:
                     S**i if i == j else -C*S**i if i < j else 0)


def exact_matrix(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[i, j]) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


def canonical_text(matrix: sp.Matrix) -> str:
    payload = {"rows": matrix.rows, "columns": matrix.cols,
               "entries": exact_matrix(matrix)}
    return json.dumps(payload, sort_keys=True, ensure_ascii=True,
                      separators=(",", ":"))


def matrix_hash(matrix: sp.Matrix) -> str:
    return hashlib.sha256(canonical_text(matrix).encode("utf-8")).hexdigest()


def fingerprint(k: int) -> dict:
    matrix = kahan_core(k)
    reversal = sp.zeros(k)
    for i in range(k):
        reversal[i, k-1-i] = 1
    transforms = {
        "identity": matrix,
        "transpose": matrix.T,
        "row_reversal": reversal*matrix,
        "column_reversal": matrix*reversal,
        "row_and_column_reversal": reversal*matrix*reversal,
        "anti_diagonal_transpose": reversal*matrix.T*reversal,
    }
    poly = matrix.charpoly(X).as_poly()
    entry_counts = {}
    for value in sorted(set(matrix)):
        entry_counts[str(value)] = sum(value == entry for entry in matrix)
    return {
        "order": k,
        "matrix": exact_matrix(matrix),
        "trace": str(sp.trace(matrix)),
        "determinant": str(matrix.det()),
        "characteristic_polynomial": {
            "variable": "x",
            "expanded": str(poly.as_expr()),
            "factorised": str(sp.factor(poly.as_expr())),
            "coefficients_descending": [str(value) for value in poly.all_coeffs()],
        },
        "eigenvalues": [
            {"value": str(value), "algebraic_multiplicity": int(multiplicity)}
            for value, multiplicity in sorted(matrix.eigenvals().items())
        ],
        "row_sums": [str(sum(matrix.row(i))) for i in range(k)],
        "column_sums": [str(sum(matrix.col(j))) for j in range(k)],
        "sorted_entries_with_repetitions": [str(value) for value in sorted(matrix)],
        "entry_multiplicities_in_increasing_order": entry_counts,
        "rank": int(matrix.rank()),
        "nonzero_pattern": [[int(matrix[i, j] != 0) for j in range(k)]
                            for i in range(k)],
        "zero_positions_one_based": [[i+1, j+1] for i in range(k)
                                     for j in range(k) if matrix[i, j] == 0],
        "canonical_matrix_serialisation": canonical_text(matrix),
        "exact_matrix_sha256": {name: matrix_hash(value)
                                for name, value in transforms.items()},
        "inverse_top_right_entry": str(matrix.inv()[0, k-1]),
    }


def make_fingerprints() -> dict:
    sequence = [{"rank": k,
                 "inverse_top_right_entry": str(kahan_core(k).inv()[0, k-1])}
                for k in range(2, 10)]
    return {
        "schema_version": 1,
        "purpose": "Known-object identification fingerprints for the novelty gate",
        "epistemic_status": (
            "The matrices and geometric inverse-entry sequence are classical "
            "Kahan ingredients. These finite fingerprints make no originality "
            "claim and do not establish novelty of the sharp coherence theorem."
        ),
        "arithmetic": "Exact SymPy rational arithmetic",
        "sympy_version": sp.__version__,
        "object_definition": {
            "family": "Untapered upper-triangular Kahan core K=D T",
            "c": str(C),
            "s": str(S),
            "relation": "c^2+s^2=1",
            "indices": "One-based mathematical indices",
            "entries": "K_ii=s^(i-1); K_ij=-c*s^(i-1) for i<j; K_ij=0 for i>j",
            "distinction": "These are the known untapered cores, not the completed Parseval examples.",
        },
        "hash_specification": {
            "algorithm": "SHA-256",
            "encoding": "UTF-8",
            "serialisation": (
                "JSON object with keys rows, columns and entries, sorted keys, "
                "ASCII strings, no whitespace, and every exact rational in "
                "reduced SymPy string form. The original serialisation is included."
            ),
            "transform_definitions": {
                "identity": "K",
                "transpose": "K^T",
                "row_reversal": "J K",
                "column_reversal": "K J",
                "row_and_column_reversal": "J K J",
                "anti_diagonal_transpose": "J K^T J",
                "J": "Order-k reversal permutation matrix J_ij=1 iff i+j=k+1",
            },
            "limitation": "Hashes certify equality of canonical finite encodings only.",
        },
        "fingerprints": [fingerprint(k) for k in (2, 3, 4)],
        "known_inverse_entry_sequence": {
            "entry": "a_k=(K_k^(-1))_(1,k), k>=2",
            "formula": "c*(1+c)^(k-2)/s^(k-1)=(3/4)*2^(k-2)",
            "values": sequence,
            "constant_ratio": "2",
            "ordinary_generating_functions": {
                "shifted_index": {
                    "definition": "sum_(n>=0) a_(n+2) x^n",
                    "rational_function": "3/(4*(1-2*x))",
                },
                "rank_index": {
                    "definition": "sum_(k>=2) a_k x^k",
                    "rational_function": "3*x^2/(4*(1-2*x))",
                },
                "interpretation": "Formal power series; analytically convergent for |x|<1/2.",
            },
            "novelty_status": "Known geometric sequence from the classical Kahan inverse identity.",
        },
    }


def main() -> None:
    default = Path(__file__).resolve().parents[1] / "fingerprints.json"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=default)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(make_fingerprints(), indent=2,
                                     ensure_ascii=True)+"\n", encoding="utf-8")
    print(f"Wrote exact known-object fingerprints to {args.output.resolve()}")


if __name__ == "__main__":
    main()
