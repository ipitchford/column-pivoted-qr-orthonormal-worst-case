# How badly can column-pivoted QR choose rows of an orthonormal basis? The sharp stretched-exponential worst case

*Greedy row selection loses at most exp(C(H)√k) against the best possible, with an exact constant, not 2^k.*

Version 0.1.0-candidate. DOI [10.5281/zenodo.24278831](https://doi.org/10.5281/zenodo.24278831) (all versions:
[10.5281/zenodo.24278830](https://doi.org/10.5281/zenodo.24278830)). Evidence Press page:
https://evidencepress.org/releases/column-pivoted-qr-orthonormal-worst-case/. Creator: Anonymous. Status: unrefereed
research candidate.

## What this is about

**The everyday problem.** Many simulations and data methods keep only a few "measurement points" out of thousands.
Suppose you have k important patterns, for example the main modes of a fluid flow sampled at N grid points, and you want
to read just k of the N values and still reconstruct the rest. The standard way to pick the points is a greedy rule built
into every linear-algebra library, column-pivoted QR (QRCP): repeatedly pick the point that the points already chosen
explain least. This is, for example, how the interpolation points of the discrete empirical interpolation method (DEIM)
are often chosen.

**How quality is measured.** One number, γ, measures a choice of points: the factor by which reconstruction from the
chosen points can be worse than the best possible reconstruction from the same patterns. Good choices with modest γ
always exist (γ at most about the square root of k times N). The question is how bad the greedy choice can be.

**What was known.** A textbook guarantee (Drmač and Gugercin, 2016) only said that γ grows at most like 2^k, which is
astronomically large already at k = 100 (about 10^30). A list of open problems from a Simons workshop (Amsel et al.,
Problem 4.3, posed by Damle and Kressner) asked whether the greedy rule does as well as the best choices, up to a
polynomial factor. Chen, Liu, He and Dong (arXiv:2606.02484) answered first, and negatively: their examples grow faster
than any polynomial, and their construction already gives growth of order exp(c√k) with a small constant c. Damle
(arXiv:2607.13532) showed that a modified greedy rule (Bischof–Stewart pivoting) does meet the existence scale.

**What is new.**

1. *An upper bound below 2^k.* The greedy rule is never worse than exp(C√k) times the existence scale, with an exact
   constant C that depends on how evenly the energy of the patterns is spread over the points (the "coherence" H): C = 2
   when all points carry equal energy, and C is about 2.25 when no point carries more than twice the average. This is
   the first guarantee for this rule below exponential growth. At k = 100 it allows about a billion (10^9 to 10^10)
   instead of about 10^30, and it beats the textbook bound from k = 13 onwards when all points carry equal energy.
2. *The constant is the right one.* Examples come arbitrarily close to it, so C cannot be lowered.
3. *No assumption needed.* Coherence can never exceed N/k, so the bound holds for every input; with N a fixed multiple
   of k the loss is at most exp(C(N/k)√k).
4. *The law in N and k.* In terms of N and k alone, the worst case is exp(Θ(√(k log N))) once N is a bit larger than
   k², with the exact constant when N is much larger than any power of k.
5. *A free check after the fact.* The pivots that the factorisation already computes give a bound on γ, with an optimal
   constant.

**The idea of the proof.** The greedy rule controls the energy of every selected column on every trailing block of
coordinates at once. Rescaling the coordinates by a carefully shaped weight turns this into a bound on a triangular
matrix with unit diagonal, whose inverse is then bounded by a classical determinant argument. The matching examples use
a classical "Kahan" triangular matrix with a two-part diagonal profile, hidden inside a large orthonormal matrix whose
extra columns are too small to be picked.

**Limits.** The results are about exact arithmetic and about worst cases; typical inputs behave far better, and
floating-point behaviour is not analysed. "Sharp" means the leading term of the exponent: upper and lower bounds still
differ by a factor polynomial in k, and the limits are approached slowly. The worst examples need N of order k³; what
happens when N is a fixed multiple of k is open. At exactly equal energies (H = 1) the matching examples rely on a stated
tie-breaking rule (smallest label first). Practical advice from the paper: compute γ after the selection (it is cheap)
and switch to a rule with a polynomial guarantee if it is large. This is an unrefereed candidate: assurance is internal
AI-assisted review and computational checks, not peer review.

## For specialists

Notation: Q ∈ ℂ^{N×k} with Q*Q = I_k; QRCP applied to Q* selects rows I; γ(Q, I) = ‖Q(I,:)^{-1}‖₂; coherence
μ(Q) = (N/k) max_j ‖Q(j,:)‖₂², so 1 ≤ μ(Q) ≤ N/k. The paper ([manuscript.pdf](manuscript.pdf), source
[manuscript.tex](manuscript.tex)) proves:

- **Theorem 2.3 (sharp growth at bounded coherence).** If μ(Q) ≤ H, every legal QRCP selection, with any tie-breaking,
  satisfies γ ≤ √(e(N−k+1)) exp(C(H)√k), where C(H) = √ρ + 1/√ρ with ρ − log ρ = 1 + log H (so C(1) = 2,
  C(2) ≈ 2.2476, C(H)²/log H → 1). The worst case G_k(H) of γ/√(N−k+1) satisfies log G_k(H) = (C(H) + o(1))√k; the
  matching examples are real, have fewer than 3k³ rows and unique pivots for H > 1, and at H = 1 use the smallest-label
  tie rule.
- **Corollary 2.6 (no coherence assumption).** γ ≤ √(e(N−k+1)) exp(min{C(N/k)√k, √((k−1) log(N−k+1))}) for every input;
  C(H)² ≤ log H + log(1 + log H) + 5.
- **Theorem 2.7 and Corollary 2.8 (the law in N and k).** Without a coherence assumption, the worst case is
  exp(Θ(√(k log(N/k)))) once N ≥ k^{2+ε} and log N = o(k), and exp((1+o(1))√(k log N)) when log N/log k → ∞ with
  log N = o(k).
- **Theorem 7.1 and Proposition 7.2 (a posteriori bound).** With L = log(d₁²/d_k²) from the diagonal of the triangular
  factor, γ ≤ (√e/d_k) exp(√((k−1)L)), and the coefficient one is optimal.

Credit, as the paper states it: Chen, Liu, He and Dong (arXiv:2606.02484) gave the first, negative answer to Problem 4.3
of Amsel et al. (arXiv:2602.05394), and their Lemma 24 already gives growth exp(c√k) for every c < (1/2)√(log H) (about
0.42 at H = 2, against C(2) ≈ 2.25). The new contributions are the matching universal upper bound, the exact constant
C(H), the endpoint H = 1 and the coherence-free law. Engler (1997) proved the pivot floor of Lemma 3.1; the determinant
estimate of Lemma 3.3 is classical, with its own short proof and no novelty claimed; the matching examples build on
Kahan's (1966) triangular matrix.

## Technical layer: files and replay

| Path | What it is |
|---|---|
| `manuscript.tex`, `manuscript.pdf` | The paper (18 pages), with the version and DOI on the title page |
| `figures/sharp-growth.pdf`, `.png`, `figures/plot-data.json` | Figure 1 and its data |
| `code/verify.py` → `results/verification.json` | Six suites: exact rational Kahan identities; 25 orthonormal inputs checked against the inequalities of the paper; eleven completed examples at k = 16 to 128, factorised by two independent QRCP implementations; the large-rank core sums up to k = 10^6 against 80-digit values; the endpoint construction at k = 4; the finite-field design for k = 4, 8 |
| `results/initial-failed-run.json` | Record of the initial run, in which a small candidate (k = 16, coherence cap 1.5) failed its intended pivot order; `verify.py` rebuilds that candidate and checks that two independent QRCP implementations reject it (a negative control) |
| `code/check_extensions.py` → `results/extensions.json` | 50-digit checks: C(H) against direct minimisation, the bound on C(H)², Table 2 and the crossovers k = 13, 15, 20, 25, the Kahan family for k = 6 to 16 with exact QRCP, inequality (7.5) on 779 grid points (smallest slack 1.4361037) |
| `code/fingerprints.py` → `results/fingerprints.json` | Formula and sequence fingerprints used in the prior-art search |
| `code/make_figure.py` | Regenerates the figure and its data |
| `refs.bib` | The 11 references, from Crossref and arXiv records |
| `AI_INDEX.md`, `RESEARCH_GATES.json`, `REVIEW_SUMMARY.md` | Claim map, research-gate record, internal review summary |
| `REPLAY_RECEIPT.md`, `ENVIRONMENT.txt`, `requirements.txt` | Release replay and environment |
| `LICENSE.md`, `CITATION.cff`, `.zenodo.json`, `MANIFEST.sha256` | Licences, citation metadata, SHA-256 of every file |

Replay (package root as working directory; each run takes a few seconds):

```sh
uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with sympy==1.14.0 --with mpmath==1.3.0 python -I -B -O code/verify.py --output results/verification.json
uv run --python 3.12 --with mpmath==1.3.0 python -I -B -O code/check_extensions.py --output results/extensions.json
uv run --python 3.12 --with numpy==2.3.5 --with sympy==1.14.0 --with mpmath==1.3.0 python -I -B -O code/fingerprints.py --output results/fingerprints.json
uv run --python 3.12 --with numpy==2.3.5 --with scipy==1.17.0 --with sympy==1.14.0 --with mpmath==1.3.0 --with matplotlib==3.10.8 python -E -s -B code/make_figure.py
```

`verify.py` prints six PASS lines and `check_extensions.py` prints "ALL CHECKS PASSED". Every check raises an error
without relying on `assert`, so the scripts stay active under `python -O`; both exit non-zero if any check fails.
Invoke the interpreter as `python` inside `uv run` and check `python --version`: on some hosts `python3` is a shim that
resolves to another version.

Cross-platform comparison policy (declared before tagging; the same policy is used in CI and in `REPLAY_RECEIPT.md`):
`results/extensions.json`, `results/fingerprints.json` and `results/initial-failed-run.json` must be byte-identical
(exact and 50-digit mpmath arithmetic; no floating-point numbers). `results/verification.json` is compared field by field with the timing, timestamp, interpreter and
platform fields ignored (`run_utc`, `elapsed_seconds`, `versions`, `vector_arithmetic`, `root_solver_evaluations`,
`root_solver_message`); exact rational strings, integers, booleans and structure must match exactly, numbers within
relative 1e-9 (absolute floor 1e-12) and decimal numeric strings within relative 1e-12, and every `relative_disagreement` field must stay below 2e-14
(the script's own criterion). Reason: the large-rank core sums use NumPy `longdouble`, which is 80-bit on x86-64 Linux and
plain double precision on arm64 macOS, so last digits can differ; the root solver may also use a different number of
evaluations. The figure is regenerated in a non-blocking CI job and not compared.

## How to cite

Anonymous (2026). *How badly can column-pivoted QR choose rows of an orthonormal basis? The sharp stretched-exponential
worst case.* Version 0.1.0-candidate. Evidence Press. https://doi.org/10.5281/zenodo.24278831

## Assurance and licence

Written proofs, with finite computational checks. Internal AI-assisted review only (see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)):
five reviewer roles, an adversarial review by a second AI model, a clean replay, one revision round with extensions, then
a re-review and a proof audit of the new proofs by Claude verification passes. This is not peer review, independent
reproduction or formal verification. Text, data and figures: CC0-1.0. Code: MIT. See [LICENSE.md](LICENSE.md).
