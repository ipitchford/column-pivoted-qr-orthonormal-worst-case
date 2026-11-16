# Replay receipt (release 0.1.0-candidate)

Replayed by the Evidence Press publisher on 16 November 2026 from fresh copies of the approved candidate, working
directory the package copy, CPython 3.12.11 (uv; `sys.version` printed inside the run; invoked as `python`), numpy
2.3.5, scipy 1.17.0, sympy 1.14.0, mpmath 1.3.0, matplotlib 3.10.8, `PYTHONDONTWRITEBYTECODE=1`, macOS arm64.

| Run | Command | Result |
|---|---|---|
| Verification | `python -I -B code/verify.py --output results/verification.json` | exit 0; six suites PASS (exact_rational, random_parseval, completed_lower_family, large_rank_core_asymptotics, coherence_one_endpoint, quantitative_endpoint_design); output identical to the recorded file apart from `run_utc` and `elapsed_seconds` |
| Verification, optimised | same with `python -O` | exit 0; same |
| Extensions | `python -I -B code/check_extensions.py --output results/extensions.json` | exit 0, "ALL CHECKS PASSED"; byte-identical (crossovers 13, 15, 20, 25; grid of 779 points, smallest slack 1.4361037; Kahan witnesses for k = 6 to 16) |
| Extensions, optimised | same with `python -O` | exit 0; byte-identical |
| Fingerprints | `python -I -B code/fingerprints.py --output results/fingerprints.json` | exit 0; byte-identical |
| Fingerprints, optimised | same with `python -O` | exit 0; byte-identical |
| Figure | `python -E -s -B code/make_figure.py` (imports `verify.py` from `code/`, hence not `-I`) | exit 0; `figures/sharp-growth.png` and `figures/plot-data.json` byte-identical; `sharp-growth.pdf` identical apart from its creation date |
| Negative control | `code/check_extensions.py` with the constant 3 in inequality (7.5) changed to 1, under `-O` | exit 1, "FAILED: 1 checks" (inequality (7.5) on grid) |
| Built-in negative control | `verify.py` rebuilds the small candidate (k = 16, coherence cap 1.5) whose failure is recorded in `results/initial-failed-run.json` | rejected by both QRCP implementations in every run ("negative_control_passed": true) |

After the replay only the title-page version/DOI line of the manuscript (and the rebuilt PDF), README and metadata files
were added or changed; the scripts and recorded outputs are byte-identical to the replayed candidate.

Cross-platform policy (declared before tagging; also in CI, README.md and ENVIRONMENT.txt): `extensions.json`,
`fingerprints.json` and `initial-failed-run.json` byte-identical; `verification.json` field by field, with structure,
booleans, integers and non-numeric strings exact, JSON numbers within relative 1e-9 (absolute floor 1e-12), decimal
numeric strings within relative 1e-12, `relative_disagreement` fields below 2e-14, and time, interpreter, platform and
root-solver bookkeeping fields ignored. The figure is regenerated on Linux in a non-blocking job and not compared.
