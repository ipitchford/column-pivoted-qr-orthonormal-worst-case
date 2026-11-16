# AI index: How badly can column-pivoted QR choose rows of an orthonormal basis? (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What this
is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is a plain
introduction: the problem, what was known, what is new, what is not claimed and a roadmap). DOI 10.5281/zenodo.24278831
(concept 10.5281/zenodo.24278830). Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

Q in C^{N x k}, Q*Q = I_k, N >= k. Exact QRCP on W = Q*: at each step select an unselected column of W whose residual
(projection onto the orthogonal complement of the selected columns) has the largest norm; any tie may be taken, and
each resulting ordered tuple I is a legal selection. gamma(Q, I) = ||Q(I,:)^{-1}||_2. Coherence mu(Q) = (N/k) max_j
||Q(j,:)||^2 in [1, N/k]. rho(H) >= 1 solves rho - log rho = 1 + log H; C(H) = sqrt(rho) + 1/sqrt(rho) (C(1) = 2,
C(1.5) = 2.1554, C(2) = 2.2476, C(4) = 2.4420, C(10) = 2.6635, C(100) = 3.1256). G_k(H) = sup gamma/sqrt(N-k+1) over
N, Q with mu(Q) <= H and legal selections. Gamma_k(N) = the same supremum at exactly N rows, no coherence assumption.
d_1 >= ... >= d_k = diagonal of the triangular factor; L = log(d_1^2/d_k^2). Natural logarithms. A26 = Amsel et al.,
arXiv:2602.05394; CLHD26 = Chen, Liu, He, Dong, arXiv:2606.02484; D26 = Damle, arXiv:2607.13532; DG16 = Drmac-Gugercin,
SIAM J. Sci. Comput. 38 (2016); E97 = Engler, Appl. Math. Lett. 10 (1997); K66 = Kahan, Canad. Math. Bull. 9 (1966).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A | Greedy trailing-energy constraints: every selected column has bounded energy on every trailing block of coordinates simultaneously; pivot floor d_i^2 >= (k-i+1)/(N-i+1) | Lemma 3.1 | written proof; pivot floor due to E97 |
| B | Weighted-similarity estimate for unit upper-triangular matrices under the trailing constraints | Lemma 3.2 | written proof |
| C | Determinant estimate for the inverse of a triangular matrix from its Frobenius norm | Lemma 3.3 | classical estimate, own short proof, no novelty claimed |
| D | Upper bound: mu(Q) <= H implies gamma <= sqrt(e(N-k+1)) exp(C(H) sqrt k) for every legal selection, real or complex | Theorem 2.3, Section 4 | written proof; first bound below exponential in k for QRCP on orthonormal rows (to the paper's knowledge) |
| E | Sharpness: lim log G_k(H)/sqrt k = C(H) for every fixed H >= 1; examples real, N < 3k^3 and unique pivots for H > 1; at H = 1, N = O(k^3 log k) under the smallest-label tie rule | Theorem 2.3, Sections 5-6, Appendix A | written proofs; generalised Kahan cores (K66) completed to orthonormal matrices; the stretched-exponential order was first reached by CLHD26 (Lemma 24, coefficients below (1/2)sqrt(log H), 0.42 at H = 2) |
| F | No coherence assumption: gamma <= sqrt(e(N-k+1)) exp(min{C(N/k) sqrt k, sqrt((k-1) log(N-k+1))}); C(H)^2 <= log H + log(1 + log H) + 5 | Corollary 2.6 | written proof; bound on C(H)^2 checked on a grid 1..1e12 |
| G | Law in N and k: log Gamma_k(N) <= 1/2 + min{...}; log Gamma_k(N) >= sqrt(kL) - L - log k - 3 with L = log(N/(66k^2)) in [1, (k-1)/2]; (1+o(1)) sqrt(k log N) when log N/log k -> infinity and log N = o(k); sqrt(alpha-2) to sqrt(alpha-1) window at N = k^alpha | Theorem 2.7, Corollary 2.8, Lemma 7.3 | written proofs; inequality (7.5) checked on 779 grid pairs (smallest slack 1.4361037); Kahan family witnesses k = 6..16 checked with exact QRCP |
| H | A posteriori bound gamma <= (sqrt e/d_k) exp(sqrt((k-1)L)); coefficient one optimal | Theorem 7.1, Proposition 7.2 | written proofs |
| I | Practical size: the new factor sqrt(e) exp(C(H) sqrt k) is below the DG16 factor sqrt(4^k + 6k - 1)/3 from k = 13, 15, 20, 25 for H = 1, 2, 10, 100; at k = 100 about 1e9 to 1e10 against about 1e30 | Table 2 | computed (check_extensions.py) |
| J | Finite-k behaviour: realised exponents 1.07, 1.30, 1.49, 1.64 (gamma) and 0.94, 1.19, 1.40, 1.56 (inverse entry) at k = 16, 32, 64, 128 for H = 2; pivot margins about 2/k^2; core coefficient 2.140 against limit 2.155 at k = 1e6 | Section 9, Figure 1 | computed (verify.py); descriptive |

Not claimed: anything about floating-point QRCP, pivoting with a tolerance or perturbed data; typical-case behaviour;
a multiplicative (rather than exponent-level) sharp constant; the worst case at a fixed aspect ratio N <= alpha k (open);
that cubic row counts are necessary. Known before this work: the existence scale and Problem 4.3 (A26); the 2^k bound
(DG16 Theorem 2.1); the pivot floor and Hadamard completions (E97); Kahan's matrix (K66); the first negative answer and
exp(c sqrt k) growth with small c (CLHD26); Bischof-Stewart pivoting meets the existence scale (D26, bound due to Osinsky
as D26 records).

## Evidence files

- [code/verify.py](code/verify.py) -> [results/verification.json](results/verification.json): A, D, E, H, J (six suites;
  includes a negative control: the small candidate whose failure is recorded in
  [results/initial-failed-run.json](results/initial-failed-run.json) is rebuilt and must be rejected).
- [code/check_extensions.py](code/check_extensions.py) -> [results/extensions.json](results/extensions.json): F, G, I and the
  C(H) table (50-digit arithmetic).
- [code/fingerprints.py](code/fingerprints.py) -> [results/fingerprints.json](results/fingerprints.json): formula and
  sequence fingerprints for the prior-art search.
- [code/make_figure.py](code/make_figure.py) -> [figures/sharp-growth.pdf](figures/sharp-growth.pdf),
  [figures/sharp-growth.png](figures/sharp-growth.png), [figures/plot-data.json](figures/plot-data.json): Figure 1.
- [requirements.txt](requirements.txt), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md):
  environment and the release replay (normal and `python -O`, negative control, cross-platform comparison policy).
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): internal review objections and dispositions. [RESEARCH_GATES.json](RESEARCH_GATES.json):
  contribution map, searches, sources read and limitation triage.
- [refs.bib](refs.bib): bibliography from Crossref and arXiv records. [LICENSE.md](LICENSE.md): CC0-1.0 text and data, MIT
  code. [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json): citation metadata. CI:
  [.github/workflows/replay.yml](.github/workflows/replay.yml). [MANIFEST.sha256](MANIFEST.sha256): SHA-256 of every file.

## Limits

No human peer review, independent human reproduction or formal verification. The general statements rest on written
proofs; the computations are finite checks. Exact arithmetic and worst cases only. Sharpness is on the logarithmic scale;
a factor polynomial in k remains at fixed H. The H = 1 lower examples depend on the smallest-label tie rule. Priority
rests on a literature search of limited scope.
