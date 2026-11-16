# Review summary (internal, AI-assisted)

This candidate was checked inside Evidence Press before release. The checks were internal and AI-assisted: five
reviewer roles (editor, methodology, domain, perspective, devil's advocate), an adversarial review by a second AI model,
clean replays of every script, one revision round with extensions, a re-review of the changes, and a proof audit of the
new proofs by Claude verification passes. A second-model audit of the new proofs was planned but did not run, because
that model's usage limit had been reached; the audit of the new proofs therefore rests on the Claude passes. This is not
peer review, independent reproduction or formal verification. No reviewer found a mathematical error in the claims; the
minor points below were applied.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| Chen, Liu, He and Dong's construction (Lemma 24) already gives growth exp(c sqrt k) for every c < (1/2) sqrt(log H) | Stated plainly in the abstract, Section 1 and Remark 2.5; the new contributions are named as the matching upper bound, the exact constant C(H), the endpoint H = 1 and the coherence-free law |
| Engler (1997) proved the pivot floor and studied QRCP on orthonormal rows; Kahan (1966) gave the triangular matrix; Count ex Machina recorded a small exact example; Amsel et al. version 3 records the updates | All credited where used (Section 1.2, Lemma 3.1, Sections 5 and 6) |
| Lemma 3.3 is a classical determinant estimate | Presented as classical, with its own short proof and no novelty claimed |
| The Drmac-Gugercin bound is Theorem 2.1 of their paper | Corrected |
| Attribution of the existence scale and of Damle's bound | The existence result is cited through Amsel et al.; Damle's Theorem 3.1 and Corollary 3.2 cited, with the bound on gamma attributed to Osinsky as Damle records |
| Significance would rise with a result that needs no coherence assumption | Extension: Corollary 2.6 (no assumption), Theorem 2.7 and Corollary 2.8 (the law in N and k, with explicit lower examples) |
| The scripts never tested C(H) itself | Extension: code/check_extensions.py checks C(H) against direct minimisation in 50-digit arithmetic, the bound on C(H)^2, Table 2, the Kahan family with exact QRCP and inequality (7.5) on a grid; a negative control (a constant changed from 3 to 1) makes it fail |
| Finite sizes were not reported | Tables 1 and 2, Figure 1 and Section 9 report crossovers with the classical bound, realised exponents and slow convergence |
| "Sharp" could be misread as a multiplicative constant | Remark 2.4 says that sharpness concerns the leading term of the exponent and that a polynomial factor remains |
| Readability for non-specialists | New title, context-first abstract and introduction, a purpose sentence at the start of each section, numbered theorem environments, practical guidance (Section 8) |
| Platform dependence of extended precision | Section 9 states that NumPy longdouble is 80-bit on x86-64 and double precision elsewhere, and that the tolerances hold in both cases |
| Re-review and proof audit (minor): a binding condition in the diagonal-sequence paragraph, the explanation of a 2 log k loss, a step N >= 66k, the floor in N = k^alpha, and "at most e^{2 eta} + o(1)" in the comparison with Chen et al. | All applied |

Recorded as open (not attempted or not landed): the worst case at a fixed aspect ratio N <= alpha k, the polynomial factor
between the bounds at fixed H, the sqrt(alpha - 2) to sqrt(alpha - 1) window at N = k^alpha, floating-point QRCP, and
which other pivot rules give trailing-energy constraints.
