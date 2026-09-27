# Naphthalene is the only golden benzenoid

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22996285.svg)](https://doi.org/10.5281/zenodo.22996285)

Krystal Guo, Korteweg-de Vries Institute for Mathematics, University of Amsterdam
(k.guo@uva.nl)

This repository contains the paper and its computational supplement.

## The result

Benzene, naphthalene and triphenylene are the only polyhexes whose graphs have no
eigenvalue in the open interval (−θ, θ), where θ = (√5 − 1)/2 is the reciprocal of
the golden ratio (Theorem 1.1). In particular, naphthalene is the only polyhex of
even order whose median eigenvalues are ±θ (Corollary 1.2); this answers a question
of Pisanski.

The proof rests on a locality theorem (Theorem 1.3): if A is the adjacency matrix of
a polyhex with at least seven hexagons, then there is a connected set C of seven
hexagons, one of which is the first hexagon of the polyhex in a fixed lexicographic
order, such that the principal submatrix of A² − θ²I on the vertices of C is not
positive semidefinite. Seven is best possible. The locality theorem reduces to a
finite computation, and the certificate for it, deposited here, reduces every
spectral claim to an inequality between integers.

## Files

| file | role |
|---|---|
| `golden-benzenoid.tex`, `golden-benzenoid.pdf` | the paper |
| `Lmacs.tex` | LaTeX macros used by the paper |
| `seven_hexagon_obstructions.txt` | the atlas: 85 cells of W7, 3,652 rooted cores, 22,009 robust obstruction patterns |
| `seven_hexagon_cover.bin` | the covering DAG, 345,554 nodes |
| `seven_hexagon_exact_verify.cpp` | the verifier; checks (i)–(v) of §7 |
| `small_polyhex_witnesses.txt` | integer witnesses for the 113 polyhexes with h ≤ 6 outside the gap class, and characteristic-polynomial certificates for the three in it |
| `seven_hexagon_reverify.py` | second verifier, written from the specification in Appendix A rather than from the C++ source; adds end-to-end tests against exact Q(√5) spectra |
| `seven_hexagon_optimality.py` | the exact checks behind Proposition 4.2 (seven is best possible) and the h = 7, 8, 9 statement scan |
| `small_polyhex_verify.py` | third verifier: regenerates the 116 polyhexes with h ≤ 6 and checks every witness in integer arithmetic |
| `*.log` | transcripts of the runs of the verifiers |
| `DEPOSIT-MANIFEST.md` | the full deposit manifest |
| `LICENSE` | MIT License |

SHA-256 digests of the two certificate files are printed in §7 of the paper, and
Appendix A specifies both file formats.

## Building and running

```
g++ -O2 -std=c++20 seven_hexagon_exact_verify.cpp -o seven_hexagon_exact_verify
./seven_hexagon_exact_verify seven_hexagon_obstructions.txt seven_hexagon_cover.bin

python3 seven_hexagon_reverify.py seven_hexagon_obstructions.txt seven_hexagon_cover.bin 2000
python3 seven_hexagon_optimality.py optimality
python3 seven_hexagon_optimality.py scan 9        # slow
python3 small_polyhex_verify.py
```

The C++ verifier needs arbitrary-precision integers, supplied by Boost's header-only
multiprecision library or, with a two-line shim aliasing `boost::multiprecision::cpp_int`
to `mpz_class`, by GMP. The Python programs need only the standard library.

## Declaration on the use of AI tools

*Reproduced from the paper.*

The strategy of the proof is the author's: to prove the stronger statement,
Theorem 1.1, rather than Corollary 1.2, by exhibiting a local block of
A² − θ²I that is not positive semidefinite, following the approach of Guo and
Royle [1]. The computations showing that six hexagons do not suffice
(Section 4.2) were carried out by the author with the assistance of Claude
(Anthropic).

The argument for seven hexagons in Section 5 is not the author's. The author's
plan for this step rested on a much larger case analysis. A draft of the paper
containing that plan and a proof sketch of a much more computationally-dependent
proof with more case analysis was given to ChatGPT 5.6 Pro (OpenAI), which
completed it with a different argument, including Lemmas 5.1, 5.2, 5.3 and 5.6,
none of which the author had drafted. This argument is shorter than the author's
original one and a forthcoming edition will reconcile the two approaches to give a
more intuitive proof. Claude was also used in editing the text of the manuscript.

[1] K. Guo, G. F. Royle, Cubic graphs with no eigenvalues in the interval (−1, 1),
J. Combin. Theory Ser. B **176** (2026) 561–583.

## Archive

This repository is archived at Zenodo. The DOI
[10.5281/zenodo.22996285](https://doi.org/10.5281/zenodo.22996285) always resolves to
the latest version; version 1.0, released on 27 September 2026, is
[10.5281/zenodo.22996286](https://doi.org/10.5281/zenodo.22996286).

## Licence

The code (`seven_hexagon_exact_verify.cpp` and the three Python programs) is released
under the MIT License; see `LICENSE`.

Copyright (c) 2026 Krystal Guo, University of Amsterdam.
