# Deposit manifest — *Naphthalene is the only golden benzenoid*

The files below, and only these, constitute the supplement to the paper.
Everything else in the working directory is a working file and is not deposited.

## The proof (required)

| file | role | size |
|---|---|---|
| `seven_hexagon_obstructions.txt` | the atlas: 85 cells of W7, 3,652 rooted cores, 22,009 robust obstruction patterns | 3.4 MB |
| `seven_hexagon_cover.bin` | the covering DAG, 345,554 nodes | 7.2 MB |
| `seven_hexagon_exact_verify.cpp` | the verifier, 350 lines; checks (i)–(v) of §7.1 | 17 KB |
| `small_polyhex_witnesses.txt` | integer witnesses for the 113 polyhexes with h ≤ 6 outside the gap class, and characteristic-polynomial certificates for the three in it | 17 KB |

SHA-256 digests of the first two are printed in §7.1 of the paper.

## Corroboration (deposited)

| file | role | size |
|---|---|---|
| `seven_hexagon_reverify.py` | second verifier, written from the specification rather than from the C++ source; adds 2,000 end-to-end tests against exact Q(√5) spectra | 15 KB |
| `seven_hexagon_optimality.py` | the exact checks behind Proposition 4.2 (seven is best possible) and the h = 7,8,9 statement scan | 8 KB |
| `small_polyhex_verify.py` | third verifier: regenerates the 116 polyhexes and checks every witness in integer arithmetic | 25 KB |

## Transcripts (deposited)

| file | of |
|---|---|
| `seven_hexagon_verification_gmp_rerun.log` | the C++ verifier |
| `seven_hexagon_reverify_shipped_2000.log` | the Python re-verifier, 22,009 patterns, 2,000 end-to-end tests |
| `seven_hexagon_statement_scan_h9.log` | the scan of all 8,353 polyhexes with 7, 8 or 9 cells |
| `small_polyhex_verify.log` | the h ≤ 6 witness verifier |

## Licence

| file | of |
|---|---|
| `README.md` | overview, instructions for running, and the AI declaration |
| `LICENSE` | the MIT License, covering the three Python programs and the C++ verifier |

Copyright (c) 2026 Krystal Guo, University of Amsterdam.

## Totals

About 1,500 lines of code and 10.7 MB, dominated by the two certificate files.

## Deliberately not deposited

- `seven_hexagon_compress.py`, `seven_hexagon_buildcover.py`, `atlas_min.txt`,
  `cover_min.bin`, `cover_min2.bin` — these supported a remark on compressing the
  atlas which the paper no longer makes.
- `seven_hexagon_bare_heptahexes.txt` — the 333 free heptahexes, for inspection only;
  they are not proof objects, as §4.1 explains.
- The search that produced the certificate. It is not part of the proof: the verifier
  reconstructs and re-derives everything and trusts nothing from it (§7.3). Appendix A
  specifies both file formats, the computation of q̂, and a method for searching afresh,
  so the certificate can be rebuilt from the paper alone.

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
to `mpz_class`, by GMP; it has been built and run with both, with identical output.
The Python programs need only the standard library.
