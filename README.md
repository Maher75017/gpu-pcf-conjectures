# gpu-pcf-conjectures

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24485588-blue.svg)](https://doi.org/10.5281/zenodo.24485588)

DOI for all versions (always the latest): [10.5281/zenodo.24485588](https://doi.org/10.5281/zenodo.24485588). Version 1.1.0: [10.5281/zenodo.24633640](https://doi.org/10.5281/zenodo.24633640). Version 1.0.0: [10.5281/zenodo.24485589](https://doi.org/10.5281/zenodo.24485589).

Three polynomial continued fractions found by a GPU search, each matching a closed form to about 220 digits, and not found in the references and database we checked.

**Update (28 November 2026).** B and B' are now **proved**: they are the cases k = 2 and k = 3 of an infinite family that comes from Euler's continued fraction of a simple series (see [below](#update-b-and-b-are-proved-an-infinite-family)). **A, for 1/(π − 3), remains a numerical conjecture**: it does not have this structure. Question on MathOverflow: [Are these continued fractions known?](https://mathoverflow.net/questions/518279)

Author: Maher El Hmissi ([Maher75017](https://github.com/Maher75017)), November 2026. Search code and analysis written with the help of Claude (Anthropic). Français : [README.fr.md](README.fr.md).

Convention (same as the Ramanujan Machine and LIReC): V = a(0) + b(1)/(a(1) + b(2)/(a(2) + b(3)/(a(3) + ...))).

## The three formulas

| | a(n) | b(n) | Value | Status |
| --- | --- | --- | --- | --- |
| **A** | 7n² + 17n + 7 | 2n²(2n − 1)² | 1 / (π − 3) | conjecture |
| **B** | 10n² + 11n + 4 | −4n²(2n − 1)² | 2 / (3 ln 3 − 4 ln 2) = 2 / ln(27/16) | proved (family, k = 2) |
| **B'** | 20n² + 26n + 9 | −9n²(2n − 1)² | 3 / (10 ln 2 − 6 ln 3) = 3 / ln(1024/729) | proved (family, k = 3) |

`python3 verifier_formules.py` (needs `mpmath`) evaluates each fraction with exact integers to 1,500 and 3,000 terms: all three agree with the closed form to 219–220 digits, and the two depths agree with each other.

**Correction.** We first wrote that none of them is a classical series in disguise, because the denominators q(n) of the convergents have large prime factors. That test was incomplete: for B and B' it is the numerators p(n) that have only small prime factors, the signature of Euler's continued fraction of a series, and this is how they are proved below. For A, both p(n) and q(n) have large prime factors (p(3) = 3² · 5 · 7 · 6,311; q(3) = 3 · 93,827).

## Update: B and B' are proved (an infinite family)

For every real k > 1,

- a(n) = 2(k² + 1)n² + (3k² − 1)n + k², b(n) = −k²n²(2n − 1)², a(0) = k²,
- V = k / L_k, with L_k = (k − 1) ln(k − 1) + (k + 1) ln(k + 1) − 2k ln k = Σ_{m≥1} 1 / (m(2m − 1) k^(2m−1)).

k = 2 gives B (L₂ = ln(27/16)), k = 3 gives B' (L₃ = ln(1024/729)), k = 4 gives 4 / ln(84375/65536) with a(n) = 34n² + 47n + 16, and so on.

**Proof.** The terms t_m = 1 / (m(2m − 1)k^(2m−1)) of the series have ratio r_m = t_(m+1)/t_m = m(2m − 1) / (k²(m + 1)(2m + 1)). Euler's continued fraction gives S = t₁ / (1 − r₁/(1 + r₁ − r₂/(1 + r₂ − ⋯))) with t₁ = 1/k, so V = k/S = k²(1 − r₁/(1 + r₁ − ⋯)). The equivalence transformation with c_n = k²(n + 1)(2n + 1) turns this into a(n) = c_n(1 + r_n) and b(n) = −r_n c_(n−1) c_n, which are the polynomials above, with a(0) = c₀ = k². The series converges for k > 1, and the continued fraction gains about log₁₀(k²) digits per term (for k = 1 it still tends to 1/ln 4, slowly).

`python3 verifier_famille.py` (< 1 s) checks the value to about 128 digits for k = 2 … 100 and for k = 5/2, √2 and π, the Euler structure in exact arithmetic for k = 2 … 20 and n ≤ 200, and the series.

**What this means.** B and B' are elementary consequences of a classical method (Euler, 1748): they are not new in a deep sense, and we had not spotted the family behind them before this update. A does not fit: it has no polynomial tails of degree 2, the numerators p(n) of its convergents have large prime factors (largest: 73, 347, 6,311, 27,953, … against 3, 5, 7, 7, … for B), and no combination u·p(n) + v·q(n) with |u|, |v| ≤ 12 has only small prime factors. A is the open question. The family and its identification with Euler's continued fraction were found with the help of Claude.

## Novelty checks (and their limits)

- **Henri Cohen, "Parametric continued fractions for π², ζ(3), and other constants" (arXiv 2304.11727)**, the reference list of low-degree families: A, B and B' are not in it. The paper gives no explicit continued fraction for π (only a mention of parametric families) and none for log(a/b).
- **The Ramanujan Machine paper (Raayoni et al., Nature 2021; arXiv 1907.00205)**: not among its conjectures.
- **LIReC, the Ramanujan Machine public database** (827 continued fractions, queried on 24 November 2026 with `search/check_lirec.py`):
  - same canonical form: A, B, B' **absent**; positive controls (their 8/π² conjecture and Apéry's ζ(3) fraction) **present**;
  - any stored fraction whose limit is (p + qV)/(r + sV), integers up to 1,000: **none** for B and B'. For A, 24 stored fractions qualify, but only because they are other formulas for π and 1/(π − 3) is a Möbius transform of π.
- LIReC canonical forms (top, bottom): A `[[8,-8,2,0,0],[49,140,79,-30,-21]]`; B `[[-16,16,-4,0,0],[100,20,-29,-3,12]]`; B' `[[-36,36,-9,0,0],[400,240,-124,-48,27]]`.

**Limits.** The public LIReC database is a selection, and the literature on continued fractions is large. Absence from these sources is evidence, not proof, that the formulas are new. Corrections and references are welcome. Update: B and B' were later proved through Euler's continued fraction (see above); no proof is known to us for A.

## How they were found

| Run | Space | Hardware, time | Outcome |
| --- | --- | --- | --- |
| 1 `cf.cu` | 6,933,454,848 fractions, a(n) of degree ≤ 3, 8 shapes of b(n) | RTX 5090, 36.8 s | 11 exact formulas, all known (Apéry, Ramanujan Machine, Zagier) or rescaled copies |
| 2 `cf2.cu` | 2,147,483,648 "balanced" degree-5 fractions aimed at ζ(5) | RTX 5090, 9.3 s | no formula for ζ(5), Catalan, γ or ζ(2) in this family; ζ(3) formulas known |
| 3 `cf3.cu` | 5,618,466,816 fractions with an alternating term E(−1)ⁿn² | RTX 5090, 34.3 s | only the known E = 0 formulas |
| 4 `cf4.cu` | run 1's space against 20 constants | RTX 5090, 36.8 s | 48 exact formulas: mostly classical (a common factor (n+1)² hiding Euler/Gauss fractions, Apéry, Ramanujan Machine, Zagier, rescaled copies), plus **A**, first seen here |
| 5 `combi.py` | 168,001 degree-2 fractions, limits tested against pairs of constants by PSLQ | 16 CPU threads, 47 min | known formulas (8/π² of the Ramanujan Machine, 18/π² of Cohen §7, 24/π²), **A** again, and **B** |
| 6 `famille.py` | 65,945,600 fractions of B's shape, limits tested for exp(m/V) rational | CPU, minutes | **B'**. We first wrote "no other member of this shape": **wrong**. k = 4, 5, … of the family are in this space, but the detection step (a double-precision test that exp(m/V) is a simple fraction, m ∈ {1, 2, 3, 4, 6}) missed them: precision for k = 4, m = 5 not tried for k = 5, fractions too large for k ≥ 6 |

Every run had a positive control (known formulas planted in the search space, e.g. Apéry's 6/ζ(3) = 5 − 1⁶/(117 − 2⁶/(535 − …))) and a negative control (random constants in the same identification table): positive controls always found, negative controls never confirmed. A 10-digit match only counts after an exact-integer recomputation agrees to at least 50 digits (80 for run 6, 200+ for the final check).

Negative results are part of the record: no "Apéry-like" balanced degree-5 formula for ζ(5) with coefficients up to (255, 511, 64) and |c| ≤ 32 (run 2).

## Reproduce

```sh
python3 verifier_formules.py                                  # the three formulas, ~1 s
python3 verifier_famille.py                                   # the family behind B and B', < 1 s
cd search
nvcc -O3 -arch=native -o cf cf.cu && ./cf && python3 verif.py hits.txt   # run 1 (GPU)
python3 combi.py 10 20 10 8                                   # run 5 (multi-core CPU)
python3 famille.py 80 160 80 64                               # run 6 (~8 GB RAM)
pip install "psycopg[binary]" psycopg2-binary sqlalchemy gmpy2 sympy mpmath
pip install --no-deps git+https://github.com/RamanujanMachine/LIReC.git
python3 check_lirec.py                                        # database check
```
Run 3 needs `-DA3N=41 -DA2N=64 -DA1N=41 -DA0N=12` at compile time. Script messages are in French.

## See also

[exact-exp32](https://github.com/Maher75017/exact-exp32) (correctly rounded exponentials checked on all 2^32 inputs) and [tt-erf-fp32](https://github.com/Maher75017/tt-erf-fp32) (FP32 erf for the Tenstorrent SFPU): the same approach of exhaustive search with explicit controls.

## License

MIT (see `LICENSE`).
