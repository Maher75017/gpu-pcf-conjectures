"""B et B' sont les cas k = 2, 3 d'une famille infinie, démontrée par la fraction continue d'Euler.

Famille (convention V = a(0) + b(1)/(a(1) + b(2)/(a(2) + ...))) :
    a(n) = 2(k^2+1) n^2 + (3k^2-1) n + k^2,   b(n) = -k^2 n^2 (2n-1)^2
    V_k  = k / L_k,   L_k = (k-1) ln(k-1) + (k+1) ln(k+1) - 2k ln k = somme_{m>=1} 1/(m (2m-1) k^(2m-1))
Preuve : c'est la fraction continue d'Euler de cette série. Les « queues » c_n = k^2 (n+1)(2n+1) vérifient
c_0 = a(0) et c_n = a(n) + b(n)/c_{n-1}, avec a(n) = c_n (1 + r_n) et b(n) = -r_n c_{n-1} c_n,
r_n = n(2n-1) / (k^2 (n+1)(2n+1)) = rapport de deux termes consécutifs de la série.
Ce script vérifie : (1) la valeur numérique pour des k entiers et réels, (2) la structure d'Euler en arithmétique exacte,
(3) la série, (4) que A (1/(pi-3)) n'a pas cette structure (numérateurs des réduites à grands facteurs premiers).
Usage : python3 verifier_famille.py      (besoin de mpmath ; sympy facultatif pour le point 4)
"""
from fractions import Fraction as Fr
import mpmath as mp

mp.mp.dps = 130

def cf(k, N):
    A2, A1, A0 = 2 * (k * k + 1), 3 * k * k - 1, k * k
    p0, p1, q0, q1 = mp.mpf(1), mp.mpf(A0), mp.mpf(0), mp.mpf(1)
    for n in range(1, N + 1):
        a = A2 * n * n + A1 * n + A0
        b = -k * k * n * n * (2 * n - 1) ** 2
        p0, p1, q0, q1 = p1, a * p1 + b * p0, q1, a * q1 + b * q0
    return p1 / q1

def chiffres(x, y):
    d = abs(x - y) / abs(y)
    return 120 if d == 0 else int(-mp.log10(d))

print("(1) valeur de la fraction continue (800 termes) contre k / L_k")
for k in [2, 3, 4, 5, 6, 8, 13, 100, mp.mpf(5) / 2, mp.sqrt(2), mp.pi]:
    k = mp.mpf(k)
    v4, v8 = cf(k, 400), cf(k, 800)
    L = (k - 1) * mp.log(k - 1) + (k + 1) * mp.log(k + 1) - 2 * k * mp.log(k)
    print(f"  k = {mp.nstr(k, 8):>10} : convergence {chiffres(v4, v8):3d} chiffres, accord avec k/L_k {chiffres(v8, k / L):3d} chiffres")

print("(2) structure d'Euler en arithmétique exacte : c_n = k^2 (n+1)(2n+1), c_0 = a(0), c_n = a(n) + b(n)/c_(n-1)")
ok = True
for k in range(2, 21):
    c_prev = Fr(k * k)
    for n in range(1, 201):
        a = 2 * (k * k + 1) * n * n + (3 * k * k - 1) * n + k * k
        b = -k * k * n * n * (2 * n - 1) ** 2
        c = a + Fr(b) / c_prev
        ok &= c == k * k * (n + 1) * (2 * n + 1)
        c_prev = c
print(f"  k = 2..20, n = 1..200 : {'OK' if ok else 'ECHEC'}")

print("(3) série somme 1/(m (2m-1) k^(2m-1)) contre L_k")
for k in [2, 3, 4, 7]:
    k = mp.mpf(k)
    S = mp.nsum(lambda m: 1 / (m * (2 * m - 1) * k ** (2 * m - 1)), [1, mp.inf])
    L = (k - 1) * mp.log(k - 1) + (k + 1) * mp.log(k + 1) - 2 * k * mp.log(k)
    print(f"  k = {int(k)} : accord {chiffres(S, L)} chiffres")

print("(4) A n'a pas cette structure : plus grand facteur premier des numérateurs p(n) des réduites")
try:
    import sympy as sp
    def numer(a, b, a0, N):
        p0, p1 = 1, a0
        out = []
        for n in range(1, N + 1):
            p0, p1 = p1, a(n) * p1 + b(n) * p0
            out.append(p1)
        return out
    for nom, a, b, a0 in [("A  1/(pi-3)", lambda n: 7 * n * n + 17 * n + 7, lambda n: 2 * n * n * (2 * n - 1) ** 2, 7),
                          ("B  k = 2   ", lambda n: 10 * n * n + 11 * n + 4, lambda n: -4 * n * n * (2 * n - 1) ** 2, 4)]:
        print(f"  {nom} : " + ", ".join(str(max(sp.factorint(p))) for p in numer(a, b, a0, 8)))
except ImportError:
    print("  (sympy absent : point 4 sauté)")
