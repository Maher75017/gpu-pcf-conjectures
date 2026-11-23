"""Famille autour de B : a_n = A2 n^2 + A1 n + A0, b_n = -c n^2 (2n-1)^2.
On cherche les V tels que exp(m/V) soit une fraction simple R (m = 1, 2, 3, 4, 6), puis on vérifie à 80 chiffres.
Contrôle positif : B (10, 11, 4, c = 4) doit donner m = 2, R = 27/16.
Usage : python3 famille.py <A2max> <A1max> <A0max> <cmax>
"""
import sys, math
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp

A2m, A1m, A0m, cm = map(int, sys.argv[1:5])
g = np.array([(A2, A1, A0, c) for A2 in range(1, A2m + 1) for A1 in range(0, A1m + 1)
              for A0 in range(1, A0m + 1) for c in range(1, cm + 1)], dtype=np.float64)
A2, A1, A0, C = g.T
N = 90
p0, p1, q0, q1 = np.ones(len(g)), A0.copy(), np.zeros(len(g)), np.ones(len(g))
ok = np.ones(len(g), bool)
vmid = None
with np.errstate(all='ignore'):
    for n in range(1, N + 1):
        a = A2 * n * n + A1 * n + A0
        b = -C * n * n * (2 * n - 1) ** 2
        p0, p1, q0, q1 = p1, a * p1 + b * p0, q1, a * q1 + b * q0
        m = np.maximum(np.abs(p1), np.abs(q1))
        ok &= np.isfinite(m) & (m > 0)
        s = np.where(m > 1e100, 1 / m, 1.0)
        p0, p1, q0, q1 = p0 * s, p1 * s, q0 * s, q1 * s
        if n == N - 20:
            vmid = p1 / q1
    V = p1 / q1
ok &= np.isfinite(V) & (np.abs(V - vmid) < 1e-12 * np.abs(V)) & (np.abs(V) > 1e-3)
print(f'{len(g)} fractions, {ok.sum()} convergees', flush=True)

def rationnel_vect(R, prof=14, eps=1e-9):
    """Développement en fraction continue vectorisé : True si R est (numériquement) une fraction simple."""
    x = R.copy(); trouve = np.zeros(len(R), bool); vivant = np.isfinite(x) & (x > 0) & (x < 1e6) & (np.abs(x - 1) > 1e-6)
    for _ in range(prof):
        a = np.floor(x); fr = x - a
        fini = vivant & (fr < eps * np.maximum(1.0, x))
        trouve |= fini
        vivant &= ~fini & (fr > 0)
        with np.errstate(all='ignore'):
            x = np.where(vivant, 1.0 / np.where(fr > 0, fr, 1.0), 1.0)
    return trouve

hits = []
idx = np.nonzero(ok)[0]
for mm in (1, 2, 3, 4, 6):
    with np.errstate(all='ignore'):
        R = np.exp(np.clip(mm / V[idx], -700, 700))
    for j in idx[rationnel_vect(R)]:
        f = Fr(float(math.exp(mm / V[j]))).limit_denominator(10 ** 9)
        hits.append((int(A2[j]), int(A1[j]), int(A0[j]), int(C[j]), mm, f))
print(f'candidats (10 chiffres) : {len(hits)}', flush=True)

def cf_exact(A2, A1, A0, c, N):
    p0, p1, q0, q1 = 1, A0, 0, 1
    for n in range(1, N + 1):
        a = A2 * n * n + A1 * n + A0
        b = -c * n * n * (2 * n - 1) ** 2
        p0, p1, q0, q1 = p1, a * p1 + b * p0, q1, a * q1 + b * q0
    return p1, q1

mp.mp.dps = 100
bons = []
for A2_, A1_, A0_, c_, mm, f in hits:
    if f <= 0 or f == 1 or f.denominator > 10 ** 8:
        continue
    p, q = cf_exact(A2_, A1_, A0_, c_, 400)
    V2 = mp.mpf(p) / q
    lg = mp.log(mp.mpf(f.numerator) / f.denominator)
    if lg == 0:
        continue
    cible = mm / lg
    if cible == 0:
        continue
    d = abs(V2 - cible)
    ch = 100 if d == 0 else int(-mp.log10(d / abs(cible)))
    if ch >= 80:
        bons.append((A2_, A1_, A0_, c_, mm, f, ch))
print(f'verifiees a >= 80 chiffres : {len(bons)}')
for A2_, A1_, A0_, c_, mm, f, ch in sorted(bons, key=lambda x: (x[3], x[0], x[1], x[2])):
    tag = '  <- B (controle)' if (A2_, A1_, A0_, c_) == (10, 11, 4, 4) else ''
    print(f'  a_n={A2_}n^2+{A1_}n+{A0_}  b_n=-{c_} n^2(2n-1)^2   V = {mm}/ln({f}){tag}')
