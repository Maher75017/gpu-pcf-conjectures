"""Vérifie les 3 formules conjecturées à 200 chiffres (convention V = a0 + b1/(a1 + b2/(a2 + ...)), a0 = a(0)).
Usage : python3 verifier_formules.py      (besoin de mpmath)"""
import mpmath as mp
mp.mp.dps = 220
F = {
    "A  : 1/(pi - 3)":             (lambda n: 7*n*n + 17*n + 7,  lambda n: 2*n*n*(2*n - 1)**2,  lambda: 1/(mp.pi - 3)),
    "B  : 2/(3 ln 3 - 4 ln 2)":    (lambda n: 10*n*n + 11*n + 4, lambda n: -4*n*n*(2*n - 1)**2, lambda: 2/(3*mp.log(3) - 4*mp.log(2))),
    "B' : 3/(10 ln 2 - 6 ln 3)":   (lambda n: 20*n*n + 26*n + 9, lambda n: -9*n*n*(2*n - 1)**2, lambda: 3/(10*mp.log(2) - 6*mp.log(3))),
}
def cf(a, b, N):
    p0, p1, q0, q1 = 1, a(0), 0, 1
    for n in range(1, N + 1):
        p0, p1, q0, q1 = p1, a(n)*p1 + b(n)*p0, q1, a(n)*q1 + b(n)*q0
    return mp.mpf(p1) / q1
for nom, (a, b, val) in F.items():
    v1, v2, c = cf(a, b, 1500), cf(a, b, 3000), val()
    d = abs(v2 - c)
    print(f"{nom:28s} chiffres corrects : {200 if d == 0 else int(-mp.log10(d/abs(c)))}   (stabilite 1500->3000 termes : {int(-mp.log10(abs(v2-v1)/abs(c))) if v2 != v1 else 200})")
