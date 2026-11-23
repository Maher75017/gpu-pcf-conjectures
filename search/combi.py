"""Fractions continues de degré 2 dont la limite mélange DEUX constantes.
a_n = A2 n^2 + A1 n + A0 ; b_n = c * forme(n), formes : n^4, n^2(2n-1)^2, (2n-1)^4, n^3(n+1)... (voir FORMES).
Filtres : pas de facteur commun a/b, b_n de signe constant, >= 30 chiffres de convergence, limite ni rationnelle ni quadratique.
Identification : PSLQ sur [V, V c1, V c2, 1, c1, c2] pour chaque paire (c1, c2), puis vérification à >= 80 chiffres.
Contrôles : Apéry zeta(2) et Catalan doivent sortir ; paires de constantes aléatoires : 0 attendu.
Usage : python3 combi.py <A2max> <A1max> <A0max> <cmax>
"""
import sys, math, os, itertools, time
from fractions import Fraction as Fr
from multiprocessing import Pool
import mpmath as mp

FORMES = ['n^4', 'n^2(2n-1)^2', '(2n-1)^4', 'n^2(n+1)^2', 'n^3(2n-1)']
def bf(f, n):
    t = 2 * n - 1
    return [n ** 4, n * n * t * t, t ** 4, n * n * (n + 1) ** 2, n ** 3 * t][f]

def base():
    mp.mp.dps = 120
    K = {'pi': mp.pi, 'pi2': mp.pi ** 2, 'ln2': mp.log(2), 'zeta3': mp.zeta(3), 'G': mp.catalan, 'ln3': mp.log(3),
         'sqrt2': mp.sqrt(2), 'sqrt3': mp.sqrt(3)}
    R = {f'alea{i}': mp.mpf(i + 1) / 7 + mp.sqrt(mp.mpf(10 + i)) / 3 + mp.log(mp.mpf(5 + i)) / 11 for i in range(4)}
    return K, R

def cf(A, f, c, N):
    p0, p1, q0, q1 = 1, A[0], 0, 1
    for n in range(1, N + 1):
        a = A[2] * n * n + A[1] * n + A[0]
        b = c * bf(f, n)
        p0, p1, q0, q1 = p1, a * p1 + b * p0, q1, a * q1 + b * q0
    return p1, q1

def candidat(args):
    A2, A1, A0, f, c = args
    A = (A0, A1, A2)
    # facteur commun a(r) = b(r) = 0 pour r rationnel simple -> formule déguisée
    for r in [Fr(k) for k in range(-6, 7)] + [Fr(2 * k + 1, 2) for k in range(-6, 6)]:
        if A2 * r * r + A1 * r + A0 == 0 and c * bf(f, r) == 0:
            return None
    if any(A2 * n * n + A1 * n + A0 <= 0 for n in range(1, 50)):
        return None
    p1, q1 = cf(A, f, c, 120)
    p2, q2 = cf(A, f, c, 240)
    if q1 == 0 or q2 == 0:
        return None
    num = abs(p1 * q2 - p2 * q1)
    if num == 0:
        return None
    if math.log10(abs(q1)) + math.log10(abs(q2)) - math.log10(num) < 30:   # < 30 chiffres de convergence
        return None
    mp.mp.dps = 50
    V = mp.mpf(p2) / q2
    if not (mp.mpf('1e-6') < abs(V) < 10 ** 6):
        return None
    if mp.pslq([V, 1], maxcoeff=10 ** 6, maxsteps=2000) or mp.pslq([V * V, V, 1], maxcoeff=10 ** 5, maxsteps=4000):
        return None
    return (A2, A1, A0, f, c, str(V))

def identifier(item):
    A2, A1, A0, f, c, Vs = item
    K, R = base()
    out = []
    for groupe, tag in ((K, 'vrai'), (R, 'alea')):
        noms = list(groupe)
        for n1, n2 in itertools.combinations(noms, 2):
            mp.mp.dps = 45
            V = mp.mpf(Vs); c1, c2 = groupe[n1], groupe[n2]
            rel = mp.pslq([V, V * c1, V * c2, 1, c1, c2], maxcoeff=500, maxsteps=20000)
            if not rel:
                continue
            x = rel
            den = x[0] + x[1] * c1 + x[2] * c2
            if den == 0 or (x[3] == 0 and x[4] == 0 and x[5] == 0):
                continue
            # vérification profonde
            p, q = cf((A0, A1, A2), f, c, 1500)
            mp.mp.dps = 120
            V2 = mp.mpf(p) / q
            c1b, c2b = mp.mpf(groupe[n1]), mp.mpf(groupe[n2])
            cible = -(x[3] + x[4] * c1b + x[5] * c2b) / (x[0] + x[1] * c1b + x[2] * c2b)
            d = abs(V2 - cible)
            ch = 120 if d == 0 else int(-mp.log10(d / abs(cible)))
            if ch >= 80:
                out.append((tag, n1, n2, list(x), ch))
                break
    return (A2, A1, A0, f, c, out)

if __name__ == '__main__':
    A2m, A1m, A0m, cm = map(int, sys.argv[1:5])
    t0 = time.time()
    taches = [(A2, A1, A0, f, c) for A2 in range(1, A2m + 1) for A1 in range(0, A1m + 1) for A0 in range(1, A0m + 1)
              for f in range(len(FORMES)) for c in list(range(-cm, 0)) + list(range(1, cm + 1))]
    ctrl = [(11, 11, 3, 0, 1), (3, 3, 1, 0, -2)]      # Apéry zeta(2), Catalan
    for t in ctrl:
        if t not in taches: taches.append(t)
    print(f'{len(taches)} fractions a tester', flush=True)
    with Pool(os.cpu_count()) as P:
        surv = [x for x in P.imap_unordered(candidat, taches, chunksize=200) if x]
        print(f'survivantes (convergentes, non triviales) : {len(surv)}  ({time.time()-t0:.0f} s)', flush=True)
        res = [r for r in P.imap_unordered(identifier, surv, chunksize=4) if r[5]]
    nv = na = 0
    for A2, A1, A0, f, c, out in sorted(res):
        for tag, n1, n2, x, ch in out:
            if tag == 'alea':
                na += 1; continue
            nv += 1
            ctl = ' <- CONTROLE' if (A2, A1, A0, f, c) in ctrl else ''
            print(f'{ch:3d} ch.  a_n={A2}n^2+{A1}n+{A0}  b_n={c}*{FORMES[f]}  constantes ({n1},{n2})  '
                  f'relation [V,Vc1,Vc2,1,c1,c2]={x}{ctl}')
    print(f'identifiees : {nv} ; controle negatif (aleatoires) : {na} ; temps total {time.time()-t0:.0f} s')
