"""Vérifie en haute précision chaque correspondance trouvée par cf.cu.
Usage : python3 verif.py hits.txt        (besoin de : pip install mpmath)
Une formule est retenue si elle colle à >= 50 chiffres avec 500 ET 1000 termes."""
import sys, collections
import mpmath as mp
mp.mp.dps = 120
C = {'zeta3': mp.zeta(3), 'zeta5': mp.zeta(5), 'zeta7': mp.zeta(7), 'catalan': mp.catalan, 'gamma': mp.euler,
     'zeta2': mp.zeta(2), 'pi': mp.pi, 'e': mp.e, 'ln2': mp.log(2), 'ln3': mp.log(3), 'ln5': mp.log(5),
     'pi3': mp.pi ** 3, 'pi_ln2': mp.pi * mp.log(2), 'ln2_2': mp.log(2) ** 2, 'pi_sqrt3': mp.pi / mp.sqrt(3),
     'pi_sqrt2': mp.pi / mp.sqrt(2), 'ln_1sqrt2': mp.log(1 + mp.sqrt(2)),
     'lemniscate': mp.gamma(mp.mpf(1) / 4) ** 2 / (2 * mp.sqrt(2 * mp.pi)), 'gauss': 1 / mp.agm(1, mp.sqrt(2)), 'pi_e': mp.pi * mp.e}
FORMES = ['n^6', 'n^4', 'n^2(2n-1)^2', 'n^4(2n-1)^2', '(2n-1)^4', 'n^3(2n-1)^3', '(2n-1)^6', 'n^4(n+1)^2']

def bval(f, c, n):
    t = 2 * n - 1
    return c * [n ** 6, n ** 4, n ** 2 * t * t, n ** 4 * t * t, t ** 4, n ** 3 * t ** 3, t ** 6, n ** 4 * (n + 1) ** 2][f]

def cf(A3, A2, A1, A0, f, c, N):
    p0, p1, q0, q1 = 1, A0, 0, 1          # entiers exacts
    for n in range(1, N + 1):
        a = ((A3 * n + A2) * n + A1) * n + A0
        b = bval(f, c, n)
        p0, p1, q0, q1 = p1, a * p1 + b * p0, q1, a * q1 + b * q0
    return mp.mpf(p1) / q1

def chiffres(x, y):
    d = abs(x - y)
    return 120 if d == 0 else int(-mp.log10(d / abs(y)))

lignes = [l.split() for l in open(sys.argv[1]) if l.strip()]
stats, retenues, vus = collections.Counter(), [], set()
for L in lignes:
    A3, A2, A1, A0, f, c = map(int, L[:6]); nom = L[6]; p, q, r, s = map(int, L[7:11])
    if nom.startswith('alea'):
        stats['alea (controle negatif) a 10 chiffres'] += 1
        continue
    cible = (p + q * C[nom]) / (r + s * C[nom])
    if chiffres(cf(A3, A2, A1, A0, f, c, 300), cible) < 25:      # tri rapide des faux succès
        stats['rejetees (< 50 chiffres)'] += 1
        continue
    d = min(chiffres(cf(A3, A2, A1, A0, f, c, 500), cible), chiffres(cf(A3, A2, A1, A0, f, c, 1000), cible))
    if d >= 50:
        stats[f'{nom} >= 50 chiffres'] += 1
        cle = (A3, A2, A1, A0, f, c)
        if cle not in vus:
            vus.add(cle); retenues.append((nom, d, A3, A2, A1, A0, f, c, p, q, r, s))
    else:
        stats['rejetees (< 50 chiffres)'] += 1
print(dict(stats))
ctrl = {(34, 51, 27, 5, 0, -1): 'Apery zeta(3)', (0, 11, 11, 3, 1, 1): 'Apery zeta(2)'}
print('controles positifs retrouves :', [v for k, v in ctrl.items() if k in vus])
for nom, d, A3, A2, A1, A0, f, c, p, q, r, s in sorted(retenues):
    tag = ctrl.get((A3, A2, A1, A0, f, c), '')
    print(f'{nom:8s} {d:3d} ch.  a_n={A3}n^3+{A2}n^2+{A1}n+{A0}  b_n={c}*{FORMES[f]}  '
          f'V=({p}+{q}C)/({r}+{s}C) {tag}')
