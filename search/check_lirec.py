"""Vérifie si nos 2 formules candidates sont déjà dans la base publique de la Ramanujan Machine (LIReC).
Deux tests : (1) même forme canonique (même fraction à équivalence près) ;
             (2) une fraction de la base a une limite liée à la nôtre par (p + qV)/(r + sV), entiers <= 1000.
Contrôles positifs : 8/pi^2 (conjecture de la Ramanujan Machine) et Apéry zeta(3) doivent être trouvés.
Installation (une fois) :
  pip install "psycopg[binary]" sqlalchemy gmpy2 sympy mpmath
  pip install --no-deps git+https://github.com/RamanujanMachine/LIReC.git
"""
import os, warnings
from multiprocessing import Pool
warnings.filterwarnings('ignore')
import mpmath as mp
from LIReC.db.access import db
from LIReC.lib.pcf import PCF

mp.mp.dps = 60
CAND = {   # nom : (a_n de la plus haute puissance à la plus basse, b_n idem, valeur exacte)
    'A : 1/(pi-3)':         ([7, 17, 7], [8, -8, 2, 0, 0], 1 / (mp.pi - 3)),
    'B : 2/ln(27/16)':      ([10, 11, 4], [-16, 16, -4, 0, 0], 2 / (3 * mp.log(3) - 4 * mp.log(2))),
    "B' : 3/ln(1024/729)":  ([20, 26, 9], [-36, 36, -9, 0, 0], 3 / (10 * mp.log(2) - 6 * mp.log(3))),
    'CONTROLE 8/pi^2 (RM)': ([3, 3, 1], [-2, 1, 0, 0, 0], 8 / mp.pi ** 2),
    'CONTROLE Apery z3':    ([34, 51, 27, 5], [-1, 0, 0, 0, 0, 0, 0], 6 / mp.zeta(3)),
}

def canon(a, b):
    t, u = PCF(a, b).canonical_form()
    return [[int(c) for c in t.all_coeffs()], [int(c) for c in u.all_coeffs()]]

def lien(args):
    v_str, prec, forme = args
    out = []
    try:
        mp.mp.dps = 40
        v = mp.mpf(v_str)
        if prec is not None and prec < 25:
            return out
        for nom, (a, b, V) in CAND.items():
            mp.mp.dps = 40
            rel = mp.pslq([v, v * V, 1, V], maxcoeff=1000, maxsteps=5000)
            if rel and rel[0] * rel[3] - rel[1] * rel[2] != 0:
                out.append((nom, rel, forme))
    except Exception:
        pass
    return out

if __name__ == '__main__':
    print('connexion a la base LIReC...', flush=True)
    rows = db.cfs
    print(f'fractions continues dans la base : {len(rows)}', flush=True)
    formes = [[[int(c) for c in r.P], [int(c) for c in r.Q]] for r in rows]
    print('\n(1) meme forme canonique :')
    for nom, (a, b, V) in CAND.items():
        print(f'   {nom:24s} -> {"PRESENTE dans la base" if canon(a, b) in formes else "absente"}')
    print('\n(2) limites liees a la notre (cela peut prendre quelques minutes)...', flush=True)
    taches = []
    for r, f in zip(rows, formes):
        try:
            taches.append((str(r.base.value), r.base.precision, f))
        except Exception:
            pass
    with Pool(os.cpu_count()) as P:
        res = [x for lst in P.imap_unordered(lien, taches, chunksize=200) for x in lst]
    par = {}
    for nom, rel, f in res:
        par.setdefault(nom, []).append((rel, f))
    for nom in CAND:
        L = par.get(nom, [])
        print(f'   {nom:24s} -> {len(L)} fraction(s) de la base liee(s)')
        for rel, f in L[:5]:
            print(f'        forme canonique {f}  relation {rel}')
