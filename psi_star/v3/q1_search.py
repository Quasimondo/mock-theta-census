"""Q1: which finite sum at the dual root of unity gives kappa(rho) = L1(rho)/Cbar ?
Candidates  T = sum_j eps^j xi^((a j^2 + b j)/2) P_j(xi),  xi a root of unity built from rho;
periodic sums are summed as a geometric series in the period (as in the odd-root rule)."""
import mpmath as mp, sys, json, itertools
from fractions import Fraction as Fr
mp.mp.dps = 30
D = json.load(open("q1_limits.json"))
w = mp.expjpi(-mp.mpf(2)/3); C = (1+w)*mp.sqrt(mp.pi*w/(1+2*w)); Cb = mp.conj(C)
rhos = sorted((Fr(k) for k in D), key=lambda r: (r.denominator, r))
kap1 = {Fr(k): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))/Cb for k, v in D.items()}
kap0 = {Fr(k): mp.mpc(mp.mpf(v[2]), mp.mpf(v[3]))/Cb for k, v in D.items()}
def order(fr):          # order of e^(2 pi i fr)
    return Fr(fr).denominator
POCH = {
 "(xi;xi^2)":   lambda xi, j: 1 - xi**(2*j-1),
 "(-xi;xi^2)":  lambda xi, j: 1 + xi**(2*j-1),
 "(xi;xi)":     lambda xi, j: 1 - xi**j,
 "(-xi;xi)":    lambda xi, j: 1 + xi**j,
 "(xi^2;xi^2)": lambda xi, j: 1 - xi**(2*j),
 "(-xi^2;xi^2)":lambda xi, j: 1 + xi**(2*j),
 "1":           lambda xi, j: 1,
}
XI = {"e(-rho/2)": lambda r: -r/2, "e(rho/2)": lambda r: r/2, "e(-rho)": lambda r: -r, "e(rho)": lambda r: r,
      "e(-(rho+1)/2)": lambda r: -(r+1)/2, "e((rho+1)/2)": lambda r: (r+1)/2, "e(-rho/4)": lambda r: -r/4, "e(rho/4)": lambda r: r/4}
def value(fr, a2, b2, eps, fac, inv):
    """exponent (a2 j^2 + b2 j)/2 ; xi = e^(2 pi i fr)"""
    n = order(fr); xi = mp.expjpi(2*mp.mpf(fr.numerator)/fr.denominator)
    xh = mp.expjpi(mp.mpf(fr.numerator)/fr.denominator)       # xi^(1/2)
    per = 2*n
    terms = []; P = mp.mpc(1)
    for j in range(2*per+1):
        if j:
            f = fac(xi, j)
            if abs(f) < mp.mpf(10)**-20:
                if inv: return None
                P = 0
            else:
                P = P/f if inv else P*f
        terms.append(eps**j * xh**(a2*j*j + b2*j) * P)
        if P == 0: break
    if P == 0: return sum(terms)
    lam = terms[per]/terms[0]
    if abs(terms[2*per]/terms[per] - lam) > mp.mpf(10)**-15: return None
    if abs(1-lam) < mp.mpf(10)**-12: return None
    return sum(terms[:per])/(1-lam)
hits = []
for (xn, xf), (pn, pf), a2, b2, eps, inv in itertools.product(XI.items(), POCH.items(), (1, 2, 3, 4), range(-4, 5), (1, -1), (False, True)):
    if (a2 + b2) % 2: continue
    if pn == "1" and inv: continue
    for tn, tgt in (("kappa1", kap1), ("kappa0", kap0)):
        best = None; ok = True; ratio = None
        for r in rhos[:10]:
            v = value(xf(r), a2, b2, eps, pf, inv)
            if v is None or abs(v) < mp.mpf(10)**-15: ok = False; break
            q = tgt[r]/v
            if ratio is None: ratio = q
            elif abs(q - ratio) > mp.mpf(10)**-8: ok = False; break
        if ok:
            # full check
            dev = max(abs(tgt[r]/value(xf(r), a2, b2, eps, pf, inv) - ratio) for r in rhos if value(xf(r), a2, b2, eps, pf, inv) is not None)
            hits.append((tn, xn, pn, inv, a2, b2, eps, ratio, dev))
            print(tn, "= const *", ("1/" if inv else "") + pn, "xi =", xn, " exponent (%d j^2 %+d j)/2" % (a2, b2), "eps", eps, " const =", mp.nstr(ratio, 10), " max deviation", mp.nstr(dev, 3), flush=True)
print(len(hits), "hits")
