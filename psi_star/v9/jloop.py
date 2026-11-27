"""Closing the loop: JJ(r) = sum_m r^(2m^2+m)/(-r^2;r^2)_m near r = 1, inside the disk.
r = e^(2 pi i/Y),  Y = Y0 + 960 k - i b,  tau = -2 pi i/Y,  dual nome q' = e^(-pi i Y/2).
   L(Y0) = lim_k (JJ(r) - Phi(tau)) sqrt(tau),      Phi = perturbative series of JJ at r = e^(-tau) (48 terms).
Step 1: for fixed Y0', the values L(Y0' + 4n), n = 0..5, are fitted by two exponentials lambda_j^n (Prony): the phases.
Step 2: the two amplitudes, as functions of q', are fitted on F0*, q F1*, G0, G1 at q'.
Usage: python3 jloop.py          (limits are cached in jloop.json)"""
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, arb, acb_series, ctx
from mpmath import mp, mpf, mpc, nstr, exp, pi, matrix, qr_solve, lu_solve, polyroots, arg, norm, sqrt
import horo
from e1_limits import extrap
HERE = os.path.dirname(os.path.abspath(__file__))
KS = [5, 6, 7, 8, 9, 10, 11, 12]
NT = 48
B = Fr(4, 5)
def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(50, radius=False)), mpf(z.imag.str(50, radius=False)))

def pert_series():
    path = os.path.join(HERE, "jloop_pert.json")
    if os.path.exists(path): return json.load(open(path))
    ctx.prec = 700; ctx.cap = NT + 1
    tot = acb_series([0]); P = acb_series([1])
    for m in range(1200):
        if m: P = P/(1 + acb_series([0, -2*m]).exp())
        tot = tot + P*acb_series([0, -(2*m*m + m)]).exp()
    co = [tot[k].real.str(150, radius=False) for k in range(NT + 1)]
    json.dump(co, open(path, "w")); return co

def JJ(Yp, prec, target=110):
    ctx.prec = prec
    Y = horo.cx(*Yp)
    r = (2*acb.pi()*acb(0, 1)/Y).exp().mid()
    K = float(abs(Y).mid()); b = -float(Y.imag.mid())
    M = int(K*math.sqrt((0.0808*K + 12 + target*0.6932)/(4*math.pi*b))) + 10
    one = acb(1); s = acb(0); den = one; e = one; r2m = one; r2 = (r*r).mid()
    for m in range(M):
        s += e/den
        e = (e*r2m*r2m*r2*r).mid()          # r^(2m^2+m) -> r^(2(m+1)^2+(m+1)): factor r^(4m+3)
        r2m = (r2m*r2).mid()
        den = (den*(one + r2m)).mid()
    return s

def limit(Y0, co):
    rows = []
    for k in KS:
        Yp = (Y0 + 960*k, B); K = float(Y0 + 960*k)
        v = JJ(Yp, int(0.0808*K/0.6932*1.15) + 220)
        tau = -2*acb.pi()*acb(0, 1)/horo.cx(*Yp)
        ph = acb(0)
        for c in reversed(co): ph = ph*tau + arb(c)
        rows.append((tomp(tau), [tomp((v - ph)*tau.sqrt())]))
    return extrap(rows, ["L"])["L"]

def duals(qd):
    F0 = F1 = G0 = G1 = mpc(0); po = mpc(1); n = 0
    while abs(qd)**(n*n) > mpf(10)**(-40):
        pn = po*(1 - qd**(2*n + 1))
        F0 += (-1)**n*qd**(2*n*n)/po; F1 += (-1)**n*qd**(2*n*n + 2*n)/pn
        if n >= 1: G0 += qd**(n*n)*po
        G1 += qd**((n + 1)**2)*po
        po = pn; n += 1
    return {"F0*": F0, "q F1*": qd*F1, "G0": G0, "G1": G1}

if __name__ == "__main__":
    co = pert_series()
    path = os.path.join(HERE, "jloop.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    bases = [Fr(0), Fr(1, 2), Fr(11, 10), Fr(17, 10), Fr(12, 5), Fr(33, 10)]
    for Y0p in bases:
        for n in range(6):
            key = str(Y0p + 4*n)
            if key not in res:
                res[key] = limit(Y0p + 4*n, co); json.dump(res, open(path, "w"), indent=0)
    mp.dps = 30
    print("perturbative series of JJ, first coefficients:", [c[:14] for c in co[:4]])
    print("largest extrapolation uncertainty:", max(v[2] for v in res.values()))
    val = lambda Y0: mpc(mpf(res[str(Y0)][0]), mpf(res[str(Y0)][1]))
    thetas = []
    for Y0p in bases:
        s = [val(Y0p + 4*n) for n in range(6)]
        A = matrix(4, 2); y = matrix(4, 1)
        for n in range(4): A[n, 0] = s[n + 1]; A[n, 1] = s[n]; y[n] = s[n + 2]
        ab, rs = qr_solve(A, y)
        lam = polyroots([1, -ab[0], -ab[1]])
        lam.sort(key=lambda z: arg(z))
        thetas.append([arg(z)/4 for z in lam])
        print(f"Y0' = {str(Y0p):6s} two-exponential fit residual {nstr(rs, 3)};  |lambda| = {nstr(abs(lam[0]), 12)}, {nstr(abs(lam[1]), 12)};  -2 theta/pi mod 1 = {nstr((-2*thetas[-1][0]/pi) % 1, 14)}, {nstr((-2*thetas[-1][1]/pi) % 1, 14)}")
    expo = [sum((-2*t[j]/pi) % 1 for t in thetas)/len(thetas) for j in range(2)]
    rat = [Fr(int(round(float(e*240))), 240) for e in expo]
    print("exponents of q' (mod 1):", [nstr(e, 14) for e in expo], " nearest k/240:", [str(x) for x in rat])
    # step 2: amplitudes with exact exponents a_j: L(Y) = sum_j q'^(a_j) A_j(q'),  q'^a := exp(-pi i Y a/2)
    amps = []
    for Y0p in bases:
        A = matrix(6, 2); y = matrix(6, 1)
        for n in range(6):
            Y = mpc(mpf((Y0p + 4*n).numerator)/(Y0p + 4*n).denominator, -mpf(B.numerator)/B.denominator)
            for j in range(2): A[n, j] = exp(-1j*pi*Y*mpf(rat[j].numerator)/rat[j].denominator/2)
            y[n] = val(Y0p + 4*n)
        x, rs = qr_solve(A, y)
        Y = mpc(mpf(Y0p.numerator)/Y0p.denominator, -mpf(B.numerator)/B.denominator)
        amps.append((exp(-1j*pi*Y/2), x, rs))
    print("residual of the fit with the rational exponents:", nstr(max(a[2] for a in amps), 3))
    for j in range(2):
        print(f"sector q'^({rat[j]}):")
        for shift in (-1, 0, 1):
            names = ["F0*", "q F1*", "G0", "G1"]
            A = matrix(6, 4); y = matrix(6, 1)
            for i, (qd, x, _) in enumerate(amps):
                d = duals(qd)
                for c, nm in enumerate(names): A[i, c] = d[nm]*qd**shift
                y[i] = x[j]
            sol, rs = qr_solve(A, y)
            print(f"   A_j = q'^({shift}) * [" + ",  ".join(f"({nstr(c, 12)}) {nm}" for c, nm in zip(sol, names)) + f"]   relative residual {nstr(rs/norm(y), 3)}")
