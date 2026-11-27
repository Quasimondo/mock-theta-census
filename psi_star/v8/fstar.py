"""The law for F0*, F1* (sign-twisted fifth-order F0, F1) at q -> 1, inside the disk.
   F0* = sum (-1)^n q^(2n^2)/(q;q^2)_n          G0 = sum_{n>=1} q^(n^2) (q;q^2)_n          F0* + G0 = theta4(q^2)(-q;q)_inf X(q^2)   (III)
   F1* = sum (-1)^n q^(2n^2+2n)/(q;q^2)_(n+1)   G1 = sum_{n>=0} q^((n+1)^2) (q;q^2)_n      q F1* - G1 = theta4(q^2)(-q;q)_inf Y(q^2) (IV)
Horocycle: XX = XX0 + 4 m - i b,  q = e^(2 pi i/XX),  t = -2 pi i/XX,  dual nome r = e^(-pi i XX/2) fixed.
Prediction from addenda v6/v7 for the bilateral parts (growing sector):
   F0* + G0   ~ Nrm0(t) EE(r),   q F1* - G1 ~ Nrm1(t) EE(r),      EE(r) = psi(r) (r^2;r^4)_inf E+'(r),
   Nrm0 = sqrt(pi/t) e^(-pi^2/(24 t)) (c/2) e^(-i Cl2/(2t)),  Nrm1 = e^(2 pi i/3) Nrm0,  c = (4/3)^(1/4) e^(-i pi/12).
The script measures the limits of F0*/Nrm0, G0/Nrm0, q F1*/Nrm1, -G1/Nrm1 on a circle in r and prints their Taylor coefficients."""
import sys, os, json, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import acb, ctx
from mpmath import mp, mpf, mpc, nstr, exp, pi
import horo
from e1_limits import extrap
HERE = os.path.dirname(os.path.abspath(__file__))
MS = [200, 225, 255, 290, 330, 375, 425, 480, 500]
NAMES = ["F0", "G0", "qF1", "mG1"]
def tomp(z):
    z = z.mid(); return mpc(mpf(z.real.str(50, radius=False)), mpf(z.imag.str(50, radius=False)))

def values(Xp, prec, target=110):
    ctx.prec = prec
    Xc = horo.cx(*Xp)
    q = (2*acb.pi()*acb(0, 1)/Xc).exp().mid()
    K = float(abs(Xc).mid()); b = -float(Xc.imag.mid())
    J = int(K*math.sqrt((0.0808*K + 12 + target*0.6932)/(2*math.pi*b))) + 10
    one = acb(1)
    F0 = acb(0); F1 = acb(0); G0 = acb(0); G1 = acb(0)
    po = one; sq = one; q2n = one            # (q;q^2)_n, q^(n^2), q^(2n)
    for n in range(J):
        pn = (po*(one - q2n*q)).mid()        # (q;q^2)_(n+1)
        a = sq*sq/po; c = sq*sq*q2n/pn       # q^(2n^2)/(q;q^2)_n,  q^(2n^2+2n)/(q;q^2)_(n+1)
        if n & 1: F0 -= a; F1 -= c
        else: F0 += a; F1 += c
        if n >= 1: G0 += sq*po
        G1 += sq*q2n*q*po                    # q^((n+1)^2) (q;q^2)_n
        sq = (sq*q2n*q).mid(); q2n = (q2n*q*q).mid(); po = pn
    return F0, G0, q*F1, -G1, q

def point(X0, b):
    rows = []; I = acb(0, 1)
    for m in MS:
        Xp = (X0 + 4*m, b); K = float(X0 + 4*m)
        F0, G0, qF1, mG1, q = values(Xp, 300 + int(0.2*K))
        pi_ = acb.pi(); t = -2*pi_*I/horo.cx(*Xp)
        cl2 = acb((pi_*I/3).exp().polylog(2).imag)
        c = (acb(4)/3).sqrt().sqrt()*(-pi_*I/12).exp()
        n0 = (pi_/t).sqrt()*(-pi_*pi_/(24*t)).exp()*(c/2)*(-I*cl2/(2*t)).exp()
        n1 = n0*(2*pi_*I/3).exp()
        rows.append((tomp(t), [tomp(F0/n0), tomp(G0/n0), tomp(qF1/n1), tomp(mG1/n1)]))
    return extrap(rows, NAMES)

def EE_series(N=24):
    """integer coefficients of EE(r) = psi(r) (r^2;r^4)_inf E+'(r)"""
    def mul(a, b):
        c = [0]*N
        for i, x in enumerate(a):
            if x:
                for j in range(N - i): c[i + j] += x*b[j]
        return c
    psi = [0]*N; k = 0
    while k*(k + 1)//2 < N: psi[k*(k + 1)//2] = 1; k += 1
    pr = [1] + [0]*(N - 1)
    for e in range(2, N, 4):
        pr = [pr[i] - (pr[i - e] if i >= e else 0) for i in range(N)]
    E = [0]*N; den = [1] + [0]*(N - 1); m = 0
    while m*(m + 1)//2 < N:
        for i in range(N - m*(m + 1)//2): E[i + m*(m + 1)//2] += den[i]
        m += 1
        new = [0]*N                                   # den / (1 - r^(2m))
        for i in range(N): new[i] = den[i] + (new[i - 2*m] if i >= 2*m else 0)
        den = new
    return mul(mul(psi, pr), E)

if __name__ == "__main__":
    b = Fr(4, 5)
    path = os.path.join(HERE, "fstar_circle.json")
    res = json.load(open(path)) if os.path.exists(path) else {}
    for j in range(20):
        key = str(Fr(j, 5))
        if key not in res:
            res[key] = point(Fr(j, 5), b); json.dump(res, open(path, "w"), indent=0)
    mp.dps = 30
    print("20 points on |r| = e^(-2 pi/5); largest extrapolation uncertainty:", max(x[2] for v in res.values() for x in v.values()))
    r0 = exp(-pi*mpf(2)/5)
    co = {}
    for nm in NAMES:
        vals = [mpc(mpf(res[str(Fr(j, 5))][nm][0]), mpf(res[str(Fr(j, 5))][nm][1])) for j in range(20)]
        co[nm] = [sum(v*exp(1j*pi*j*n/10) for j, v in enumerate(vals))/20/r0**n for n in range(14)]
        print(f"{nm:4s} Taylor coefficients in r:", "  ".join(nstr(c, 9) for c in co[nm]))
    print("F0 + G0   :", "  ".join(nstr(co["F0"][n] + co["G0"][n], 9) for n in range(14)))
    print("qF1 - G1  :", "  ".join(nstr(co["qF1"][n] + co["mG1"][n], 9) for n in range(14)))
    print("EE(r)     :", EE_series()[:14])
