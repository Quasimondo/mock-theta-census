"""The 2x2 block of Omega_gamma on both half-planes by series, and its Fourier modes on circles around points of the negative real axis.
   Im x > 0: I(u; x) by Theorem C.   Im x < 0: I(u; x) = I(u; 1/x) (Phi_b = Phi_{1/b}), Theorem C at 1/x."""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, matrix, nstr
from common import Setup, I
from a1_formula import Iser
from frames import Pelem, mnorm
mp.dps = 25

def Iany(x, u):
    return Iser(Setup(x), u) if x.imag > 0 else Iser(Setup(1/x), u)

def block(a, x):
    b = sqrt(x); P = Pelem(Setup(x))
    def col(a):
        f = I*P*exp(-pi*I*a)
        return f*Iany(x, a/b), f*Iany(x, (a - 1)/b)
    c0 = col(a); c1 = col(a + x)
    return matrix([[c0[0], c1[0]], [c0[1], c1[1]]])

def Omblock(a, x):
    return block(a, x)**-1*block(a/(1 + 2*x), x/(1 + 2*x))

if __name__ == "__main__":
    a = mpf("0.3")
    for x0, r, n in ((mpf("0.5"), mpf("0.2"), 32), (mpf("-0.25"), mpf("0.1"), 32), (mpf("-0.2"), mpf("0.15"), 64), (mpf("-0.7"), mpf("0.1"), 32)):
        th = [2*pi*(k + mpf(1)/2)/n for k in range(n)]
        vals = [Omblock(a, x0 + r*exp(I*t)) for t in th]
        neg = [mnorm(sum((v*exp(I*k*t) for v, t in zip(vals, th)), matrix(2, 2))/n) for k in (1, 2, 3)]
        pos = [mnorm(sum((v*exp(-I*k*t) for v, t in zip(vals, th)), matrix(2, 2))/n) for k in (0, 1, 2, 3)]
        jump = mnorm(Omblock(a, x0 + mpc(0, "0.002")) - Omblock(a, x0 - mpc(0, "0.002"))) if x0 < 0 else 0
        print(f"x0 = {nstr(x0, 4)}, r = {nstr(r, 3)}: negative modes k=1..3: " + ", ".join(nstr(v, 3) for v in neg) + ";  modes k=0..3: " + ", ".join(nstr(v, 3) for v in pos) + (f";  |Omega(x0+0.002i) - Omega(x0-0.002i)| = {nstr(jump, 3)}" if x0 < 0 else ""))
    print("jump of the block across the negative real axis, by series on both sides: max|Omega(x0 + i eps) - Omega(x0 - i eps)| / max|Omega|")
    for x0 in ("-0.1", "-0.2", "-0.3", "-0.35", "-0.4", "-0.45", "-0.6"):
        x0 = mpf(x0); out = []
        for e in ("0.004", "0.002", "0.001"):
            e = mpc(0, e)
            up = Omblock(a, x0 + e); lo = Omblock(a, x0 - e)
            out.append(nstr(mnorm(up - lo), 3) + " / " + nstr(mnorm(up), 3))
        print(f"   x0 = {nstr(x0, 3)}:  eps = 0.004, 0.002, 0.001:  " + ";  ".join(out))
