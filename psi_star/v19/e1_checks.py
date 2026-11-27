"""Addendum v19: checks for the asymptotics of the state integrals.
Potential (u = a'/b):  V(w) = -w^2/2 - 2 pi i a' w + Li2(-e^w);  size of e^(V/t) at arg t = -90 deg + eps:  exp(cos(eps) H(w)/|t|),
   H(w) = -Im V(w) + tau Re V(w),  tau = tan(eps)  (eps = arg x).
 a  growing integral a' = -1/2 (B4, B3): contour = ray {i pi - rho e^(-i 30deg)} + segment [i pi, 0] + positive real axis; saddle 2 pi i/3
 b  decaying integrals: horizontal line through the saddle w1 = Log W (lower root); a' = a = 1/2 (B2, B1) and other a
 c  the elementary inequality log(2 cos(eta/2)) >= (log 2/pi) (pi - 3 eta/2) on [0, 2 pi/3]
 d  B4, B2, B3, B1 computed exactly (Theorem C of v12) against (4/3)^(1/4) e^(-+ i pi/12) e^(-+ i Cl2/t)
 e  the coefficients of Theorem 8:  P I(1/(2b)) = +A e^(-5 pi i/12) t^(-1/2) e^(-S/t),  P I(-1/(2b)) = -A e^(5 pi i/12) t^(-1/2) e^(-conj(S)/t),
    q^(1/2) P I(b +- 1/(2b)) = e^(-+ 2 pi i/3) P I(+-1/(2b)),   A = sqrt(pi/sqrt 3),  S = pi^2/3 - i Cl2(pi/3)
Usage: python3 e1_checks.py [a b c d e]"""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr, radians, clsin
I = mpc(0, 1)

def li2_unit_disc(y):                      # numpy, |y| <= 1 (series with 400 terms plus direct mp for safety near 1 not needed here)
    return np.array([complex(polylog(2, complex(v))) for v in y])

def V_strip(w, ap):                        # principal Li2, |Im w| <= pi
    return -w*w/2 - 2j*np.pi*ap*w + li2_unit_or_any(-np.exp(w))
def li2_unit_or_any(y):
    return np.array([complex(polylog(2, complex(v))) for v in np.atleast_1d(y)])
def V_left(w, ap):                         # Re w < 0, any Im w: Li2(-e^w) with |e^w| < 1
    return -w*w/2 - 2j*np.pi*ap*w + li2_unit_or_any(-np.exp(w))

def Hfun(V, tau): return -V.imag + tau*V.real

def report(name, w, V, w0, V0, tau):
    d = np.abs(w - w0); H = Hfun(V, tau) - Hfun(np.array([V0]), tau)[0]
    m = d > 1e-9
    ratio = H[m]/np.minimum(d[m]**2, d[m])
    return f"{name}: max(H - H0) = {H[m].max():+.4f}, sup (H - H0)/min(d^2, d) = {ratio.max():+.4f}"

def part_a():
    mp.dps = 20
    print("(a) growing integral, a' = -1/2, saddle w0 = 2 pi i/3;  H0 = Cl2(pi/3) - tau 5 pi^2/12")
    ap = -0.5; w0 = 2j*np.pi/3; V0 = V_strip(np.array([w0]), ap)[0]
    print(f"    V(w0) = {V0:.12f}   (-5 pi^2/12 = {-5*np.pi**2/12:.12f}, -Cl2(pi/3) = {-float(clsin(2, pi/3)):.12f})")
    rho = np.concatenate([np.linspace(1e-6, 12, 2400)]); wl = 1j*np.pi - rho*np.exp(-1j*np.pi/6)
    eta = np.linspace(0, np.pi, 3001)[1:-1]; wa = 1j*eta
    r = np.linspace(1e-6, 12, 2400); wr = r + 0j
    Vl = V_left(wl, ap); Va = V_strip(wa, ap); Vr = V_strip(wr, ap)
    for tau in (0.0, 0.1, 0.2):
        print(f"    tau = {tau}:  " + report("left ray", wl, Vl, w0, V0, tau) + ";  " + report("segment", wa, Va, w0, V0, tau) + ";  " + report("real axis", wr, Vr, w0, V0, tau))
    print(f"    on the left ray: min Im V = {Vl.imag.min():+.5f} (claim >= 0), max Re V = {Vl.real.max():+.5f} (claim <= -pi^2/3 = {-np.pi**2/3:.5f})")
    print(f"    on the real axis: max |Im V - pi r| = {np.abs(Vr.imag - np.pi*r).max():.1e}, max Re V = {Vr.real.max():+.5f} (claim <= -pi^2/12 = {-np.pi**2/12:.5f})")
    print(f"    on the segment: max |Im V + Cl2(pi - eta)| = {max(abs(Va.imag[i] + float(clsin(2, np.pi - eta[i]))) for i in range(0, len(eta), 50)):.1e}, max |Re V - (3 eta^2/4 - pi eta - pi^2/12)| = {np.abs(Va.real - (0.75*eta**2 - np.pi*eta - np.pi**2/12)).max():.1e}")

def part_b():
    mp.dps = 20
    print("(b) decaying integrals: horizontal line through w1 = Log W, W the root of W^2 + W = e^(-2 pi i a) in the lower half-plane")
    for a, taus in ((0.5, (0.0, 0.2, 1.0, 1.7)), (0.25, (0.0, 0.05)), (1/3, (0.0, 0.05)), (5/6, (0.0, 0.05)), (0.03, (0.0, 0.01)), (0.97, (0.0, 0.01))):
        c = np.exp(-2j*np.pi*a); d = np.sqrt(1 + 4*c + 0j); W = (-1 + d)/2
        if W.imag > 0: W = (-1 - d)/2
        w1 = np.log(W); V1 = V_strip(np.array([w1]), a)[0]
        defect = abs(w1 + 2j*np.pi*a + np.log(1 + W))
        r = np.concatenate([np.linspace(-14, 14, 5601)]); w = r + 1j*w1.imag
        Vv = V_strip(w, a)
        out = "; ".join(f"tau = {tau}: " + report("line", w, Vv, w1, V1, tau).split(": ", 1)[1] for tau in taus)
        print(f"    a = {a:.4f}: Im w1/pi = {w1.imag/np.pi:+.4f}, critical-point defect {defect:.1e}, -Im V(w1) = {-V1.imag:+.6f};  " + out)

def part_c():
    eta = np.linspace(0, 2*np.pi/3, 200001)[:-1]
    R = np.log(2*np.cos(eta/2))/(np.pi - 1.5*eta)
    print(f"(c) min over [0, 2 pi/3) of log(2 cos(eta/2))/(pi - 3 eta/2) = {R.min():.6f} at eta = {eta[R.argmin()]:.4f}  (log 2/pi = {np.log(2)/np.pi:.6f}); limit at 2 pi/3: {R[-1]:.6f} (1/sqrt 3 = {1/np.sqrt(3):.6f})")
    eta2 = np.linspace(2*np.pi/3, np.pi, 200001)[1:-1]
    R2 = -np.log(2*np.cos(eta2/2))/(1.5*eta2 - np.pi)
    print(f"    min over (2 pi/3, pi) of -log(2 cos(eta/2))/(3 eta/2 - pi) = {R2.min():.6f}")

def part_d():
    from common import Setup
    mp.dps = 110
    Cl2 = clsin(2, pi/3); A = (mpf(4)/3)**(mpf(1)/4)
    print("(d) exact B_j (Theorem C of v12) divided by the asymptotic value, minus 1;  B4, B3 growing, B2, B1 decaying")
    for deg in ("0.5", "5", "11"):
        for r in (mpf(1)/100, mpf(1)/200, mpf(1)/400):
            x = r*exp(I*radians(mpf(deg))); S = Setup(x); b = S.b; t = -2*pi*I*x
            cp = S.cconst()*S.ppow(mpf(1)/8); q4 = S.qpow(mpf(1)/4)
            Iu = lambda u: S.cconst()*exp(-pi*I*u*u/2)*S.D(u)
            B2 = Iu(1/(2*b))/cp; B4 = Iu(-1/(2*b))/cp
            B1 = Iu(b + 1/(2*b))/(-I*cp/q4); B3 = Iu(b - 1/(2*b))/(I*cp/q4)
            e = exp(I*Cl2/t)
            vals = [B4/(A*exp(-I*pi/12)/e) - 1, B3/(A*exp(I*pi/12)/e) - 1, B2/(A*exp(I*pi/12)*e) - 1, B1/(A*exp(-I*pi/12)*e) - 1]
            print(f"    arg x = {deg:>4s} deg, |x| = 1/{int(1/r)}: " + "  ".join(f"{n}: {nstr(abs(v), 3)}" for n, v in zip(("B4", "B3", "B2", "B1"), vals)) + f"   (|t| = {nstr(abs(t), 3)}, |B4| = {nstr(abs(B4), 3)}, |B2| = {nstr(abs(B2), 3)})")

def part_e():
    from common import Setup
    mp.dps = 110
    Cl2 = clsin(2, pi/3); A = sqrt(pi/sqrt(mpf(3))); S_ = pi**2/3 - I*Cl2; Sb = pi**2/3 + I*Cl2
    print("(e) coefficients of Theorem 8: exact value divided by the asymptotic value, minus 1")
    for deg in ("1", "10"):
        for r in (mpf(1)/100, mpf(1)/400):
            x = r*exp(I*radians(mpf(deg))); S = Setup(x); b = S.b; t = -2*pi*I*x
            P = exp(pi**2/(12*t) + t/6)/(sqrt(2)*b)
            Iu = lambda u: S.cconst()*exp(-pi*I*u*u/2)*S.D(u)
            dec = P*Iu(1/(2*b)); gro = P*Iu(-1/(2*b))
            dec0 = S.qpow(mpf(1)/2)*P*Iu(b + 1/(2*b)); gro0 = S.qpow(mpf(1)/2)*P*Iu(b - 1/(2*b))
            vals = [dec/(A*exp(-5*pi*I/12)/sqrt(t)*exp(-S_/t)) - 1, gro/(-A*exp(5*pi*I/12)/sqrt(t)*exp(-Sb/t)) - 1,
                    dec0/(exp(-2*pi*I/3)*dec) - 1, gro0/(exp(2*pi*I/3)*gro) - 1]
            print(f"    arg x = {deg:>2s} deg, |x| = 1/{int(1/r)}: " + "  ".join(f"{n}: {nstr(abs(v), 3)}" for n, v in zip(("P I(1/(2b))", "P I(-1/(2b))", "psi0 decaying", "psi0 growing"), vals)))

if __name__ == "__main__":
    for p_ in (sys.argv[1:] or list("abcde")): {"a": part_a, "b": part_b, "c": part_c, "d": part_d, "e": part_e}[p_]()
