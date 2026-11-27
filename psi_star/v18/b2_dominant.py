"""Addendum v18, task 2: which term of the exact formula carries the volume.
For 0 < a < 1 and x -> 0 close to the real axis (x = r e^(i phi), phi = 1 degree; |q| < 1):
  (i)  the kernel-free integral I((a-1)/b), computed exactly from Theorem C of v12, against its saddle-point value
         I_sp = i sqrt(2) b e^(pi i a) (1 + W) e^(w/2) sqrt(pi/(t (1 + 2W))) e^(V(w)/t)       (W, w, V by the rule of b1_roots.py);
  (ii) F(z) against  -i P e^(-pi i a) I((a-1)/b) + 1/(1 - 2z)   (the m = 1 term of the dual sum plus the leading perturbative term).
Usage: python3 b2_dominant.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr, radians
from common import Setup, I
from b1_roots import rule

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    for name, a in (("-1", mpf(1)/2), ("i", mpf(1)/4), ("omega", mpf(1)/3), ("-omega", mpf(5)/6)):
        W, w, k, defect, V = rule(a); z = exp(2*pi*I*a)
        print(f"z = {name}")
        for r in (mpf(1)/100, mpf(1)/200, mpf(1)/400):
            x = r*exp(I*radians(1)); S = Setup(x); b = S.b; t = -2*pi*I*x
            u = (a - 1)/b
            Iex = S.cconst()*exp(-pi*I*u*u/2)*S.D(u)
            Isp = I*sqrt(2)*b*exp(pi*I*a)*(1 + W)*exp(w/2)*sqrt(pi/(t*(1 + 2*W)))*exp(V/t)
            P = exp(pi**2/(12*t) + t/6)/(sqrt(2)*b)
            F = S.Psi(z); dom = -I*P*exp(-pi*I*a)*Iex
            # the decaying partner I(a/b): root in the lower half-plane, w = Log W, generic one-loop formula
            Wl = -1 - W; wl = mp.log(Wl); kl = (-wl - mp.log(1 + Wl))/(2*pi*I) - a
            Vl = -wl*wl/2 - 2*pi*I*(a + kl)*wl + mp.polylog(2, -Wl)
            I0 = S.cconst()*exp(-pi*I*(a/b)**2/2)*S.D(a/b)
            I0sp = sqrt(2*pi*t*(1 + Wl)/(1 + 2*Wl))/(2*pi*b)*exp(Vl/t)
            extra = f";  I(a/b)/saddle(lower root, k = {nstr(kl, 3)}) = {nstr(I0/I0sp, 6)}, |I(a/b)/I((a-1)/b)| = {nstr(abs(I0/Iex), 3)}"
            print(f"    |x| = 1/{int(1/r)}: I/I_sp - 1 = {nstr(Iex/Isp - 1, 5)};  |F| = {nstr(abs(F), 5)},  |F - dominant|  = {nstr(abs(F - dom), 5)},  |F - dominant - 1/(1-2z)| = {nstr(abs(F - dom - 1/(1 - 2*z)), 3)}" + extra)
