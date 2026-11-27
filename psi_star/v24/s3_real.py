"""Addendum v24: the Mordell identification at real x = 1/X (q on the unit circle), double precision.
   M1 = 2 i q^(11/60) P J(1/b),  M0 = 2 i q^(29/60) P J(b + 1/b)  against  sqrt(60 t/pi) L(2/5, 10 t), sqrt(60 t/pi) L(1/5, 10 t),  t = -2 pi i/X.
   P J over the line L (10 degrees) with Faddeev's integral for Phi_b (n1_zroots.py of v22, fad.py of v16)."""
import numpy as np
from mpmath import mp, mpf, mpc, exp, pi, sqrt, nstr
from n1_zroots import PJ
from s2_law import L
mp.dps = 20
I = mpc(0, 1)
for X in (5, 8, 20, 40):
    t = -2*pi*I/X; q = exp(-t)
    M1 = 2*I*exp(-t*mpf(11)/60)*mpc(complex(PJ(1.0, X)))
    M0 = 2*I*exp(-t*mpf(29)/60)*mpc(complex(PJ(1.0 + 1.0/X, X)))
    m1 = sqrt(60*t/pi)*L(mpf(2)/5, 10*t); m0 = sqrt(60*t/pi)*L(mpf(1)/5, 10*t)
    print(f"X = {X:3d}: M1 = {nstr(M1, 12)}, |M1 - Mordell| = {nstr(abs(M1 - m1), 3)};   M0 = {nstr(M0, 12)}, |M0 - Mordell| = {nstr(abs(M0 - m0), 3)}")
