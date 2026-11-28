"""Addendum v30: the matrix of state integrals for real x and for Im x < 0 (numpy double precision, npint.py).
   (a) against the series of Theorem 20 for Im x > 0;  (b) Fourier modes on circles around points of the positive real axis;
   (c) Im x < 0: the q-difference equation of the three rows and the symmetry I(u; x) = I(u; 1/x)."""
import sys
import numpy as np
from mpmath import mp, mpf, mpc, nstr, matrix
from npint import NP
import frames
from common import Setup
mp.dps = 25
np.set_printoptions(linewidth=200, precision=6)
def tonp(M): return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])
a = 0.3
print("(a) numpy integrals against Theorem 20 (series) in the upper half-plane")
for x in (0.5 + 0.15j, 0.35 + 0.05j, 0.62 + 0.03j):
    S = Setup(mpc(x)); am = mpf(a)
    ser = np.linalg.solve(tonp(frames.frame_dual(S, am)), tonp(frames.Delta(S, am)*frames.frame_q(S, am)))
    num = NP(x).Smat(a)
    print(f"   x = {x}: max|S_integral - S_series| = {np.abs(num - ser).max():.2e}   (max entry {np.abs(ser).max():.3g})")
print("(b) holomorphy across the positive real axis: Fourier modes of S(a, x0 + r e^(i theta)) on a circle around a real point")
for x0, r, n in ((0.5, 0.15, 32), (0.3, 0.1, 32), (1.3, 0.3, 64)):
    th = 2*np.pi*np.arange(n)/n
    vals = np.array([NP(x0 + r*np.exp(1j*t)).Smat(a) for t in th])
    c0 = NP(x0).Smat(a)
    mean = vals.mean(axis=0)
    neg = [np.abs((vals*np.exp(1j*k*th)[:, None, None]).mean(axis=0)).max() for k in (1, 2, 3, 4)]   # coefficients of (x-x0)^(-k) r^(-k)
    pos = [np.abs((vals*np.exp(-1j*k*th)[:, None, None]).mean(axis=0)).max() for k in (1, 2, 3, 4)]
    print(f"   x0 = {x0}, r = {r}: max|mean - S(x0)| = {np.abs(mean - c0).max():.2e};  negative modes k=1..4: " + ", ".join(f"{v:.1e}" for v in neg) + ";  positive modes: " + ", ".join(f"{v:.1e}" for v in pos))
print("(c) lower half-plane: the q-difference equation of the first row, and I(u; x) = I(u; 1/x)")
for x in (0.5 - 0.15j, 0.3 - 0.1j, 0.8 - 0.3j):
    N = NP(x); z = np.exp(2j*np.pi*a); q = np.exp(2j*np.pi*x)
    c = [N.col(a + k*x) for k in range(3)]
    L = [c[0][i] - z*q*c[1][i] - z*q*q*c[2][i] for i in range(3)]
    N2 = NP(1/x); u = a/N.b
    print(f"   x = {x} (|q| = {abs(q):.3g}): L Jc - 1 = {abs(L[0] - 1):.2e}, L i0 = {abs(L[1]):.2e}, L i1 = {abs(L[2]):.2e};  |I(u;x) - I(u;1/x)|/|I| = {abs(N.Iu(u) - N2.Iu(u))/abs(N.Iu(u)):.2e}")
