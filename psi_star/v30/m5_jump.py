"""Does Omega_gamma(a, x) = S(a,x)^-1 S(a',x') continue across the negative real axis?
   Above the axis: Omega = J(a;x)^-1 Lambda^-1 J(a';x') by series.  Below: state integrals (lowint.py), cross-checked against npint.py and the series of the block."""
import numpy as np
from mpmath import mp, mpf, mpc, nstr
from lowint import Low
from npint import NP
import frames, m3_modular, m4_block
from common import Setup
mp.dps = 20
def tonp(M): return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])
a = 0.3

def Om_up(x):
    x = mpc(x); am = mpf(a); xp = x/(1 + 2*x); ap = am/(1 + 2*x)
    return tonp(frames.frame_q(Setup(x), am)**-1*m3_modular.Lam(am, x)**-1*frames.frame_q(Setup(xp), ap))

def Om_low(x):
    xp = x/(1 + 2*x); ap = a/(1 + 2*x)
    um = max(abs(a/np.sqrt(x)), abs((a - 1)/np.sqrt(x)), abs((a + x)/np.sqrt(x)), abs((a + x - 1)/np.sqrt(x))) + 0.2
    ump = max(abs(ap/np.sqrt(xp)), abs((ap - 1)/np.sqrt(xp)), abs((ap + xp)/np.sqrt(xp)), abs((ap + xp - 1)/np.sqrt(xp))) + 0.2
    return np.linalg.solve(Low(x, umax=um).Smat(a), Low(xp, umax=ump).Smat(ap))

if __name__ == "__main__":
    print("(a) lowint.py against npint.py (integral representation of Phi_b) and against the series of the block, Im x < 0")
    for x in (0.5 - 0.15j, 0.2 - 0.4j, -0.2 - 0.3j):
        S1 = Low(x).Smat(a)
        d1 = np.abs(S1 - NP(x).Smat(a)).max() if x.real > 0 else float("nan")
        blk = tonp(m4_block.block(mpf(a), mpc(x)))
        print(f"   x = {x}: max|S_low - S_np| = {d1:.2e};  max|block - series| = {np.abs(S1[1:, 1:] - blk).max():.2e}")
    print("(b) jump of Omega_gamma across the negative real axis: max over the first row | max over the block, at x0 +- i eps")
    for x0 in (-0.15, -0.25, -0.35, -0.7):
        out = []
        for eps in (0.08, 0.04, 0.02):
            U = Om_up(x0 + 1j*eps); L = Om_low(x0 - 1j*eps); D = np.abs(U - L)
            out.append(f"{D[0].max():.3g} | {D[1:].max():.3g}  (sizes {np.abs(U[0]).max():.3g} | {np.abs(U[1:]).max():.3g})")
        print(f"   x0 = {x0}: eps = 0.08: {out[0]};   0.04: {out[1]};   0.02: {out[2]}")
