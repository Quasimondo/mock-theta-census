# Independent check of addendum-v32 (own code from check_v31.py / check_v31_np.py: state integrals by own quadrature,
# Phi_b by the product formula at b or 1/b, own q- and p-series).  Omega^(3)(a,x) := S(a,x+1)^-1 S^c(a',x'-1).
# (1) the word gamma = T S T^2 S T in SL2(Z) (up to sign);  (2) Theorem 23 above the axis against J^-1 Lambda^-1 J' (22 digits);
# (3) Omega^(3) against Omega^(1) below the axis and against Omega^(2) above it at Re x < 0 (double precision);
# Uses S, Sc, frames from ../v31/check_v31.py and ../v31/check_v31_np.py (own code of the v31 check).
# (4) Omega_gamma(a, 0) by Cauchy's mean value of Omega^(3) on circles |x| = r, never evaluated on the real axis.
import mpmath as mp, numpy as np, warnings, os
V31 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'v31')
warnings.filterwarnings('ignore')
src = open(os.path.join(V31, 'check_v31.py')).read()
exec(src.split("def mx(M)")[0]); mxm = lambda M: max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))
print('(1) the word')
Tm = np.array([[1, 1], [0, 1]]); Sm = np.array([[0, -1], [1, 0]]); g = np.array([[1, 0], [2, 1]])
w = Tm@Sm@Tm@Tm@Sm@Tm; print(f'  T S T^2 S T = {w.tolist()},  gamma = {g.tolist()}  (equal up to sign: {bool((w == g).all() or (w == -g).all())})')
print('(2) Theorem 23, Im x > 0, 22 digits: max|S(a,x+1)^-1 S^c(a\',x\'-1) - J^-1 Lambda^-1 J\'|')
for a, x in ((mp.mpf(0.25), mp.mpc(0.3, 0.4)), (mp.mpc(0.3, 0.1), mp.mpc(-0.15, 0.25))):
    ap, xp = a/(1 + 2*x), x/(1 + 2*x)
    Jq, _, dl, _, _ = frames(a, x); Jq1, _, dl1, _, _ = frames(ap, xp); lam = dl/dl1
    ser = Jq**-1*mp.diag([1, lam, -lam])**-1*Jq1
    O3 = S(a, x + 1)**-1*Sc(ap, xp - 1)
    print(f'  a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: {mp.nstr(mxm(O3 - ser),2)}  (|Om| {mp.nstr(mxm(ser),3)})', flush=True)
del S, Sc
srcn = open(os.path.join(V31, 'check_v31_np.py')).read()
exec(srcn.split("print('(3)")[0])
def O1(a, x): ap, xp = gam(a, x); return np.linalg.inv(S(a, x))@S(ap, xp)
def O2(a, x): ap, xp = gam(a, x); return np.linalg.inv(Sc(a, x))@Sc(ap, xp)
def O3(a, x): ap, xp = gam(a, x); return np.linalg.inv(S(a, x + 1))@Sc(ap, xp - 1)
print('(3) double precision: Omega^(3) against Omega^(1) (Im x < 0) and Omega^(2) (Im x > 0, Re x < 0)')
for x in (0.5-0.15j, -0.25-0.05j, -0.45-0.05j, 0.15-0.1j):
    print(f'  x = {x}: |O3 - O1| = {mx(O3(0.3, x) - O1(0.3, x)):.1e}', flush=True)
for x in (-0.25+0.05j, -0.4+0.1j):
    print(f'  x = {x}: |O3 - O2| = {mx(O3(0.3, x) - O2(0.3, x)):.1e}', flush=True)
print('(4) Omega_gamma(0.3, 0) = mean of Omega^(3) over |x| = r, N points at angles (k + 1/2) 2 pi/N')
for r, N in ((0.1, 16), (0.2, 24)):
    vals = [O3(0.3, r*np.exp(2j*pi*(k + 0.5)/N)) for k in range(N)]
    m = sum(vals)/N; neg = max(mx(sum(v*np.exp(2j*pi*kk*(k + 0.5)/N) for k, v in enumerate(vals))/N) for kk in (1, 2))
    print(f'  r = {r}: first row {np.round(m[0], 9)};  negative Fourier modes 1, 2: {neg:.1e}', flush=True)
print('  addendum §2: (1, 0.4612731569 - 0.3395605587i, 0.0197242516 + 0.4838433601i)')
