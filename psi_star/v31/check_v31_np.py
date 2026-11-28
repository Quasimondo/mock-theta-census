# Parts (3), (4) of check_v31.py: Theorem 22 (i), R(a',x') = R(a,x) for Im x < 0; Omega^(1) = Omega^(2) below the axis; Omega_gamma on (-1/2, 0) by Cauchy's mean value.
# Note: circles that come close to x = 0 lose accuracy (the line integrals at |x| < 0.07 are not reliable here) in double precision (numpy trapezoid on the line e^{-i d}R, same conventions), for points near the axis.
import numpy as np, warnings
warnings.filterwarnings('ignore')
pi = np.pi
def logprod(a, q):           # log prod_{k>=0} (1 + a q^k), a an array (summed in logs: no overflow)
    s = np.zeros_like(a); f = a.copy()
    while np.max(np.abs(f)) > 1e-18:
        s += np.log1p(f); f = f*q
    return s
def S(a, x, L=14.0, n=120001):
    b = np.sqrt(complex(x)); bb = b if x.imag > 0 else 1/b; xx = bb*bb
    qq = np.exp(2j*pi*xx); pp = np.exp(-1j*pi/xx)
    d = (pi/2 - abs(np.angle(b)))/2; e = np.exp(-1j*d)
    r = np.linspace(-L, L, n); y = r*e; h = r[1] - r[0]
    lPhi = logprod(np.exp(1j*pi*xx + 2*pi*bb*y), qq) - logprod(pp*np.exp(2*pi*y/bb), pp**2)
    lbase = -1j*pi*y*y - lPhi; ker = 1/(1 + np.exp(2*pi*y/b)); eh = e*h
    P = 2**-0.5*np.exp(-1j*pi*x/3)*np.exp(1j*pi/(24*x))/b
    def col(A):
        f = 1j*P*np.exp(-1j*pi*A)*eh; w = np.exp(lbase + 2*pi*A/b*y)
        return f*np.sum(w*ker), f*np.sum(w), f*np.sum(np.exp(lbase + 2*pi*(A - 1)/b*y))
    j0, i0, k0 = col(a); j1, i1, k1 = col(a + x)
    return np.array([[1, j0, j1], [0, i0, i1], [0, k0, k1]])
def Sc(a, x): return np.conj(S(-np.conj(a), -np.conj(x)))
def gam(a, x): return a/(1 + 2*x), x/(1 + 2*x)
mx = lambda M: np.max(np.abs(M))
print('(3) Theorem 22 (i), Im x < 0: R(a\',x\') = R(a,x);  (4) Omega^(1) = Omega^(2)')
for a, x in ((0.3, 0.5-0.15j), (0.25, -0.3-0.1j), (0.3, -0.25-0.05j), (0.3, -0.45-0.05j), (0.3+0.1j, 0.3-0.3j)):
    ap, xp = gam(a, x); S0, S1, C0, C1 = S(a, x), S(ap, xp), Sc(a, x), Sc(ap, xp)
    R0 = C0@np.linalg.inv(S0); R1 = C1@np.linalg.inv(S1); O1 = np.linalg.inv(S0)@S1; O2 = np.linalg.inv(C0)@C1
    print(f'  a = {a}, x = {x}: max|R\' - R| = {mx(R1-R0):.1e} (|R| {mx(R0):.3g});  max|Om1 - Om2| = {mx(O1-O2):.1e} (|Om| {mx(O1):.3g})')
print('    on the interval: Omega_gamma(x0) = mean of Omega over the circle |x - x0| = r (Omega^(2) above, Omega^(1) below), a = 0.3, N = 20')
for x0, rad in ((-0.2, 0.06), (-0.2, 0.1), (-0.3, 0.06), (-0.3, 0.1)):
    N = 20; acc = 0
    for k in range(N):
        x = x0 + rad*np.exp(2j*pi*(k + 0.5)/N); ap, xp = gam(0.3, x)
        acc = acc + (np.linalg.inv(Sc(0.3, x))@Sc(ap, xp) if x.imag > 0 else np.linalg.inv(S(0.3, x))@S(ap, xp))
    print(f'  x0 = {x0}, r = {rad}: first row {np.round(acc[0]/N, 7)}', flush=True)
