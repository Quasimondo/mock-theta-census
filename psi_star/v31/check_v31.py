# Independent check of addendum-v31 (own code: state integrals by quadrature on e^{-i d}R with d = (90deg - |arg b|)/2,
# Phi_b from the product formula at b if Im x > 0 and at 1/b if Im x < 0 (Phi_b = Phi_{1/b}); own q-, p-series for the frames).
# (1) Proposition 1: Jt(-a;x) S^c(a,x) = Delta(a,x) J(a;x), Im x > 0.   (2) the dual shift S(a+1) = C S(a), both half-planes.
# (3), (4): see check_v31_np.py.  Here also: Omega^(2) against the series J^-1 Lambda^-1 J' above the axis at Re x < 0.  (5) the example V != 1 of Lemma U for |q| < 1.
import mpmath as mp
mp.mp.dps = 22
def pinf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-3): return s
        k += 1
def S(a, x):
    b = mp.sqrt(x); bb = b if mp.im(x) > 0 else 1/b; xx = bb*bb
    qq = mp.expjpi(2*xx); pp = mp.expjpi(-1/xx)
    Phi = lambda y: pinf(mp.expjpi(xx)*mp.exp(2*mp.pi*bb*y), qq)/pinf(pp*mp.exp(2*mp.pi*y/bb), pp**2)
    d = (mp.pi/2 - abs(mp.arg(b)))/2; e = mp.expj(-d)
    def integ(u, ker):
        f = lambda r: mp.exp(-mp.pi*1j*(r*e)**2 + 2*mp.pi*u*r*e)/Phi(r*e)*(1/(1 + mp.exp(2*mp.pi*r*e/b)) if ker else 1)*e
        return mp.quad(f, mp.linspace(-9, 9, 37))
    P = 2**mp.mpf(-0.5)*mp.expjpi(-x/3)*mp.expjpi(1/(24*x))/b
    def col(A):
        f = 1j*P*mp.expjpi(-A); u = A/b
        return f*integ(u, True), f*integ(u, False), f*integ((A - 1)/b, False)
    j0, i0, k0 = col(a); j1, i1, k1 = col(a + x)
    return mp.matrix([[1, j0, j1], [0, i0, i1], [0, k0, k1]])
def Sc(a, x): return S(-mp.conj(a), -mp.conj(x)).apply(mp.conj)
def H(y, r, r14, part=None):
    s = 0; rm = mp.mpc(1)
    for m in range(120):
        if m: rm *= 1 - r**m
        if (part == 'e' and m % 2) or (part == 'o' and not m % 2): continue
        s += y**m*r14**(m*m)/rm
    return s
def pm(q, n):
    r = mp.mpc(1)
    for j in range(1, n+1): r *= 1 + q**j
    return r
def frames(a, x):
    q = mp.expjpi(2*x); p = mp.expjpi(-1/x); ze = mp.expjpi(a); z = ze*ze; zt = mp.expjpi(a/x); T = p*zt*zt
    F = lambda Z: sum(Z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(200))
    q14 = mp.expjpi(x/2); mu = mp.expjpi(-3*x/2)/ze; qh = mp.expjpi(x)
    Jq = mp.matrix([[1, F(z), F(z*q)], [0, H(ze, q, q14), mu*H(ze*qh, q, q14)], [0, H(-ze, q, q14), -mu*H(-ze*qh, q, q14)]])
    def W(TT):
        s = 0; po = mp.mpc(1)
        for k in range(80): s += (-TT)**k*p**(k*k+k)*po; po *= 1 - p**(2*k+1)
        return s
    p2 = p*p; p12 = mp.expjpi(-1/(2*x))
    Mt = mp.matrix([[H(zt*p, p2, p12, 'e'), -H(zt*p, p2, p12, 'o')], [-zt*p12*H(zt*p2, p2, p12, 'o'), zt*p12*H(zt*p2, p2, p12, 'e')]])**-1
    Jt = mp.matrix([[1, p*T*W(p2*T), -W(T)], [0, Mt[0, 0], Mt[0, 1]], [0, Mt[1, 0], Mt[1, 1]]])
    dl = mp.expjpi(mp.mpf(1)/4)/(2*mp.sqrt(x))*mp.expjpi(-5*x/12)*mp.expjpi(-1/(24*x))*mp.expjpi(-a - a*a/(2*x))
    return Jq, Jt, dl, p, T
def mx(M): return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))
def gam(a, x): return a/(1 + 2*x), x/(1 + 2*x)
print('(1) Proposition 1')
for a, x in ((mp.mpf(1)/4, mp.mpc(0.3, 0.2)), (mp.mpc(0.3, 0.1), mp.mpc(-0.2, 0.6))):
    Jq, _, dl, _, _ = frames(a, x); _, Jtm, _, _, _ = frames(-a, x)
    print(f'  a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: {mp.nstr(mx(Jtm*Sc(a, x) - mp.diag([1, dl, dl])*Jq),2)}  (max entry {mp.nstr(mx(Jq),3)})', flush=True)
print('(2) dual shift S(a+1) = C S(a), C = [1,-1,0; 0,p,-p/T; 0,-1,0]')
for a, x in ((mp.mpf(1)/4, mp.mpc(0.3, 0.2)), (mp.mpf(0.3), mp.mpc(0.5, -0.15)), (mp.mpc(0.3, 0.1), mp.mpc(-0.2, -0.3))):
    p = mp.expjpi(-1/x); T = p*mp.expjpi(2*a/x)
    C = mp.matrix([[1, -1, 0], [0, p, -p/T], [0, -1, 0]]); S0 = S(a, x)
    print(f'  a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: {mp.nstr(mx(S(a + 1, x) - C*S0),2)}  (max entry {mp.nstr(mx(S0),3)})', flush=True)
print('(4\') Omega^(2) above the axis at Re x < 0 against the series J^-1 Lambda^-1 J\'')
for a, x in ((mp.mpf(0.3), mp.mpc(-0.25, 0.05)), (mp.mpf(0.3), mp.mpc(-0.4, 0.1))):
    ap, xp = gam(a, x); Jq, _, dl, _, _ = frames(a, x); Jq1, _, dl1, _, _ = frames(ap, xp)
    lam = dl/dl1; ser = Jq**-1*mp.diag([1, lam, -lam])**-1*Jq1
    print(f'  a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: {mp.nstr(mx(Sc(a, x)**-1*Sc(ap, xp) - ser),2)}  (|Om| {mp.nstr(mx(ser),3)})', flush=True)
# parts (3) and (4) are in check_v31_np.py (double precision; the points near the axis are too slow at 22 digits)
print('(5) Lemma U fails for |q| < 1: V = [1,-F(z),-F(zq); 0,2,0; 0,0,2]')
q = mp.mpc(0.2, 0.3); z = mp.mpc(0.7, -0.4)
F = lambda Z: sum(Z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(120))
A = lambda z: mp.matrix([[1, 0, -1/(z*q*q)], [0, 0, 1/(z*q*q)], [0, 1, -1/q]])
V = lambda z: mp.matrix([[1, -F(z), -F(z*q)], [0, 2, 0], [0, 0, 2]])
print(f'  max|A(z) V(zq) - V(z) A(z)| = {mp.nstr(mx(A(z)*V(z*q) - V(z)*A(z)),2)}')
