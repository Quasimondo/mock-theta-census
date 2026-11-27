# Independent check of addendum-v7 sec. 1:
#   I(u) = int_C Phi_b(y) exp(-pi i y^2/2 + 2 pi u y) dy  =  -i b (q;q)/(qt;qt) exp(-pi i c_b^2/2 + 2 pi u c_b) * B_j
# written from the formulas in the addendum (own Phi_b, own contour, own series).
import mpmath as mp
mp.mp.dps = 30
def poch_inf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a * q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps - 5): return s
        k += 1
def qq(q): return poch_inf(-1, q) if False else _qq(q)
def _qq(q):                         # (q;q)_inf
    s = mp.mpc(1); k = 1
    while True:
        f = q**k; s *= 1 - f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps - 5): return s
        k += 1
def setup(x):
    b = mp.sqrt(x); q = mp.expjpi(2*x); qt = mp.expjpi(-2/x)
    global QH, QTH
    QH = mp.expjpi(x); QTH = mp.expjpi(-1/x)       # q^{1/2}, qt^{1/2} on the right branch
    return b, q, qt
def Phi(y, b, q, qt):
    num = poch_inf(QH * mp.exp(2*mp.pi*b*y), q)            # (-q^{1/2} e^{2 pi b y}; q)_inf
    den = poch_inf(QTH * mp.exp(2*mp.pi*y/b), qt)          # (-qt^{1/2} e^{2 pi y/b}; qt)_inf
    return num / den
def I(u, x, delta=mp.mpf('0.4'), R=9):
    b, q, qt = setup(x)
    f = lambda y: Phi(y, b, q, qt) * mp.exp(-mp.pi*1j*y**2/2 + 2*mp.pi*u*y)
    e1 = mp.expj(mp.pi - delta); e2 = mp.expj(delta)
    pts = [R*k/18 for k in range(19)]
    left = -mp.quad(lambda r: f(r*e1) * e1, pts)            # from infinity e^{i(pi-delta)} to 0
    right = mp.quad(lambda r: f(r*e2) * e2, pts)
    return left + right
def sums(q, p):
    # N0, Y at q; E+, E+', E-, E-' at p (E's: sum sign^m p^{m(m+a)/2}/(p^2;p^2)_m)
    n0 = y = 0; P = mp.mpc(1)
    for n in range(200):
        n0 += (-1)**n * q**(n*n+n) / P
        P1 = P * (1 - q**(2*n+1)); y += (-1)**n * q**((n+1)**2) / P1
        P = P1 * (1 - q**(2*n+2))
    def E(sign, a):
        s = 0; Q = mp.mpc(1)
        for m in range(200):
            s += sign**m * p**(mp.mpf(m*(m+a))/2) / Q; Q *= 1 - p**(2*m+2)
        return s
    return n0, y, E(-1, 1), E(1, 1), E(-1, 3), E(1, 3)
for x in [1/mp.mpc(3, -0.4), mp.mpc(0.3, 0.2), mp.mpc(1.2, 0.5), mp.mpc(0.11, 0.03)]:
    b, q, qt = setup(x); p = mp.expjpi(-1/x)
    cb = 1j*(b + 1/b)/2
    n0, y, Ep, Epp, Em, Emp = sums(q, p)
    qm14 = mp.expjpi(-x/2)                                     # q^{-1/4}
    B = {1: n0*Epp + 1j*qm14*y*Ep, 3: n0*Emp - 1j*qm14*y*Em}
    U = {1: b/4, 3: b/4 - 1/(2*b)}
    for j in (1, 3):
        if True:
            rhs = -1j*b*_qq(q)/_qq(qt) * mp.exp(-mp.pi*1j*cb**2/2 + 2*mp.pi*U[j]*cb) * B[j]
            lhs = I(U[j], x)
            print(f"x = {mp.nstr(x, 4):>22}  B{j}:  |I/rhs - 1| = {mp.nstr(abs(lhs/rhs - 1), 3)}", flush=True)
