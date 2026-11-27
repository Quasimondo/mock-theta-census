# Independent check of addendum-v12, Theorems C and A(up), own Phi_b, own quadrature, own series.
#  I(u) = int_L e^{-pi i y^2 + 2 pi u y}/Phi_b(y) dy,  J(u) = int_L (same)/(1 + e^{2 pi y/b}) dy,  L = e^{-i delta} R
#  Thm C:  I(u) = c e^{-pi i u^2/2} D(u)
#  Thm A:  J(u) = -i b q^{1/8} e^{pi i u b} (p;p^2)/(-q;q) Psi(e^{2 pi i u b}) - sum_{m>=1} (-1)^m I(u - m/b)
import mpmath as mp
mp.mp.dps = 30
def poch_inf(a, q):                 # prod (1 + a q^k)
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-5): return s
        k += 1
def setup(x):
    b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x)
    return b, q, p
def Phi(y, x):
    b, q, p = setup(x)
    return poch_inf(mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q) / poch_inf(p*mp.exp(2*mp.pi*y/b), p**2)
def quadL(fun, x, R=6, pieces=24):
    b = mp.sqrt(x); beta = mp.arg(b); delta = (mp.pi/2 - beta)/2
    e = mp.expj(-delta)
    pts = [-R + 2*R*k/pieces for k in range(pieces+1)]
    return mp.quad(lambda r: fun(r*e)*e, pts)
def I_quad(u, x):
    b = mp.sqrt(x)
    return quadL(lambda y: mp.exp(-mp.pi*1j*y**2 + 2*mp.pi*u*y)/Phi(y, x), x)
def J_quad(u, x):
    b = mp.sqrt(x)
    return quadL(lambda y: mp.exp(-mp.pi*1j*y**2 + 2*mp.pi*u*y)/(Phi(y, x)*(1 + mp.exp(2*mp.pi*y/b))), x)
def I_series(u, x, K=60):
    b, q, p = setup(x)
    c = 1/mp.sqrt(2) * mp.expjpi(-mp.mpf(1)/4) * mp.expjpi(-x/12) * mp.expjpi(-1/(12*x))
    D = 0; qk = mp.mpc(1)
    for k in range(K):
        if k: qk *= 1 - q**k
        pn = mp.mpc(1)
        for n in range(K):
            if n: pn *= 1 - p**(2*n)
            D += (-1)**(k*n+n) * mp.exp(mp.pi*1j*u*(b*k + n/b)) * mp.expjpi(x*k*k/2) * mp.expjpi(-(n*n/mp.mpf(2) + n)/x) / (qk*pn)
    return c * mp.exp(-mp.pi*1j*u**2/2) * D
def Psi(z, q, N=200):
    s = 0; P = mp.mpc(1)
    for n in range(N):
        s += q**(n*(n+1)//2) * z**n * P; P *= 1 + q**(n+1)
    return s
for x in [mp.mpc(0.3, 0.2)]:
    b, q, p = setup(x)
    for u in [mp.mpc(0.3, -0.2), 1/(2*b), b + 1/(2*b)]:
        Iq = I_quad(u, x); Is = I_series(u, x)
        print(f'x={mp.nstr(x,3)} u={mp.nstr(u,4)}:  Thm C |I_quad/I_series - 1| = {mp.nstr(abs(Iq/Is - 1), 3)}', flush=True)
    for u in [mp.mpc(0.3, -0.2), 1/(2*b)]:
        pp2 = mp.mpc(1)
        for k in range(300): pp2 *= 1 - p**(2*k+1)
        mq = poch_inf(q, q)
        res = -1j*b*mp.expjpi(x/4) * mp.exp(mp.pi*1j*u*b) * pp2/mq * Psi(mp.exp(2j*mp.pi*u*b), q)
        dual = mp.fsum((-1)**m * I_series(u - m/b, x) for m in range(1, 40))
        Jq = J_quad(u, x)
        print(f'   u={mp.nstr(u,4)}:  Thm A(up)  J = {mp.nstr(Jq, 12)}   rhs = {mp.nstr(res - dual, 12)}   |diff| = {mp.nstr(abs(Jq - res + dual), 3)}', flush=True)
