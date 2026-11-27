# Independent check of addendum-v11: T(x) = (-q;q)_inf / (b q^{1/8} (p;p^2)_inf) * int_L e^{-pi i y^2 + pi y/b} / (Phi_b(y)(1 + e^{2 pi y/b})) dy
import mpmath as mp
from tfun import run
mp.mp.dps = 40
def poch_inf(a, q):                 # prod_{k>=0} (1 + a q^k)
    s = mp.mpc(1); k = 0
    while True:
        f = a * q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps - 5): return s
        k += 1
def kernel_side(x, delta=mp.mpf('0.15'), R=7, pieces=28):
    b = mp.sqrt(x); q = mp.expjpi(2*x); qt = mp.expjpi(-2/x)
    QH, QTH = mp.expjpi(x), mp.expjpi(-1/x)              # q^{1/2}, qt^{1/2} on the right branch
    Phi = lambda y: poch_inf(QH*mp.exp(2*mp.pi*b*y), q) / poch_inf(QTH*mp.exp(2*mp.pi*y/b), qt)
    f = lambda y: mp.exp(-mp.pi*1j*y**2 + mp.pi*y/b) / (Phi(y) * (1 + mp.exp(2*mp.pi*y/b)))
    e = mp.expj(-delta)
    pts = [-R + 2*R*k/pieces for k in range(pieces+1)]
    I = mp.quad(lambda r: f(r*e)*e, pts)
    p = mp.expjpi(-1/x)
    mq = poch_inf(q, q); pp2 = mp.mpc(1)
    for k in range(300): pp2 *= 1 - p**(2*k+1)
    return mq / (b * mp.expjpi(x/4) * pp2) * I, abs(f(R*e)), abs(f(-R*e))
for x in [mp.mpc(0.3, 0.2), 1/mp.mpc(3, -0.4), mp.mpc(0.11, 0.03)]:
    with mp.workdps(90):
        T, psi1 = run(1/x)
    K, e1, e2 = kernel_side(x)
    print(f'x = {mp.nstr(x,4):>22}:  T = {mp.nstr(T,15)}   kernel side = {mp.nstr(K,15)}   |diff| = {mp.nstr(abs(T-K),3)}   (psi1* = {mp.nstr(psi1,6)}; integrand at ends {mp.nstr(e1,2)}, {mp.nstr(e2,2)})', flush=True)
