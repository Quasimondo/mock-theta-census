# Independent check of addendum-v18 Theorem Z (strip form) at z = i and z = -omega:
#  F(z) = i P e^{-pi i a} [ J(a/b) - I((a-1)/b) W(T) + p T I(a/b) W(p^2 T) ],  T = e^{pi i (2a-1)/x},
#  P = (-q;q)/(b q^{1/8} (p;p^2)),  W(t) = sum (-t)^k p^{k^2+k} (p;p^2)_k.   Own I (Theorem C series), own J (quadrature).
import sys
import mpmath as mp
sys.path.insert(0, '../v12')
import check_thmAC as C
mp.mp.dps = 30
def Wd(tt, p, K=80):
    s = 0; P = mp.mpc(1)
    for k in range(K):
        s += (-tt)**k * p**(k*k+k) * P; P *= 1 - p**(2*k+1)
    return s
for x in [mp.mpc(0.3, 0.2), mp.mpc(0.45, 0.3)]:
    b, q, p = C.setup(x)
    mq = C.poch_inf(q, q); pp2 = mp.mpc(1)
    for k in range(300): pp2 *= 1 - p**(2*k+1)
    P = mq/(b*mp.expjpi(x/4)*pp2)
    for a, name in [(mp.mpf(1)/4, 'z = i'), (mp.mpf(5)/6, 'z = -omega')]:
        z = mp.expjpi(2*a); T = mp.expjpi((2*a - 1)/x)
        F = C.Psi(z, q)
        rhs = 1j*P*mp.expjpi(-a)*(C.J_quad(a/b, x) - C.I_series((a-1)/b, x)*Wd(T, p) + p*T*C.I_series(a/b, x)*Wd(p**2*T, p))
        print(f'x = {mp.nstr(x,3)}, {name}: F = {mp.nstr(F, 14)}   rhs = {mp.nstr(rhs, 14)}   |diff| = {mp.nstr(abs(F - rhs), 3)}', flush=True)
