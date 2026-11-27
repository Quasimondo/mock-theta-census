# psi(q) = q^{-1/8} sqrt(pi/(2 tau)) theta4(q'),  q = e^{2 pi i/Y}, tau = -2 pi i/Y, q' = e^{-pi i Y}, at complex Y
import mpmath as mp
mp.mp.dps = 40
for Y in [mp.mpc(7.3, -1.2), mp.mpc(-3.1, -0.7), mp.mpc(25, -2)]:
    q = mp.exp(2j*mp.pi/Y); tau = -2j*mp.pi/Y; qp = mp.exp(-1j*mp.pi*Y)
    psi = mp.nsum(lambda n: q**(n*(n+1)/2), [0, mp.inf])
    th4 = 1 + 2*mp.nsum(lambda n: (-1)**n * qp**(n*n), [1, mp.inf])
    rhs = mp.exp(-2j*mp.pi/(8*Y)) * mp.sqrt(mp.pi/(2*tau)) * th4
    print(mp.nstr(abs(psi/rhs - 1), 3))
