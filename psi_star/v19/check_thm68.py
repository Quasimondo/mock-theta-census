# Independent check of addendum-v19 (own code, no common.py).
# (1) Theorem 8 against the q-series themselves: psi1*(q), psi0*(q) summed directly, divided by the growing
#     term A t^{-1/2} e^{5 pi i/12} e^{-Sbar/t} W(p)  (psi0*: -A e^{i pi/12} ...); trivial and decaying sectors
#     are exponentially smaller here.  Expect ratio - 1 = O(t), i.e. (ratio - 1)/t roughly constant.
# (2) Theorem 6 for the growing integral I(-1/(2b)) by direct quadrature (own Phi_b, own contour) against
#     (t/2 pi)^{1/2} 3^{-1/4} e^{-i pi/12} b^{-1} exp((-5 pi^2/12 - i Cl2)/t).
import mpmath as mp, sys
mp.mp.dps = 40
Cl2 = mp.clsin(2, mp.pi/3)
A = mp.sqrt(mp.pi/mp.sqrt(3))
Sb = mp.pi**2/3 + 1j*Cl2

def psistar(x, shift):              # shift 0: psi1* = sum (-1)^n q^{n(n+1)/2}(-q)_n ; shift 1: psi0* = sum (-1)^n q^{(n+1)(n+2)/2}(-q)_n
    q = mp.expjpi(2*x); s = 0; P = mp.mpc(1); n = 0
    while True:
        e = (n+shift)*(n+shift+1)//2
        term = (-1)**n * q**e * P
        s += term
        if abs(q)**e * abs(P) < mp.mpf(10)**(-mp.mp.dps) and n > 10: return s
        n += 1; P *= 1 + q**n

def Wdual(p):
    s = 0; po = mp.mpc(1); j = 0
    while abs(p)**(j*j) > mp.mpf(10)**(-mp.mp.dps):
        s += (-1)**j * p**(j*j+j) * po; po *= 1 - p**(2*j+1); j += 1
    return s

print('(1) Theorem 8, growing sector, against the q-series')
for deg in (1, 5, 10):
    for N in (50, 100, 200, 400):
        x = mp.expj(mp.radians(deg))/N
        with mp.workdps(40 + N//5):
            t = -2j*mp.pi*x; p = mp.expjpi(-1/x)
            G = A*t**(-0.5)*mp.exp(-Sb/t)*Wdual(p)
            r1 = psistar(x, 0)/(mp.expjpi(mp.mpf(5)/12)*G)
            r0 = psistar(x, 1)/(-mp.expjpi(mp.mpf(1)/12)*G)
            print(f'  arg x = {deg:2d} deg, |x| = 1/{N}:  psi1*: (ratio-1)/t = {mp.nstr((r1-1)/t, 6):>24}   psi0*: (ratio-1)/t = {mp.nstr((r0-1)/t, 6):>24}   |psi1*| = {mp.nstr(abs(G), 3)}', flush=True)

if len(sys.argv) > 1 and sys.argv[1] == 'quick': sys.exit()

print('(2) Theorem 6, growing integral I(-1/(2b)), direct quadrature')
mp.mp.dps = 20
def poch(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-3): return s
        k += 1
for deg, N in ((10, 6), (10, 12), (10, 24)):
    x = mp.expj(mp.radians(deg))/N; b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x); qh = mp.expjpi(x)
    Phi = lambda y: poch(qh*mp.exp(2*mp.pi*b*y), q)/poch(p*mp.exp(2*mp.pi*y/b), p**2)
    u = -1/(2*b)
    f = lambda y: mp.exp(-mp.pi*1j*y**2 + 2*mp.pi*u*y)/Phi(y)
    delta = (mp.pi/2 - mp.arg(b))/2; e = mp.expj(-delta); R = 4 + N/3
    I = mp.quad(lambda r: f(r*e)*e, mp.linspace(-R, R, 40))
    t = -2j*mp.pi*x
    asy = mp.sqrt(t/(2*mp.pi))*mp.mpf(3)**(-0.25)*mp.expjpi(-mp.mpf(1)/12)/b*mp.exp((-5*mp.pi**2/12 - 1j*Cl2)/t)
    print(f'  arg x = {deg} deg, |x| = 1/{N}:  |I| = {mp.nstr(abs(I), 4)}   (I/asy - 1)/t = {mp.nstr((I/asy-1)/t, 6)}', flush=True)

print('(3) the Gamma_2 inequality of v19 sec. 5: g(eta) - g(2pi/3) <= -0.015 (eta - 2pi/3)^2, tau in [0, 0.2]')
mp.mp.dps = 20
worst = -1
for i in range(41):
    tau = mp.mpf(i)/200
    g = lambda e: mp.clsin(2, mp.pi - e) + tau*(0.75*e**2 - mp.pi*e - mp.pi**2/12)
    g0 = g(2*mp.pi/3)
    for k in range(1, 800):
        e = mp.pi*k/800
        if abs(e - 2*mp.pi/3) < 1e-9: continue
        worst = max(worst, (g(e) - g0)/(e - 2*mp.pi/3)**2)
print(f'  sup over tau, eta of (g - g0)/(eta - 2pi/3)^2 = {mp.nstr(worst, 5)}  (claim <= -0.015)')
