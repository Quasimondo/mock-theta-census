# Independent check of addendum-v15 Proposition 1 (own code, mpmath):
# psi(rho) = Re(e^{-i theta} G(w)) - pi rho cos(gamma), w = +-rho e^{i alpha}; exact formula (i), bound (iii), kappa (iv).
import random
import mpmath as mp
mp.mp.dps = 20
def G(w):
    if mp.re(w) <= 0: return -w**2/2 + mp.polylog(2, -mp.exp(w)) + mp.pi**2/12
    return -w**2 - mp.polylog(2, -mp.exp(-w)) - mp.pi**2/12
def lam(al): return 0 if al == 0 else mp.log(1 - mp.exp(-(mp.pi/2)/mp.tan(abs(al))))
random.seed(1); worst_i = 0; viol = 0
for _ in range(400):
    beta = mp.radians(random.uniform(-85, 85)); delta = mp.radians(random.uniform(0.01, 89.99 - abs(float(mp.degrees(beta)))))
    al, ga = beta - delta, beta + delta; theta = 2*beta - mp.pi/2; rho = mp.mpf(random.uniform(0.01, 30))
    for side in (+1, -1):
        w = side*rho*mp.expj(al)
        psi = mp.re(mp.expj(-theta)*G(w)) - mp.pi*rho*mp.cos(ga)
        ell = mp.quad(lambda l: mp.log(1 + mp.exp(-l*w)) if side > 0 else mp.log(1 + mp.exp(l*w)), mp.linspace(0, 1, int(4 + 3*rho)))
        X, Y = mp.re(ell), mp.im(ell)
        f = (-(mp.pi - Y)*mp.cos(ga) - X*mp.sin(ga) - rho*mp.sin(2*delta)) if side > 0 else (-(mp.pi + Y)*mp.cos(ga) + X*mp.sin(ga) - rho/2*mp.sin(2*delta))
        worst_i = max(worst_i, abs(psi/rho - f))
        kap = (mp.pi/2)*mp.cos(ga) - max(mp.log(2), -lam(al))*abs(mp.sin(ga))
        bound = -kap*rho - (1 if side > 0 else mp.mpf(1)/2)*rho**2*mp.sin(2*delta)
        if psi > bound + 1e-12: viol += 1
print('(i) exact formula, 800 random points: max |psi/rho - formula| =', mp.nstr(worst_i, 3))
print('(iii) bound psi <= -kappa rho - c rho^2 sin 2delta: violations =', viol)
d = mp.radians(10); kmin = min(((mp.pi/2)*mp.cos(b + d) - max(mp.log(2), -lam(b - d))*abs(mp.sin(b + d)), float(mp.degrees(b))) for b in [mp.radians(-60 + 0.05*k) for k in range(2201)])
print('(iv) delta = 10 deg, -60 <= beta <= 50: min kappa =', mp.nstr(kmin[0], 5), 'at beta =', round(kmin[1], 2))
