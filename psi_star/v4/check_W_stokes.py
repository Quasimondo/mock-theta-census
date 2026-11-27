# psi(q) P+(q) sqrt(t) e^{-pi^2/(30 t)} -> pi c0 = sqrt(2 pi/sqrt5)/2 at q = e^{-t}; first correction s1 = (20+3 sqrt5)/120
import mpmath as mp
mp.mp.dps = 60
def val(t):
    q = mp.e**(-t)
    psi = mp.nsum(lambda n: q**(n*(n+1)/2), [0, mp.inf])
    s = 0; P = mp.mpf(1); k = 0
    while True:
        term = q**(2*k*k+k) / P
        s += term
        if k > 5 and term < mp.mpf(10)**-70 * s: break
        k += 1; P *= 1 - q**(2*k)
    return psi * s * mp.sqrt(t) * mp.e**(-mp.pi**2 / (30*t))
target = mp.sqrt(2*mp.pi/mp.sqrt(5)) / 2
ts = [mp.mpf(1)/n for n in (200, 300, 400, 600, 800, 1200)]
vs = [val(t)/target for t in ts]
def nev(xs, ys):
    P = list(ys)
    for k in range(1, len(P)):
        for i in range(len(P)-k): P[i] = (xs[i+k]*P[i] - xs[i]*P[i+1])/(xs[i+k]-xs[i])
    return P[0]
print('ratio -> 1 ?', mp.nstr(nev(ts, vs), 20))
d = [(v-1)/t for v, t in zip(vs, ts)]
print('s1 estimate', mp.nstr(nev(ts, d), 15), ' claimed (20+3 sqrt5)/120 =', mp.nstr((20+3*mp.sqrt(5))/120, 15))
