import mpmath as mp
mp.mp.dps = 40
def exact_log_inv_poch(w, h):        # log 1/(w;q)_inf, q = e^{-h}
    q = mp.exp(-h); s = 0; k = 0
    while True:
        t = w*q**k; s -= mp.log(1 - t)
        if abs(t) < mp.mpf(10)**-40: return s
        k += 1
h = mp.mpc(0.01, 0.005); z = -mp.exp(-mp.mpc(0.4, 2.0))
ex = exact_log_inv_poch(z*mp.exp(-h), h)        # w = z q, as for (-q^{n+1};q)_inf with z = -q^n
for s in [-1, 0, 1]:
    approx = sum(mp.bernpoly(k, s)/mp.factorial(k) * mp.polylog(2-k, z) * h**(k-1) for k in range(4))
    print(f'B_m({s}): |error| = {mp.nstr(abs(ex - approx), 3)}')
