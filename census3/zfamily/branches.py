import mpmath as mp
mp.mp.dps = 30
meas = {mp.mpf(1)/2: 0.16196603, mp.mpf(1)/4: 0.11131682, mp.mpf(5)/6: 0.077880743, mp.mpf(1)/3: 0.13812295}
def V(a, x): return 2j*mp.pi*a*x - x**2/2 + mp.polylog(2, -mp.exp(-x)) + mp.pi**2/12
for a, m in meas.items():
    cands = []
    for sgn in (1, -1):
        w = (-1 + sgn*mp.sqrt(1 + 4*mp.expjpi(-2*a)))/2
        for k in range(-3, 4):
            x = -mp.log(w) + 2j*mp.pi*k
            # check stationarity with the branch-consistent derivative: f'(x) = 2 pi i a - x + log(1 + e^{-x}) (+ 2 pi i j)
            v = V(a, x)
            cands.append((abs(abs(mp.im(v))/(2*mp.pi) - m), sgn, k, v))
    cands.sort(key=lambda c: c[0])
    d, sgn, k, v = cands[0]
    print(f'a = {mp.nstr(a,4):6s} measured {m:.8f}  best saddle: root {sgn:+d}, branch k = {k:+d}, V = {mp.nstr(v, 10)},  |Im V|/2pi = {mp.nstr(abs(mp.im(v))/(2*mp.pi), 8)}  (diff {mp.nstr(d, 2)})')
