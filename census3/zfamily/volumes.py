# Is |Im V| at the selected saddle the Bloch-Wigner dilogarithm (volume of an ideal tetrahedron) of the saddle shape?
import mpmath as mp
mp.mp.dps = 30
def D(x): return mp.im(mp.polylog(2, x)) + mp.arg(1 - x) * mp.log(abs(x))
sel = {mp.mpf(1)/2: (-1, -3), mp.mpf(1)/4: (1, 0), mp.mpf(5)/6: (-1, 0), mp.mpf(1)/3: (-1, 1), mp.mpf(0): (1, 0)}
for a, (sgn, k) in sel.items():
    w = (-1 + sgn*mp.sqrt(1 + 4*mp.expjpi(-2*a)))/2
    x = -mp.log(w) + 2j*mp.pi*k
    V = 2j*mp.pi*a*x - x**2/2 + mp.polylog(2, -mp.exp(-x)) + mp.pi**2/12
    cands = {'-w': -w, '1+w': 1+w, '-1/w': -1/w, 'w/(1+w)': w/(1+w), '1/(1+w)': 1/(1+w), 'w': w}
    best = min(cands.items(), key=lambda kv: abs(abs(D(kv[1])) - abs(mp.im(V))))
    print(f'a = {mp.nstr(a,4):6s} |Im V| = {mp.nstr(abs(mp.im(V)), 15)}   |D({best[0]})| = {mp.nstr(abs(D(best[1])), 15)}   shape {best[0]} = {mp.nstr(best[1], 8)}')
print('regular ideal tetrahedron volume D(e^{i pi/3}) =', mp.nstr(D(mp.expjpi(mp.mpf(1)/3)), 15))
