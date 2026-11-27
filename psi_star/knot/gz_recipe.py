# Garoufalidis-Zagier sec. 6.3 recipe applied to psi1*(q) = sum_n (-1)^n q^{n(n+1)/2} (-q;q)_n, q = e^{-h}:
#  (-q;q)_n = (-q;q)_inf / (-q^{n+1};q)_inf,   1/(z e^{s h};q)_inf ~ exp(sum_m B_m(s)/m! Li_{2-m}(z) h^{m-1})   (their eq. 6.8)
#  here z = -e^{-x}, x = n h, s = 0 (B_m(0); see check_term2.py);  (-q;q)_inf = 2^{-1/2} e^{pi^2/(12h) + h/24} (exact up to exponentially small terms)
#  log T(n) = f(x)/h + l(x) + h m(x) + O(h^2), sum -> (1/h) int dx (Poisson, k = 0 term), Laplace expansion at the stationary point.
# Compare the first correction with addendum-v3/v5: psi1* ~ C t^{-1/2} e^{-Sbar/t} W (1 + (1/4 + i/(18 sqrt3)) t), t = h.
import mpmath as mp
mp.mp.dps = 30
B1, B2 = mp.bernpoly(1, 0), mp.bernpoly(2, 0)            # -1/2, 1/6  (checked numerically: check_term2.py)
f  = lambda x: 1j*mp.pi*x - x**2/2 + mp.polylog(2, -mp.exp(-x)) + mp.pi**2/12
f1 = lambda x: 1j*mp.pi - x + mp.log(1 + mp.exp(-x))
f2 = lambda x: -1 - mp.exp(-x)/(1 + mp.exp(-x))
f3 = lambda x: mp.exp(-x)/(1 + mp.exp(-x))**2
f4 = lambda x: mp.exp(x)*(1 - mp.exp(x))/(mp.exp(x) + 1)**3
Li1 = lambda z: -mp.log(1 - z); Li0 = lambda z: z/(1 - z)
l  = lambda x: -x/2 + B1*Li1(-mp.exp(-x)) - mp.log(2)/2          # O(1): from q^{n(n+1)/2}, the B_1 term and 2^{-1/2}
L1 = lambda x: -mp.mpf(1)/2 + B1 * mp.exp(-x)/(1 + mp.exp(-x))   # l'(x)   [d/dx Li1(-e^{-x}) = e^{-x}/(1+e^{-x})]
L2 = lambda x: -B1 * mp.exp(-x)/(1 + mp.exp(-x))**2              # l''(x)
m  = lambda x: B2/2 * Li0(-mp.exp(-x)) + mp.mpf(1)/24          # O(h): B_2 term and the h/24 of (-q;q)_inf
Sbar = mp.pi**2/3 + 1j*mp.clsin(2, mp.pi/3)
for x0 in [2j*mp.pi/3, 4j*mp.pi/3]:
    x = mp.findroot(f1, x0)
    F2, F3, F4 = f2(x), f3(x), f4(x)
    G0 = 1; G1 = L1(x); G2 = L2(x) + L1(x)**2                    # G = e^{l}, normalised by G0
    K = m(x) + (-G2/(2*F2) + G1*F3/(2*F2**2) + F4/(8*F2**2) - 5*F3**2/(24*F2**3))
    print(f'stationary point x = {mp.nstr(x, 8)}:  f(x) = {mp.nstr(f(x), 12)}  (-Sbar = {mp.nstr(-Sbar, 12)})')
    print(f'   first correction = {mp.nstr(K, 15)}    measured (1/4 + i/(18 sqrt3)) = {mp.nstr(mp.mpf(1)/4 + 1j/(18*mp.sqrt(3)), 15)}')
