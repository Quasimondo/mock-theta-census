# Independent check of addendum-v18 sec. 2.2: F(z; e^{2 pi i/X}) / prediction -> 1, with
#  prediction = (1+W) e^{w/2} (pi/(t(1+2W)))^{1/2} exp((V_k(w) + pi^2/12)/t), t = -2 pi i/X,
#  W the root of W^2 + W = 1/z with Im W > 0, w = Log W, k from the critical-point equation,
#  V_k(w) = -w^2/2 - 2 pi i (a+k) w + Li2(-e^w).  Exact F: terminating sum at even X (own code).
import mpmath as mp
def Fexact(a, X):
    with mp.workdps(40 + int(0.08*X)):
        q = mp.expjpi(mp.mpf(2)/X); z = mp.expjpi(2*a); s = 0; P = mp.mpc(1)
        for n in range(X//2 + 1):
            s += z**n * q**(n*(n+1)//2) * P; P *= 1 + q**(n+1)
        return +s
def pred(a, X):
    with mp.workdps(40 + int(0.08*X)):
        z = mp.expjpi(2*a); rts = [(-1 + sg*mp.sqrt(1 + 4/z))/2 for sg in (1, -1)]
        W = max(rts, key=lambda r: mp.im(r)); w = mp.log(W)
        kf = (-w - mp.log(1 + W))/(2j*mp.pi) - a; k = mp.nint(mp.re(kf))
        assert abs(kf - k) < 1e-20, kf
        t = -2j*mp.pi/X
        V = -w**2/2 - 2j*mp.pi*(a + k)*w + mp.polylog(2, -mp.exp(w))
        return +((1 + W)*mp.exp(w/2)*mp.sqrt(mp.pi/(t*(1 + 2*W)))*mp.exp((V + mp.pi**2/12)/t)), int(k), W
mp.mp.dps = 30
for a, name in [(mp.mpf(1)/4, 'z = i'), (mp.mpf(1)/7, 'z = e^{2 pi i/7}'), (mp.mpf(5)/6, 'z = -omega'), (1/mp.sqrt(2), 'z = e^{2 pi i/sqrt2}')]:
    Xs = (802, 1602, 3202); rs = []
    for X in Xs:
        P, k, W = pred(a, X); rs.append(Fexact(a, X)/P)
    xs = [mp.mpf(1)/X for X in Xs]
    # Neville in 1/X
    Pn = list(rs)
    for j in range(1, 3):
        for i in range(3 - j): Pn[i] = (xs[i+j]*Pn[i] - xs[i]*Pn[i+1])/(xs[i+j] - xs[i])
    print(f'{name:22s} W = {mp.nstr(W, 6):>20} k = {k:+d}   |ratio - 1| at X = 802, 1602, 3202: {", ".join(mp.nstr(abs(r-1), 3) for r in rs)}   extrapolated: {mp.nstr(abs(Pn[0]-1), 3)}', flush=True)
