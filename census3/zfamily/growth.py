# F(z) = sum z^n q^{n(n+1)/2} (-q;q)_n at q = e^{2 pi i/X}, z = e^{2 pi i a}: exact values
# (even order: the sum terminates; odd order m: A_z/(1 - 2 z^m), A_z = first m terms; addendum-v3 sec. 2).
# Prediction (saddle of addendum-v3 sec. 2): |F| ~ exp(|Im V(a)| X/(2 pi)), V = value at w^2 + w = e^{-2 pi i a}.
import mpmath as mp
mp.mp.dps = 30
def Fz(a, X):
    with mp.workdps(30 + X // 8):
        q = mp.expjpi(mp.mpf(2)/X); z = mp.expjpi(2*a); s = 0; P = mp.mpc(1); terms = []
        for n in range(X):
            terms.append(z**n * q**(n*(n+1)//2) * P); P *= 1 + q**(n+1)
        A = mp.fsum(terms)
        return +(A if X % 2 == 0 else A / (1 - 2*z**X))
def Vsaddles(a):
    out = []
    for sgn in (1, -1):
        c = mp.expjpi(-2*a); w = (-1 + sgn*mp.sqrt(1 + 4*c))/2; x = -mp.log(w)
        out.append(2j*mp.pi*a*x - x**2/2 + mp.polylog(2, -mp.exp(-x)) + mp.pi**2/12)
    return out
for a, name in [(mp.mpf(1)/2, 'z = -1 (psi1*)'), (mp.mpf(1)/4, 'z = i'), (mp.mpf(5)/6, 'z = -omega'), (mp.mpf(1)/3, 'z = omega'), (mp.mpf(0), 'z = 1 (psi1, mock)')]:
    Vs = Vsaddles(a)
    pred = max(abs(mp.im(V)) for V in Vs) / (2*mp.pi)
    rates = []
    for X in (400, 800, 1600):
        v = Fz(a, X); rates.append((X, mp.log(abs(v))/X))
    meas = (1600*rates[2][1] - 800*rates[1][1]) / 800        # slope of log|F| between X = 800 and 1600
    print(f'{name:22s} predicted rate {mp.nstr(pred, 8):>11}   measured slope {mp.nstr(meas, 8):>11}   (log|F|/X: {", ".join(mp.nstr(r, 5) for _, r in rates)})', flush=True)
