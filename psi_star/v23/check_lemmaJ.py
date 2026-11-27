# Independent check of addendum-v23 (own code, no fad.py).
# (1) Lemma U: |Phi_b(y)| = 1 for real b, y, with Faddeev's integral computed here (trapezoid on Im xi = c, numpy).
# (2) Lemma J: P J(a/b) at x = 1/X over the real line, |P J| against 1/(2 sqrt2 sin(pi a)); compare with
#     Fable's values (run_log.txt) and with 1/(1 - 2z).
# (3) Theorem 15: X (F/prediction - 1) from the exact finite sums, Richardson-extrapolated.
import numpy as np, mpmath as mp
def logPhi(y, b, n=20001):
    c = np.pi*min(b, 1/b)/2
    L = 40/(b + 1/b); r = np.linspace(-L, L, n); xi = r + 1j*c; h = r[1] - r[0]
    y = np.atleast_1d(y)[:, None]
    f = np.exp(-2j*y*xi)/(4*xi*np.sinh(b*xi)*np.sinh(xi/b))
    return f.sum(axis=1)*h
print('(1) Lemma U: max |Re log Phi_b(y)| over 50 real y in [-3, 3]')
for X in (8, 40, 200):
    b = 1/np.sqrt(X); print(f'   X = {X}: {np.abs(logPhi(np.linspace(-3, 3, 50), b).real).max():.1e}')
print('(2) Lemma J at x = 1/X: |P J(a/b)| over the real line  (bound; Fable; |1/(1-2z)|)')
fable = {(0.25, 8): 0.468003, (0.25, 40): 0.470865, (0.25, 200): 0.452311, (0.5, 8): 0.326445, (0.5, 40): 0.333029, (0.5, 200): 0.333321, (0.7, 8): 0.345924, (0.05, 40): 1.818911, (0.95, 40): 0.731923}
for (a, X), fv in fable.items():
    b = 1/np.sqrt(X); t = -2j*np.pi/X
    S = 14/min(a, 1 - a); s = np.linspace(-S, S, 8001); h = s[1] - s[0]
    lp = np.concatenate([logPhi(s[i:i+500]*b, b) for i in range(0, len(s), 500)])
    ker = np.where(s < 0, np.exp(2*np.pi*a*np.minimum(s, 0))/(1 + np.exp(2*np.pi*np.minimum(s, 0))), np.exp(2*np.pi*(a - 1)*np.maximum(s, 0))/(1 + np.exp(-2*np.pi*np.maximum(s, 0))))
    g = np.exp(-np.pi*1j*(b*s)**2 - lp)*ker*b
    J = g.sum()*h; P = 2**-0.5*np.exp(np.pi**2/(12*t) + t/6)/b; z = np.exp(2j*np.pi*a)
    print(f'   a = {a:<5} X = {X:3d}: {abs(P*J):.6f}  (bound {1/(2*np.sqrt(2)*np.sin(np.pi*a)):.6f}; Fable {fv}; |1/(1-2z)| = {abs(1/(1-2*z)):.5f})', flush=True)
print('(3) Theorem 15: X (F/prediction - 1), extrapolated')
for a in (mp.mpf(1)/4, mp.mpf(1)/3, mp.mpf(5)/6, mp.mpf(1)/20):
    rows = []
    for X in (400, 800, 1600, 3200):
        with mp.workdps(30 + X//10):
            z = mp.expj(2*mp.pi*a); q = mp.expjpi(mp.mpf(2)/X); t = -2j*mp.pi/X
            W = [r for r in mp.polyroots([1, 1, -1/z]) if mp.im(r) > 0][0]; w0 = mp.log(W)
            V = -w0**2/2 - 2j*mp.pi*(a - 1)*w0 + mp.polylog(2, -W)
            pred = (1 + W)*mp.exp(w0/2)*mp.sqrt(mp.pi/(t*(1 + 2*W)))*mp.exp((V + mp.pi**2/12)/t)
            s = 0; P_ = mp.mpc(1)
            for n in range(X//2):
                s += z**n*q**(n*(n+1)//2)*P_; P_ *= 1 + q**(n+1)
            rows.append(X*(s/pred - 1))
    xs = [mp.mpf(1)/X for X in (400, 800, 1600, 3200)]; v = rows[:]
    for j in range(1, 4):
        for i in range(4 - j): v[i] = (xs[i+j]*v[i] - xs[i]*v[i+1])/(xs[i+j] - xs[i])
    print(f'   a = {mp.nstr(a,4)}: X=3200: {mp.nstr(rows[-1], 8)}   extrapolated: {mp.nstr(v[0], 10)}', flush=True)
