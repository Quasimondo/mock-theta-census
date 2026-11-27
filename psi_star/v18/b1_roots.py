"""Addendum v18, task 2: F(z) at q = e^(2 pi i/X), X even (the sum is finite, n < X/2), against the saddle-point formula
     F ~ (1 + W) e^(w/2) sqrt(pi/(t (1 + 2W))) exp( (V(w) + pi^2/12)/t ),   t = -2 pi i/X,
   with the RULE:  W = the root of W^2 + W = 1/z in the upper half-plane,  w = Log W (principal, 0 < Im w < pi),
     V(w) = -w^2/2 - 2 pi i (a + k) w + Li2(-W),  k the integer with  w + 2 pi i (a + k) + Log(1 + W) = 0.
   Height: -Im V = -D(-W) = Vol(z) (Bloch-Wigner), the same for every branch k and with the opposite sign for the other root.
Usage: python3 b1_roots.py [digits]"""
import sys
from mpmath import mp, mpf, mpc, exp, pi, sqrt, log, polylog, nstr, arg, im, re
I = mpc(0, 1)

def F_root(a, X):
    # terms reach e^(0.1615 X) while the sum can be much smaller: 0.0702 X extra digits
    old = mp.dps; mp.dps = old + int(0.0702*X) + 10
    try:
        z = exp(2*pi*I*a); q = exp(2*pi*I/X); s = mpc(0); term = mpc(1); qn = mpc(1)
        for n in range(X//2):
            s += term; qn *= q; term *= qn*z*(1 + qn)
        return +s
    finally:
        mp.dps = old

def bloch_wigner(y):
    return im(polylog(2, y)) + arg(1 - y)*log(abs(y))

def rule(a):
    z = exp(2*pi*I*a); d = sqrt(1 + 4/z)
    W = (-1 + d)/2
    if im(W) < 0: W = (-1 - d)/2
    w = log(W)
    kk = (-w - log(1 + W))/(2*pi*I) - a
    k = int(round(float(re(kk))))
    V = -w*w/2 - 2*pi*I*(a + k)*w + polylog(2, -W)
    return W, w, k, abs(kk - k), V

def predict(W, w, V, X):
    t = -2*pi*I/X
    return (1 + W)*exp(w/2)*sqrt(pi/(t*(1 + 2*W)))*exp((V + pi**2/12)/t)

if __name__ == "__main__":
    mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 45
    cases = [("-1", mpf(1)/2), ("i", mpf(1)/4), ("omega", mpf(1)/3), ("-omega", mpf(5)/6), ("-i", mpf(3)/4), ("e^(2 pi i/7)", mpf(1)/7),
             ("e^(2 pi i 0.1)", mpf(1)/10), ("e^(2 pi i 0.9)", mpf(9)/10), ("e^(2 pi i 0.03)", mpf(3)/100), ("e^(2 pi i sqrt(2)/2)", sqrt(mpf(2))/2)]
    Xs = [802 + 800*j for j in range(12)]
    print("X =", Xs[0], "...", Xs[-1], "(12 values); ratio = F/prediction, extrapolated polynomially in 1/X")
    for name, a in cases:
        W, w, k, defect, V = rule(a)
        Fs = [F_root(a, X) for X in Xs]
        ratios = [F/predict(W, w, V, X) for F, X in zip(Fs, Xs)]
        xs = [mpf(1)/X for X in Xs]; tab = ratios[:]
        for m in range(1, len(tab)):
            tab = [(xs[i]*tab[i + 1] - xs[i + m]*tab[i])/(xs[i] - xs[i + m]) for i in range(len(tab) - 1)]
        vol = -bloch_wigner(-W)
        print(f"z = {name}: W = {nstr(W, 12)}, Im w/pi = {nstr(im(w)/pi, 6)}, k = {k} (defect {nstr(defect, 2)}), a + k = {nstr(a + k, 6)}")
        print(f"    Vol = -D(-W) = {nstr(vol, 18)},  Vol/(2 pi) = {nstr(vol/(2*pi), 12)},  -Im V = {nstr(-im(V), 18)},  Re V + pi^2/12 = {nstr(re(V) + pi**2/12, 15)}")
        print(f"    ratio at X = {Xs[0]}: {nstr(ratios[0], 8)};  at X = {Xs[-1]}: {nstr(ratios[-1], 8)};  |extrapolated ratio - 1| = {nstr(abs(tab[0] - 1), 3)}")
    worst = mpf(0)
    for name, a in cases:
        z = exp(2*pi*I*a); d = sqrt(1 + 4/z)
        for W in ((-1 + d)/2, (-1 - d)/2):
            for k in range(-3, 4):
                w = -log(1 + W) - 2*pi*I*(a + k)
                V = -w*w/2 - 2*pi*I*(a + k)*w + polylog(2, -W)
                worst = max(worst, abs(im(V) - bloch_wigner(-W)), abs(exp(w) - W))
    print("identity Im V_k(w_k) = D(-W) and e^(w_k) = W, both roots, k = -3..3, all ten z: largest deviation " + nstr(worst, 3))
    print("near z = 1: Vol(e^(2 pi i a))/a for a = 0.03, 0.01, 0.001, 0.0001  (2 pi log((1+sqrt 5)/2) = " + nstr(2*pi*log((1 + sqrt(mpf(5)))/2), 12) + ")")
    print("    " + ", ".join(nstr(-bloch_wigner(-rule(mpf(s))[0])/mpf(s), 10) for s in ("0.03", "0.01", "0.001", "0.0001")))
