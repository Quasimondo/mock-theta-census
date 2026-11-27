# Independent check of addendum-v9 sec. 1:  (J(r) - Phi_J(tau)) sqrt(tau) -> c [ q'^{-1/120}(F0*(q') - 1) + i q'^{71/120} F1*(q') ]
# along Y = Y0 + 960k - i b (q' and the phases fixed), r = e^{2 pi i/Y}, tau = -2 pi i/Y, q' = e^{-pi i Y/2}; Neville in tau.
import mpmath as mp
from phiJ import PHI
def J(r, dps):
    s = 0; P = mp.mpc(1); m = 0
    while True:
        term = r**(2*m*m+m) / P
        s += term
        if m > 10 and abs(term) < mp.mpf(10)**(-dps+15) * abs(s): return s
        m += 1; P *= 1 + r**(2*m)
def F0s(q):
    s = 0; P = mp.mpc(1)
    for n in range(60):
        s += (-1)**n * q**(2*n*n) / P; P *= 1 - q**(2*n+1)
    return s
def F1s(q):
    s = 0; P = 1 - q
    for n in range(60):
        s += (-1)**n * q**(2*n*n+2*n) / P; P *= 1 - q**(2*n+3)
    return s
def neville(xs, ys):
    P = list(ys)
    for k in range(1, len(P)):
        for i in range(len(P)-k): P[i] = (xs[i+k]*P[i] - xs[i]*P[i+1])/(xs[i+k]-xs[i])
    return P[0]
for Y0, b in [(3, mp.mpf('0.8')), (1.5, mp.mpf('0.5'))]:
    xs, ys = [], []
    for k in (1, 2, 3, 4):
        dps = int(60 + (Y0 + 960*k)/8)
        with mp.workdps(dps):
            Y = mp.mpc(Y0 + 960*k, -b); r = mp.exp(2j*mp.pi/Y); tau = -2j*mp.pi/Y
            terms = [mp.mpf(PHI[n].numerator)/PHI[n].denominator * tau**n for n in range(len(PHI))]
            nopt = min(range(1, len(terms)), key=lambda n: abs(terms[n]))
            val = (J(r, dps) - sum(terms[:nopt])) * mp.sqrt(tau)
            xs.append(tau); ys.append(val)
            print(f'  Y0={Y0} k={k}: smallest Phi term {mp.nstr(abs(terms[nopt]),2)} at n={nopt}', flush=True)
    with mp.workdps(40):
        Y = mp.mpc(Y0, -b); qp = mp.exp(-1j*mp.pi*Y/2)
        c = mp.sqrt(2*mp.pi/mp.sqrt(5))
        pred = c * (mp.exp(1j*mp.pi*Y/240)*(F0s(qp) - 1) + 1j*mp.exp(-71j*mp.pi*Y/240)*F1s(qp))
        lim = neville(xs, ys)
        print(f'Y0={Y0} b={mp.nstr(b,2)}: limit {mp.nstr(lim,12)}  prediction {mp.nstr(pred,12)}  |diff|/|pred| = {mp.nstr(abs(lim-pred)/abs(pred),3)}', flush=True)
