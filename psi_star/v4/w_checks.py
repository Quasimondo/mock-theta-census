"""Further checks for W:
 (a) even-root rule: at xi of even order 2s, the geometric-series value of W equals -F_a(xi), F_a = sum_{n>=1} xi^(2n^2-n)/(xi;xi^2)_n
     (ratio -1/2 per s steps), and both equal the radial limit of W;  identity (V) at complex q
 (b) perturbative series of W at xi = e^(-t) (terminating) and its Borel singularity (prediction: -pi^2/30)
 (c) the term psi(q) R_+(q^2) at q = e^(-t): constant and first correction against the large-order data
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from math import gcd
from wval import W
mp.mp.dps = 40
half = mp.mpf(1)/2
print("(a) even-root rule")
worst = 0
for n in range(2, 41, 2):
    s = n//2
    for p in range(1, n):
        if gcd(p, n) != 1: continue
        xi = mp.expjpi(mp.mpf(2*p)/n)
        Fa = 0; P = mp.mpc(1)
        for k in range(1, s+1):
            P *= 1 - xi**(2*k-1)
            Fa += xi**(2*k*k-k)/P
        Fa *= mp.mpf(2)/3
        worst = max(worst, abs(W(p, n) + Fa))
print("    |W_geometric(xi) + F_a(xi)| over all primitive roots of even order <= 40:", mp.nstr(worst, 3))
def Wser(q, dps):
    mp.mp.dps = dps
    s = mp.mpc(0); P = mp.mpc(1); j = 0
    while True:
        if j: P *= 1 - q**(2*j-1)
        t = (-1)**j*q**(j*j+j)*P
        s += t
        if j > 30 and abs(q)**(j*j) < mp.mpf(10)**(-dps): break
        j += 1
    return s
def Faser(q, dps):
    mp.mp.dps = dps
    s = mp.mpc(0); P = mp.mpc(1); k = 1; big = mp.mpf(0)
    while True:
        P *= 1 - q**(2*k-1)
        term = q**(2*k*k-k)/P
        s += term; big = max(big, abs(term))
        if k > 30 and abs(term) < mp.mpf(10)**(-dps)*big: break
        k += 1
    return s
def Rp(q, dps):          # sum q^(2k^2+k)/(q^2;q^2)_k
    mp.mp.dps = dps
    s = mp.mpc(0); P = mp.mpc(1); k = 0; big = mp.mpf(0)
    while True:
        if k: P *= 1 - q**(2*k)
        term = q**(2*k*k+k)/P
        s += term; big = max(big, abs(term))
        if k > 30 and abs(term) < mp.mpf(10)**(-dps)*big: break
        k += 1
    return s
def psi_theta(q, dps):
    mp.mp.dps = dps
    s = mp.mpc(0); n = 0
    while abs(q)**(n*(n+1)//2) > mp.mpf(10)**(-dps-5):
        s += q**(n*(n+1)//2); n += 1
    return s
for (p, n) in ((1, 2), (1, 4), (1, 6), (3, 8), (3, 10)):
    xi = mp.expjpi(mp.mpf(2*p)/n); s = n//2
    for u in (mp.mpf(1)/8, mp.mpf(1)/32):
        t = u/s**2
        q = xi*mp.exp(-t)
        w = Wser(q, 80); fa = Faser(q, 80); pr = psi_theta(q, 80)*Rp(q, 80)
        mp.mp.dps = 40
        print("    root %d/%d, s^2 t = %s: |W + F_a - psi R_+| = %s,  |psi R_+| = %s,  |W - W_geometric(xi)| = %s" % (p, n, mp.nstr(u, 3), mp.nstr(abs(w+fa-pr), 3), mp.nstr(abs(pr), 3), mp.nstr(abs(w - W(p, n)), 3)))
print("(b) perturbative series of W at xi = e^(-t)")
from flint import acb_series, acb, arb, ctx
K = 120
ctx.prec = 3000; ctx.cap = K+1
tot = acb_series([0]); P = acb_series([1]); j = 0
while True:
    if j: P = P*(1 - acb_series([0, -(2*j-1)]).exp())
    term = P*acb_series([0, -(j*j+j)]).exp()
    tot = tot + (term if j % 2 == 0 else -term)
    if j > K + 2: break
    j += 1
mp.mp.dps = 100
a = [mp.mpf(tot[k].real.mid().str(100, radius=False)) for k in range(K+1)]
print("    k! a_k for k = 0..7:", [int(mp.nint(mp.factorial(k)*a[k])) for k in range(8)], " (integers within", mp.nstr(max(abs(mp.factorial(k)*a[k] - mp.nint(mp.factorial(k)*a[k])) for k in range(31)), 3), "for k <= 30)")
est = [-(k+half)*a[k]/a[k+1] for k in range(K)]
r2 = [((k+2)**2*est[k+2] - 2*(k+1)**2*est[k+1] + k**2*est[k])/2 for k in range(len(est)-2)]
print("    singularity estimate A (ratio test, Richardson):", mp.nstr(r2[110], 15), "  pi^2/30 =", mp.nstr(mp.pi**2/30, 15))
A = mp.pi**2/30
rows, rhs = [], []
J = 12
for k in range(60, K+1):
    sc = (-1)**k*mp.gamma(k+half)*A**(-(k+half))
    rows.append([mp.gamma(k+half-j)*A**(-(k+half-j))*(-1)**k/sc for j in range(J)]); rhs.append(a[k]/sc)
sol, res = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
c0 = sol[0]; s1 = -sol[1]/sol[0]
print("    a_k = (-1)^k sum_j c_j Gamma(k+1/2-j) A^-(k+1/2-j):  c0 =", mp.nstr(c0, 20), " residual", mp.nstr(res, 3))
print("    pi c0^2 - 1/(2 sqrt5) =", mp.nstr(mp.pi*c0**2 - 1/(2*mp.sqrt(5)), 3), "   first coefficient of the attached series s1 =", mp.nstr(s1, 18))
print("    PSLQ for s1 in Q(sqrt5):", mp.pslq([s1, 1, mp.sqrt(5)], maxcoeff=5000, maxsteps=100000, tol=mp.mpf(10)**-20))
print("(c) psi(q) R_+(q^2) at q = e^(-t): limit of psi R_+ sqrt(t) e^(-pi^2/(30 t)) and first correction")
vals = []; ts = [mp.mpf(1)/(40*2**i) for i in range(7)]
for t in ts:
    q = mp.exp(-t)
    R = Rp(q, 80)
    mp.mp.dps = 80
    psi = mp.sqrt(mp.pi/(2*t))*mp.exp(t/8)*(1 + 2*sum((-1)**n*mp.exp(-2*mp.pi**2*n*n/t) for n in range(1, 5)))   # psi(e^-t) by the theta transformation
    vals.append(psi*R*mp.sqrt(t)*mp.exp(-mp.pi**2/(30*t)))
def rich(v):
    tab = [v[:]]
    for j in range(1, len(v)):
        p_ = tab[-1]; tab.append([(2**j*p_[i+1]-p_[i])/(2**j-1) for i in range(len(p_)-1)])
    return tab[-1][0]
lim = rich(vals); corr = rich([(vals[i]/lim - 1)/ts[i] for i in range(len(ts))])
print("    limit =", mp.nstr(lim, 18), "  pi*c0 =", mp.nstr(mp.pi*c0, 18), "  ratio - 1 =", mp.nstr(lim/(mp.pi*c0) - 1, 3))
print("    first correction =", mp.nstr(corr, 15), "  s1 from large order =", mp.nstr(s1, 15))
