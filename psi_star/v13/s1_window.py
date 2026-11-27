"""Addendum v13, step 1: the perturbative series from the 1/cosh window of the kernel integral, in exact rational arithmetic.

With y = b s, q = e^(-t), sigma = i s, and the p-product dropped, the integrand of P J(1/(2b)) is
      e^(-t sigma^2/2 + t/8) / (2 cosh(pi s)) * (-q;q)_inf / (-q^(1/2 - sigma); q)_inf .
Euler-Maclaurin:  log (-q^a;q)_inf ~ pi^2/(12 t) - sum_{k>=1} B_k(a) g^(k-1)(0) t^(k-1)/k!,   g(u) = log(1 + e^(-u)),
so the ratio of products is  2^(-1/2 - sigma) exp( - sum_{k>=2} [B_k(1) - B_k(1/2 - sigma)] g^(k-1)(0) t^(k-1)/k! ),
a power series in t whose coefficients are polynomials in sigma, and
      int sigma^j 2^(-1/2 - i s) / (2 cosh(pi s)) ds = m_j := j! [u^j] e^(u/2)/(2 + e^u)        (m_0 = 1/3).
Output: the coefficients a_0..a_N of the resulting series (exact rationals), for psi1* (u = 1/(2b)) and psi0* (u = b + 1/(2b),
extra factor e^(t(sigma - 1/2))), compared with the 151 coefficients of Phi from v3 (q2_series_v3.json, 120 digits);
written to phi_exact.json.          Usage: python3 s1_window.py [N]   (default 150)"""
import sys, json, time
from fractions import Fraction
from flint import fmpq, fmpq_poly
N = int(sys.argv[1]) if len(sys.argv) > 1 else 150

def fact(n):
    r = 1
    for i in range(2, n + 1): r *= i
    return r

def series_div(num, den, n):        # coefficient lists of fmpq, den[0] != 0
    out = []
    for k in range(n + 1):
        s = num[k] if k < len(num) else fmpq(0)
        for j in range(1, k + 1):
            if j < len(den): s -= den[j]*out[k - j]
        out.append(s/den[0])
    return out

M = 2*N + 2
invfact = [fmpq(1, fact(k)) for k in range(M + 3)]
# Bernoulli numbers: t/(e^t - 1)
bern = [c*fact(k) for k, c in enumerate(series_div([fmpq(1)], [invfact[k + 1] for k in range(N + 3)], N + 2))]
# g^(m)(0), m >= 1:  g'(u) = -1/(1 + e^u)
inv1 = series_div([fmpq(1)], [fmpq(2)] + invfact[1:N + 3], N + 2)
gder = [None] + [-inv1[m - 1]*fact(m - 1) for m in range(1, N + 2)]
# moments m_j = j! [u^j] e^(u/2)/(2 + e^u)
h = series_div([invfact[k]/fmpq(2)**k for k in range(M + 1)], [fmpq(3)] + invfact[1:M + 1], M)
mom = [h[j]*fact(j) for j in range(M + 1)]
# P_k(sigma) = B_k(1/2 - sigma):  P_k' = -k P_(k-1),  P_k(0) = B_k(1/2) = (2^(1-k) - 1) B_k
P = [fmpq_poly([1])]
for k in range(1, N + 2):
    P.append(fmpq_poly([(fmpq(2)/fmpq(2)**k - 1)*bern[k]]) - k*P[-1].integral())
Bk1 = [bern[k] if k != 1 else fmpq(1, 2) for k in range(N + 2)]

def series(case):
    F = [None]*(N + 1)
    for m in range(1, N + 1):
        F[m] = (P[m + 1] - fmpq_poly([Bk1[m + 1]]))*(gder[m]*invfact[m + 1])
    F[1] = F[1] + fmpq_poly([fmpq(1, 8), 0, fmpq(-1, 2)])
    if case == "psi0": F[1] = F[1] + fmpq_poly([fmpq(-1, 2), 1])
    E = [fmpq_poly([1])]
    for m in range(1, N + 1):
        acc = fmpq_poly([0])
        for j in range(1, m + 1): acc = acc + (j*F[j])*E[m - j]
        E.append(acc*fmpq(1, m))
    a = []
    for m in range(N + 1):
        c = E[m].coeffs(); a.append(sum((c[j]*mom[j] for j in range(len(c))), fmpq(0)))
    return a

if __name__ == "__main__":
    ref = json.load(open("q2_series_v3.json")); out = {}
    for case in ("psi1", "psi0"):
        t0 = time.time(); a = series(case); out[case] = [f"{int(c.p)}/{int(c.q)}" for c in a]
        print(f"{case}*: a_0..a_5 =", ", ".join(out[case][:6]), f"  ({time.time() - t0:.0f} s)")
        from mpmath import mp, mpf
        mp.dps = 130; worst = 0; nref = len(ref[case + "@1"])
        for k in range(min(N + 1, nref)):
            x = mpf(int(a[k].p))/mpf(int(a[k].q)); r = mpf(ref[case + "@1"][k][0])
            d = abs(x - r)/max(abs(x), mpf(10)**(-100)) if a[k] != 0 else abs(r)
            worst = max(worst, d)
        nz = [k for k in range(N + 1) if a[k] == 0]
        print(f"   compared with v3: {min(N + 1, nref)} coefficients, worst relative deviation {float(worst):.2e};  exact zeros at k = {nz[:8]}")
        from mpmath import nstr
        print(f"   a_{N} = {nstr(mpf(int(a[N].p))/mpf(int(a[N].q)), 7)},  denominator of a_{N} has {len(str(int(a[N].q)))} digits")
        sys.stdout.flush()
    json.dump({"N": N, **out}, open("phi_exact.json", "w"))
