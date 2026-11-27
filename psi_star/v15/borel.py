"""Borel-Pade-Laplace for the series Phi (exact coefficients from phi_exact.json), and the kernel integral P J on both sides of the real x-axis."""
import json, os
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, quad, pade, polyval, polyroots, factorial, nstr, arg
I = mpc(0, 1)
HERE = os.path.dirname(os.path.abspath(__file__))
_EX = json.load(open(os.path.join(HERE, "phi_exact.json")))
_CACHE = {}

def coeff(case, k):
    n, d = _EX[case][k].split("/"); return mpf(int(n))/mpf(int(d))

def pade_borel(case, n, dps=None):
    """Pade [n/n] of the Borel transform sum a_k xi^k/k!; returns (P, Q) highest degree first."""
    key = (case, n)
    if key not in _CACHE:
        old = mp.dps; mp.dps = dps or max(150, 3*n)
        c = [coeff(case, k)/factorial(k) for k in range(2*n + 1)]
        P, Q = pade(c, n, n)
        _CACHE[key] = (P[::-1], Q[::-1]); mp.dps = old
    return _CACHE[key]

def borel_sum(case, t, theta, n):
    """S_theta[Phi](t) = (e^(i theta)/t) int_0^inf e^(-r e^(i theta)/t) B(r e^(i theta)) dr."""
    P, Q = pade_borel(case, n); e = exp(I*theta); s = abs(t)
    f = lambda r: exp(-r*e/t)*polyval(P, r*e)/polyval(Q, r*e)
    return e/t*quad(f, [0] + [s*k for k in (2, 5, 10, 20, 40, 80, 160, 320)])

def qpg(a, q):
    """(a;q)_inf for |q| < 1, and its continuation 1/(a/q; 1/q)_inf for |q| > 1."""
    return qp(a, q, maxterms=10**6) if abs(q) < 1 else 1/qp(a/q, 1/q, maxterms=10**6)

def PJ(x, delta=mpf(2)/5, case="psi1", smax=None):
    """P J(1/(2b)) (psi1) or P q^(1/2) J(b + 1/(2b)) (psi0), with P = 2^(-1/2) e^(pi^2/(12t) + t/6)/b  (t = -2 pi i x),
    which is (-q;q)_inf/(b q^(1/8) (p;p^2)_inf) for Im x > 0.  Integration variable s = y/b on the line e^(-i delta) R / b."""
    x = mpc(x); b = sqrt(x); t = -2*pi*I*x; q = exp(-t); p = exp(-pi*I/x); qh = exp(-t/2)
    def integrand(s):
        y = b*s
        Phi = qpg(-qh*exp(2*pi*x*s), q)/qpg(-p*exp(2*pi*s), p*p)
        extra = exp(2*pi*x*s) if case == "psi0" else 1
        return exp(-pi*I*x*s*s + pi*s)*extra/(Phi*(1 + exp(2*pi*s)))
    e = exp(-I*delta)/b; e = e/abs(e)                      # direction of L/b
    R = smax or (sqrt(2*mp.dps*mpf("2.31")/(abs(t)*abs((e*e*t/abs(t)).real))) + 5)
    pts = [-R] + [-R*mpf(k)/8 for k in range(7, 0, -1)] + [-3, -1.5, -0.7, 0, 0.7, 1.5, 3] + [R*mpf(k)/8 for k in range(1, 8)] + [R]
    pts = sorted(set(pts))
    Jb = quad(lambda r: integrand(r*e)*e, pts)             # = J/b
    pref = exp(pi**2/(12*t) + t/6)/sqrt(2)
    if case == "psi0": pref *= exp(-t/2)
    return pref*Jb
