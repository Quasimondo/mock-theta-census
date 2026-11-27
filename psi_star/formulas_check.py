"""Numerically checks every equivalent form listed in docs/formulas.md.
Each check evaluates both sides at a few points |q| < 1 (and x, z where relevant) at 50 digits.
Run:  python3 formulas_check.py      (prints one line per identity; all should say OK)"""
import mpmath as mp
mp.mp.dps = 50
N = 90                                    # terms; q^(n^2/2) is far below 1e-50 for |q| <= 0.65

def poch(a, q, n):                        # (a;q)_n
    p = mp.mpf(1)
    for k in range(n):
        p *= (1 - a * q ** k)
    return p
def poch_inf(a, q):
    p = mp.mpf(1); k = 0
    while True:
        t = a * q ** k
        if abs(t) < mp.mpf(10) ** -60: return p
        p *= (1 - t); k += 1
def S(term, n0=0, n1=N):
    return mp.fsum(term(n) for n in range(n0, n1))
def binom2(n): return n * (n - 1) // 2

# --- the objects -------------------------------------------------------------
psi1  = lambda q: S(lambda n: q ** (n * (n + 1) // 2) * poch(-q, q, n))
psi0  = lambda q: S(lambda n: q ** ((n + 1) * (n + 2) // 2) * poch(-q, q, n))
psi1s = lambda q: S(lambda n: (-1) ** n * q ** (n * (n + 1) // 2) * poch(-q, q, n))
psi0s = lambda q: S(lambda n: (-1) ** n * q ** ((n + 1) * (n + 2) // 2) * poch(-q, q, n))
f0s   = lambda q: S(lambda n: (-1) ** n * q ** (n * n) / poch(-q, q, n))
f1s   = lambda q: S(lambda n: (-1) ** n * q ** (n * n + n) / poch(-q, q, n))
R0    = lambda q: S(lambda n: (-1) ** n * q ** (n * n) / poch(q, q, n))
R1    = lambda q: S(lambda n: (-1) ** n * q ** (n * n + n) / poch(q, q, n))
th4   = lambda q: 1 + 2 * S(lambda n: (-1) ** n * q ** (n * n), 1)
Aq    = lambda a, q: S(lambda n: (-a) ** n * q ** (n * n) / poch(q, q, n))      # Ramanujan's entire function
F     = lambda z, q: S(lambda n: z ** n * q ** (n * (n + 1) // 2) * poch(-q, q, n))
def phi_rs(a, b, q, z, nmax=N):           # basic hypergeometric r phi s (Gasper-Rahman convention)
    r, s = len(a), len(b)
    def t(n):
        num = mp.mpf(1); den = poch(q, q, n)
        for x in a: num *= poch(x, q, n)
        for y in b: den *= poch(y, q, n)
        return num / den * ((-1) ** n * q ** binom2(n)) ** (1 + s - r) * z ** n
    return S(t, 0, nmax)
H  = lambda x, q: S(lambda n: q ** binom2(n) * poch(-q, q, n) / (poch(x * q, q, n) * poch(q / x, q, n)))
g2 = lambda x, q: S(lambda n: poch(-q, q, n) * q ** (n * (n + 1) // 2) / (poch(x, q, n + 1) * poch(q / x, q, n + 1)))

checks = []
def check(name, lhs, rhs, pts):
    err = max(abs(lhs(*p) - rhs(*p)) for p in pts)
    checks.append((name, err))
Q = [(mp.mpf('0.3'),), (mp.mpc('-0.45', '0.2'),), (mp.mpc('0.1', '0.6'),)]

# psi1* rewrites
check("psi1*: (q^2;q^2)_n/(q;q)_n form", psi1s,
      lambda q: S(lambda n: (-1) ** n * q ** (n * (n + 1) // 2) * poch(q * q, q * q, n) / poch(q, q, n)), Q)
check("psi1*: binomial exponent q^binom(n+1,2)", psi1s,
      lambda q: S(lambda n: (-1) ** n * q ** ((n + 1) * n // 2) * mp.fprod([1 + q ** k for k in range(1, n + 1)])), Q)
check("psi1* = F(-1)", psi1s, lambda q: F(-1, q), Q)
check("psi1* = 2phi2(q, -q; 0, 0; q, q)", psi1s, lambda q: phi_rs([q, -q], [0, 0], q, q), Q)
# psi0* rewrites
check("psi0* = q*sum (-1)^n q^{n(n+3)/2}(-q;q)_n", psi0s,
      lambda q: q * S(lambda n: (-1) ** n * q ** (n * (n + 3) // 2) * poch(-q, q, n)), Q)
check("psi0* = -sum_{n>=1} (-1)^n q^{n(n+1)/2}(-q;q)_{n-1}", psi0s,
      lambda q: -S(lambda n: (-1) ** n * q ** (n * (n + 1) // 2) * poch(-q, q, n - 1), 1), Q)
check("psi0* = q*F(-q)", psi0s, lambda q: q * F(-q, q), Q)
check("1 - psi1* - psi0* = sum (-1)^n q^{n(n-1)/2}(-q;q)_n", lambda q: 1 - psi1s(q) - psi0s(q),
      lambda q: S(lambda n: (-1) ** n * q ** binom2(n) * poch(-q, q, n)), Q)
# halves
check("sum q^{2n^2+n}(-q;q)_{2n} = (psi1 + psi1*)/2",
      lambda q: S(lambda n: q ** (2 * n * n + n) * poch(-q, q, 2 * n), 0, N // 2), lambda q: (psi1(q) + psi1s(q)) / 2, Q)
check("sum q^{2n^2+3n}(-q;q)_{2n} = (psi0 + psi0*)/(2q)",
      lambda q: S(lambda n: q ** (2 * n * n + 3 * n) * poch(-q, q, 2 * n), 0, N // 2), lambda q: (psi0(q) + psi0s(q)) / (2 * q), Q)
# the identities
check("(I)  2 psi0* = f0* - theta4 R0", lambda q: 2 * psi0s(q), lambda q: f0s(q) - th4(q) * R0(q), Q)
check("(II) 2 psi1* = f1* + theta4 R1", lambda q: 2 * psi1s(q), lambda q: f1s(q) + th4(q) * R1(q), Q)
check("theta4 = (q;q)_inf/(-q;q)_inf", th4, lambda q: poch_inf(q, q) / poch_inf(-q, q), Q)
check("R0 = A_q(1), R1 = A_q(q)", lambda q: R0(q) + 3 * R1(q), lambda q: Aq(1, q) + 3 * Aq(q, q), Q)
check("G = A_q(-1), H = A_q(-q)  (Rogers-Ramanujan)", lambda q: Aq(-1, q),
      lambda q: 1 / (poch_inf(q, q ** 5) * poch_inf(q ** 4, q ** 5)), Q)
# double sums
def ds1(q):
    return S(lambda k: S(lambda n: (-1) ** n * q ** (k * (k + 1) // 2 + n * n + n + k * n), 0, 50) / poch(q, q, k) if k % 2 == 0 else 0, 0, 60)
def ds0(q):
    return S(lambda k: S(lambda n: (-1) ** n * q ** (k * (k + 1) // 2 + n * n + k * n), 0, 50) / poch(q, q, k) if k % 2 == 1 else 0, 0, 60)
check("(-q;q)_inf psi1* = double sum over k even", lambda q: poch_inf(-q, q) * psi1s(q), ds1, Q)
check("(-q;q)_inf psi0* = double sum over k odd", lambda q: poch_inf(-q, q) * psi0s(q), ds0, Q)
# q-difference equation
XZ = [(mp.mpf('0.3'), mp.mpf('0.2')), (mp.mpc('-0.7', '0.4'), mp.mpc('0.5', '0.3')), (mp.mpf(-1), mp.mpf('0.6'))]
check("F(z) = 1 + zq F(zq) + zq^2 F(zq^2)", lambda z, q: F(z, q), lambda z, q: 1 + z * q * F(z * q, q) + z * q * q * F(z * q * q, q), XZ)
# H(x) identity
XQ = [(mp.mpf(3), mp.mpf('0.3')), (mp.mpc('0.4', '0.7'), mp.mpc('0.2', '-0.4')), (mp.mpf('-0.5'), mp.mpf('0.55'))]
check("H(x) = (1 - x^2) g2(x,q) + 1 - x", lambda x, q: H(x, q), lambda x, q: (1 - x * x) * g2(x, q) + 1 - x, XQ)
w = mp.exp(2j * mp.pi / 3)
check("sum q^binom(n,2) (-q;q)_n^2/(-q^3;q^3)_n = H(-w)",
      lambda q: S(lambda n: q ** binom2(n) * poch(-q, q, n) ** 2 / poch(-q ** 3, q ** 3, n)), lambda q: H(-w, q), Q)
check("sum q^binom(n,2) (q^2;q^2)_n/(q^3;q^3)_n = H(w)",
      lambda q: S(lambda n: q ** binom2(n) * poch(q * q, q * q, n) / poch(q ** 3, q ** 3, n)), lambda q: H(w, q), Q)

# --- section 9: identities (III)-(VII) ----------------------------------------
F0s = lambda q: S(lambda n: (-1) ** n * q ** (2 * n * n) / poch(q, q * q, n), 0, 60)
F1s = lambda q: S(lambda n: (-1) ** n * q ** (2 * n * n + 2 * n) / poch(q, q * q, n + 1), 0, 60)
Xs  = lambda q: S(lambda n: (-1) ** n * q ** (n * n) / poch(q, q, 2 * n), 0, 60)
Ys  = lambda q: S(lambda n: (-1) ** n * q ** ((n + 1) ** 2) / poch(q, q, 2 * n + 1), 0, 60)
check("(III) F0* + sum q^{n^2}(q;q^2)_n = theta4(q^2)(-q;q)_inf X(q^2)",
      lambda q: F0s(q) + S(lambda n: q ** (n * n) * poch(q, q * q, n), 1), lambda q: th4(q * q) * poch_inf(-q, q) * Xs(q * q), Q)
check("(IV) qF1* - sum q^{(n+1)^2}(q;q^2)_n = theta4(q^2)(-q;q)_inf Y(q^2)",
      lambda q: q * F1s(q) - S(lambda n: q ** ((n + 1) ** 2) * poch(q, q * q, n)), lambda q: th4(q * q) * poch_inf(-q, q) * Ys(q * q), Q)
Wd  = lambda q: S(lambda j: (-1) ** j * q ** (j * j + j) * poch(q, q * q, j), 0, 60)
W3d = lambda q: S(lambda j: (-1) ** j * q ** (j * j + 3 * j) * poch(q, q * q, j), 0, 60)
Fa  = lambda q: S(lambda n: q ** (2 * n * n - n) / poch(q, q * q, n), 1, 60)
Fb  = lambda q: S(lambda n: q ** (2 * n * n - 3 * n + 1) / poch(q, q * q, n), 1, 60)
psiR = lambda q: S(lambda n: q ** (n * (n + 1) // 2))
Pp  = lambda q: S(lambda k: q ** (2 * k * k + k) / poch(q * q, q * q, k), 0, 60)
Pm  = lambda q: S(lambda k: q ** (2 * k * k - k) / poch(q * q, q * q, k), 0, 60)
check("(V)  W + F_a = psi P_+", lambda q: Wd(q) + Fa(q), lambda q: psiR(q) * Pp(q), Q)
check("(VI) q W3 + F_b = psi P_-", lambda q: q * W3d(q) + Fb(q), lambda q: psiR(q) * Pm(q), Q)
Jc = lambda r: S(lambda m: r ** (2 * m * m + m) / poch(-r * r, r * r, m), 0, 60)
Kc = lambda r: 2 * S(lambda n: r ** ((n + 1) ** 2) * poch(-r * r, r * r, n), 0, 60)
check("(VII) J + K = psi(r)(r^2;r^4)_inf sum r^{m(m+1)/2}/(r^2;r^2)_m (not proved)", lambda r: Jc(r) + Kc(r),
      lambda r: psiR(r) * poch_inf(r * r, r ** 4) * S(lambda m: r ** (m * (m + 1) // 2) / poch(r * r, r * r, m)), Q)
# --- section 13: the dual series at z = 1, and the four fifth-order laws with the Mordell integral ------------
calW = lambda T, p: S(lambda k: (-T) ** k * p ** (k * k + k) * poch(p, p * p, k), 0, 60)
phi0 = lambda q: S(lambda n: q ** (n * n) * poch(-q, q * q, n), 0, 60)
phi1 = lambda q: S(lambda n: q ** ((n + 1) ** 2) * poch(-q, q * q, n), 0, 60)
check("calW(1/p; p) = phi0(-p)", lambda p: calW(1 / p, p), lambda p: phi0(-p), Q)
check("p calW(p; p) = -phi1(-p)", lambda p: p * calW(p, p), lambda p: -phi1(-p), Q)
def mordell(r, al):                        # L(r, alpha), ray turned so that the Gaussian decays
    ph = mp.expj(-mp.arg(al) / 2)
    g = lambda s: (lambda y: mp.exp(-3 * al * y * y / 2) * (mp.cosh((3 * r - 2) * al * y) + mp.cosh((3 * r - 1) * al * y)) / mp.cosh(3 * al * y / 2))(s * ph) * ph
    return mp.quad(g, [c / mp.sqrt(abs(al)) for c in (0, 1, 3, 8, 20)])
f1u = lambda q: S(lambda n: q ** (n * n + n) / poch(-q, q, n))
f0u = lambda q: S(lambda n: q ** (n * n) / poch(-q, q, n))
F0u = lambda q: S(lambda n: q ** (2 * n * n) / poch(q, q * q, n), 0, 60)
F1u = lambda q: S(lambda n: q ** (2 * n * n + 2 * n) / poch(q, q * q, n + 1), 0, 60)
def laws(x):
    q, p, t = mp.expjpi(2 * x), mp.expjpi(-1 / x), -2j * mp.pi * x
    qa = lambda al: mp.expjpi(2 * x * al); pa = lambda al: mp.expjpi(-al / x)
    cp = mp.sqrt(mp.pi * (5 + mp.sqrt(5)) / (10 * t)); cm = mp.sqrt(mp.pi * (5 - mp.sqrt(5)) / (10 * t))
    M1 = mp.sqrt(60 * t / mp.pi) * mordell(mp.mpf(2) / 5, 10 * t); M0 = mp.sqrt(60 * t / mp.pi) * mordell(mp.mpf(1) / 5, 10 * t)
    F0, F1, f0, f1 = F0u(p) - 1, F1u(p), phi0(-p), phi1(-p)
    return [qa(mp.mpf(11) / 60) * f1u(q) - (M1 + 2 * cp * pa(-mp.mpf(1) / 120) * F0 - 2 * cm * pa(mp.mpf(71) / 120) * F1),
            2 * qa(mp.mpf(11) / 60) * psi1(q) - (-M1 + 2 * cp * pa(-mp.mpf(1) / 120) * f0 + 2 * cm * pa(-mp.mpf(49) / 120) * f1),
            qa(-mp.mpf(1) / 60) * f0u(q) - (M0 + 2 * cm * pa(-mp.mpf(1) / 120) * F0 + 2 * cp * pa(mp.mpf(71) / 120) * F1),
            2 * qa(-mp.mpf(1) / 60) * psi0(q) - (-M0 + 2 * cm * pa(-mp.mpf(1) / 120) * f0 - 2 * cp * pa(-mp.mpf(49) / 120) * f1)]
XL = [(mp.mpc('0.3', '0.25'),), (mp.mpc('-0.1', '0.3'),)]
for k, name in enumerate(("f1", "psi1", "f0", "psi0")):
    check(f"fifth-order law for {name} with the Mordell integral (sec. 13)", lambda x, k=k: laws(x)[k], lambda x: 0, XL)

w_ = max(len(n) for n, _ in checks)
for name, err in checks:
    print(f"{name:<{w_}}  max|diff| = {mp.nstr(err, 3):>9}  {'OK' if err < mp.mpf(10) ** -30 else 'FAIL'}")
