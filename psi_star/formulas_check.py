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

w_ = max(len(n) for n, _ in checks)
for name, err in checks:
    print(f"{name:<{w_}}  max|diff| = {mp.nstr(err, 3):>9}  {'OK' if err < mp.mpf(10) ** -30 else 'FAIL'}")
