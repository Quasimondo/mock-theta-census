# Independent check of addendum-v10 sec. 1 at finite X:
#   psi1*(q) - (kappa/2) [B2 p W3(p) - B4 W(p)]  =  Borel sum of Phi  ~  Phi optimally truncated (error ~ e^{-|S|/|t|})
# The test that matters: no extra multiple lambda W(p) in the decaying term (lambda = 0).
import mpmath as mp
import sys; sys.path.insert(0,'../v10'); from phi_psi1 import PHI
mp.mp.dps = 90
def ser(f, nmax): return mp.fsum(f(n) for n in range(nmax))
def run(Xc):
    X = mp.mpc(Xc); q = mp.exp(2j*mp.pi/X); t = -2j*mp.pi/X; p = mp.exp(-1j*mp.pi*X)
    nq = int(40 * abs(X)) + 50
    def poch(a, base, n):               # prod_{k<n} (1 - a base^k)
        P = mp.mpc(1)
        for k in range(n): P *= 1 - a * base**k
        return P
    # psi1*(q) directly
    s = 0; P = mp.mpc(1)
    for n in range(nq):
        s += (-1)**n * q**(n*(n+1)//2) * P; P *= 1 + q**(n+1)
    psi1 = s
    # X(q), N1(q)
    Xq = 0; N1 = 0; P = mp.mpc(1)
    for n in range(nq):
        Xq += (-1)**n * q**(n*n) / P
        P1 = P * (1 - q**(2*n+1)); N1 += (-1)**n * q**(n*n+n) / P1
        P = P1 * (1 - q**(2*n+2))
    th4 = 1 + 2 * mp.fsum((-1)**n * q**(n*n) for n in range(1, nq))
    mq = mp.mpc(1)
    for k in range(1, 40 * nq):
        f = q**k; mq *= 1 + f
        if abs(f) < mp.mpf(10)**-95: break
    def E(sign, a):
        s = 0; Q = mp.mpc(1)
        for m in range(120):
            s += sign**m * p**(m*(m+a)//2) / Q; Q *= 1 - p**(2*m+2)
        return s
    Ep, Epp, Em, Emp = E(-1, 1), E(1, 1), E(-1, 3), E(1, 3)
    def G(a):
        s = 0; P = mp.mpc(1)
        for j in range(120):
            s += (-1)**j * p**(j*j+a*j) * P; P *= 1 - p**(2*j+1)
        return s
    W, W3 = G(1), G(3)
    pp2 = poch(p, p**2, 200); psip = mp.fsum(p**(n*(n+1)//2) for n in range(200))
    q14 = mp.exp(-t/4)
    kappa = th4 * pp2 / (2j * q14 * mq * psip)
    B2 = Xq * Ep + 1j * q14 * N1 * Epp
    B4 = Xq * Em - 1j * q14 * N1 * Emp
    grow = -(kappa/2) * B4 * W; dec = (kappa/2) * B2 * p * W3; decW = (kappa/2) * B2 * W
    T = psi1 - grow - dec
    terms = [mp.mpf(PHI[n].numerator)/PHI[n].denominator * t**n for n in range(len(PHI))]
    nopt = min(range(1, len(terms)), key=lambda n: abs(terms[n]))
    phi_trunc = mp.fsum(terms[:nopt])
    return T, psi1
