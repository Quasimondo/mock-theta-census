"""Shared definitions for addendum v12.   b^2 = x (Im x > 0),  q = e^(2 pi i x),  p = e^(-pi i/x),  qt = p^2.
    Phi_b(y) = (-q^(1/2) e^(2 pi b y); q)_inf / (-p e^(2 pi y/b); p^2)_inf          (Garoufalidis-Kashaev convention)
    I(u)  = int_L e^(-pi i y^2 + 2 pi u y) / Phi_b(y) dy                                 (kernel-free)
    J(u)  = int_L e^(-pi i y^2 + 2 pi u y) / ( Phi_b(y) (1 + e^(2 pi y/b)) ) dy          (with kernel)
    L = e^(-i delta) R + y0  (y0 = 0 unless stated)."""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, qp, quad, nstr
I = mpc(0, 1)


class Setup:
    def __init__(self, x):
        self.x = x = mpc(x)
        self.b = sqrt(x)
        self.q = exp(2*pi*I*x); self.p = exp(-pi*I/x)
        self.qh = exp(pi*I*x)
        self.tiny = mpf(10)**(-mp.dps - 10)

    def qpow(self, a):           # q^a, principal in x
        return exp(2*pi*I*self.x*a)

    def ppow(self, a):           # p^a
        return exp(-pi*I/self.x*a)

    def Phi(self, y):
        b = self.b
        return qp(-self.qh*exp(2*pi*b*y), self.q)/qp(-self.p*exp(2*pi*y/b), self.p**2)

    # ---- integrals
    def line(self, f, delta=mpf(2)/5, y0=0, Lmax=10):
        e = exp(-I*delta)
        pts = [-Lmax, -7, -5, -3.5, -2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5, 3.5, 5, 7, Lmax]
        return quad(lambda s: f(y0 + s*e)*e, pts)

    def Iu(self, u, **kw):
        return self.line(lambda y: exp(-pi*I*y*y + 2*pi*u*y)/self.Phi(y), **kw)

    def Ju(self, u, **kw):
        b = self.b
        return self.line(lambda y: exp(-pi*I*y*y + 2*pi*u*y)/(self.Phi(y)*(1 + exp(2*pi*y/b))), **kw)

    # ---- q-series
    def Psi(self, z):            # sum_{n>=0} q^(n(n+1)/2) z^n (-q;q)_n
        q = self.q; s = mpc(0); term = mpc(1); n = 0
        while True:
            s += term
            n += 1
            term *= q**n*z*(1 + q**n)
            if abs(term) < self.tiny and abs(q)**(n*(n + 1)//2)*abs(z)**n < self.tiny: break
        return s

    def fstar(self, z):          # sum_{n>=0} q^(n^2+n) z^n / (-q;q)_n
        q = self.q; s = mpc(0); term = mpc(1); n = 0
        while abs(term) > self.tiny or n < 5:
            s += term; n += 1; term *= q**(2*n)*z/(1 + q**n)
        return s

    def qside(self):             # X, N0, N1, Y
        q = self.q; X = N0 = N1 = Y = mpc(0); n = 0; pe = mpc(1)
        while True:
            po = pe*(1 - q**(2*n + 1)); s = (-1)**n
            a = q**(n*n)/pe; c = q**(n*n)/po
            X += s*a; N0 += s*a*q**n; N1 += s*c*q**n; Y += s*c*q**(2*n + 1)
            if abs(q)**(n*n) < self.tiny*mpf(10)**(-15): break
            n += 1; pe = po*(1 - q**(2*n))
        return X, N0, N1, Y

    def pside(self):             # E+, E-, E+', E-', W, W3 at the dual nome p
        p = self.p; Ep = Em = Epp = Emp = mpc(0); m = 0; den = mpc(1)
        while True:
            a = p**(m*(m + 1)//2)/den
            Epp += a; Emp += a*p**m; Ep += (-1)**m*a; Em += (-1)**m*a*p**m
            if abs(a) < self.tiny: break
            m += 1; den *= 1 - p**(2*m)
        W = W3 = mpc(0); po = mpc(1); j = 0
        while abs(p)**(j*j) > self.tiny:
            W += (-1)**j*p**(j*j + j)*po; W3 += (-1)**j*p**(j*j + 3*j)*po; po *= 1 - p**(2*j + 1); j += 1
        return Ep, Em, Epp, Emp, W, W3

    def Fab(self):               # F_a(p), F_b(p)
        p = self.p; Fa = Fb = mpc(0); po = mpc(1); n = 1
        while abs(p)**(2*n*n - 3*n + 1) > self.tiny:
            po *= 1 - p**(2*n - 1)
            Fa += p**(2*n*n - n)/po; Fb += p**(2*n*n - 3*n + 1)/po; n += 1
        return Fa, Fb

    def B(self):                 # the four bilinear combinations of v6
        X, N0, N1, Y = self.qside(); Ep, Em, Epp, Emp, W, W3 = self.pside()
        q4 = self.qpow(mpf(1)/4)
        return [N0*Epp + I/q4*Y*Ep, X*Ep + I*q4*N1*Epp, N0*Emp - I/q4*Y*Em, X*Em - I*q4*N1*Emp]

    def D(self, u):              # sum_{k,n>=0} (-1)^(kn+n) e^(pi i u (b k + n/b)) q^(k^2/4) p^(n^2/2+n) / ((q;q)_k (p^2;p^2)_n)
        b = self.b; q = self.q; p = self.p
        s = exp(pi*I*u*b); t = exp(pi*I*u/b)
        ke = ko = mpc(0); k = 0; den = mpc(1)
        while True:
            a = s**k*self.qpow(mpf(k*k)/4)/den
            if k % 2 == 0: ke += a
            else: ko += a
            if abs(a) < self.tiny and k > 4: break
            k += 1; den *= 1 - q**k
        ne = no = mpc(0); n = 0; den = mpc(1)
        while True:
            a = t**n*self.ppow(mpf(n*n)/2 + n)/den
            ne += (-1)**n*a; no += a
            if abs(a) < self.tiny and n > 4: break
            n += 1; den *= 1 - p**(2*n)
        return ke*ne + ko*no

    def cconst(self):            # c = 2^(-1/2) e^(-i pi/4) q^(-1/24) qt^(1/24)
        return exp(-I*pi/4)/sqrt(2)*self.qpow(-mpf(1)/24)*self.ppow(mpf(1)/12)

    def P(self):                 # prefactor of v11:  (-q;q)_inf / (b q^(1/8) (p;p^2)_inf)
        return qp(-self.q, self.q)/(self.b*self.qpow(mpf(1)/8)*qp(self.p, self.p**2))

    def kappa(self):             # kappa of v10
        q = self.q; p = self.p
        psip = sum(p**(k*(k + 1)//2) for k in range(200))
        return qp(q, q)/qp(-q, q)*qp(p, p*p)/(2*I*self.qpow(mpf(1)/4)*qp(-q, q)*psip)


POINTS = [mpc("0.3", "0.2"), 1/mpc(3, -mpf(2)/5), mpc("0.11", "0.03"), mpc("0.45", "0.3")]
