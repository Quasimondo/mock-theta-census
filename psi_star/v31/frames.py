"""Addendum v30: the q-side frame, the dual frame, the automorphy factor and the matrix of state integrals.
   x in the upper half-plane, b = sqrt(x), q = e^(2 pi i x), p = e^(-pi i/x);  zeta = e^(pi i a), zetat = e^(pi i a/x), z = zeta^2, T = p zetat^2.
   H(y; r) = sum y^m r^(m^2/4)/(r;r)_m,   F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n,   Wd(T) = sum (-T)^k p^(k^2+k) (p;p^2)_k."""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, matrix, nstr
from common import Setup, I
from a1_formula import Wt, Hfun


def Hq(S, y):
    return Hfun(y, S.q, S.qpow(mpf(1)/4))


def Hp(S, y):
    return Hfun(y, S.p**2, S.ppow(mpf(1)/2))


def frame_q(S, a):
    """J(zeta; q): rows (1, F(z), F(zq)), (0, H(zeta), q^(-3/4) zeta^-1 H(zeta q^(1/2))), (0, H(-zeta), -q^(-3/4) zeta^-1 H(-zeta q^(1/2)))"""
    ze = exp(pi*I*a); z = ze*ze; q = S.q; qh = S.qpow(mpf(1)/2); m = S.qpow(-mpf(3)/4)/ze
    return matrix([[1, S.Psi(z), S.Psi(z*q)],
                   [0, Hq(S, ze), m*Hq(S, ze*qh)],
                   [0, Hq(S, -ze), -m*Hq(S, -ze*qh)]])


def Mdual(S, a):
    zt = exp(pi*I*a/S.x); p = S.p; ph = S.ppow(mpf(1)/2)
    e1 = (Hp(S, zt*p) + Hp(S, -zt*p))/2; o1 = (Hp(S, zt*p) - Hp(S, -zt*p))/2
    e2 = (Hp(S, zt*p*p) + Hp(S, -zt*p*p))/2; o2 = (Hp(S, zt*p*p) - Hp(S, -zt*p*p))/2
    return matrix([[e1, -o1], [-zt*ph*o2, zt*ph*e2]])


def frame_dual(S, a):
    """Jt(zetat; p) = [[1, pT Wd(p^2 T), -Wd(T)], [0, Mt^-1]]"""
    zt = exp(pi*I*a/S.x); p = S.p; T = p*zt*zt
    Mi = Mdual(S, a)**-1
    return matrix([[1, p*T*Wt(S, p*p*T), -Wt(S, T)], [0, Mi[0, 0], Mi[0, 1]], [0, Mi[1, 0], Mi[1, 1]]])


def delta(S, a):
    """i P c e^(-pi i a - pi i a^2/(2x)) = e^(i pi/4)/(2b) q^(-5/24) p^(1/24) e^(-pi i a - pi i a^2/(2x))"""
    return exp(I*pi/4)/(2*S.b)*S.qpow(-mpf(5)/24)*S.ppow(mpf(1)/24)*exp(-pi*I*a - pi*I*a*a/(2*S.x))


def Delta(S, a):
    d = delta(S, a)
    return matrix([[1, 0, 0], [0, d, 0], [0, 0, d]])


def Pelem(S):
    return S.qpow(-mpf(1)/6)*S.ppow(-mpf(1)/24)/(sqrt(2)*S.b)


def state_matrix(S, a, **kw):
    """S(a, x) by quadrature: rows (1, Jc(a), Jc(a+x)), (0, i0(a), i0(a+x)), (0, i1(a), i1(a+x))
       Jc(a) = i P e^(-pi i a) J(a/b),  i0(a) = i P e^(-pi i a) I(a/b),  i1(a) = i P e^(-pi i a) I((a-1)/b)"""
    b = S.b; x = S.x; P = Pelem(S)
    from mpmath import arg
    if "delta" not in kw:        # the line e^(-i delta) R must avoid the rays of zeros and poles: |arg b| + delta < 90 degrees
        kw["delta"] = min(mpf(2)/5, (pi/2 - abs(arg(b)))/2)
        if kw["delta"] < mpf(2)/5: kw.setdefault("Lmax", 16)
    def col(a):
        f = I*P*exp(-pi*I*a)
        return f*S.Ju(a/b, **kw), f*S.Iu(a/b, **kw), f*S.Iu((a - 1)/b, **kw)
    j0, i00, i10 = col(a); j1, i01, i11 = col(a + x)
    return matrix([[1, j0, j1], [0, i00, i01], [0, i10, i11]])


def mnorm(A):
    return max(abs(A[i, j]) for i in range(A.rows) for j in range(A.cols))
