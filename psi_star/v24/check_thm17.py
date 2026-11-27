# Independent check of addendum-v24 (own code: product Phi_b, own quadrature, own q-series, own Mordell integral).
# (1) Proposition M: I(0) by quadrature and I(1/b), I(b), I(b+1/b) by quadrature, against the closed forms.
# (2) Theorem 17: M1 = 2i q^{11/60} P J(1/b), M0 = 2i q^{29/60} P J(b+1/b) by quadrature, against
#     (60t/pi)^{1/2} L(2/5, 10t), L(1/5, 10t).
# (3) Theorem 16 with M_j replaced by the Mordell integrals: the four classical laws for f1, psi1, f0, psi0,
#     using only Ramanujan's series (no state integrals at all).
import mpmath as mp
mp.mp.dps = 25
def pinf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-5): return s
        k += 1
def L(r, al):                      # Mordell integral, ray turned to make the Gaussian decay: x = s e^{-i arg(al)/2}
    ph = mp.expj(-mp.arg(al)/2)
    f = lambda s: (lambda x: mp.exp(-3*al*x**2/2)*(mp.cosh((3*r-2)*al*x) + mp.cosh((3*r-1)*al*x))/mp.cosh(3*al*x/2))(s*ph)*ph
    return mp.quad(f, [c/mp.sqrt(abs(al)) for c in (0, 1, 3, 8, 20)])
def ser(name, q, N=400):
    out = 0; A = mp.mpc(1); B = mp.mpc(1)
    for n in range(N):
        if name == 'f1':   out += q**(n*n+n)/A
        if name == 'f0':   out += q**(n*n)/A
        if name == 'psi1': out += q**(n*(n+1)//2)*A
        if name == 'psi0': out += q**((n+1)*(n+2)//2)*A
        if name == 'F0':   out += q**(2*n*n)/B
        if name == 'F1':   out += q**(2*n*n+2*n)/(B*(1 - q**(2*n+1)))
        if name == 'phi0': out += q**(n*n)*B
        if name == 'phi1': out += q**((n+1)**2)*B
        A *= 1 + q**(n+1)
        B *= (1 + q**(2*n+1)) if name in ('phi0', 'phi1') else (1 - q**(2*n+1))
    return out
for x in (mp.mpc(0.3, 0.2), mp.mpc(0, 0.3), mp.mpc(0.11, 0.08)):
    b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x); t = -2j*mp.pi*x
    qa = lambda al: mp.expjpi(2*x*al); pa = lambda al: mp.expjpi(-al/x)
    Phi = lambda y: pinf(mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q)/pinf(p*mp.exp(2*mp.pi*y/b), p**2)
    delta = (mp.pi/2 - mp.arg(b))/2; e = mp.expj(-delta)
    def integ(u, kern):
        g = lambda r: (lambda y: mp.exp(-mp.pi*1j*y**2 + 2*mp.pi*u*y)/Phi(y)*(1/(1 + mp.exp(2*mp.pi*y/b)) if kern else 1))(r*e)*e
        return mp.quad(g, mp.linspace(-7, 7, 29))
    eps = mp.expjpi(-mp.mpf(1)/4); s1 = mp.sin(mp.pi/5); s2 = mp.sin(2*mp.pi/5); k = 2/mp.sqrt(5)
    closed = {'I(0)': (0, eps*k*s2*qa(-mp.mpf(1)/60)*pa(mp.mpf(1)/30)), 'I(1/b)': (1/b, eps*k*s1*qa(-mp.mpf(1)/60)*pa(mp.mpf(19)/30)),
              'I(b)': (b, eps*k*s1*qa(-mp.mpf(19)/60)*pa(mp.mpf(1)/30)), 'I(b+1/b)': (b + 1/b, -eps*k*s2*qa(-mp.mpf(19)/60)*pa(mp.mpf(19)/30))}
    print(f'x = {mp.nstr(x, 3)}')
    print('  (1) Prop. M, |quadrature/closed - 1|: ' + '  '.join(f'{n}: {mp.nstr(abs(integ(u, False)/c - 1), 2)}' for n, (u, c) in closed.items()))
    P = 2**mp.mpf(-0.5)*mp.exp(mp.pi**2/(12*t) + t/6)/b
    M1 = 2j*qa(mp.mpf(11)/60)*P*integ(1/b, True); M0 = 2j*qa(mp.mpf(29)/60)*P*integ(b + 1/b, True)
    L1 = mp.sqrt(60*t/mp.pi)*L(mp.mpf(2)/5, 10*t); L0 = mp.sqrt(60*t/mp.pi)*L(mp.mpf(1)/5, 10*t)
    print(f'  (2) Thm 17: |M1 - Mordell| = {mp.nstr(abs(M1 - L1), 2)},  |M0 - Mordell| = {mp.nstr(abs(M0 - L0), 2)}   (M1 = {mp.nstr(M1, 10)})')
    cp = mp.sqrt(mp.pi*(5 + mp.sqrt(5))/(10*t)); cm = mp.sqrt(mp.pi*(5 - mp.sqrt(5))/(10*t))
    F0, F1, ph0, ph1 = ser('F0', p) - 1, ser('F1', p), ser('phi0', -p), ser('phi1', -p)
    laws = [('f1', qa(mp.mpf(11)/60)*ser('f1', q), L1 + 2*cp*pa(-mp.mpf(1)/120)*F0 - 2*cm*pa(mp.mpf(71)/120)*F1),
            ('psi1', 2*qa(mp.mpf(11)/60)*ser('psi1', q), -L1 + 2*cp*pa(-mp.mpf(1)/120)*ph0 + 2*cm*pa(-mp.mpf(49)/120)*ph1),
            ('f0', qa(-mp.mpf(1)/60)*ser('f0', q), L0 + 2*cm*pa(-mp.mpf(1)/120)*F0 + 2*cp*pa(mp.mpf(71)/120)*F1),
            ('psi0', 2*qa(-mp.mpf(1)/60)*ser('psi0', q), -L0 + 2*cm*pa(-mp.mpf(1)/120)*ph0 - 2*cp*pa(-mp.mpf(49)/120)*ph1)]
    print('  (3) classical laws with Mordell integrals, |lhs - rhs|: ' + '  '.join(f'{n}: {mp.nstr(abs(l - r), 2)}' for n, l, r in laws), flush=True)
