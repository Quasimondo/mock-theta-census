# Independent check of addendum-v30 (own code: product Phi_b, own quadrature for J, own double series for I
# (Theorem C of v12), own q- and p-series for the frames).
# (1) Theorem 20: Jt.S = Delta.J, det S formula, and the Corollary (calW(T), pT calW(p^2T) from F and S).
# (2) Theorem 21 (T) and (gamma), with lambda as stated.
# (3) The relations (A0), (A1), (R1)-(R4) in the form of Andrews-Garvan (2.12)-(2.21), read with q -> -q.
import mpmath as mp
mp.mp.dps = 30
def pinf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-mp.mp.dps-5): return s
        k += 1
def pm(q, n):
    p = mp.mpc(1)
    for j in range(1, n+1): p *= 1 + q**j
    return p
def setup(a, x):
    b = mp.sqrt(x); q = mp.expjpi(2*x); p = mp.expjpi(-1/x)
    qa = lambda al: mp.expjpi(2*x*al); pa = lambda al: mp.expjpi(-al/x)
    zeta = mp.expjpi(a); z = zeta**2; zt = mp.expjpi(a/x); T = p*zt**2
    c = 1/mp.sqrt(2)*mp.expjpi(-mp.mpf(1)/4)*qa(-mp.mpf(1)/24)*pa(mp.mpf(1)/12)
    def I(u, K=50):
        D = 0; qk = mp.mpc(1)
        for k in range(K):
            if k: qk *= 1 - q**k
            pn = mp.mpc(1)
            for n in range(K):
                if n: pn *= 1 - p**(2*n)
                D += (-1)**(k*n+n)*mp.exp(mp.pi*1j*u*(b*k + n/b))*mp.expjpi(x*k*k/2)*mp.expjpi(-(n*n/mp.mpf(2) + n)/x)/(qk*pn)
        return c*mp.exp(-mp.pi*1j*u**2/2)*D
    Phi = lambda y: pinf(mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q)/pinf(p*mp.exp(2*mp.pi*y/b), p**2)
    e = mp.expj(-(mp.pi/2 - mp.arg(b))/2)
    J = lambda u: mp.quad(lambda r: mp.exp(-mp.pi*1j*(r*e)**2 + 2*mp.pi*u*r*e)/(Phi(r*e)*(1 + mp.exp(2*mp.pi*r*e/b)))*e, mp.linspace(-7, 7, 29))
    P = 2**mp.mpf(-0.5)*qa(-mp.mpf(1)/6)*pa(-mp.mpf(1)/24)/b
    pre = lambda A: 1j*P*mp.expjpi(-A)
    S = mp.matrix([[1, pre(a)*J(a/b), pre(a+x)*J((a+x)/b)], [0, pre(a)*I(a/b), pre(a+x)*I((a+x)/b)], [0, pre(a)*I((a-1)/b), pre(a+x)*I((a+x-1)/b)]])
    def H(y, r, r14, part=None):
        s = 0; rm = mp.mpc(1)
        for m in range(80):
            if m: rm *= 1 - r**m
            if part == 'e' and m % 2: continue
            if part == 'o' and not m % 2: continue
            s += y**m*r14**(m*m)/rm
        return s
    q14 = mp.expjpi(x/2); p12 = mp.expjpi(-1/(2*x))
    F = lambda Z: sum(Z**n*q**(n*(n+1)//2)*pm(q, n) for n in range(150))
    mu = qa(-mp.mpf(3)/4)/zeta; qh = mp.expjpi(x)
    Jq = mp.matrix([[1, F(z), F(z*q)], [0, H(zeta, q, q14), mu*H(zeta*qh, q, q14)], [0, H(-zeta, q, q14), -mu*H(-zeta*qh, q, q14)]])
    def W(TT):
        s = 0; po = mp.mpc(1)
        for k in range(60): s += (-TT)**k*p**(k*k+k)*po; po *= 1 - p**(2*k+1)
        return s
    p2 = p**2
    Ht = lambda y, part: H(y, p2, p12, part)
    # lower block of the dual frame: Mt^-1, Mt = [[He(zt p), -Ho(zt p)], [-zt p^(1/2) Ho(zt p^2), zt p^(1/2) He(zt p^2)]], H(.; p^2)
    Mt = mp.matrix([[Ht(zt*p, 'e'), -Ht(zt*p, 'o')], [-zt*p12*Ht(zt*p2, 'o'), zt*p12*Ht(zt*p2, 'e')]])
    Mi = Mt**-1
    Jt = mp.matrix([[1, p*T*W(p2*T), -W(T)], [0, Mi[0, 0], Mi[0, 1]], [0, Mi[1, 0], Mi[1, 1]]])
    delta = mp.expjpi(mp.mpf(1)/4)/(2*b)*qa(-mp.mpf(5)/24)*pa(mp.mpf(1)/24)*mp.expjpi(-a - a*a/(2*x))
    return dict(S=S, Jq=Jq, Jt=Jt, delta=delta, q=q, p=p, z=z, zeta=zeta, zt=zt, T=T, F=F, W=W, b=b)
def mx(M): return max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))
print('(1) Theorem 20')
for a, x in ((mp.mpf(1)/4, mp.mpc(0.3, 0.2)), (mp.mpc(0.3, 0.1), mp.mpc(-0.2, 0.35)), (mp.mpf(0.8), mp.mpc(0.45, 0.3))):
    d = setup(a, x)
    D = mp.diag([1, d['delta'], d['delta']])
    res = d['Jt']*d['S'] - D*d['Jq']
    detS = d['S'][1, 1]*d['S'][2, 2] - d['S'][1, 2]*d['S'][2, 1]
    qa = lambda al: mp.expjpi(2*x*al); pa = lambda al: mp.expjpi(-al/x)
    detF = -(1j/(2*x))*qa(-mp.mpf(7)/6)*pa(mp.mpf(7)/12)*mp.expjpi(-3*a + a/x - a*a/x)
    detF2 = d['delta']**2*(-2*qa(-mp.mpf(3)/4)/d['zeta'])*d['zt']*mp.expjpi(-1/(2*x))
    Si = d['S']**-1; row = mp.matrix([[1, d['F'](d['z']), d['F'](d['z']*d['q'])]])*Si
    cor = max(abs(row[0, 1] - d['p']*d['T']*d['W'](d['p']**2*d['T'])), abs(row[0, 2] + d['W'](d['T'])))
    print(f'  a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: max|Jt S - Delta J| = {mp.nstr(mx(res),2)} (max entry {mp.nstr(mx(D*d["Jq"]),3)});  det S (2x2 block) vs closed form {mp.nstr(abs(detS/detF - 1),2)}, vs delta^2 detJ/detJt {mp.nstr(abs(detS/detF2 - 1),2)};  Corollary {mp.nstr(cor,2)}', flush=True)
print('(2) Theorem 21')
CT = mp.matrix([[1, 0, 0], [0, (1+1j)/2, (1-1j)/2], [0, (1-1j)/2, (1+1j)/2]])
for a, x in ((mp.mpf(1)/4, mp.mpc(0.3, 0.2)), (mp.mpc(0.3, 0.1), mp.mpc(0.45, 0.3))):
    d = setup(a, x); d1 = setup(a, x + 1)
    print(f'  (T) a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: {mp.nstr(mx(d1["Jq"] - CT*d["Jq"]),2)}')
    xp = x/(1 + 2*x); ap = a/(1 + 2*x); dp = setup(ap, xp)
    lam = (1 + 2*x)**mp.mpf(-0.5)*mp.expjpi(mp.mpf(1)/12)*mp.expjpi(-mp.mpf(5)*x*xp/6)*mp.expjpi(-(a - ap))*mp.expjpi(-a*ap)
    Om = d['S']**-1*dp['S']; Lam = mp.diag([1, lam, -lam])
    print(f'  (gamma) a = {mp.nstr(a,3)}, x = {mp.nstr(x,3)}: |lam/(delta/delta\') - 1| = {mp.nstr(abs(lam/(d["delta"]/dp["delta"]) - 1),2)};  max|J(a\';x\') - Lambda J Omega| = {mp.nstr(mx(dp["Jq"] - Lam*d["Jq"]*Om),2)}', flush=True)
print('(3) Andrews-Garvan (2.12)-(2.21) with q -> -q give (A0), (A1), (R1)-(R4)')
mp.mp.dps = 40
for q in (mp.mpf('0.3'), mp.mpc('0.2', '0.45')):
    S = lambda t, N=120: sum(t(n) for n in range(N))
    def poch(a, qq, n):
        r = mp.mpc(1)
        for j in range(n): r *= 1 - a*qq**j
        return r
    f0 = lambda q: S(lambda n: q**(n*n)/poch(-q, q, n)); f1 = lambda q: S(lambda n: q**(n*n+n)/poch(-q, q, n))
    phi0 = lambda q: S(lambda n: q**(n*n)*poch(-q, q*q, n)); phi1 = lambda q: S(lambda n: q**((n+1)**2)*poch(-q, q*q, n))
    psi0 = lambda q: S(lambda n: q**((n+1)*(n+2)//2)*poch(-q, q, n)); psi1 = lambda q: S(lambda n: q**(n*(n+1)//2)*poch(-q, q, n))
    F0 = lambda q: S(lambda n: q**(2*n*n)/poch(q, q*q, n)); F1 = lambda q: S(lambda n: q**(2*n*n+2*n)/poch(q, q*q, n+1))
    chi0 = lambda q: S(lambda n: q**n/poch(q**(n+1), q, n)); chi1 = lambda q: S(lambda n: q**n/poch(q**(n+1), q, n+1))
    th4 = lambda q: 1 + 2*sum((-1)**n*q**(n*n) for n in range(1, 40))
    G = lambda q: 1/(poch(q, q**5, 200)*poch(q**4, q**5, 200)); Hh = lambda q: 1/(poch(q**2, q**5, 200)*poch(q**3, q**5, 200))
    psiR = lambda q: S(lambda n: q**(n*(n+1)//2))
    rows = [('(2.12) phi0(-q) + chi0(q) = 2F0(q)', phi0(-q) + chi0(q) - 2*F0(q)),
            ('(2.13) f0(-q) + 2F0(q^2) - 2 = th4(q)G(q)', f0(-q) + 2*F0(q*q) - 2 - th4(q)*G(q)),
            ('(2.14) phi0(-q^2) + psi0(-q) = th4(q)G(q)', phi0(-q*q) + psi0(-q) - th4(q)*G(q)),
            ('(2.15) 2phi0(-q^2) - f0(q) = th4(q)G(q)', 2*phi0(-q*q) - f0(q) - th4(q)*G(q)),
            ('(2.16) psi0(q) - F0(q^2) + 1 = q psi(q^2) H(q^4)', psi0(q) - F0(q*q) + 1 - q*psiR(q*q)*Hh(q**4)),
            ('(2.17) chi1(q) - phi1(-q)/q = 2F1(q)', chi1(q) - phi1(-q)/q - 2*F1(q)),
            ('(2.18) f1(-q) - 2qF1(q^2) = th4(q)H(q)', f1(-q) - 2*q*F1(q*q) - th4(q)*Hh(q)),
            ('(2.19) phi1(-q^2)/q + psi1(-q) = th4(q)H(q)', phi1(-q*q)/q + psi1(-q) - th4(q)*Hh(q)),
            ('(2.20) 2phi1(-q^2)/q + f1(q) = th4(q)H(q)', 2*phi1(-q*q)/q + f1(q) - th4(q)*Hh(q)),
            ('(2.21) psi1(q) - qF1(q^2) = psi(q^2)G(q^4)', psi1(q) - q*F1(q*q) - psiR(q*q)*G(q**4)),
            ('(R3) psi1 - f1 = q chi1(q^2)', psi1(q) - f1(q) - q*chi1(q*q)),
            ('(R4) psi0 - f0 + 2 = chi0(q^2)', psi0(q) - f0(q) + 2 - chi0(q*q))]
    print(f'  q = {mp.nstr(q,3)}: ' + '  '.join(f'{n.split()[0]} {mp.nstr(abs(v),2)}' for n, v in rows))
