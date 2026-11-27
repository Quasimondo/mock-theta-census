# First-principles derivation of the asymptotics of B_j (addendum-v6 sec. 2) from the state integrals of addendum-v7:
#   I(u) = int Phi_b(y) e^{-pi i y^2/2 + 2 pi u y} dy = -i b (q;q)/(qt;qt) e^{-pi i c_b^2/2 + 2 pi u c_b} B_j,  u = s1 b + s2/b
# With y = z/(2 pi b), hbar = 2 pi i b^2 = -t:  integrand = exp(f(z)/hbar) g(z) (1 - (hbar/24) Li0(-e^z) + O(hbar^3)),
#   f(z) = Li2(-e^z) + z^2/4 + 2 pi i s2 z,  g(z) = e^{s1 z}.
# Laplace: int e^{f/hbar} G dz ~ e^{f0/hbar} sqrt(2 pi hbar/(-f2)) [G0 + hbar(-G2/(2f2) + G1 f3/(2f2^2) + G0 f4/(8 f2^2) - 5 G0 f3^2/(24 f2^3))].
# Prefactor (exact): (q;q)/(qt;qt) = e^{-hbar/24} e^{pi^2/(6 hbar)} sqrt(-2 pi/hbar);
#   -pi i c_b^2/2 + 2 pi u c_b = hbar(1/16 + s1/2) + pi i(1/4 + s1 + s2) - (2 pi^2/hbar)(1/8 + s2).
# Result: B = e^{A/hbar} C (1 + kappa1 hbar + ...), so the coefficient of t is -kappa1.
import mpmath as mp
mp.mp.dps = 30
Cl2 = mp.clsin(2, mp.pi/3); r3 = mp.sqrt(3)
claims = {1: (1/4, 0, '+', -mp.mpf(1)/24 - 1j/(18*r3)), 2: (-1/4, 0, '+', mp.mpf(1)/24 - 1j/(18*r3)),
          3: (1/4, -1/2, '-', -mp.mpf(1)/24 + 1j/(18*r3)), 4: (-1/4, -1/2, '-', mp.mpf(1)/24 + 1j/(18*r3))}
c_claim = mp.mpf(4/3)**0.25 * mp.expjpi(-mp.mpf(1)/12)
for j, (s1, s2, sgn, corr_claim) in claims.items():
    s1 = mp.mpf(s1); s2 = mp.mpf(s2)
    f  = lambda z: mp.polylog(2, -mp.exp(z)) + z**2/4 + 2j*mp.pi*s2*z
    f1 = lambda z: -mp.log(1 + mp.exp(z)) + z/2 + 2j*mp.pi*s2
    f2 = lambda z: -mp.exp(z)/(1 + mp.exp(z)) + mp.mpf(1)/2
    f3 = lambda z: -mp.exp(z)/(1 + mp.exp(z))**2
    f4 = lambda z: -mp.exp(z)*(1 - mp.exp(z))/(1 + mp.exp(z))**3
    target = (1j if sgn == '+' else -1j) * Cl2          # B ~ e^{target/t} = e^{-target/hbar}
    found = None
    for z0 in [2j*mp.pi/3, -2j*mp.pi/3, 4j*mp.pi/3, -4j*mp.pi/3, 2j*mp.pi, 0.5j]:
        try: z = mp.findroot(f1, z0)
        except Exception: continue
        A = f(z) - mp.pi**2/6 + 2*mp.pi**2*(mp.mpf(1)/8 + s2)
        # A is defined up to 2 pi i * (integer) * (something) from branch choices; compare mod 4 pi^2 (Li2 branch shifts)
        if abs(mp.re(A)) < 1e-20 and abs(abs(mp.im(A)) - Cl2) < 1e-20 and mp.sign(mp.im(-A)) == mp.sign(mp.im(target)):
            found = (z, A); break
    if not found:
        print(f'B{j}: no matching saddle found'); continue
    z, A = found
    G0 = mp.exp(s1*z); G1 = s1*G0; G2 = s1**2*G0
    F2, F3, F4 = f2(z), f3(z), f4(z)
    K = -(mp.mpf(1)/24)*G0*mp.polylog(0, -mp.exp(z)) + (-G2/(2*F2) + G1*F3/(2*F2**2) + G0*F4/(8*F2**2) - 5*G0*F3**2/(24*F2**3))
    kappa1 = K/G0 + mp.mpf(1)/24 - (mp.mpf(1)/16 + s1/2)
    # constant: (1/(2 pi b)) sqrt(2 pi/(-f2)) sqrt(hbar) G0 / (-i b sqrt(-2 pi/hbar) e^{pi i (1/4+s1+s2)}) = i G0 / (sqrt(-f2) sqrt(-1)) e^{-pi i(1/4+s1+s2)} (up to sign)
    C = 1j*G0/mp.sqrt(F2) * mp.expjpi(-(mp.mpf(1)/4 + s1 + s2))
    cc = c_claim if j in (1, 4) else mp.conj(c_claim)
    print(f'B{j}: saddle z = {mp.nstr(z, 6)};  |C| = {mp.nstr(abs(C), 12)} (claimed {mp.nstr(abs(cc), 12)});  C/c_claimed = {mp.nstr(C/cc, 6)}')
    print(f'     coefficient of t: derived {mp.nstr(-kappa1, 15)}   claimed {mp.nstr(corr_claim, 15)}   |diff| = {mp.nstr(abs(-kappa1 - corr_claim), 2)}')
