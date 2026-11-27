# Independent check of addendum-v22 (own code).
# (1) Theorem 13 at real x (tau = 0): saddle w0 = Log W (Im W > 0) of V(w) = -w^2/2 - 2 pi i (a-1) w + Li2(-e^w);
#     V'(w0) = 0, Im V(w0) = -Vol(z) (Bloch-Wigner), and the height H = -Im V on the five pieces C1..C5 of the
#     contour: sup of (H - H(w0))/min(d^2, d) must be < 0 (condition (L1)). Li2 continued as stated in v22 sec. 1
#     (principal in the strip, boundary values from inside; for Re w <= 0 the 2 pi i-periodic power series).
# (2) The claim of sec. 5: along the horizontal ray w0 - s, -Im V rises (+7.4 at s = 8 for a = 1/2).
# (3) Theorem 14's limits: on the path x = 1/X + i eps, the dual series W(T), W(p^2 T) -> 1 and p T -> z^X.
import mpmath as mp
mp.mp.dps = 25
def Li2(w, where):
    z = -mp.exp(w)
    if where == 'left':                      # Re w <= 0: power series, 2 pi i-periodic
        return mp.polylog(2, z)
    if where == 'strip_top':                 # Im w = pi from inside: -e^w = e^r - i0
        return mp.polylog(2, mp.mpc(z.real, -mp.mpf(10)**-30))
    return mp.polylog(2, z)                  # principal (strip interior / real axis)
def V(w, a, where='strip'):
    return -w**2/2 - 2j*mp.pi*(a-1)*w + Li2(w, where)
def D(z): return mp.im(mp.polylog(2, z)) + mp.arg(1 - z)*mp.log(abs(z))   # Bloch-Wigner
worst_all = []
for a in (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99):
    a = mp.mpf(a); z = mp.expj(2*mp.pi*a)
    W = [r for r in mp.polyroots([1, 1, -1/z]) if mp.im(r) > 0][0]
    w0 = mp.log(W); r0, e0 = w0.real, w0.imag
    Vp = -w0 - 2j*mp.pi*(a-1) - mp.log(1 + mp.exp(w0))
    H0 = -mp.im(V(w0, a)); vol = abs(D(-W))
    pieces = {
      'C1': [(mp.mpc(r, 2*mp.pi), 'left') for r in mp.linspace(-40, -1e-6, 400)],
      'C2': [(mp.mpc(0, e), 'left') for e in mp.linspace(2*mp.pi, mp.pi, 300)],
      'C3': [(mp.mpc(r, mp.pi), 'strip_top') for r in mp.linspace(0, r0, 200)[1:-1]],
      'C4': [(mp.mpc(r0, e), 'strip') for e in mp.linspace(mp.pi - 1e-9, 1e-9, 600)],
      'C5': [(mp.mpc(r, 0), 'strip') for r in mp.linspace(r0, 40, 400)]}
    row = []
    for name, pts in pieces.items():
        s = -mp.inf
        for w, where in pts:
            d = abs(w - w0)
            if d < 1e-3: continue
            s = max(s, (-mp.im(V(w, a, where)) - H0)/min(d*d, d))
        row.append(s)
    worst_all.append((a, max(row)))
    print(f'a = {mp.nstr(a,3):5s}: |V\'(w0)| = {mp.nstr(abs(Vp),2):8s} Im V(w0) + Vol = {mp.nstr(mp.im(V(w0,a)) + vol, 2):8s}  sup per piece C1..C5: ' + '  '.join(mp.nstr(s, 3) for s in row))
print('max over a of the worst sup:', mp.nstr(max(s for _, s in worst_all), 4), ' (must be < 0)')
print('(2) horizontal ray from the growing saddle, a = 1/2 and 0.2: rise of -Im V at s = 2, 8')
for a in (mp.mpf(1)/2, mp.mpf('0.2')):
    z = mp.expj(2*mp.pi*a); W = [r for r in mp.polyroots([1, 1, -1/z]) if mp.im(r) > 0][0]; w0 = mp.log(W)
    print('  a =', mp.nstr(a, 2), ' '.join(f's={s}: {mp.nstr(-mp.im(V(w0 - s, a)) + mp.im(V(w0, a)), 4)}' for s in (2, 8)))
print('(3) limits in Theorem 14, X = 8, z = e^{2 pi i a}')
X = 8
for a in (mp.mpf('0.25'), mp.mpf('0.3')):
    for eps in (mp.mpf('1e-3'), mp.mpf('1e-4')):
        with mp.workdps(60):
            x = mp.mpf(1)/X + 1j*eps; p = mp.expjpi(-1/x); T = p*mp.expjpi(2*a/x)
            def Wc(T):
                s = 0; po = mp.mpc(1); k = 0
                while abs(p)**(k*k) > mp.mpf(10)**-40 or k < 5:
                    s += (-T)**k*p**(k*k+k)*po; po *= 1 - p**(2*k+1); k += 1
                return s
            print(f'  a = {mp.nstr(a,2)}, eps = {mp.nstr(eps,1)}: |W(T) - 1| = {mp.nstr(abs(Wc(T)-1),3)}  |W(p^2 T) - 1| = {mp.nstr(abs(Wc(p*p*T)-1),3)}  |pT - z^X| = {mp.nstr(abs(p*T - mp.expj(2*mp.pi*a*X)),3)}')
