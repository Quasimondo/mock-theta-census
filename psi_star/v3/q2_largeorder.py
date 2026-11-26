"""Q2: large-order behaviour of the perturbative series of psi1*, psi0* (from q2_series.json) and
comparison with the asymptotics of the terms -+ (1/2) theta4 R of identities (I), (II).

Model at an even root (A > 0):   a_k = (-1)^k sum_j c_j Gamma(k+1/2-j) A^-(k+1/2-j),  attached series  sum (-1)^j (c_j/c_0) t^j
Model at q -> 1 (S complex):     a_k = 2 Re sum_j c_j Gamma(k+1/2-j) S^-(k+1/2-j),    attached series  sum (c_j/c_0) t^j
"""
import json, mpmath as mp
mp.mp.dps = 110
D = json.load(open("q2_series.json")); K = D["K"]
def coef(key): return [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i, _ in D[key]]
half = mp.mpf(1)/2
Cl2 = mp.clsin(2, mp.pi/3)

def fit_even(a, A, J=12, k0=80):
    rows, rhs = [], []
    for k in range(k0, K+1):
        sc = (-1)**k * mp.gamma(k+half) * A**(-(k+half))
        rows.append([mp.gamma(k+half-j) * A**(-(k+half-j)) * (-1)**k / sc for j in range(J)])
        rhs.append(a[k]/sc)
    M = mp.matrix(rows)
    cre = mp.qr_solve(M, mp.matrix([x.real for x in rhs]))[0]
    cim = mp.qr_solve(M, mp.matrix([x.imag for x in rhs]))[0]
    c = [mp.mpc(cre[j], cim[j]) for j in range(J)]
    res = max(abs(sum(rows[i][j]*c[j] for j in range(J)) - rhs[i]) for i in range(len(rhs)))
    return c, res
def fit_one(a, S, J=14, k0=75):
    rows, rhs = [], []
    for k in range(k0, K+1):
        sc = mp.gamma(k+half) * abs(S)**(-(k+half))
        row = []
        for j in range(J):
            g = mp.gamma(k+half-j) * S**(-(k+half-j)) / sc
            row += [2*g.real, -2*g.imag]          # 2 Re((x+iy) g) = 2x Re g - 2y Im g
        rows.append(row); rhs.append(a[k].real/sc)
    sol, _ = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
    c = [mp.mpc(sol[2*j], sol[2*j+1]) for j in range(J)]
    res = max(abs(sum(rows[i][n]*sol[n] for n in range(2*J)) - rhs[i]) for i in range(len(rhs)))
    return c, res

def Rsum(q, shift, dps):                 # R1 (shift 0): sum (-1)^n q^(n^2+n)/(q)_n ; R0 (shift 1): q^(n^2)
    mp.mp.dps = dps
    s = mp.mpc(0); P = mp.mpc(1); qn = mp.mpc(1); n = 0; lim = mp.mpf(10)**(-dps)
    while True:
        if n: qn *= q; P /= (1 - qn)
        e = n*n + n if shift == 0 else n*n
        term = (-1)**n * q**e * P
        s += term
        if n > 20 and abs(q)**e < lim * min(1, abs(s)) and abs(term) < lim*abs(s): break
        n += 1
    return s
def theta4(q):
    s = mp.mpc(1); n = 1
    while abs(q)**(n*n) > mp.mpf(10)**(-mp.mp.dps-5):
        s += 2*(-1)**n*q**(n*n); n += 1
    return s
def theta4_near_one(t):            # theta4(e^-t) = sqrt(pi/t) * 2 sum_{n>=0} e^(-pi^2 (n+1/2)^2/t)
    return mp.sqrt(mp.pi/t)*2*sum(mp.exp(-mp.pi**2*(n+half)**2/t) for n in range(6))
def rich(v, ratio=2):
    tab = [v[:]]
    for j in range(1, len(v)):
        p = tab[-1]; tab.append([(ratio**j*p[i+1]-p[i])/(ratio**j-1) for i in range(len(p)-1)])
    return tab[-1][0]

out = {}
print("=== q = -1 (control: addendum-v2 table) ===")
A = mp.pi**2/60
for key, K2pi, s1 in (("psi1@-1", (3*mp.sqrt(5)-5)/40, -(23+6*mp.sqrt(5))/60), ("psi0@-1", (5+3*mp.sqrt(5))/40, (13-6*mp.sqrt(5))/60)):
    c, res = fit_even(coef(key), A)
    print(" %s: K = c0 = %s   K^2 pi - table = %s   s1 - table = %s   fit residual %s" % (key, mp.nstr(c[0].real, 20),
          mp.nstr(c[0].real**2*mp.pi - K2pi, 3), mp.nstr((-c[1]/c[0]).real*(-1) - s1, 3) , mp.nstr(res, 3)))

print("=== q = i : singularity at -pi^2/240 ===")
A = mp.pi**2/240
fits = {}
for key in ("psi1@i", "psi0@i"):
    c, res = fit_even(coef(key), A)
    c2, res2 = fit_even(coef(key), A, J=10, k0=95)
    fits[key] = c
    print(" %s: c0 = %s  (change with other window: %s), fit residual %s" % (key, mp.nstr(c[0], 22), mp.nstr(abs(c[0]-c2[0]), 3), mp.nstr(res, 3)))
    print("      attached series 1 + s1 t + s2 t^2: s1 = %s, s2 = %s" % (mp.nstr(-c[1]/c[0], 18), mp.nstr(c[2]/c[0], 15)))
# direct asymptotics of (1/2) theta4 R1 and -(1/2) theta4 R0 at q = i e^(-t)
print(" direct: (1/2) theta4 R at q = i e^(-t), times sqrt(t) e^(-pi^2/(240 t)), Richardson to t = 0:")
for key, shift, sgn in (("psi1@i", 0, +1), ("psi0@i", 1, -1)):
    vals = []; ts = [mp.mpf(1)/(1500*2**j) for j in range(6)]
    for t in ts:
        q = 1j*mp.exp(-t)
        v = sgn*half*theta4(q)*Rsum(q, shift, 70)
        mp.mp.dps = 110
        vals.append(v*mp.sqrt(t)*mp.exp(-mp.pi**2/(240*t)))
    c0d = rich(vals)
    s1d = rich([(vals[i]/c0d - 1)/ts[i] for i in range(len(ts))])
    c = fits[key]
    print("   %s: limit = %s ;  pi * c0(large order) = %s ; ratio = %s" % (key, mp.nstr(c0d, 18), mp.nstr(mp.pi*c[0], 18), mp.nstr(mp.pi*c[0]/c0d, 15)))
    print("        first correction of theta4*R: %s ;  s1 from large order: %s" % (mp.nstr(s1d, 14), mp.nstr(-c[1]/c[0], 14)))
    out[key] = dict(c0=str(c[0]), direct=str(c0d))

print("=== q -> 1 : singularities at pi^2/3 -+ i Cl2(pi/3) ===")
S = mp.pi**2/3 - 1j*Cl2
w = mp.expjpi(-mp.mpf(2)/3)                       # saddle of F(-1): w^2 + w + 1 = 0
for key in ("psi1@1", "psi0@1"):
    c, res = fit_one(coef(key), S)
    c2, res2 = fit_one(coef(key), S, J=12, k0=95)
    cb, resb = fit_one(coef(key), S*(1+mp.mpf(10)**-12))
    print(" %s: c0 = %s  (change with other window: %s), fit residual %s; with S moved by 1e-12 relative: residual %s" % (key, mp.nstr(c[0], 20), mp.nstr(abs(c[0]-c2[0]), 3), mp.nstr(res, 3), mp.nstr(resb, 3)))
    print("      attached series: s1 = %s, s2 = %s" % (mp.nstr(c[1]/c[0], 16), mp.nstr(c[2]/c[0], 14)))
    fits[key] = c
# saddle prediction for the non-perturbative term of psi1* = F(-1) at q = e^(-t):  2 Re[ C t^(-1/2) e^(-S/t) ],
C_pred = (1+w)*mp.sqrt(mp.pi*w/(1+2*w))
print(" saddle prefactor C = (1+w) sqrt(pi w/(1+2w)), w = e^(-2 pi i/3):", mp.nstr(C_pred, 18))
for key in ("psi1@1",):
    for cand, nm in ((C_pred, "C"), (mp.conj(C_pred), "conj C")):
        print("   c0 * pi /", nm, "=", mp.nstr(fits[key][0]*mp.pi/cand, 15))
# direct: (1/2) theta4 R1 = psi1* - f1*/2 at q = e^(-t): fit 2 Re[ C t^-1/2 e^(-S/t) (1 + s1 t + ...) ]
print(" direct: (1/2) theta4 R1 at q = e^(-t), least squares for C, s1, s2 with S fixed:")
ts = [mp.mpf(1)/n for n in range(60, 141, 4)]
rows, rhs = [], []
for t in ts:
    q = mp.exp(-t)
    R = Rsum(q, 0, int(1.6/float(t)/2.302) + 60)
    mp.mp.dps = 110
    v = half*theta4_near_one(t)*R
    sc = mp.exp(-S.real/t)/mp.sqrt(t)
    row = []
    for j in range(6):
        g = t**j*mp.exp(-S/t)/mp.sqrt(t)/sc
        row += [2*g.real, -2*g.imag]
    rows.append(row); rhs.append(v.real/sc)
sol, _ = mp.qr_solve(mp.matrix(rows), mp.matrix(rhs))
Cd = mp.mpc(sol[0], sol[1]); s1d = mp.mpc(sol[2], sol[3])/Cd
print("   C(direct) =", mp.nstr(Cd, 14), "  saddle C =", mp.nstr(C_pred, 14), "  s1(direct) =", mp.nstr(s1d, 10), "  s1(large order) =", mp.nstr(fits["psi1@1"][1]/fits["psi1@1"][0], 10))
print("   2 pi c0(large order) / (i C(direct)) =", mp.nstr(2*mp.pi*fits["psi1@1"][0]/(1j*Cd), 12))
# closed forms
s5, s3 = mp.sqrt(5), mp.sqrt(3)
print("=== closed forms (PSLQ found them; differences shown) ===")
c = fits["psi1@i"]; print(" q=i  psi1*: pi c0^2 - (5-sqrt5)(1-2i)/80 =", mp.nstr(abs(mp.pi*c[0]**2 - (5-s5)*(1-2j)/80), 3), "  s1 - [(175+84 sqrt5)/300 - i(5+sqrt5)/25] =", mp.nstr(abs(-c[1]/c[0] - ((175+84*s5)/300 - 1j*(5+s5)/25)), 3))
c = fits["psi0@i"]; print(" q=i  psi0*: pi c0^2 + (5+sqrt5)(1+2i)/80 =", mp.nstr(abs(mp.pi*c[0]**2 + (5+s5)*(1+2j)/80), 3), "  s1 - [(84 sqrt5-125)/300 - i(5-sqrt5)/25] =", mp.nstr(abs(-c[1]/c[0] - ((84*s5-125)/300 - 1j*(5-s5)/25)), 3))
c = fits["psi1@1"]; print(" q->1 psi1*: c0 - e^(i pi/12)/(2 3^(1/4) sqrt(pi)) =", mp.nstr(abs(c[0] - mp.expjpi(mp.mpf(1)/12)/(2*mp.root(3, 4)*mp.sqrt(mp.pi))), 3), "  s1 - (1/4 - i sqrt3/54) =", mp.nstr(abs(c[1]/c[0] - (mp.mpf(1)/4 - 1j*s3/54)), 3))
c = fits["psi0@1"]; print(" q->1 psi0*: c0 - e^(-7 i pi/12)/(2 3^(1/4) sqrt(pi)) =", mp.nstr(abs(c[0] - mp.expjpi(-mp.mpf(7)/12)/(2*mp.root(3, 4)*mp.sqrt(mp.pi))), 3), "  s1 - (-1/12 - i sqrt3/54) =", mp.nstr(abs(c[1]/c[0] - (-mp.mpf(1)/12 - 1j*s3/54)), 3))
