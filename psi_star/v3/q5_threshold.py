"""Q5: growth of F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n at roots of odd order, from a saddle point.
 V(x) = 2 pi i a x - x^2/2 + Li2(-e^(-x)) + pi^2/12 ,   saddle: w^2 + w = e^(-2 pi i a), w = e^(-x)
 leading term:  F(e^(2 pi i a); e^(-t)) ~ (1+w) sqrt(pi w / ((1+2w) t)) e^(V/t)
 At a root of odd order m the same holds with a = arg(z^m)/(2 pi), t -> m^2 t (kappa = Re V / (4 pi^2)).
"""
import json, sys, math
import mpmath as mp
mp.mp.dps = 30
def saddle(a):
    a = mp.mpf(a)
    if a == mp.mpf(1)/2: a -= mp.mpf(10)**(-mp.mp.dps+5)      # stay on the branch continued from a = 0
    w = (-1 + mp.sqrt(1 + 4*mp.expjpi(-2*a)))/2
    x = -mp.log(w)
    V = 2j*mp.pi*a*x - x*x/2 + mp.polylog(2, -w) + mp.pi**2/12
    return w, x, V
ReV = lambda a: mp.re(saddle(a)[2])
astar = mp.findroot(ReV, 0.107)
print("threshold: Re V = 0 at a* =", mp.nstr(astar, 20), "= 1/%s of a turn" % mp.nstr(1/astar, 10))
print("Re V(0) - pi^2/60 =", mp.nstr(ReV(0) - mp.pi**2/60, 3))

# ---- comparison with census3/zfam_results.jsonl (path as argument)
if len(sys.argv) > 1:
    rows = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
    def cls(r):
        if "error" in r: return "?"
        if r["resid"] < 1e-25 and r["kappa_err"] < 1e-20: return "E"
        if abs(r["kappa"]) < 1e-4: return "Q"
        return "X"
    from fractions import Fraction
    stat = {}
    worst = {}
    for r in rows:
        if not r["name"].startswith("F["): continue
        k, Nn = map(int, r["name"][2:-1].split("/"))
        m = r["m"]
        if m % 2 == 0: continue
        a = Fraction(m*k, Nn) % 1; a = min(a, 1-a)
        c = cls(r)
        key = str(a)
        stat.setdefault(key, {"a": a, "cls": {}, "kap": []})
        stat[key]["cls"][c] = stat[key]["cls"].get(c, 0) + 1
        if c in "EX" and "kappa" in r: stat[key]["kap"].append(r["kappa"])
    print("\nF(z) at roots of odd order m in zfam_results.jsonl, by a = distance of arg(z^m)/2pi from an integer")
    print("  a        predicted kappa = max(Re V,0)/(4 pi^2)   classes found        fitted kappa (range over rows with growth)")
    for key in sorted(stat, key=lambda s: stat[s]["a"]):
        a = stat[key]["a"]; pv = float(ReV(mp.mpf(a.numerator)/a.denominator))/(4*math.pi**2)
        kap = stat[key]["kap"]
        print("  %-7s  %-38s  %-20s %s" % (key, ("%.10f" % pv) if pv > 0 else "bounded (Re V = %.4f)" % (pv*4*math.pi**2),
              " ".join(f"{c}:{n}" for c, n in sorted(stat[key]["cls"].items())),
              ("%.10f .. %.10f" % (min(kap), max(kap))) if kap else ""))

# ---- high-precision check of the leading term at q -> 1
def F(z, t, dps):
    mp.mp.dps = dps
    q = mp.exp(-t); s = mp.mpc(0); P = mp.mpf(1); qn = mp.mpf(1); zn = mp.mpc(1); n = 0
    lim = mp.mpf(10)**(-dps)
    while True:
        if n: qn *= q; P *= 1 + qn; zn *= z
        term = zn * q**(n*(n+1)//2) * P
        s += term
        if n*t > 3 and abs(term) < lim: break
        n += 1
    return s
print("\ncheck of V and of the prefactor at q = e^(-t)  (Richardson in t)")
for an, ad in ((0, 1), (1, 12)):
    a = mp.mpf(an)/ad
    vals = []; ts = [mp.mpf(1)/(800*2**j) for j in range(5)]
    for t in ts:
        dps = int(0.1645/float(t)/2.302 - float(ReV(a))/float(t)/2.302) + 45
        Fv = F(mp.expjpi(2*a), t, dps)
        mp.mp.dps = 40
        w, x, V = saddle(a)
        pred = (1+w)*mp.sqrt(mp.pi*w/((1+2*w)*t))*mp.exp(V/t)
        vals.append(Fv/pred)
    mp.mp.dps = 40
    tab = [vals[:]]
    for j in range(1, len(vals)):
        p = tab[-1]; tab.append([(2**j*p[i+1]-p[i])/(2**j-1) for i in range(len(p)-1)])
    print("  a = %d/%d: F/prediction at t = 1/800 .. 1/12800:" % (an, ad), [mp.nstr(v, 8) for v in vals])
    print("           extrapolated to t = 0:", mp.nstr(tab[-1][0], 18), "  |limit - 1| =", mp.nstr(abs(tab[-1][0]-1), 3))

# ---- bounded case: radial limit = (sum of the first m terms) / (1 - 2 z^m)   (m odd; for z = -1 this is the 1/3 rule)
print("\nlimit rule in the bounded case, F(z; zeta e^(-t)) -> A_z/(1 - 2 z^m):")
def Froot(z, zeta, t, dps):
    mp.mp.dps = dps
    q = zeta*mp.exp(-t); s = mp.mpc(0); P = mp.mpc(1); qn = mp.mpc(1); zn = mp.mpc(1); n = 0
    lim = mp.mpf(10)**(-dps)
    while True:
        if n: qn *= q; P *= 1 + qn; zn *= z
        term = zn * q**(n*(n+1)//2) * P
        s += term
        if n*n*t > 40 and abs(term) < lim: break
        n += 1
    return s
for (zn, zd, h, m) in ((1, 4, 0, 1), (1, 5, 0, 1), (1, 5, 1, 3), (1, 6, 1, 5), (1, 2, 2, 7)):
    z = mp.expjpi(mp.mpf(2*zn)/zd); zeta = mp.expjpi(mp.mpf(2*h)/m)
    A = 0; P = mp.mpc(1)
    for n in range(m):
        if n: P *= 1 + zeta**n
        A += z**n * zeta**(n*(n+1)//2) * P
    pred = A/(1 - 2*z**m)
    ts = [mp.mpf(1)/(4*2**j)/m**2 for j in range(8)]
    vals = [Froot(z, zeta, t, 260) for t in ts]
    mp.mp.dps = 40
    tab = [vals[:]]
    for j in range(1, len(vals)):
        p_ = tab[-1]; tab.append([(2**j*p_[i+1]-p_[i])/(2**j-1) for i in range(len(p_)-1)])
    print("  z = e^(2 pi i %d/%d), root %d/%d: |extrapolated limit - prediction| = %s" % (zn, zd, h, m, mp.nstr(abs(tab[-1][0]-pred), 3)))
