"""The A = 1/2 family.  H(z;q) = sum_m z^m q^(m^2/4)/(q;q)_m.
q side (even/odd parts of H at z = i q^(1/2) and z = i):
    N0 = sum (-1)^n q^(n^2+n)/(q)_(2n)       N1 = sum (-1)^n q^(n^2+n)/(q)_(2n+1)
    X  = sum (-1)^n q^(n^2)/(q)_(2n)         Y  = sum (-1)^n q^((n+1)^2)/(q)_(2n+1)
dual side (H at z = -+ q^(1/4), -+ q^(3/4), written in the nome p with p^2 = q):
    Ep  = sum (-1)^m p^(m(m+1)/2)/(p^2;p^2)_m     Em  = sum (-1)^m p^(m(m+3)/2)/(p^2;p^2)_m
    Ep' = sum        p^(m(m+1)/2)/(p^2;p^2)_m     Em' = sum        p^(m(m+3)/2)/(p^2;p^2)_m
"""
import sys, os, math; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import acb, ctx
import horo

def qside(Xp, prec, target=110):
    """N0, N1, X, Y at q = e^(2 pi i/X), X = X0 - i b; midpoint arithmetic"""
    ctx.prec = prec
    Xc = horo.cx(*Xp)
    q = (2*acb.pi()*acb(0, 1)/Xc).exp().mid()
    K = float(abs(Xc).mid()); b = -float(Xc.imag.mid())
    J = int(K*math.sqrt((0.1615*K + 12 + target*0.6932)/(2*math.pi*b))) + 10
    one = acb(1)
    N0 = acb(0); N1 = acb(0); X = acb(0); Y = acb(0)
    sq = one       # q^(n^2)
    qn = one       # q^n
    pe = one       # (q)_(2n)
    q2n = one      # q^(2n)
    for n in range(J):
        po = (pe*(one - q2n*q)).mid()          # (q)_(2n+1)
        a = sq/pe; c = sq/po
        if n & 1:
            X -= a; N0 -= a*qn; N1 -= c*qn; Y -= c*q2n*q
        else:
            X += a; N0 += a*qn; N1 += c*qn; Y += c*q2n*q
        sq = (sq*q2n*q).mid(); qn = (qn*q).mid(); q2n = (q2n*q*q).mid()
        pe = (po*(one - q2n)).mid()
    return dict(N0=N0, N1=N1, X=X, Y=Y, q=q)

def dside(p, M=400):
    """Ep, Em, Ep', Em' at the nome p (|p| < 1)"""
    one = acb(1)
    Ep = acb(0); Em = acb(0); Epp = acb(0); Emp = acb(0)
    tri = one; pm = one; den = one
    for m in range(M):
        a = tri/den
        Epp += a; Emp += a*pm
        if m & 1: Ep -= a; Em -= a*pm
        else: Ep += a; Em += a*pm
        pm = pm*p; tri = tri*pm; den = den*(one - pm*pm)
    return dict(Ep=Ep, Em=Em, Epp=Epp, Emp=Emp)

if __name__ == "__main__":
    from fractions import Fraction as Fr
    I = acb(0, 1)
    ctx.prec = 300
    p = acb(3)/10 + I/7
    d = dside(p)
    print("Casoratian  Ep' Em + Em' Ep - 2 at p = 0.3 + i/7:", abs((d["Epp"]*d["Em"] + d["Emp"]*d["Ep"] - 2).mid()).str(5, radius=False))
    for X0 in (Fr(3, 10), Fr(11, 10)):
        for m in (100, 200):
            Xp = (X0 + 2*m, Fr(2, 5)); K = float(X0 + 2*m)
            prec = int(0.9*K) + 400
            r = qside(Xp, prec, target=int(0.5*K) + 150)
            Xc = horo.cx(*Xp)
            p = (-acb.pi()*I*Xc).exp()
            d = dside(p)
            q4 = (2*acb.pi()*I/Xc/4).exp()
            f = lambda z: float(abs(z.mid()).str(5, radius=False))
            print(f"X0 = {X0}, K = {K}:  |N0| = {f(r['N0']):.2e}")
            for s, lab in ((1, "+1"), (-1, "-1"), (I, "+i"), (-I, "-i")):
                T1 = r["N0"]*d["Epp"] + s*I/q4*r["Y"]*d["Ep"]
                T2 = r["X"]*d["Ep"] + s*I*q4*r["N1"]*d["Epp"]
                T3 = r["N0"]*d["Emp"] + s*I/q4*r["Y"]*d["Em"]
                T4 = r["X"]*d["Em"] + s*I*q4*r["N1"]*d["Emp"]
                print(f"    s = {lab}:  |N0 Ep' + s i q^(-1/4) Y Ep| = {f(T1):.2e}   |X Ep + s i q^(1/4) N1 Ep'| = {f(T2):.2e}   |N0 Em' + s.. Y Em| = {f(T3):.2e}   |X Em + s.. N1 Em'| = {f(T4):.2e}")
