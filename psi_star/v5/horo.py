"""Evaluation inside the unit disk, close to q = 1 along a horocycle.
q = exp(2 pi i / Y) with Im Y = -b < 0, so |q| < 1 and the dual nome q' = exp(-pi i Y) has |q'| = exp(-pi b) fixed.
near(Y, prec)  -> dict with W, W3, Fa, Fb, Pp, Pm at q   (flint balls)
dual(Y, prec)  -> dict with psi1, psi0, f1, f0, R1, R0, th4 at q'
"""
from flint import acb, arb, fmpq, ctx
import math

def cx(Y0, b):
    return acb(arb(fmpq(Y0.numerator, Y0.denominator)), -arb(fmpq(b.numerator, b.denominator)))

def near(Yp, prec, target=110):
    ctx.prec = prec
    Y = cx(*Yp)
    q = (2*acb.pi()*acb(0, 1)/Y).exp()
    absY = float(abs(Y).mid()); b = -float(Y.imag.mid())
    J = int(absY*math.sqrt((0.0808*absY + 12 + target*0.6932)/(2*math.pi*b))) + 10
    return near_q(q, J)

def near_q(q, J):
    """W, W3, Fa, Fb, Pp, Pm at the nome q (|q| < 1), J terms"""
    q = q.mid(); q2 = (q*q).mid()
    one = acb(1)
    W = acb(0); W3 = acb(0); Fa = acb(0); Fb = acb(0); Pp = acb(0); Pm = acb(0)
    poch = one              # (q;q^2)_j
    poch2 = one             # (q^2;q^2)_j
    e1 = one                # q^(j^2+j)
    qj2 = one               # q^(2j)
    g = one                 # q^(2j^2-j)
    qodd = q                # q^(2j+1)
    for j in range(J):
        t = e1*poch
        if j & 1:
            W -= t; W3 -= t*qj2
        else:
            W += t; W3 += t*qj2
        if j >= 1:
            u = g/poch
            Fa += u; Fb += u*q/qj2
        v = g/poch2
        Pm += v; Pp += v*qj2
        # advance j -> j+1
        # midpoints only: ball radii would blow up through the wrapping effect; accuracy is checked by rerunning at higher precision
        poch = (poch*(one - qodd)).mid()
        e1 = (e1*qj2*q2).mid()
        g = (g*qj2*qj2*q).mid()
        qj2 = (qj2*q2).mid()
        poch2 = (poch2*(one - qj2)).mid()
        qodd = (qodd*q2).mid()
    return dict(W=W, W3=W3, Fa=Fa, Fb=Fb, Pp=Pp, Pm=Pm, q=q)

def fwd(Xp, prec, target=110):
    """psi1*, psi0*, f1*, f0* at q = e^(2 pi i/X), X = X0 - i b (midpoint arithmetic)"""
    ctx.prec = prec
    X = cx(*Xp)
    q = (2*acb.pi()*acb(0, 1)/X).exp().mid()
    K = float(abs(X).mid()); b = -float(X.imag.mid())
    J = int(K*math.sqrt((0.1615*K + 12 + target*0.6932)/(math.pi*b))) + 10
    one = acb(1)
    psi1 = acb(0); psi0 = acb(0); f1 = acb(0); f0 = acb(0)
    pp = one; qn = one; tri = one; sq = one
    for n in range(J):
        a = tri*pp; c = sq/pp
        if n & 1:
            psi1 -= a; psi0 -= a*qn*q; f0 -= c; f1 -= c*qn
        else:
            psi1 += a; psi0 += a*qn*q; f0 += c; f1 += c*qn
        sq = (sq*qn*qn*q).mid()
        qn = (qn*q).mid()
        tri = (tri*qn).mid()
        pp = (pp*(one + qn)).mid()
    return dict(psi1=psi1, psi0=psi0, f1=f1, f0=f0, q=q)

def mfwd(Zp, prec, target=110):
    """psi1*, psi0* at q = -e^(2 pi i/Z), Z = Z0 - i b (midpoint arithmetic)"""
    ctx.prec = prec
    Z = cx(*Zp)
    q = (-(2*acb.pi()*acb(0, 1)/Z).exp()).mid()
    K = float(abs(Z).mid()); b = -float(Z.imag.mid())
    J = int(K*math.sqrt((0.1*K + 20 + target*0.6932)/(math.pi*b))) + 10
    one = acb(1)
    psi1 = acb(0); psi0 = acb(0)
    pp = one; qn = one; tri = one
    for n in range(J):
        a = tri*pp
        if n & 1:
            psi1 -= a; psi0 -= a*qn*q
        else:
            psi1 += a; psi0 += a*qn*q
        qn = (qn*q).mid()
        tri = (tri*qn).mid()
        pp = (pp*(one + qn)).mid()
    return dict(psi1=psi1, psi0=psi0, q=q)

def dual(Yp, prec=400):
    ctx.prec = prec
    Y = cx(*Yp)
    qd = (-acb.pi()*acb(0, 1)*Y).exp()
    return dual_q(qd, prec)

def dual_q(qd, prec=400):
    """psi1, psi0, f1, f0, R1, R0, th4 at the nome qd (|qd| < 1)"""
    one = acb(1)
    psi1 = acb(0); psi0 = acb(0); f1 = acb(0); f0 = acb(0); R1 = acb(0); R0 = acb(0); th4 = one
    pp = one   # (-q;q)_n
    pm = one   # (q;q)_n
    qn = one   # q^n
    tri = one  # q^(n(n+1)/2)
    sq = one   # q^(n^2)
    n = 0
    while True:
        s = -1 if n & 1 else 1
        psi1 += s*tri*pp
        psi0 += s*tri*qn*qd*pp          # q^((n+1)(n+2)/2) = q^(n(n+1)/2) q^(n+1)
        f0 += s*sq/pp; f1 += s*sq*qn/pp
        R0 += s*sq/pm; R1 += s*sq*qn/pm
        if n >= 1: th4 += 2*s*sq
        if float(abs(tri).mid()) < 2.0**(-prec + 20) and n > 5: break
        n += 1
        sq = sq*qn*qn*qd
        qn = qn*qd
        tri = tri*qn
        pp = pp*(one + qn); pm = pm*(one - qn)
    return dict(psi1=psi1, psi0=psi0, f1=f1, f0=f0, R1=R1, R0=R0, th4=th4, qd=qd)

if __name__ == "__main__":
    from fractions import Fraction as Fr
    Yp = (Fr(240) + Fr(3, 10), Fr(1, 2))
    for prec in (300, 400, 600):
        r = near(Yp, prec)
        d = dual(Yp, prec)
        tau = -2*acb.pi()*acb(0, 1)/cx(*Yp)
        psi = (tau/8).exp()*(acb.pi()/(2*tau)).sqrt()*d["th4"]     # psi(q) = q^(-1/8) sqrt(pi/(2 tau)) theta4(q')
        print(prec, "W =", r["W"].mid().str(40, radius=False), " Fa =", r["Fa"].mid().str(40, radius=False))
        print("  (V):", abs((r["W"] + r["Fa"] - psi*r["Pp"]).mid()).str(5, radius=False), " (VI):", abs((r["q"]*r["W3"] + r["Fb"] - psi*r["Pm"]).mid()).str(5, radius=False))
    dd = dual((Fr(3, 10), Fr(1, 2)))
    print(" (I):", abs(2*dd["psi0"] - dd["f0"] + dd["th4"]*dd["R0"]), " (II):", abs(2*dd["psi1"] - dd["f1"] - dd["th4"]*dd["R1"]))
