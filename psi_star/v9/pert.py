"""Perturbative series at q = e^(-t), exact rational coefficients, cached in pert.json.
  W(e^-t)  = sum a_k t^k,   W3(e^-t) = sum a3_k t^k      (W = G(-q), W3 = G(-q^3))"""
import json, os
from flint import fmpq_series, fmpq, ctx
HERE = os.path.dirname(os.path.abspath(__file__))

def series(N=60):
    path = os.path.join(HERE, "pert.json")
    if os.path.exists(path):
        d = json.load(open(path))
        if len(d["W"]) >= N + 1: return d
    ctx.cap = N + 1
    t = fmpq_series([0, 1])
    def E(m): return (-m*t).exp()
    W = fmpq_series([0]); W3 = fmpq_series([0]); poch = fmpq_series([1])
    for j in range(N + 1):
        s = -1 if j & 1 else 1
        W += s*E(j*j + j)*poch; W3 += s*E(j*j + 3*j)*poch
        poch = poch*(1 - E(2*j + 1))
    d = {"W": [str(W[k]) for k in range(N + 1)], "W3": [str(W3[k]) for k in range(N + 1)]}
    json.dump(d, open(path, "w"))
    return d

def series_m1(N=60):
    """psi1*(-e^-t), psi0*(-e^-t) as power series in t (exact rationals), cached in pert_m1.json"""
    path = os.path.join(HERE, "pert_m1.json")
    if os.path.exists(path):
        d = json.load(open(path))
        if len(d["psi1"]) >= N + 1: return d
    ctx.cap = N + 1
    t = fmpq_series([0, 1])
    def E(m): return (-m*t).exp()
    p1 = fmpq_series([0]); p0 = fmpq_series([0]); poch = fmpq_series([1])
    for n in range(2*N + 3):
        T = n*(n + 1)//2; T0 = (n + 1)*(n + 2)//2
        p1 += (-1)**(n + T)*E(T)*poch; p0 += (-1)**(n + T0)*E(T0)*poch
        i = n + 1
        poch = poch*(1 + (-1)**i*E(i))
    d = {"psi1": [str(p1[k]) for k in range(N + 1)], "psi0": [str(p0[k]) for k in range(N + 1)]}
    json.dump(d, open(path, "w"))
    return d

if __name__ == "__main__":
    dm = series_m1()
    print("psi1*(-e^-t):", dm["psi1"][:6]); print("psi0*(-e^-t):", dm["psi0"][:6])
    import math
    from fractions import Fraction as Fr
    d = series()
    print("k! a_k  (W):", [int(Fr(d["W"][k])*math.factorial(k)) for k in range(8)])
    print("k! a3_k (W3):", [int(Fr(d["W3"][k])*math.factorial(k)) for k in range(8)])
