"""Exact formal expansion of psi1*, psi0* at q = -e^{-t} (the series terminate order by order).
psi(-e^{-t}) = sum_k c_k t^k / k!,  c_k rational.
"""
from flint import fmpq_series, fmpq, ctx
import sys, json

K = int(sys.argv[1]) if len(sys.argv) > 1 else 60
ctx.cap = K + 1


def E(a):
    # exp(-a t) as series
    x = fmpq_series([0, -a], prec=K + 1)
    return x.exp()


def expand(shift):
    """sum_n (-1)^n q^{(n+shift)(n+shift+1)/2} (-q;q)_n at q=-e^{-t}; shift=0: psi1*, shift=1: psi0*"""
    tot = fmpq_series([0], prec=K + 1)
    P = fmpq_series([1], prec=K + 1)
    for n in range(0, 2 * K + 3):
        if n > 0:
            P = P * (1 + (-1) ** n * E(n))
        e = (n + shift) * (n + shift + 1) // 2
        tot = tot + ((-1) ** n) * ((-1) ** e) * E(e) * P
    co = tot.coeffs()
    co = co + [fmpq(0)] * (K + 1 - len(co))
    f = fmpq(1)
    out = []
    for k in range(K + 1):
        if k:
            f *= k
        out.append(co[k] * f)
    return out


c1 = expand(0)
c0 = expand(1)
print("psi1*(-e^-t): c_k =", [str(x) for x in c1[:12]])
print("psi0*(-e^-t): c_k =", [str(x) for x in c0[:12]])
json.dump({"c1": [str(x) for x in c1], "c0": [str(x) for x in c0]}, open("exp_m1.json", "w"))
import math
for name, c in (("psi1*", c1), ("psi0*", c0)):
    print(name, "ratio c_{k+1}/c_k / ((2k+1)(2k+2)) :")
    print([round(float(c[k + 1] / c[k]) / ((2 * k + 1) * (2 * k + 2)), 6) if c[k] != 0 else None for k in range(K - 12, K)])
