"""Exact values of psi1*, psi0* at roots of unity e^(2 pi i p/n):
 n even: the series terminates;  n odd: one third of the first n terms (odd-root rule, section Q3)."""
import mpmath as mp
from math import gcd
def J(p, n, dps=None):
    """returns (psi1*(zeta), psi0*(zeta)), zeta = e^(2 pi i p/n), p/n in lowest terms"""
    g = gcd(p, n); p //= g; n //= g
    if n < 0: p, n = -p, -n
    z = mp.expjpi(mp.mpf(2*p)/n)
    s1 = s0 = mp.mpc(0); P = mp.mpc(1); zn = mp.mpc(1)
    # z^(k(k+1)/2) by recurrence: e_k = e_{k-1} * z^k
    e = mp.mpc(1)            # z^(k(k+1)/2) at k = 0
    last = n if n % 2 else n//2
    for k in range(last):
        if k:
            zn *= z; P *= 1 + zn; e *= zn
        sg = -1 if k % 2 else 1
        s1 += sg*e*P
        s0 += sg*e*zn*z*P        # z^((k+1)(k+2)/2) = z^(k(k+1)/2) * z^(k+1)
    if n % 2: s1 /= 3; s0 /= 3
    return s1, s0
if __name__ == "__main__":
    mp.mp.dps = 30
    print("psi1*(-1), psi1*(i), psi1*(e^(i pi/3)):", J(1, 2)[0], J(1, 4)[0], J(1, 6)[0])
    print("psi1*(1) =", J(0, 1)[0], " psi1*(e^(2 pi i/5)) =", J(1, 5)[0])
