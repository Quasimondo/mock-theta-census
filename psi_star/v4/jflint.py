"""psi1*, psi0* at a root of unity e^(2 pi i p/n) with ball arithmetic (python-flint).
 n even: terminating sum;  n odd: one third of the first n terms (odd-root rule, addendum v3 section 1)."""
from flint import acb, ctx
from math import gcd
def J_both(p, n, prec):
    g = gcd(p, n); p //= g; n //= g
    if n < 0: p, n = -p, -n
    ctx.prec = prec
    z = (acb(2*p)/n).exp_pi_i()
    s1 = acb(0); s0 = acb(0); P = acb(1); zn = acb(1); e = acb(1)
    last = n if n % 2 else n//2
    for k in range(last):
        if k:
            zn *= z; P *= 1 + zn; e *= zn
        t = e*P if k % 2 == 0 else -(e*P)
        s1 += t
        s0 += t*zn*z
    if n % 2:
        s1 /= 3; s0 /= 3
    return s1, s0
