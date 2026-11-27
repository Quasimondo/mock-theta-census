"""W(xi) = sum_j (-1)^j xi^(j^2+j) (xi;xi^2)_j at xi = e^(2 pi i p/n) with ball arithmetic (python-flint).
 n odd: terminating sum.  n = 2s even: the terms repeat with factor -2 after s steps, value = (1/3) * (first s terms).
 Also the partner W3 (exponent j^2+3j)."""
from flint import acb, arb, ctx
from math import gcd
def W_both(p, n, prec):
    g = gcd(p, n); p //= g; n //= g
    ctx.prec = prec
    xi = (acb(2*p)/n).exp_pi_i()
    x2 = xi*xi
    last = (n+1)//2 if n % 2 else n//2          # number of terms summed
    s1 = acb(0); s3 = acb(0); P = acb(1)
    xo = xi                                     # xi^(2j-1), starting at j = 1
    e = acb(1)                                  # xi^(j^2+j)
    x2j = acb(1)                                # xi^(2j)
    for j in range(last):
        if j:
            P *= (1 - xo); xo *= x2
            x2j *= x2; e *= x2j
        t = e*P if j % 2 == 0 else -(e*P)
        s1 += t
        s3 += t*x2j                             # exponent j^2+3j = (j^2+j) + 2j
    if n % 2 == 0:
        s1 /= 3; s3 /= 3
    return s1, s3
if __name__ == "__main__":
    import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import mpmath as mp
    from wval import W
    mp.mp.dps = 30
    for (p, n) in ((1, 3), (1, 4), (1, 6), (3, 10), (2, 7)):
        a, b = W_both(p, n, 200)
        print(p, n, a.real.mid().str(15, radius=False), a.imag.mid().str(15, radius=False), "| mp:", mp.nstr(W(p, n), 15))
