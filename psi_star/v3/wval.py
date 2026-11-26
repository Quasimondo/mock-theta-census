"""W(xi) = sum_j (-1)^j xi^(j^2+j) (xi;xi^2)_j at a root of unity xi = e^(2 pi i p/n).
 n odd: the sum terminates.  n even: the terms repeat with a constant factor after n/2 steps... (period found numerically)
 and the sum is taken as the geometric series (as in the odd-root rule for psi*)."""
import mpmath as mp
from math import gcd
def W(p, n, shift=0):
    """shift 0: exponent j^2+j ; other exponents via (a2, b2) = shift tuple meaning (a2 j^2 + b2 j)/2 with xi^(1/2)"""
    g = gcd(p, n); p //= g; n //= g
    xi = mp.expjpi(mp.mpf(2*p)/n)
    terms = []; P = mp.mpc(1)
    per = n if n % 2 == 0 else None
    jmax = (n+1)//2 if n % 2 else 2*n
    xo = xi                     # xi^(2j-1)
    e = mp.mpc(1)               # xi^(j^2+j)
    x2 = xi*xi
    for j in range(jmax+1):
        if j:
            P *= (1 - xo); xo *= x2
            e *= xi**(2*j)      # (j^2+j) - ((j-1)^2+(j-1)) = 2j
        terms.append((-1)**j * e * P)
    if n % 2: return sum(terms)
    # even order: find the period
    for per in (n//2, n):
        lam = terms[per]/terms[0]
        if all(abs(terms[per+i] - lam*terms[i]) < mp.mpf(10)**(-mp.mp.dps+8)*(abs(terms[i])+1) for i in range(1, min(per, 6))):
            return sum(terms[:per])/(1-lam)
    raise ValueError("no period")
if __name__ == "__main__":
    mp.mp.dps = 30
    for (p, n) in ((0, 1), (1, 2), (1, 4), (3, 4), (1, 3), (1, 6)):
        print(p, n, mp.nstr(W(p, n), 15))
