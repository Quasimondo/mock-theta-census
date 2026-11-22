import json, math
from fractions import Fraction as Fr
from mpmath import mp, mpf, pi, gamma, matrix, lu_solve, sqrt, pslq, log, taylor
mp.dps = 120
d = json.load(open("exp_m1.json"))
A = pi**2/60
res = {}
for name in ("c0","c1"):
    c = [Fr(x) for x in d[name]]
    K = len(c)-1
    J = 24
    ks = list(range(K-J, K+1))
    M = matrix(len(ks), J+1); b = matrix(len(ks),1)
    for i,k in enumerate(ks):
        ak = mpf(c[k].numerator)/mpf(c[k].denominator)/mp.factorial(k)
        b[i] = ak * (-1)**k * A**(k+mpf(1)/2)/gamma(k+mpf(1)/2)
        for j in range(J+1):
            M[i,j] = A**j * gamma(k+mpf(1)/2-j)/gamma(k+mpf(1)/2)
    x = lu_solve(M,b)
    s = [x[j]/x[0] for j in range(J+1)]
    res[name] = (x[0], s)
    print(name, "K0 =", mp.nstr(x[0], 40))
    print("   s1 =", mp.nstr(s[1], 40))
    print("   s2 =", mp.nstr(s[2], 35))
    print("   s3 =", mp.nstr(s[3], 30))
r5 = sqrt(5)
for name in ("c0","c1"):
    K0, s = res[name]
    mp.dps = 30
    print(name)
    print("  pslq s1 [1, sqrt5]:", pslq([s[1], 1, r5], maxcoeff=10**6, maxsteps=10**5))
    print("  pslq s2 [1, sqrt5]:", pslq([s[2], 1, r5], maxcoeff=10**7, maxsteps=10**5))
    print("  pslq s3 [1, sqrt5]:", pslq([s[3], 1, r5], maxcoeff=10**9, maxsteps=10**5))
    print("  pslq K0^2 [1, sqrt5]:", pslq([K0**2, 1, r5], maxcoeff=10**6, maxsteps=10**5))
    print("  pslq K0^2*pi [1, sqrt5]:", pslq([K0**2*pi, 1, r5], maxcoeff=10**6, maxsteps=10**5))
    print("  pslq K0^2/pi [1, sqrt5]:", pslq([K0**2/pi, 1, r5], maxcoeff=10**6, maxsteps=10**5))
    mp.dps = 120
