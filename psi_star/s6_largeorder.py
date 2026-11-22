import json, math
from fractions import Fraction as Fr
from mpmath import mp, mpf, pi, gamma, matrix, lu_solve, identify, sqrt, exp, taylor
mp.dps = 80
d = json.load(open("exp_m1.json"))
A = pi**2/60   # |A|; series alternates
for name in ("c0","c1"):
    c = [Fr(x) for x in d[name]]
    K = len(c)-1
    J = 10
    ks = list(range(K-J, K+1))
    # model: a_k (-1)^k A^{k+1/2}/Gamma(k+1/2) = sum_j x_j * (A)^j * Gamma(k+1/2-j)/Gamma(k+1/2)   (signs absorbed in x_j)
    M = matrix(len(ks), J+1); b = matrix(len(ks),1)
    for i,k in enumerate(ks):
        ak = mpf(c[k].numerator)/mpf(c[k].denominator)/mp.factorial(k)
        b[i] = ak * (-1)**k * A**(k+mpf(1)/2)/gamma(k+mpf(1)/2)
        for j in range(J+1):
            M[i,j] = A**j * gamma(k+mpf(1)/2-j)/gamma(k+mpf(1)/2)
    x = lu_solve(M,b)
    print(name, "K0 =", mp.nstr(x[0], 25))
    s = [x[j]/x[0] for j in range(J+1)]
    print("   s_j:", [mp.nstr(v, 18) for v in s[:7]])
    # compare with pure exponential e^{s1 t}
    print("   s2/(s1^2/2) =", mp.nstr(s[2]/(s[1]**2/2), 15), "  s3/(s1^3/6) =", mp.nstr(s[3]/(s[1]**3/6), 12))
    print("   s1*240 =", mp.nstr(s[1]*240, 15), " s1*144=", mp.nstr(s[1]*144,15))
    print("   log-series of S: ", [mp.nstr(v,15) for v in taylor(lambda t: mp.log(sum(s[j]*t**j for j in range(7))), 0, 5)])
    print("   K0^2 =", mp.nstr(x[0]**2, 20), " K0*sqrt(pi)=", mp.nstr(x[0]*sqrt(pi),20), " K0/sqrt(pi)=", mp.nstr(x[0]/sqrt(pi),20))
