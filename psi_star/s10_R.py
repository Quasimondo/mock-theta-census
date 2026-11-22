from mpmath import mp, mpf, pi, exp, sqrt, matrix, lu_solve, log
mp.dps = 700
r5 = sqrt(5)
def R(t, lin):
    q = -exp(-t)
    tot = mpf(0); P = mpf(1); n = 0
    while True:
        if n > 0: P *= (1 - q**n)
        term = (-1)**n * q**(n*n + lin*n) / P
        tot += term
        if n > 20 and abs(term) < mpf(10)**(-680) : break
        n += 1
    return tot
A = pi**2/60
ts = [mpf(1)/k for k in range(60, 84)]
for lin, name, s1, s2, K0sq in ((0,"R0",(13-6*r5)/60,(781-300*r5)/7200,(5+3*r5)/40),(1,"R1",(-23-6*r5)/60,(1141+420*r5)/7200,(-5+3*r5)/40)):
    vals = [R(t, lin)*exp(-A/t) for t in ts]
    J = len(ts)
    M = matrix(J,J); b = matrix(J,1)
    for i,t in enumerate(ts):
        b[i] = vals[i]
        for j in range(J): M[i,j] = t**j
    x = lu_solve(M,b)
    print(name, "rho =", mp.nstr(x[0],25), " rho^2 =", mp.nstr(x[0]**2, 25))
    print("   x1/x0 =", mp.nstr(x[1]/x[0],18), "  -s1 =", mp.nstr(-s1,18))
    print("   x2/x0 =", mp.nstr(x[2]/x[0],15), "   s2 =", mp.nstr(s2,15))
    K0 = sqrt(K0sq/pi)
    print("   rho/|K0| =", mp.nstr(x[0]/K0, 20), "   (rho/K0)^2 =", mp.nstr((x[0]/K0)**2, 20), "  /pi:", mp.nstr((x[0]/K0)**2/pi, 20))
