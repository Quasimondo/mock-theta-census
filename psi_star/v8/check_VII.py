# Independent check of addendum-v8: identity (VII) and the integer coefficients quoted in secs. 1, 2 (variable r)
exec(open('/home/claude/v4check/check_V_VI.py').read().split("W  = tot(")[0])
J  = tot([sh(poch(2, 2, +1, m, True), 2*m*m+m) for m in range(14) if 2*m*m+m < N])          # sum r^{2m^2+m}/(-r^2;r^2)_m
K1 = tot([sc(sh(poch(2, 2, +1, n), (n+1)**2), 2) for n in range(19) if (n+1)**2 < N])          # 2 sum r^{(n+1)^2}(-r^2;r^2)_n
def m1_r2(n):  # (-1;r^2)_n = 2 (-r^2;r^2)_{n-1} for n >= 1
    return sc(poch(2, 2, +1, n-1), 2)
K2 = tot([sh(m1_r2(n), n*n) for n in range(1, 19) if n*n < N])                                  # sum_{n>=1} r^{n^2}(-1;r^2)_n
psi = [0]*N
for n in range(30):
    if n*(n+1)//2 < N: psi[n*(n+1)//2] = 1
Epp = tot([sh(poch(2, 2, -1, m, True), m*(m+1)//2) for m in range(30) if m*(m+1)//2 < N])      # sum r^{m(m+1)/2}/(r^2;r^2)_m
calE = mul(mul(psi, poch(2, 4, -1)), Epp)                                                          # psi(r) (r^2;r^4)_inf E+'(r)
print('two forms of K agree:', K1 == K2)
print('(VII) J + K = psi (r^2;r^4) E+\' to r^%d:' % (N-1), [a+b for a, b in zip(J, K1)] == calE)
print('E(r) coefficients:', calE[:14], ' claimed 1, 2, 0, 1, 2, -1, 2, 1, 0, 1, 1, 3, -1, 1')
print('J coefficients:   ', J[:14], ' claimed 1, 0, 0, 1, 0, -1, 0, 1, 0, -1, 1, 1, -1, -1')
print('K coefficients:   ', K1[:14], ' claimed 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 2')
