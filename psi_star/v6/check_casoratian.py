# Independent check of addendum-v6 sec. 1: the reformulations N0, N1, E+-, and the two Casoratian identities,
# exact integer q-series to q^(N-1).
exec(open('/home/claude/v4check/check_V_VI.py').read().split("W  = tot(")[0])
R0 = tot([sc(sh(poch(1, 1, -1, n, True), n*n), (-1)**n) for n in range(20) if n*n < N])
R1 = tot([sc(sh(poch(1, 1, -1, n, True), n*n+n), (-1)**n) for n in range(20) if n*n+n < N])
N0 = tot([sc(sh(poch(1, 1, -1, 2*n, True), n*n+n), (-1)**n) for n in range(20) if n*n+n < N])
N1 = tot([sc(sh(poch(1, 1, -1, 2*n+1, True), n*n+n), (-1)**n) for n in range(20) if n*n+n < N])
X  = tot([sc(sh(poch(1, 1, -1, 2*n, True), n*n), (-1)**n) for n in range(20) if n*n < N])
Y  = tot([sc(sh(poch(1, 1, -1, 2*n+1, True), (n+1)**2), (-1)**n) for n in range(20) if (n+1)**2 < N])
Pp = tot([sh(poch(2, 2, -1, k, True), 2*k*k+k) for k in range(14) if 2*k*k+k < N])
Pm = tot([sh(poch(2, 2, -1, k, True), 2*k*k-k) for k in range(14) if 2*k*k-k < N])
def E(sign, a):   # sum sign^m q^{m(m+a)/2}/(q^2;q^2)_m
    return tot([sc(sh(poch(2, 2, -1, m, True), m*(m+a)//2), sign**m) for m in range(30) if m*(m+a)//2 < N])
Ep, Em, Epp, Emp = E(-1, 1), E(-1, 3), E(1, 1), E(1, 3)
mq = poch(1, 1, +1)            # (-q;q)_inf
q_q2 = poch(1, 2, -1)          # (q;q^2)_inf
one = [1] + [0]*(N-1); two = [2] + [0]*(N-1)
print('N0 = (-q;q)R0:', mul(mq, R0) == N0, '  N1 = (-q;q)R1:', mul(mq, R1) == N1)
print('E+ = (q;q^2)P+:', mul(q_q2, Pp) == Ep, '  E- = (q;q^2)P-:', mul(q_q2, Pm) == Em)
print('X N0 + N1 Y = 1:', [a+b for a, b in zip(mul(X, N0), mul(N1, Y))] == one)
print("E+' E- + E-' E+ = 2:", [a+b for a, b in zip(mul(Epp, Em), mul(Emp, Ep))] == two)
