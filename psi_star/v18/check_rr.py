# Independent check of addendum-v18 sec. 3.1: four Rogers-Ramanujan-type identities, exact to q^(N-1)
exec(open('../v4/check_V_VI.py').read().split("W  = tot(")[0])
def RRG(sub, sign=1):   # G(sign*q^sub) = 1/((q;q^5)(q^4;q^5)) at q -> sign*q^sub, as a series in q
    s = [1] + [0]*(N-1)
    for k in range(0, N):
        for r in (1, 4):
            e = sub*(5*k + r)
            if e >= N: continue
            c = sign**(5*k + r)
            # divide by (1 - c q^e)
            for i in range(e, N): s[i] += c*s[i-e]
    return s
def RRH(sub, sign=1):
    s = [1] + [0]*(N-1)
    for k in range(0, N):
        for r in (2, 3):
            e = sub*(5*k + r)
            if e >= N: continue
            c = sign**(5*k + r)
            for i in range(e, N): s[i] += c*s[i-e]
    return s
q_q2 = poch(1, 2, -1)                                   # (q;q^2)_inf
L1 = tot([sh(poch(1, 1, -1, 2*m, True), m*m) for m in range(20) if m*m < N])
L2 = tot([sh(poch(1, 1, -1, 2*m+1, True), m*m+m) for m in range(20) if m*m+m < N])
print('sum q^{m^2}/(q)_{2m}        = G(q^4)/(q;q^2):', mul(L1, q_q2) == RRG(4))
print('sum q^{m^2+m}/(q)_{2m+1}    = H(-q)/(q;q^2): ', mul(L2, q_q2) == RRH(1, -1))
mr2 = poch(2, 2, +1)                                    # (-r^2;r^2)_inf
L3 = tot([sh(poch(4, 4, -1, n, True), n*n) for n in range(20) if n*n < N])
L4 = tot([sh(poch(4, 4, -1, n, True), n*n+2*n) for n in range(20) if n*n+2*n < N])
print('sum r^{n^2}/(r^4;r^4)_n     = G(r)/(-r^2;r^2):', mul(L3, mr2) == RRG(1))
print('sum r^{n^2+2n}/(r^4;r^4)_n  = H(r)/(-r^2;r^2):', mul(L4, mr2) == RRH(1))
