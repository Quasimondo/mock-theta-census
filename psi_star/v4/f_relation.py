# Exact check to q^319 of the q-difference relation of F(z) = sum z^n q^(n(n+1)/2) (-q;q)_n:
#   F(z) = 1 + z q F(zq) + z q^2 F(zq^2)
# at z = -1:  psi1* = 1 - psi0* - q^2 F(-q^2),  with psi1* = F(-1), psi0* = q F(-q).
N = 320
def F(shift):                       # F(-q^shift) mod q^N
    out = [0]*N; poch = [1] + [0]*(N-1); n = 0
    while n*(n+1)//2 + shift*n < N:
        e = n*(n+1)//2 + shift*n; s = -1 if n % 2 else 1
        for i in range(N-e): out[i+e] += s*poch[i]
        n += 1
        new = poch[:]
        for i in range(N-n): new[i+n] += poch[i]
        poch = new
    return out
F0, F1, F2 = F(0), F(1), F(2)
lhs = F0
rhs = [0]*N; rhs[0] = 1
for i in range(N-1): rhs[i+1] -= F1[i]
for i in range(N-2): rhs[i+2] -= F2[i]
print("psi1* = 1 - psi0* - q^2 F(-q^2):", lhs == rhs)
print("psi1* =", F0[:12])
