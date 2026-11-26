"""Q2: arithmetic of the perturbative coefficients in q2_series.json (numerical check with 120 digits).
 q = i   : k! a_k is a Gaussian integer
 q -> 1  : 3^(2k+1) k! a_k is an integer
"""
import json, mpmath as mp
from fractions import Fraction
mp.mp.dps = 130
D = json.load(open("q2_series.json")); K = D["K"]
for key in ("psi1@i", "psi0@i", "psi1@-1", "psi0@-1"):
    worst = 0; kbad = None
    for k in range(K+1):
        r, i, _ = D[key][k]
        f = mp.factorial(k)
        x, y = f*mp.mpf(r), f*mp.mpf(i)
        mag = max(abs(x), abs(y), 1)
        if mag > mp.mpf(10)**110: break          # beyond the stored digits
        d = max(abs(x - mp.nint(x)), abs(y - mp.nint(y)))
        worst = max(worst, d); kbad = k
    print(f"{key}: k! a_k within {mp.nstr(worst, 3)} of (Gaussian) integers for k <= {kbad}")
    if key.endswith("@i"):
        print("    first values of k! a_k:", [f"{int(mp.nint(mp.factorial(k)*mp.mpf(D[key][k][0])))}{int(mp.nint(mp.factorial(k)*mp.mpf(D[key][k][1]))):+d}i" for k in range(6)])
for key in ("psi1@1", "psi0@1"):
    worst = 0; kbad = None
    for k in range(K+1):
        x = mp.factorial(k)*mp.mpf(3)**(2*k+1)*mp.mpf(D[key][k][0])
        if abs(x) > mp.mpf(10)**110: break
        worst = max(worst, abs(x - mp.nint(x))); kbad = k
    print(f"{key}: 3^(2k+1) k! a_k within {mp.nstr(worst, 3)} of integers for k <= {kbad}")
    print("    a_0..a_6 =", [str(Fraction(int(mp.nint(mp.factorial(k)*mp.mpf(3)**(2*k+1)*mp.mpf(D[key][k][0]))), 3**(2*k+1)*int(mp.factorial(k)))) for k in range(7)])
