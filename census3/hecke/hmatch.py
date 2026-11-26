"""Step 2: find f_{a,b,c}(x,y,q) (Hickerson-Mortenson notation, x,y = +-q^alpha, +-q^beta)
equal to a target series up to q^shift and a constant."""
import sys, math
from sparse import T, N
M = 120
def fabc(a, b, c, sx, al, sy, be):
    s = [0] * M
    # f = sum_{sg r = sg s} sg(r) (-x)^r (-y)^s q^{a C(r,2) + b r s + c C(s,2)},  x = sx q^al
    for sgn in (1, -1):
        r = 0 if sgn == 1 else -1
        while True:
            hit = False
            t = 0 if sgn == 1 else -1
            while True:
                e = a * r * (r - 1) // 2 + b * r * t + c * t * (t - 1) // 2 + al * r + be * t
                if e < M:
                    if e < -50: return None  # not bounded below in this cone
                    if e >= 0:
                        s[e] += sgn * (-sx) ** r * (-sy) ** t if r >= 0 and t >= 0 else sgn * (-sx) ** (-r) * (-sy) ** (-t)
                    hit = True
                elif t * sgn > 40:
                    break
                t += sgn
                if abs(t) > 400: break
            if not hit and abs(r) > 40: break
            r += sgn
            if abs(r) > 400: break
    return s
def norm(s):
    for i, x in enumerate(s):
        if x:
            return i, tuple(y * (1 if x > 0 else -1) for y in s[i:i + 60])
    return None
targets = {nm: norm(T[nm][:M]) for nm in sys.argv[1:] or ['psi3*']}
hits = []
for a in range(1, 7):
    for c in range(a, 7):
        for b in range(-int(math.isqrt(a * c)), 9):
            if b * b == a * c and b < 0: continue
            for al in range(-6, 7):
                for be in range(-6, 7):
                    for sx in (1, -1):
                        for sy in (1, -1):
                            s = fabc(a, b, c, sx, al, sy, be)
                            if s is None: continue
                            n = norm(s)
                            if n is None or len(n[1]) < 60: continue
                            for nm, (i0, v) in targets.items():
                                if n[1] == v:
                                    print(nm, f"= const q^{i0 - n[0]} f_{{{a},{b},{c}}}({'+' if sx > 0 else '-'}q^{al}, {'+' if sy > 0 else '-'}q^{be}; q)", flush=True)
