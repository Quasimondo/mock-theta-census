"""Addendum v13, step 2: singularities of the Borel transform of Phi from Pade approximants [n/n].
Prints the poles nearest the origin and a count of all poles by direction.  Poles that sit on a zero of the numerator
(distance below 1e-10) are spurious and are listed separately.
Reference points: S = pi^2/3 - i Cl2(pi/3) and its conjugate (the actions of v5/v7).   Usage: python3 s2_poles.py [n ...]"""
import sys
from mpmath import mp, mpf, mpc, pi, polyroots, nstr, arg, clsin
from borel import pade_borel
mp.dps = 30
Cl2 = clsin(2, pi/3); S = mpc(pi**2/3, -Cl2)
print("S = pi^2/3 - i Cl2 =", nstr(S, 15), "  |S| =", nstr(abs(S), 8), "  arg S (deg) =", nstr(arg(S)*180/pi, 6))
for case in ("psi1", "psi0"):
    for n in ([int(a) for a in sys.argv[1:]] or [40, 60]):
        P, Q = pade_borel(case, n)
        old = mp.dps; mp.dps = 200
        roots = polyroots(Q, maxsteps=500, extraprec=600); zeros = polyroots(P, maxsteps=500, extraprec=600)
        mp.dps = old
        info = sorted(((abs(r), r, min(abs(r - z) for z in zeros)) for r in roots), key=lambda v: v[0])
        genuine = [v for v in info if v[2] > mpf(10)**(-10)]; spurious = [v for v in info if v[2] <= mpf(10)**(-10)]
        print(f"{case}, Pade [{n}/{n}]: {len(genuine)} poles, {len(spurious)} spurious pole-zero pairs")
        print("   nearest six (pole, |pole|, arg in degrees):")
        for a, r, dz in genuine[:6]:
            print("      ", nstr(r, 12), "  ", nstr(a, 8), "  ", nstr(arg(r)*180/pi, 6))
        print("   distance of the nearest pole to S:", nstr(min(abs(r - S) for a, r, dz in genuine), 3), "  to conj(S):", nstr(min(abs(r - S.conjugate()) for a, r, dz in genuine), 3))
        bins = {}
        for a, r, dz in genuine:
            d = float(arg(r)*180/pi); key = "|arg| < 20" if abs(d) < 20 else ("20 <= |arg| < 60" if abs(d) < 60 else ("60 <= |arg| < 120" if abs(d) < 120 else "|arg| >= 120"))
            bins.setdefault(key, []).append((a, r))
        for key in ("|arg| < 20", "20 <= |arg| < 60", "60 <= |arg| < 120", "|arg| >= 120"):
            lst = bins.get(key, [])
            print(f"   {key} degrees: {len(lst)} poles" + (";  nearest: " + nstr(lst[0][1], 8) + "  |.| = " + nstr(lst[0][0], 6) if lst else ""))
        sys.stdout.flush()
