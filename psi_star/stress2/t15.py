"""T15: Theorem 15 of v23, the volume rule at even roots of unity.

 K(a, X) := X (F(z; e^{2 pi i/X}) / prediction - 1) for even X up to 6400 and a from 1e-3 to 1 - 1e-3
 (exact sums in mpmath with 30 + 0.08 X digits).  Questions: does K converge (O(1/X)), how does the
 constant behave as a -> 0, 1, and is the remainder consistent with |E| <= 1/(2 sqrt2 sin(pi a)) + small?
Run:  python3 t15.py   (about 8 min)
"""
import time
import numpy as np
import mpmath as mp
from math import pi, sin, sqrt
import common2 as c

XS = [50, 100, 200, 400, 800, 1600, 3200, 6400]
AS = [0.001, 0.003, 0.01, 0.03, 0.05, 0.1, 0.25, 1 / 3, 0.5, 2 / 3, 0.75, 5 / 6, 0.9, 0.95, 0.97, 0.99, 0.997, 0.999]


def main():
    print("K(a, X) = X (F/prediction - 1); last column: Richardson limit 2K(2X) - K(X) from the two largest X")
    print("and X Vol(a)/2pi at the first X from which |K| stays within 1% of its limit")
    rows = []
    for a in AS:
        Ks = []
        for X in XS:
            F = c.F_finite(a, X)
            P = c.F_prediction(a, X)
            K = X * (F / P - 1)
            Ks.append(complex(K))
        rich = 2 * Ks[-1] - Ks[-2]
        rows.append((a, Ks, rich))
        ons = [X for i, X in enumerate(XS) if all(abs(k / rich - 1) < 0.01 for k in Ks[i:])]
        onset = "%.1f" % (ons[0] * c.vol(a) / (2 * pi)) if ons else "> %.1f" % (XS[-1] * c.vol(a) / (2 * pi))
        print("  a=%.4f Vol=%.4f: " % (a, c.vol(a)) + " ".join("%9.3f" % abs(k) for k in Ks)
              + "   lim %.4f%+.4fi  |lim| %.3f   onset X Vol/2pi %s" % (rich.real, rich.imag, abs(rich), onset))
    print("  (|K| printed for X = " + ", ".join(str(x) for x in XS) + ")")
    print("\n  growth of the constant at the ends: |lim K| * Vol(a) and |lim K| * a(1-a)")
    for a, Ks, rich in rows:
        print("   a=%.4f: |K| Vol = %.4f, |K| a(1-a) = %.4f, Vol X/2pi at X=6400: %.2f"
              % (a, abs(rich) * c.vol(a), abs(rich) * a * (1 - a), c.vol(a) * 6400 / (2 * pi)))
    print("\n  comparison with the Richardson values quoted in v23: a = 1/4, 1/3, 5/6, 1/20")
    quoted = {0.25: -0.21657 - 1.34534j, 1 / 3: -0.18821 - 1.50019j, 5 / 6: 0.09094 - 1.15995j, 0.05: -0.05388 - 1.15907j}
    for a, q in quoted.items():
        r = [rr for aa, kk, rr in rows if abs(aa - a) < 1e-12]
        if r:
            print("   a=%.4f: here %.5f%+.5fi, quoted %.5f%+.5fi" % (a, r[0].real, r[0].imag, q.real, q.imag))


if __name__ == "__main__":
    T = time.time()
    main()
    print("\nrun time %.0f s" % (time.time() - T))
