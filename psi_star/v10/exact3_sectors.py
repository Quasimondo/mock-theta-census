"""The same test as exact3.py away from the horocycle: different directions of t (XX = K - i b, arg t = -90 deg - arg XX)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf
import exact3
mp.dps = 60
for K, b in ((40, 10), (30, 15), (20, 20), (10, 20), (5, 20), (0, 20), (-10, 20), (-20, 20), (-40, 10)):
    exact3.run(mpf(K), mpf(b))
