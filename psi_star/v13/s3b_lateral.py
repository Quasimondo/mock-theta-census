"""Addendum v13, step 3b: at the directions where s3_sectors.py is limited by the resummation, P J (values copied from its log, 20 digits)
against lateral Borel sums S_theta for several theta and several Pade orders."""
import sys
from mpmath import mp, mpf, mpc, exp, pi, nstr
from borel import borel_sum, pade_borel, I
mp.dps = 30
vals = {-30: mpc("0.36362091262845505409", "-0.019378584639813526961"), 0: mpc("0.36964085484843653381", "0.000040277505701394339041"),
        30: mpc("0.36454904346831840384", "0.019989486366166202783"), -250: mpc("0.32005458059598281113", "0.029421332111284624254")}
T = mpf("0.3")
for phi, thetas in ((-30, (-30, -45, -60)), (0, (-25, -45, -60, -75)), (30, (-25, -45, -55)), (-250, (-170, -165, -175))):
    t = T*exp(I*phi*pi/180)
    for th in thetas:
        row = []
        for n in (40, 60, 75, 90, 120):
            s = borel_sum("psi1", t, th*pi/180, n); row.append(f"[{n}] {nstr(abs(s - vals[phi]), 2)}")
        print(f"phi = {phi:5d}  theta = {th:5d}   |P J - S_theta|, Pade order:  " + "   ".join(row)); sys.stdout.flush()
