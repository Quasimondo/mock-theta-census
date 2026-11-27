"""Closed forms for the first-correction constants, against the fitted values printed by e1_fit.py and e3_mirror.py."""
from mpmath import mp, mpf, sqrt, pi, nstr
mp.dps = 30
s5 = sqrt(5); phi = (1 + s5)/2; c = sqrt(2*pi/s5); al = mpf(5)**(-mpf(1)/4); kap = sqrt(pi/s5)
rows = [
 ("W,  sector u0", mpf("-0.373087647859312573476")/(-c),            (20 + 3*s5)/120, "(20+3 sqrt5)/120"),
 ("W,  sector u1", mpf("0.1856734710871066989069")/c,               (20 - 3*s5)/120, "(20-3 sqrt5)/120"),
 ("W3, sector u0", mpf("-2.980695269792955782095")/(-c*phi),        (152 - 9*s5)/120, "(152-9 sqrt5)/120"),
 ("W3, sector u1", mpf("-1.48600927661128422847")/(-c/phi),         (152 + 9*s5)/120, "(152+9 sqrt5)/120"),
 ("P+, sector u0", mpf("0.06524789890454809788127")/al,             (5 + 3*s5)/120, "(5+3 sqrt5)/120"),
 ("P+, sector u1", mpf("-0.009519540156512929214329")/al,           (5 - 3*s5)/120, "(5-3 sqrt5)/120"),
 ("P-, sector u0", mpf("-0.02817474287320855140331")/(al*phi),      (17 - 9*s5)/120, "(17-9 sqrt5)/120"),
 ("P-, sector u1", mpf("-0.1278646616212899208641")/(-al/phi),      (17 + 9*s5)/120, "(17+9 sqrt5)/120"),
 ("mirror psi1*, u0", mpf("-0.44462169415799582369")/(-kap/phi),    (23 + 6*s5)/60, "(23+6 sqrt5)/60"),
 ("mirror psi1*, u1", mpf("0.30633537246335717236")/(kap*phi),      (23 - 6*s5)/60, "(23-6 sqrt5)/60"),
 ("mirror psi0*, u0", mpf("0.013310297081110392547")/(kap*phi),     (-13 + 6*s5)/60, "(-13+6 sqrt5)/60"),
 ("mirror psi0*, u1", mpf("0.3225279127487263321")/(-kap/phi),      (-13 - 6*s5)/60, "(-13-6 sqrt5)/60"),
]
for lab, fitted, exact, name in rows:
    print(f"{lab:18s} fitted {nstr(fitted, 18):>22s}   {name:18s} = {nstr(exact, 18):>22s}   diff {nstr(abs(fitted - exact), 2)}")
print("leading constants: 5^(-1/4) =", nstr(al, 20), " phi 5^(-1/4) =", nstr(al*phi, 20), " 5^(-1/4)/phi =", nstr(al/phi, 20))
print("                   c = sqrt(2 pi/sqrt5) =", nstr(c, 20), " c phi =", nstr(c*phi, 20), " c/phi =", nstr(c/phi, 20))
