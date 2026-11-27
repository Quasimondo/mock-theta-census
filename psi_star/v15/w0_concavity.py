"""Addendum v15: is the leading exponent psi concave along admissible straight contours?
psi''(rho) = -[ sin(2 delta) (1 + Re sigma(w)) - cos(2 delta) Im sigma(w) ],  sigma(w) = 1/(1 + e^(-w)),  w = +- rho e^(i alpha), alpha = beta - delta.
Grid: theta = arg t from -255 to 60 degrees in steps of 5, delta in steps of 1 degree over the admissible range with |delta + beta| < 77.56.
Floating point.  Usage: python3 w0_concavity.py"""
import numpy as np
rho = np.concatenate([np.linspace(1e-3, 3, 600), np.linspace(3, 40, 800)])
g0 = np.degrees(np.arctan(np.pi/np.log(2)))
tot = 0; conc = 0; none_at = []; all_at = []
for theta in range(-255, 65, 5):
    beta = (theta + 90)/2; n = 0; c = 0
    for d in np.arange(1, 90 - abs(beta), 1.0):
        if abs(d + beta) >= g0: continue
        al = np.radians(beta - d); S, C = np.sin(2*np.radians(d)), np.cos(2*np.radians(d)); ok = True
        for sgn in (1, -1):
            sig = 1/(1 + np.exp(-sgn*rho*np.exp(1j*al)))
            if np.min(S*(1 + sig.real) - C*sig.imag) <= 0: ok = False
        n += 1; c += ok
    tot += n; conc += c
    if c == 0: none_at.append(theta)
    if c == n: all_at.append(theta)
    print(f"theta = {theta:5d}: concave on {c:3d} of {n:3d} admissible contours")
print(f"total: concave on {conc} of {tot};  no concave contour at theta = {none_at};  all concave at theta = {all_at}")
