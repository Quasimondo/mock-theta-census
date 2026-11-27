"""Addendum v14, exploration: for which directions theta = arg t does a STRAIGHT contour descend at leading order?
Leading-order size of the integrand of P J on the ray s = r e^(-i gamma)  (gamma = delta + beta, beta = arg b = (theta + 90 deg)/2):
     |integrand| ~ exp( psi(rho)/|t| ),  rho = |t| |r|,
     psi(rho) = Re( e^(-i theta) [ -w^2/2 + L(w) + pi^2/12 ] ) - pi rho cos(gamma),      w = i t s = +- i rho e^(i (theta - gamma)),
     L(w) = Li2(-e^w) for Re w <= 0,   L(w) = -pi^2/6 - w^2/2 - Li2(-e^(-w)) for Re w >= 0     (the continuation along the ray).
For each theta the script scans delta in (0, 90 - |beta|) and prints the delta-range on which  sup_rho psi(rho)/rho < 0  (both halves
of the ray), and the best value; at the end it compares the set found with the set predicted by the linear term at the origin.
This is floating-point exploration, not a proof.  Usage: python3 u2_scan.py"""
import numpy as np
from scipy.special import spence

def Li2(z): return spence(1 - z)

def L(w):
    out = np.empty_like(w)
    m = w.real <= 0
    out[m] = Li2(-np.exp(w[m]))
    out[~m] = -np.pi**2/6 - w[~m]**2/2 - Li2(-np.exp(-w[~m]))
    return out

rho = np.concatenate([np.linspace(0.002, 1, 300), np.linspace(1.01, 12, 500), np.linspace(12.1, 60, 300)])

def worst(theta, gamma):
    th, ga = np.radians(theta), np.radians(gamma); best = -1e9
    for sgn in (1, -1):
        w = sgn*1j*rho*np.exp(1j*(th - ga))
        psi = (np.exp(-1j*th)*(-w**2/2 + L(w) + np.pi**2/12)).real - np.pi*rho*np.cos(ga)
        best = max(best, np.max(psi/rho))
    return best

if __name__ == "__main__":
    print("theta = arg t | arg x | beta | delta-range with descent (degrees) | best sup psi/rho | at delta | value for delta = beta (tau real)")
    g0 = np.degrees(np.arctan(np.pi/np.log(2))); mism = 0; total = 0; margin = []
    for theta in list(range(-265, 90, 5)):
        beta = (theta + 90)/2; dmax = 90 - abs(beta)
        ds = np.arange(0.5, dmax, 0.5)
        vals = np.array([worst(theta, d + beta) for d in ds])
        ok = ds[vals < 0]
        pred = np.abs(ds + beta) < g0
        close = np.abs(np.abs(ds + beta) - g0) < 0.26          # grid points within half a step of the predicted boundary
        mism += int(np.sum((pred != (vals < 0)) & ~close)); total += len(ds)
        real_tau = f"{worst(theta, 2*beta):8.3f}" if 0 < beta < dmax else "    -   "
        if len(ok):
            j = int(np.argmin(vals))
            # report whether the admissible set is one interval
            gaps = np.sum(np.diff(ok) > 0.75)
            print(f"{theta:6d} | {theta + 90:5d} | {beta:6.1f} | {ok[0]:5.1f} .. {ok[-1]:5.1f}{' (not one interval)' if gaps else ''} | {vals[j]:8.3f} | {ds[j]:5.1f} | {real_tau}")
        else:
            j = int(np.argmin(vals))
            print(f"{theta:6d} | {theta + 90:5d} | {beta:6.1f} | none | {vals[j]:8.3f} | {ds[j]:5.1f} | {real_tau}")
    print(f"predicted set: 0 < delta < 90 - |beta| and |delta + beta| < arctan(pi/ln 2) = {g0:.3f} degrees")
    print(f"grid points (theta, delta) where 'descent' differs from the prediction, not counting points within 0.26 degrees of the boundary: {mism} of {total}")
