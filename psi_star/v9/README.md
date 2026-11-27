# Scripts for psi-star-addendum-v9

Needs Python 3 with mpmath and python-flint. Run from this folder.

    python3 jloop.py        # limits of (JJ - Phi) sqrt(tau) at 36 values of Y0, two-exponential fit, first amplitude fit
                            #   (uses jloop.json and jloop_pert.json if present; recomputing takes about 5 min)
    python3 jloop_fit.py    # final fit with exponents -1/120, 71/120 and a constant in the basis

horo.py, e1_limits.py, pert.py, pert.json are helpers from the v5 scripts. run_log.txt is the output of both commands.
The long recurrences use midpoint arithmetic at high precision (not rigorous enclosures).
