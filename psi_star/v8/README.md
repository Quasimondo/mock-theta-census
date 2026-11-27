# Scripts for psi-star-addendum-v8

Needs Python 3 with mpmath and python-flint. Run from this folder.

    python3 fstar.py          # limits of F0*, G0, q F1*, -G1 on a circle in the dual nome, Taylor coefficients
                              #   (uses fstar_circle.json if present; delete it to recompute, about 1 min)
    python3 fstar_check.py    # exact check of identity (VII) to r^299; comparison of the limits with JJ, KK

horo.py, e1_limits.py, pert.py, pert.json are helpers from the v5 scripts. run_log.txt is the output of both commands.
As before, the long recurrences use midpoint arithmetic at high precision (not rigorous enclosures).
