# Scripts for psi-star-addendum-v6

Needs Python 3 with mpmath and python-flint. Run from this folder.

    python3 exact_checks.py   # section 1: the four exact identities, to q^309 (about 1 min)
    python3 hfam.py           # section 2: Casoratian at a complex point; which combinations cancel the growing sector
    python3 e5_state.py       # section 2: limits of the four bilinear combinations at six points
                              #   (uses e5_state.json if present; delete it to recompute, about 1 min)

horo.py, e1_limits.py, pert.py, pert.json are helpers from the v5 scripts (evaluation, extrapolation).
run_log.txt is the complete output of the three commands.
As in v5, the long recurrences use midpoint arithmetic at high precision (not rigorous enclosures).
