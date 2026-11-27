# Scripts for psi-star-addendum-v5

Needs Python 3 with mpmath and python-flint. Run from this folder. All limit files are included,
so the fit scripts run in seconds; the commands marked (slow) regenerate the limit files.

    python3 horo.py          # self-test of the evaluators: identities (I), (II), (V), (VI) at complex q
    python3 pert.py          # exact perturbative series of W, W3 (q -> 1) and psi1*, psi0* (q -> -1)
    python3 wronskian.py     # exact q-series of Wr = q W3 Fa - W Fb

    python3 e1_limits.py     # (slow, 1 min)  converse direction: limits along horocycles -> e1_limits.json
    python3 e1_fit.py        # section 2: fits on functions of the dual nome, leading order and first correction
    python3 e2_forward.py    # (slow, 1 min)  forward direction, growing sector -> e2_forward.json
    python3 e2_fit.py        # section 4
    python3 e3_mirror.py     # section 3 (uses e3_mirror.json if present; delete it to recompute, 1.5 min)
    python3 e4_decay.py      # (slow, 25 s)  decaying sector at 12 points -> e4_decay.json
    python3 e4_run.py        # section 5: 20 points on a circle, Taylor coefficients (uses e4_circle_2_5.json if present)
    python3 e4_fit.py        # section 5: comparison with 2 Wr on all points
    python3 constants.py     # closed forms of the first-correction constants

q2_series_v3.json is the perturbative series of psi1*, psi0* at q -> 1 from the v3 scripts (q2_series.py).
run_log.txt is the complete output of the fast commands.

Arithmetic: the long recurrences use midpoint arithmetic at high precision (ball radii would blow up),
so the values are not rigorous enclosures. Stability was checked by rerunning at higher precision (horo.py).
