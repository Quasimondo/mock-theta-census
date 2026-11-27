# Scripts for psi-star-addendum-v4

Needs Python 3 with mpmath and python-flint. Run from this folder.

    python3 w_identities.py                    # (V), (VI), G relations, exact to q^319
    python3 f_relation.py                      # psi1* = 1 - psi0* - q^2 F(-q^2)
    python3 w_checks.py                        # even-root rule, resurgence of W at xi -> 1
    python3 conv_fit.py conv_limits_hi.json    # converse law: fits and closed forms
    python3 mirror_fit.py                      # mirror law at q = -1: fits and closed forms

The two limit files are included. To regenerate them (slow):

    python3 conv_limits.py 0,1/2,1/3,2/3,1/4,3/4,1/5,2/5 long   # -> conv_limits_hi.json
    python3 mirror_limits.py 0,1/2,1/3,2/3                      # -> mirror_limits.json

run_log.txt is the complete output of the five fast commands above.
Helpers: qs.py (exact integer q-series), jval.py / jflint.py (psi1*, psi0* at roots of unity),
wval.py / wflint.py (W, W3 at roots of unity).
