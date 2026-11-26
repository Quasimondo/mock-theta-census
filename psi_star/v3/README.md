# Scripts for psi-star-addendum-v3.md

Run from this folder with Python 3. Needs `mpmath`; `q2_series.py` also needs `python-flint`.

| Script | Section | What it does |
|---|---|---|
| `q3_oddroot.py` | 1 | finite identity A = B, convergence of f*, size of theta4*R, radial values |
| `q5_threshold.py [zfam_results.jsonl]` | 2 | saddle function V, threshold, comparison with the census data, precision check, limit rule |
| `q2_series.py 150` | 3 | perturbative series at q = -1, i and q -> 1, writes `q2_series.json` (about 2 minutes) |
| `q2_largeorder.py` | 3 | large-order fits, comparison with theta4*R, closed forms |
| `q2_integrality.py` | 3 | arithmetic of the coefficients |
| `q4_twists.py` | 4 | identities (III), (IV) to q^319 (`qs.py`: exact series helpers) |
| `jval.py`, `wval.py` | 5 | exact values of psi1*, psi0* and of W at roots of unity |
| `q1_limits.py 9` | 5 | limits L(rho), writes `q1_limits.json` |
| `q1_fit.py` | 5 | linear fits against 1, psi0*, psi1* (negative result) |
| `q1_search.py` | 5 | the search over finite sums that found W |
| `q1_law.py 9` | 5 | test of the law at higher precision, with the first correction |
| `q1_second.py` | 5 | second sector (exploratory) |

`run_log.txt` is the output of all scripts as run for the document. `q2_series.json` and `q1_limits.json` are included so that the later scripts run without recomputing.
