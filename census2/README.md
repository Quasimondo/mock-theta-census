# census2: double-sum (rank-2) census

Same classifier as `../census`, applied to double sums

  Σ_{m,n≥0} s₁^m s₂^n q^{(a m² + 2b mn + c n² + d m + e n)/2} · F₁(q)_m · F₂(q)_n · F₁₂(q)_{m+n}

where F₁ and F₂ are products like 1/(q)_m or (−q)_m, and F₁₂ is an optional coupling factor on m+n.

## Run (low priority, resumable)

```bash
pip install python-flint mpmath numpy scipy
python3 gen2.py --size small        # ~2 min: candidates, dedupe, q→1 pre-filter -> cands2.json
python3 run2.py --bench --workers 4 # times 12 probes, estimates stage A
nohup python3 run2.py --workers 4 > run2.log 2>&1 &
tail -f run2.log                    # progress and ETA
python3 report2.py                  # summary -> report2.md, report2.json
```

- The workers run at `nice 19`: they only use CPU time your own work leaves idle. Choose `--workers` to leave cores free; RAM use is small (well under 1 GB per worker).
- Stopping (Ctrl-C, `kill`, reboot) is safe. Rerun the same command and it continues where it stopped.
- Stages: A tests q → 1 only, B the roots of order 2–5, C all 19 roots. Each stage only gets the survivors of the one before.
- Sizes: `small` = 8,442 forms, 6,327 distinct functions, 4,843 after the pre-filter. `medium` (112k forms) and `large` (457k) are for later.

Rough timing for `small` on 4 workers: stage A a few hours, then roughly 1–2 days for stages B and C. The bench prints an estimate for your machine.

## Files

- `gen2.py`: candidate generator, exact coefficients (dedupe up to constant, monomial, sign and q → −q), and the q → 1 dilogarithm pre-filter for sums with positive terms. Self-check: the Andrews–Gordon sum Σ q^{N₁²+N₂²}/((q)_{N₁−N₂}(q)_{N₂}) gives the right product coefficients and κ = 1/42.
- `sum2_fast.py`: evaluator used by default. The double sum is split as A_m·B_n·C_{m+n} and computed as one FLINT polynomial product, with the truncation added to the error ball. `python3 sum2_fast.py` runs the self-test; `python3 sum2_fast.py compare 12` compares it with `sum2.py` on real candidates.
- `sum2.py`: the original evaluator (`run2.py --old-eval`).
- `run2.py`: staged runner (`stageA.jsonl`, `stageB.jsonl`, `stageC.jsonl`, one line per probe).
- `report2.py`: classification (mock, modular, quantum type, R type, neither), with matches to the single-sum census marked. R type: a rational exponent κ at some root but no exact law there (power-series corrections, as for R₀, R₁ at q = −1). Such functions are not modular or mock, so they leave the mock search, but they are kept as their own class. A κ counts as rational if it lies within 10⁻⁸ of a fraction with denominator ≤ 500, so a few irrational values may be misfiled.
- `catalog.py`, `probe.py`: copies from `../census`.
