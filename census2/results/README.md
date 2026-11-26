# Results of the first `small` run (Nov 22–25, 2026)

- `stageA.jsonl`, `stageB.jsonl`, `stageC.jsonl`: one line per probe, as written by `run2.py`.
- `report2.md`, `report2.json`: the classification by the current `report2.py`, with the R class (rational exponent, inexact law) separated from "neither".

The run used the original evaluator (`sum2.py`); 573 functions were parked with timeouts. They are being rerun with `sum2_fast.py` (`run2.py --retry`); those results will be added separately. The one mock candidate, d1151, is described in `../../docs/d1151.md`.
