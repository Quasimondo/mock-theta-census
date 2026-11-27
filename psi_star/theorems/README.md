# Source of `docs/psi-star-theorems-v1.md`

- `overview.md`: the hand-written part (title, status, §0: results, map of the proof, what is quoted, where a referee should look first, concordance).
- `assemble.py`: appends to it the proofs from `docs/addendum-v2.md`, `-v3`, `-v4`, `-v12`, `-v14` … `-v17`, `-v19` … `-v21`, in logical order, text unchanged; inherited headings are demoted and prefixed with their source ("[v12 §4.1] …").

To change the consolidated document, edit `overview.md` or the addenda and run `python3 psi_star/theorems/assemble.py`; do not edit the file in `docs/` by hand. `python3 psi_star/theorems/assemble.py --check` tells whether the file in `docs/` is current.
