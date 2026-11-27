"""Builds docs/psi-star-theorems-v1.md: the hand-written overview (overview.md, section 0) followed by the proofs,
taken unchanged from docs/addendum-v*.md in logical order. Headings of the inherited sections are demoted by one level
and prefixed with their source, e.g. "## 2. Theorem 9" in addendum v20 becomes "### [v20 §2] Theorem 9".
Usage (from anywhere):  python3 psi_star/theorems/assemble.py            writes docs/psi-star-theorems-v1.md
                        python3 psi_star/theorems/assemble.py --check    only reports whether the file in docs/ is up to date"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "..", "docs") + os.sep
OUT = os.path.join(D, "psi-star-theorems-v1.md")
def grab(n, start, end):
    lines = open(D + f"addendum-v{n}.md").read().split("\n")
    i = next(k for k, l in enumerate(lines) if re.match(start, l))
    j = next(k for k, l in enumerate(lines) if k > i and re.match(end, l))
    out = []
    for l in lines[i:j]:
        m = re.match(r"^(#{2,3}) (\d+(?:\.\d+)*)\.? (.*)$", l)
        if m: l = "#" + m.group(1) + f" [v{n} §{m.group(2)}] " + m.group(3)      # demote and prefix with the source
        elif l.startswith("### "): l = "####" + l[3:]
        elif l.startswith("## "): l = "###" + l[2:]
        out.append(l)
    while out and out[-1].strip() in ("", "---"): out.pop()
    return "\n".join(out)
def src(n, what): return f"*From `docs/addendum-v{n}.md`, {what}; text unchanged.*\n\n"
parts = [open(os.path.join(HERE, "overview.md")).read()]
def part(title, intro, blocks):
    s = f"## {title}\n\n{intro}\n\n"
    for n, what, a, b in blocks:
        s += src(n, what) + grab(n, a, b) + "\n\n"
    parts.append(s + "---\n")
part("Part 1. The kernel identity (Theorem I)", "Local notation: see §0.6. Statements in §1, proofs in §2–§5.",
     [(12, "§1–§5", r"^## 1\. ", r"^## 6\. ")])
part("Part 2. Faddeev's dilogarithm", "Lemma Φ and the strip lemmas S, F, D (v16), then the bounds in the form used for the state integrals (Lemma B) and Laplace's method (Lemma L) from v19.",
     [(16, "§1–§3", r"^## 1\. ", r"^## 4\. "), (19, "§2", r"^## 2\. ", r"^## 3\. ")])
part("Part 3. Descent of the leading exponent", "The contour angle is δ = 10° from here on.",
     [(15, "§1", r"^## 1\. ", r"^## 2\. ")])
part("Part 4. The uniform bound for the kernel integrand", "Needed later: Lemma 3′ (away from the real x-axis) and the first bullet of the proof of Theorem 3 (near it), which together give |H| ≤ C·e^{−κ|s|} on the closed sector, and the Corollary (analyticity). The expansions themselves (Theorems 1–3) are superseded by Theorem 4 of Part 5.",
     [(14, "§2.1–§2.3 (contour, Euler–Maclaurin lemma, second form of the integrand)", r"^### 2\.1 ", r"^### 2\.4 "),
      (15, "§2 up to Lemma 3′", r"^## 2\. ", r"^\*\*Lemma 4′"),
      (16, "§4", r"^## 4\. ", r"^## 5\. ")])
part("Part 5. Factorial bounds and Borel summability (Theorem II)", "The coefficients called ãₖ in the proof of Theorem 4 are the aₖ of §0.1, by definition; the last step of that proof (ãₖ = aₖ) only identifies them with the coefficients computed in v13.", [(17, "§1–§5", r"^## 1\. ", r"^## 6\. ")])
part("Part 6. The state integrals and the forward law (Theorems III, IV)", "Lemmas B and L are in Part 2.",
     [(19, "§1", r"^## 1\. ", r"^## 2\. "), (19, "§3–§6", r"^## 3\. ", r"^## 7\. ")])
part("Part 7. Roots of unity (Theorem V)", "", [(20, "§1–§3 (even X)", r"^## 1\. ", r"^## 4\. "), (21, "§1–§3 (odd X)", r"^## 1\. ", r"^## 4\. "),
      (25, "§1 (X = 1, 2, 3)", r"^## 1\. ", r"^## 2\. ")])
part("Appendix. Identities and the odd-root rule", "Used in Part 7: identity (II) (and (I)) in the odd-root rule, and identities (V), (VI) in Lemma R2′.",
     [(2, "§2", r"^## 2\. ", r"^## 3\. "), (3, "§1", r"^## 1\. ", r"^## 2\. "), (4, "§1", r"^## 1\. ", r"^## 2\. ")])
out = "\n".join(parts).rstrip() + "\n"
if "--check" in sys.argv:
    same = os.path.exists(OUT) and open(OUT).read() == out
    print("docs/psi-star-theorems-v1.md is", "up to date" if same else "NOT up to date")
    sys.exit(0 if same else 1)
open(OUT, "w").write(out)
print(len(out), "characters;", out.count("\n"), "lines written to docs/psi-star-theorems-v1.md")
