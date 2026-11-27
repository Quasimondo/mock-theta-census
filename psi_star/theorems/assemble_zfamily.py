"""Builds docs/zfamily-theorems-v1.md: the hand-written overview (overview_zfamily.md, section 0) followed by the
proofs for the z-family, taken unchanged from docs/addendum-v*.md in logical order. Same conventions as assemble.py:
headings of the inherited sections are demoted by one level and prefixed with their source, e.g. "[v22 §2]".
Usage (from anywhere):  python3 psi_star/theorems/assemble_zfamily.py            writes docs/zfamily-theorems-v1.md
                        python3 psi_star/theorems/assemble_zfamily.py --check    only reports whether it is up to date"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "..", "docs") + os.sep
OUT = os.path.join(D, "zfamily-theorems-v1.md")
def grab(n, start, end):
    lines = open(D + f"addendum-v{n}.md").read().split("\n")
    i = next(k for k, l in enumerate(lines) if re.match(start, l))
    j = next(k for k, l in enumerate(lines) if k > i and re.match(end, l))
    out = []
    for l in lines[i:j]:
        m = re.match(r"^(#{2,3}) (\d+(?:\.\d+)*)\.? (.*)$", l)
        if m: l = "#" + m.group(1) + f" [v{n} §{m.group(2)}] " + m.group(3)
        elif l.startswith("### "): l = "####" + l[3:]
        elif l.startswith("## "): l = "###" + l[2:]
        out.append(l)
    while out and out[-1].strip() in ("", "---"): out.pop()
    return "\n".join(out)
def src(n, what): return f"*From `docs/addendum-v{n}.md`, {what}; text unchanged.*\n\n"
parts = [open(os.path.join(HERE, "overview_zfamily.md")).read()]
def part(title, intro, blocks):
    s = f"## {title}\n\n{intro}\n\n"
    for n, what, a, b in blocks:
        s += src(n, what) + grab(n, a, b) + "\n\n"
    parts.append(s + "---\n")
part("Part 1. The exact formula (Theorem Z)", "Uses Theorems A–C of v12, Part 1 of `docs/psi-star-theorems-v1.md`.",
     [(18, "§1", r"^## 1\. ", r"^## 2\. ")])
part("Part 2. Saddles and the volume", "The identity for Im V is proved here; the selection of the saddle was found numerically in this section and is proved by Parts 4 and 6.",
     [(18, "§2.1–§2.3", r"^### 2\.1 ", r"^### 2\.4 ")])
part("Part 3. The decaying state integral (Theorem 7′)", "Lemmas B and L are in Part 2 of `docs/psi-star-theorems-v1.md` (v19 §2, completed in v26 §2). The statements of v19 §3 include Theorem 6 (z = −1) for reference.",
     [(19, "§3–§4", r"^## 3\. ", r"^## 5\. "), (26, "§3 (Lemma V, Theorem 7′)", r"^## 3\. ", r"^## 4\. ")])
part("Part 4. The growing state integral (Theorem 13)", "",
     [(22, "§1–§3", r"^## 1\. ", r"^## 4\. "), (22, "§5", r"^## 5\. ", r"^## 6\. ")])
part("Part 5. The kernel integral at real x (Lemmas U, J)", "",
     [(23, "§1", r"^## 1\. ", r"^## 2\. "), (27, "§3", r"^## 3\. ", r"^## 4\. ")])
part("Part 6. Even roots of unity: the exact formula and the volume rule (Theorems 14, 15, 15′)", "Lemmas R1–R3′ are in Part 7 of `docs/psi-star-theorems-v1.md`.",
     [(22, "§4", r"^## 4\. ", r"^## 5\. "), (23, "§2", r"^## 2\. ", r"^## 3\. "), (27, "§1–§2", r"^## 1\. ", r"^## 3\. ")])
part("Part 7. z = 1 and z = q: the fifth-order mock theta functions (Proposition M, Theorems 16, 17)", "",
     [(18, "§3.1", r"^### 3\.1 ", r"^### 3\.2 "), (24, "§1–§4", r"^## 1\. ", r"^## 5\. "), (27, "§5", r"^## 5\. ", r"^## 6\. ")])
part("Part 8. Remarks (inference)", "Not used by the theorems above.",
     [(18, "§3.2–§4", r"^### 3\.2 ", r"^## 5\. ")])
out = "\n".join(parts).rstrip() + "\n"
if "--check" in sys.argv:
    same = os.path.exists(OUT) and open(OUT).read() == out
    print("docs/zfamily-theorems-v1.md is", "up to date" if same else "NOT up to date")
    sys.exit(0 if same else 1)
open(OUT, "w").write(out)
print(len(out), "characters;", out.count("\n"), "lines written to docs/zfamily-theorems-v1.md")
