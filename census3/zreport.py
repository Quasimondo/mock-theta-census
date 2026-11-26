"""Fingerprints of the z-family scan over all roots h/m (both halves of the circle).
x at (m-h)/m is read off conj(x) at h/m."""
import json, sys
from fractions import Fraction
from collections import defaultdict
import probe
from zfam import ct

R = {}
for l in open('zfam_results.jsonl'):
    r = json.loads(l)
    k = (r['name'], r['h'], r['m'])
    if k not in R or ('error' in R[k] and 'error' not in r):
        R[k] = r

def parse(name):
    fam, a = name[:-1].split('[')
    return fam, Fraction(a)

def partner(name):
    fam, a = parse(name)
    b = (1 - a) % 1
    return f"{fam}[{b.numerator}/{b.denominator}]"

FULL = []
for (h, m) in probe.ROOTS:
    FULL.append((h, m, False))
for (h, m) in probe.ROOTS:
    if h != 0 and 2 * h != m:
        FULL.append((h, m, True))  # stands for root (m-h)/m

def cells(name):
    out, kap = [], []
    for (h, m, mirror) in FULL:
        n = partner(name) if mirror else name
        r = R.get((n, h, m))
        out.append('.' if r is None else ct(r))
        kap.append(None if r is None or 'error' in r else r['kappa_rat'])
    return ''.join(out), kap

names = sorted({k[0] for k in R}, key=lambda n: (parse(n)[0], parse(n)[1].denominator, parse(n)[1]))
groups = defaultdict(list)
for n in names:
    fp, kap = cells(n)
    if '.' in fp:
        g = 'incomplete'
    elif set(fp) <= {'Q'}:
        g = 'bounded everywhere'
    elif set(fp) <= {'E', 'Q'}:
        g = 'E/Q only (modular or mock type)'
    elif '?' in fp and set(fp) <= {'E', 'Q', '?'}:
        g = 'E/Q with errors'
    else:
        g = 'neither'
    groups[g].append((n, fp, kap))
nroot = len(FULL)
print(f"{len(names)} functions, {nroot} roots (first {len(probe.ROOTS)}: h<=m/2, then the mirrored ones)\n")
for g in ['bounded everywhere', 'E/Q only (modular or mock type)', 'E/Q with errors', 'neither', 'incomplete']:
    if groups[g]:
        print(f"## {g}: {len(groups[g])}")
        for n, fp, kap in groups[g]:
            ks = sorted({k for k, c in zip(kap, fp) if c == 'E' and k}, key=Fraction)
            print(f"  {n:12s} {fp[:len(probe.ROOTS)]} {fp[len(probe.ROOTS):]}  kappa {' '.join(ks)}")
        print()
