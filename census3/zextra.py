"""Test the rule 'F(z) grows at odd-order roots of order m exactly when z^m is near 1'
at roots of order 9, 11, 13 that are not in the 19-root list."""
import sys, json
from multiprocessing import Pool
import zfam
jobs = []
for name in ['F[1/2]', 'F[1/4]', 'F[1/6]', 'F[5/11]', 'F[1/8]', 'F[1/10]', 'F[1/12]', 'g3[1/2]', 'g2[1/2]']:
    for (h, m) in [(1, 9), (2, 9), (4, 9), (1, 11), (1, 13), (1, 15)]:
        jobs.append((name, h, m))
with Pool(2, initializer=zfam.init_worker) as p, open('zextra_results.jsonl', 'a') as f:
    for r in p.imap_unordered(zfam.work, jobs):
        f.write(json.dumps(r) + '\n'); f.flush()
        print(r['name'], f"{r['h']}/{r['m']}", zfam.ct(r), r.get('error') or r['kappa_rat'], r['secs'], flush=True)
