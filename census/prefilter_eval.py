import json, collections, mpmath as mp
from fractions import Fraction
from multiprocessing import Pool
from dilog import predict
reps = {r['id']: r for r in json.load(open('cands.json'))}
label = {int(k): v for k, v in json.load(open('labels.json')).items()}
def work(i):
    r = dict(reps[i]); r['pat'] = [tuple(x) for x in r['pat']]
    k = predict(r)
    if k < mp.mpf(10) ** -12: return i, 0.0, True
    f = Fraction(str(k)).limit_denominator(5000)
    return i, float(k), abs(mp.mpf(f.numerator) / f.denominator - k) < mp.mpf(10) ** -15
if __name__ == '__main__':
    plain = [i for i in reps if reps[i]['s'] > 0]
    with Pool(2) as p: out = p.map(work, plain, chunksize=20)
    T = collections.Counter()
    miss = []
    for i, k, rational in out:
        lab = label[i]
        kept = lab in ('mock', 'modular', 'known(catalog)')
        T[('pred-out' if not rational else 'pred-pass', 'keep' if kept else lab if lab in ('quantum',) else ('out@0/1' if lab == 'out@0/1' else 'out-later'))] += 1
        if not rational and kept: miss.append((i, k, lab))
    for k in sorted(T): print(k, T[k])
    print('kept functions wrongly excluded:', miss)
    json.dump([[i, k, r] for i, k, r in out], open('dilog_pred.json', 'w'))
