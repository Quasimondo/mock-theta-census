import sys
import sparse, ref
from sparse import mul, etaq, series_sum
sparse.T.clear()
sparse.T['f0(q)_inf'] = mul(series_sum(ref.f0), etaq({1: 1}))
sparse.T['psi1(q2;q2)'] = mul(series_sum(ref.psi1), etaq({2: 1}))
sparse.T['psi1(q)_inf'] = mul(series_sum(ref.psi1), etaq({1: 1}))
sys.argv = ['x'] + list(sparse.T)
exec(open('hmatch.py').read().replace('from sparse import T, N', 'from sparse import T, N'))
