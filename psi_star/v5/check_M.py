# The decaying-sector function measured in addendum-v5 sec. 5 vs the exact q-series 2(q W3 F_a - W F_b)
exec(open('/home/claude/v4check/check_V_VI.py').read().split("print('(V)")[0])
Fa, Fb = L5, L6
M = [2*(a - b) for a, b in zip(sh(mul(W3, Fa), 1), mul(W, Fb))]
print('2(q W3 F_a - W F_b) =', M[:14])
print('measured           =', [-2, -2, 2, -2, -2, 0, -8, -2, 0, -4, -4, -4, 0, -6])
