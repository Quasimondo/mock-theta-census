# Check the expansion of log T(n), T(n) = (-1)^n q^{n(n+1)/2} (-q;q)_n, q = e^{-h}, x = n h, against exact values
import mpmath as mp
mp.mp.dps = 40
exec(open('gz_recipe.py').read().split('Sbar =')[0])
def exact_logT(n, h):
    q = mp.exp(-h); P = mp.mpc(1)
    for k in range(1, n+1): P *= 1 + q**k
    return 1j*mp.pi*n - h*n*(n+1)/2 + mp.log(P)
for h in [mp.mpc(0.02, 0.01), mp.mpc(0.01, 0.005)]:
    for n in [30, 60]:
        x = n*h
        approx0 = f(x)/h + l(x)
        approx1 = approx0 + h*m(x)
        ex = exact_logT(n, h)
        d = lambda a: abs((ex - a) - 2j*mp.pi*mp.nint(mp.im(ex - a)/(2*mp.pi)))
        print(f'h={mp.nstr(h,3)} n={n}: |err| order 0: {mp.nstr(d(approx0),3)}   order 1: {mp.nstr(d(approx1),3)}')
