# Check: log Phi_b(z/(2 pi b)) ~ Li2(-e^z)/hbar + sum_{n>=1} B_{2n}(1/2) hbar^{2n-1}/(2n)! Li_{2-2n}(-e^z),  hbar = 2 pi i b^2
import mpmath as mp
mp.mp.dps = 40
def poch_inf(a, q):
    s = mp.mpc(1); k = 0
    while True:
        f = a*q**k; s *= 1 + f
        if abs(f) < mp.mpf(10)**(-45): return s
        k += 1
def Phi(y, x):
    b = mp.sqrt(x); q = mp.expjpi(2*x); qt = mp.expjpi(-2/x)
    return poch_inf(mp.expjpi(x)*mp.exp(2*mp.pi*b*y), q) / poch_inf(mp.expjpi(-1/x)*mp.exp(2*mp.pi*y/b), qt)
z = mp.mpc(0.3, 0.5)
for x in [mp.mpc(0.02, 0.01), mp.mpc(0.01, 0.005), mp.mpc(0.005, 0.0025)]:
    b = mp.sqrt(x); hb = 2j*mp.pi*x
    exact = mp.log(Phi(z/(2*mp.pi*b), x))
    a0 = mp.polylog(2, -mp.exp(z))/hb
    a1 = a0 + mp.bernpoly(2, 0.5)*hb/2*mp.polylog(0, -mp.exp(z))
    a2 = a1 + mp.bernpoly(4, 0.5)*hb**3/24*mp.polylog(-2, -mp.exp(z))
    # exact log is defined mod 2 pi i: compare differences modulo 2 pi i
    def d(a): 
        v = exact - a; v -= 2j*mp.pi*mp.nint(v.imag/(2*mp.pi)); return abs(v)
    print(f'x={mp.nstr(x,3)}:  err order0 {mp.nstr(d(a0),3)}  order1 {mp.nstr(d(a1),3)}  order3 {mp.nstr(d(a2),3)}')
