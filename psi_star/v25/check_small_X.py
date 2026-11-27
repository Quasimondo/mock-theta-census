# Independent check of addendum-v25 (own code): the root-of-unity formulas at X = 1, 2, 3, with an own
# Faddeev integral for Phi_b at real b (functional equation to reach the strip), on the line e^{-45 deg i} R (no pole or zero between it and L).
import numpy as np
def logPhi0(y, b, n=6001):            # Faddeev's integral, valid for |Im y| < (b + 1/b)/2
    c = np.pi*min(b, 1/b)/2
    L = 30/(b + 1/b - 2*np.abs(np.imag(y)).max()); s = np.linspace(-L, L, n); xi = s + 1j*c; h = s[1] - s[0]
    f = np.exp(-2j*np.asarray(y)[:, None]*xi)/(4*xi*np.sinh(b*xi)*np.sinh(xi/b))
    return f.sum(axis=1)*h
def logPhi(y, b):                      # bring Im y into |Im y| <= b/2 with Phi(Y) = Phi(Y - ib)/(1 + e^{2 pi b (Y - ib/2)})
    y = np.asarray(y); k = np.round(np.imag(y)/b).astype(int); y0 = y - 1j*k*b
    out = logPhi0(y0, b)
    for i in range(len(y)):
        Y = y[i]
        if k[i] > 0:
            for j in range(k[i]): out[i] -= np.log(1 + np.exp(2*np.pi*b*(Y - 1j*b*j - 1j*b/2)))
        else:
            for j in range(-k[i]): out[i] += np.log(1 + np.exp(2*np.pi*b*(Y + 1j*b*j + 1j*b/2)))
    return out
def integrals(X, us, kern_us):
    b = 1/np.sqrt(X); e = np.exp(-1j*np.radians(45)); r = np.linspace(-6, 6, 4801); y = r*e; h = (r[1] - r[0])*e
    lp = np.concatenate([logPhi(y[i:i+400], b) for i in range(0, len(y), 400)])
    base = np.exp(-np.pi*1j*y**2 - lp)
    I = {u: (base*np.exp(2*np.pi*u*y)).sum()*h for u in us}
    J = {u: (base*np.exp(2*np.pi*u*y)/(1 + np.exp(2*np.pi*y/b))).sum()*h for u in kern_us}
    t = -2j*np.pi/X; P = 2**-0.5*np.exp(np.pi**2/(12*t) + t/6)/b
    return b, P, I, J
# X = 2: psi1*(-1) = 1 = P[J(1/2b) + I(1/2b) - I(-1/2b)];  F(z;-1) = 1 = i P e^{-pi i a}[J(a/b) - I((a-1)/b) + z^2 I(a/b)]
b = 1/np.sqrt(2); a_s = (0.25, 0.3, 0.7)
_, P, I, J = integrals(2, [1/(2*b), -1/(2*b)] + [a/b for a in a_s] + [(a-1)/b for a in a_s], [1/(2*b)] + [a/b for a in a_s])
print('X = 2: psi1* formula =', np.round(P*(J[1/(2*b)] + I[1/(2*b)] - I[-1/(2*b)]), 10), ' (claim 1)')
for a in a_s:
    z = np.exp(2j*np.pi*a)
    print(f'       F(z) formula, a = {a}:', np.round(1j*P*np.exp(-1j*np.pi*a)*(J[a/b] - I[(a-1)/b] + z**2*I[a/b]), 10), ' (claim 1)')
# X = 1 and 3: psi1* radial limit = A/3 = P[J - I(+)/3 - I(-)/3];  X = 3, psi0*: J(b+1/2b) = (I(b+1/2b) + I(b-1/2b))/3
for X, claim in ((1, 1/3), (3, 1)):
    b = 1/np.sqrt(X)
    _, P, I, J = integrals(X, [1/(2*b), -1/(2*b), b + 1/(2*b), b - 1/(2*b)], [1/(2*b), b + 1/(2*b)])
    print(f'X = {X}: psi1* formula =', np.round(P*(J[1/(2*b)] - I[1/(2*b)]/3 - I[-1/(2*b)]/3), 10), f' (claim {claim:.6f})')
    v0 = np.exp(1j*np.pi/X)*P*(J[b + 1/(2*b)] - I[b + 1/(2*b)]/3 - I[b - 1/(2*b)]/3)
    print(f'       psi0* formula =', np.round(v0, 10), f" (claim {'1/3' if X == 1 else '0'})")
