exec(open('check_laws.py').read().split("c = mp.sqrt(2 * mp.pi")[0])
c = mp.sqrt(2 * mp.pi / mp.sqrt(5))
rho, M0 = Fr(0), 3
xs, ys = [], []
for k in (1, 2, 4, 8, 16, 32, 48):
    Y = rho + M0 + 120 * k; Yf = mp.mpf(Y.numerator) / Y.denominator
    tau = -2j * mp.pi / Yf
    q1 = -Y / 2
    e1 = mp.expjpi(Yf / 60); e2 = mp.expjpi(-11 * Yf / 60)
    lhs = Gsum(1 / Y, 1); pred = c * (-e1 * psi_star(q1, 0) + 1j * e2 * psi_star(q1, 1))
    xs.append(1 / mp.sqrt(Yf)); ys.append(lhs * mp.sqrt(tau) - pred)
    print(f"Y={int(Yf):6d}  raw diff {mp.nstr(abs(ys[-1]), 3):>9}   Neville on all points so far: {mp.nstr(abs(neville(xs, ys)), 3) if len(xs) > 1 else '-'}", flush=True)
