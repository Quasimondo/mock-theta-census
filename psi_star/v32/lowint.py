"""State integrals in the lower half-plane of x (also Re x < 0), numpy double precision, Phi_b by the product formula
   Phi_b(y) = (-p1 e^(2 pi y/b); p1^2)_inf / (-q1^(1/2) e^(2 pi b y); q1)_inf,   q1 = e^(-2 pi i x), p1 = e^(pi i/x)   (Im x < 0: |q1|, |p1| < 1),
   which is Phi_{1/b} written out.  Line e^(-i delta) R with delta = (90 degrees - |arg b|)/2."""
import numpy as np


class Low:
    def __init__(self, x, umax=2.0, h=None, tol=30.0):
        self.x = x = complex(x); assert x.imag < 0
        self.b = b = np.sqrt(x)
        self.delta = d = (np.pi/2 - abs(np.angle(b)))/2
        g = np.pi*np.sin(2*d)                              # Gaussian rate
        R = (2*np.pi*umax + np.sqrt((2*np.pi*umax)**2 + 4*g*tol))/(2*g)
        if h is None: h = min(0.004, R/20000)
        s = np.arange(-R, R + h/2, h); self.e = e = np.exp(-1j*d); self.y = y = s*e; self.h = h; self.R = R
        q1 = np.exp(-2j*np.pi*x); p1 = np.exp(1j*np.pi/x)
        lp = np.zeros_like(y)
        K = int(40/(-np.log(abs(p1)**2))) + 5
        wd = 2*np.pi*y/b
        for k in range(K): lp += np.log(1 + np.exp((2*k + 1)*1j*np.pi/x + wd))
        K = int(40/(-np.log(abs(q1)))) + 5
        w = 2*np.pi*b*y
        for k in range(K): lp -= np.log(1 + np.exp(-(k + 0.5)*2j*np.pi*x + w))
        self.logbase = -1j*np.pi*y*y - lp
        self.kern = 1.0/(1.0 + np.exp(wd))
        self.P = np.exp(-2j*np.pi*x/6 + 1j*np.pi/(24*x))/(np.sqrt(2)*b)

    def Iu(self, u):
        return np.sum(np.exp(self.logbase + 2*np.pi*u*self.y))*self.h*self.e

    def Ju(self, u):
        return np.sum(np.exp(self.logbase + 2*np.pi*u*self.y)*self.kern)*self.h*self.e

    def col(self, a):
        f = 1j*self.P*np.exp(-1j*np.pi*a); b = self.b
        return np.array([f*self.Ju(a/b), f*self.Iu(a/b), f*self.Iu((a - 1)/b)])

    def Smat(self, a):
        c0 = self.col(a); c1 = self.col(a + self.x)
        return np.array([[1, c0[0], c1[0]], [0, c0[1], c1[1]], [0, c0[2], c1[2]]])
