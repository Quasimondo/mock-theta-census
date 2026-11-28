"""State integrals for any x with Re x > 0 (also real x and Im x < 0), numpy double precision.
   I(u) = int_L e^(-pi i y^2 + 2 pi u y)/Phi_b(y) dy,  J(u): the same with the kernel 1/(1 + e^(2 pi y/b)),  L = e^(-i delta) R.
   Phi_b from fad.logPhi_np (integral representation plus functional equations), which does not use the product formula."""
import numpy as np
from fad import logPhi_np


class NP:
    def __init__(self, x, delta=0.4, R=9.0, h=0.02):
        self.x = x = complex(x); self.b = b = np.sqrt(x)
        s = np.arange(-R, R + h/2, h); self.e = e = np.exp(-1j*delta); self.y = y = s*e; self.h = h
        lp = np.array([logPhi_np(2*np.pi*b*yy, x) for yy in y])
        self.base = np.exp(-1j*np.pi*y*y - lp)                       # e^(-pi i y^2)/Phi_b(y)
        self.kern = 1.0/(1.0 + np.exp(2*np.pi*y/b))
        self.P = np.exp(-2j*np.pi*x/6 + 1j*np.pi/(24*x))/(np.sqrt(2)*b)

    def Iu(self, u):
        return np.sum(self.base*np.exp(2*np.pi*u*self.y))*self.h*self.e

    def Ju(self, u):
        return np.sum(self.base*self.kern*np.exp(2*np.pi*u*self.y))*self.h*self.e

    def col(self, a):
        f = 1j*self.P*np.exp(-1j*np.pi*a); b = self.b
        return np.array([f*self.Ju(a/b), f*self.Iu(a/b), f*self.Iu((a - 1)/b)])

    def Smat(self, a):
        c0 = self.col(a); c1 = self.col(a + self.x)
        return np.array([[1, c0[0], c1[0]], [0, c0[1], c1[1]], [0, c0[2], c1[2]]])
