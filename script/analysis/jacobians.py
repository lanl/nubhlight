# Authors: Jonah Miller (jonahm@lanl.gov) and Ben Ryan
# (brryan@lanl.gov)

# ======================================================================
# copyright 2020. Triad National Security, LLC. All rights
# reserved. This program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
# is operated by Triad National Security, LLC for the U.S. Department
# of Energy/National Nuclear Security Administration. All rights in
# the program are reserved by Triad National Security, LLC, and the
# U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others
# acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# ======================================================================

# Purpose: Computes Jacobians in python as a cross-check for those
# output by the code
import numpy as np
SMALL = 1.e-40
NDIM = 4

class Jacobians:
    def __init__(self, startx, mks_smooth,
                 hslope, poly_alpha, poly_xt,
                 derefine_poles):
        poly_norm = (0.5*np.pi
                     *1./(1. + 1./(poly_alpha + 1.)*1./(poly_xt**poly_alpha)))
        self.startx = startx
        self.mks_smooth = mks_smooth
        self.hslope = hslope
        self.poly_alpha = poly_alpha
        self.poly_norm = poly_norm
        self.poly_xt = poly_xt
        self.derefine_poles = bool(derefine_poles)

    @classmethod
    def fromhdr(cls, hdr):
        return cls(hdr['startx'],hdr['mks_smooth'],
                   hdr['hslope'],hdr['poly_alpha'],hdr['poly_xt'],
                   hdr['DEREFINE_POLES'])

    def get_r(self,X1,X2,X3):
        return np.exp(X1)

    def get_thG(self,X1,X2,X3):
        hslope = self.hslope
        return np.pi*X2 + ((1. - hslope)/2.)*np.sin(2.*np.pi*X2)

    def get_thJ(self,X1,X2,X3):
        poly_norm = self.poly_norm
        poly_xt = self.poly_xt
        poly_alpha = self.poly_alpha
        y = 2.*X2 - 1.
        thJ = poly_norm*y*(1. + pow(y/poly_xt, poly_alpha) / (poly_alpha + 1.)) + 0.5*np.pi
        return y,thJ

    def get_th(self,X1,X2,X3):
        thG = self.get_thG(X1,X2,X3)
        if self.derefine_poles:
            mks_smooth = self.mks_smooth
            startx = self.startx
            y,thJ = self.get_thJ(X1,X2,X3)
            return thG + np.exp(mks_smooth*(startx[1] - X1))*(thJ - thG)
        return thG

    def get_h2bl(self,X1,X2,X3):
        hslope = self.hslope
        poly_norm = self.poly_norm
        poly_alpha = self.poly_alpha
        poly_xt = self.poly_xt
        mks_smooth = self.mks_smooth
        startx = self.startx

        Jcov = np.zeros((NDIM,NDIM))
        Jcon = np.zeros_like(Jcov)

        drdX1 = np.exp(X1)
        dthGdX2 = np.pi + np.pi*(1-hslope)*np.cos(2*np.pi*X2)
        if self.derefine_poles:
            thG = self.get_thG(X1,X2,X3)
            y,thJ = self.get_thJ(X1,X2,X3)
            dydX2 = 2.
            dthJdy = poly_norm*(1+(y/poly_xt)**poly_alpha)
            dthJdX2 = dthJdy*dydX2
            dthdX1 = -mks_smooth*(thJ-thG)*np.exp(mks_smooth*(startx[1]-X1))
            dthdX2 = dthGdX2 + np.exp(mks_smooth*(startx[1]-X1))*(dthJdX2-dthGdX2)
        else:
            dthdX1 = 0.0
            dthdX2 = dthGdX2

        Jcon[0,0] = 1.
        Jcon[1,1] = drdX1
        Jcon[2,1] = dthdX1
        Jcon[2,2] = dthdX2
        Jcon[3,3] = 1.

        Jcov[0,0] = 1.
        Jcov[1,1] = 1./(drdX1 + SMALL)
        Jcov[2,1] = -dthdX1 / (drdX1*dthdX2 + SMALL)
        Jcov[2,2] = 1./(dthdX2 + SMALL)
        Jcov[3,3] = 1.

        return Jcov, Jcon

    def get_bl2cart(self,X1,X2,X3):
        from numpy import sin,cos
        r = self.get_r(X1,X2,X3)
        th = self.get_th(X1,X2,X3)
        ph = X3

        x = r*sin(th)*cos(ph)
        y = r*sin(th)*sin(ph)
        z = r*cos(th)

        Jcov = np.zeros((NDIM,NDIM))
        Jcon = np.zeros_like(Jcov)

        Jcon[0,0] = 1.
        Jcon[1,1] = sin(th) * cos(ph)
        Jcon[1,2] = r * cos(th) * cos(ph)
        Jcon[1,3] = -r * sin(th) * sin(ph)
        Jcon[2,1] = sin(th) * sin(ph)
        Jcon[2,2] = r * cos(th) * sin(ph)
        Jcon[2,3] = r * sin(th) * cos(ph)
        Jcon[3,1] = cos(th)
        Jcon[3,2] = -r * sin(th)

        Jcov[0,0] = 1.
        Jcov[1,1] = cos(ph) * sin(th)
        Jcov[1,2] = sin(ph) * sin(th)
        Jcov[1,3] = cos(th)
        Jcov[2,1] = cos(ph) * cos(th) / (r + SMALL)
        Jcov[2,2] = cos(th) * sin(ph) / (r + SMALL)
        Jcov[2,3] = -sin(th) / (r + SMALL)
        Jcov[3,1] = -sin(ph) / (r * sin(th) + SMALL)
        Jcov[3,2] = cos(ph) / (r * sin(th) + SMALL)

        return Jcov, Jcon

    def get_h2cart(self,X1,X2,X3):
        J_h2bl_cov,J_h2bl_con = self.get_h2bl(X1,X2,X3)
        J_bl2c_cov,J_bl2c_con = self.get_h2cart(X1,X2,X3)

        Jcov = np.zeros((NDIM,NDIM))
        Jcon = np.zeros_like(Jcov)

        for mupp in range(NDIM):
            for mu in range(NDIM):
                for mup in range(NDIM):
                    Jcon[mupp,mu] += J_bl2c_con[mupp,mup]*J_h2bl_con[mup,mu]
                    Jcov[mu,mupp] += J_h2bl_cov[mu,mup]*J_bl2c_cov[mup,mupp]

        return Jcov,Jcon
