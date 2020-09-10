################################################################################
#                                                                              #
#  Weak Equilibrium                                                            #
#                                                                              #
################################################################################

import h5py
import numpy as np
from units import cgs

class Opacity:
    def __init__(self,opacfile):
        from scipy.interpolate import RegularGridInterpolator
        from scipy.integrate import trapz
        with h5py.File(opacfile,'r') as f:
            self.Ye = f['Ye'][()]
            self.lT = f['lT'][()]
            self.lnu = f['lnu'][()]
            self.lrho = f['lrho'][()]
            self.opac = f['opac'][()]
            self.emis = f['emis'][()]
            self.dim = f['dimensions']
            self.n_rho = self.dim.attrs['numRho']
            self.n_T = self.dim.attrs['numT']
            self.n_Ye = self.dim.attrs['numYe']
            self.n_tp = self.dim.attrs['numRadTypes']
            self.n_nu = self.dim.attrs['numNu']
            self.iord = self.dim.attrs['index order']
            self.emis_nue = self.emis[...,0,:]
            self.emis_nubar = self.emis[...,1,:]
            self.emis_nue_lep_tot = 4*np.pi*trapz(
                self.emis_nue/cgs['HPL'],
                x=self.lnu,axis=-1)
            self.emis_nubar_lep_tot = 4*np.pi*trapz(
                self.emis_nubar/cgs['HPL'],
                x=self.lnu,axis=-1)
            self.nue_interp = RegularGridInterpolator(
                (self.Ye,self.lT,self.lrho),
                np.log10(self.emis_nue_lep_tot),
                bounds_error=False)
            self.nubar_interp = RegularGridInterpolator(
                (self.Ye,self.lT,self.lrho),
                np.log10(self.emis_nubar_lep_tot),
                bounds_error=False)

class WeakEquilibriumFinder:
    def __init__(self,opacfile):
        self.opac = Opacity(opacfile)

    def find_Ye(self,lT,lRho):
        from scipy.optimize import brentq
        Ye = self.opac.Ye
        def objective(ye):
            p = (ye,lT,lRho)
            return self.opac.nue_interp(p) - self.opac.nubar_interp(p)
        try:
            sol = brentq(objective,Ye.min(),Ye.max())
        except:
            if objective(Ye.max()) > 0 and objective(Ye.min()) > 0:
                sol = Ye.min()
            elif objective(Ye.min()) < 0 and objective(Ye.max()) < 0:
                sol = Ye.max()
            else:
                sol = 0.5
        return sol

    def get_Ye_eq_dmp(self,dump,hdr):
        Ye_equil = np.empty_like(dump['Ye'])
        lrho = np.log10(dump['RHO']*hdr['RHO_unit'])
        lT = np.log10(dump['TEMP'])
        for i in range(Ye_equil.shape[0]):
            for j in range(Ye_equil.shape[1]):
                for k in range(Ye_equil.shape[2]):
                    Ye_equil[i,j,k] = self.find_Ye(lT[i,j,k],lrho[i,j,k])
        return Ye_equil
