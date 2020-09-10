################################################################################
#                                                                              #
#  Weak Equilibrium                                                            #
#                                                                              #
################################################################################

import h5py
import numpy as np

class Opacity:
    def __init__(self,opacfile):
        from scipy.interpolate import RegularGridInterpolator
        with h5py.File(opacfile) as f:
            self.Ye = f['Ye'][()]
            self.lT = f['lT'][()]
            self.lnu = f['lnu'][()]
            self.lrho = f['lrho'][()]
            self.opac = f['opac'][()]
            self.emis = f['emis'][()]
            self.dim = f['dimensions']
            self.n_rho = dim.attrs['numRho']
            self.n_T = dim.attrs['numT']
            self.n_Ye = dim.attrs['numYe']
            self.n_tp = dim.attrs['numRadTypes']
            self.n_nu = dim.attrs['numNu']
            self.iord = dim.attrs['index order']
            self.nue_interp = RegularGridInterpolator((Ye,lT,lrho),
                                                      np.log10(emis_nue_lep_tot),
                                                      bounds_error=False)
            self.nubar_interp = RegularGridInterpolator((Ye,lT,lrho),
                                                        np.log10(emis_nubar_lep_tot),
                                                        bounds_error=False)

class WeakEquilibriumFinder:
    def __init__(self,opacfile):
        self.opac = Opacity(opacfile)

    def find_Ye(lT,lRho):
        from scipy.optimize import brentq
        Ye = self.opac.Ye
        def objective(ye):
            p = (ye,lT,lrho)
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
