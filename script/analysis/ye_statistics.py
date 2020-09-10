################################################################################
#                                                                              #
#  Weak Equilibrium                                                            #
#                                                                              #
################################################################################

import numpy as np
from scipy import integrate

import hdf5_to_dict as io
from ye_equil import WeakEquilibriumFinder

SMALL = 1e-20

def integrate_3d(q,hdr,geom):
    if hdr['N3'] > 1:
        return integrate.simps(
            integrate.simps(
                integrate.simps(
                    q,dx = hdr['dx'][3],axis=2),
                dx = hdr['dx'][2],axis=1),
            dx = hdr['dx'][1],axis=0)
    else:
        return 2*np.pi*integrate.simps(
            integrate.simps(
                q[:,:,0],dx = hdr['dx'][2],axis=1                
            ),
            dx = hdr['dx'][1],axis=0
        )

class YeStatistics:
    def __init__(self, finder=None):
        self.finder = finder
        
    def get_from_dump(self,dump,hdr,geom):
    
        integrand_num = dump['RHO']*dump['Ye']*geom['gdet'][:,:,np.newaxis]
        integrand_den = dump['RHO']*geom['gdet'][:,:,np.newaxis]
    
        num = integrate_3d(integrand_num,hdr,geom)
        den = integrate_3d(integrand_den,hdr,geom)
        mean = num/(den+SMALL)
    
        integrand_num = dump['RHO']*geom['gdet'][:,:,np.newaxis]*(dump['Ye']-mean)**2
        num = integrate_3d(integrand_num,hdr,geom)
        std = np.sqrt(num/(den+SMALL))
    
        if self.finder is not None:
            Ye_equil = self.finder.get_Ye_eq_dmp(dump,hdr)
            integrand_num = dump['RHO']*Ye_equil*geom['gdet'][:,:,np.newaxis]
            num = integrate_3d(integrand_num,hdr,geom)
            equil = num/(den+SMALL)
    
            integrand_num = dump['RHO']*geom['gdet'][:,:,np.newaxis]*(Ye_equil-equil)**2
            num = integrate_3d(integrand_num,hdr,geom)
            equil_std = np.sqrt(num/(den+SMALL))
    
            integrand_num = dump['RHO']*geom['gdet'][:,:,np.newaxis]*(dump['Ye']-Ye_equil)**2
            num = integrate_3d(integrand_num,hdr,geom)
            diff = np.sqrt(num/(den+SMALL))
            
            return mean,std,equil,equil_std,diff

        else:
            return mean,std

    def get_from_name(self,name,hdr,geom):
        dump = io.load_dump(name,geom=geom)
        t = dump['t']
        print("t = ",t)
        
        if self.finder is not None:
            mean,std,equil,equil_std,diff = self.get_from_dump(dump,hdr,geom)
        else:
            mean,std = self.get_from_dump(dump,hdr,geom)           
        
        print("\tt = ",t," complete")
        if self.finder is not None:
            return t,mean,std,equil,equil_std,diff
        else:
            return t,mean,std
