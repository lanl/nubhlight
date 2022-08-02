################################################################################
#                                                                              #
#  Profiles Integrated over disk                                               #
#                                                                              #
################################################################################

import numpy as np
from scipy import integrate,interpolate
import sadw
import matplotlib as mpl
#matplotlib.use('agg')
import matplotlib.pyplot as plt

def get_scale_heights(hdr,geom,dump):
    r = geom['r'][:,0,0]
    th0 = sadw.get_density_weighted_profile(dump,geom['th'] - np.pi/2)
    th0 = th0[:,np.newaxis,np.newaxis] + np.pi/2.
    thdsqr = sadw.get_density_weighted_profile(dump,(geom['th'] - th0)**2)
    thd = np.sqrt(thdsqr)
    H = r*np.tan(thd)
    zoH = np.tan(geom['th'] - th0)/np.tan(thd)[:,np.newaxis,np.newaxis]
    return th0,thd,H,zoH

class zoHProjector:
    def __init__(self,hdr,geom,dump,
                 r_range=None,
                 zoH_range=None):
        self.hdr = hdr
        self.geom = geom
        self.dump = dump
        self.r = geom['r'][:,0,0]
        self.phi = geom['phi'][0,0,:]

        self.r_range = r_range
        self.zoH_range = zoH_range
        self.irmin = None
        self.irmax = None
        self.zoH_uniform = None

        self.th0,self.thd,self.H,self.zoH = get_scale_heights(hdr,geom,dump)

        if r_range is not None:
            self.set_r_range(r_range[0],r_range[1])
        if zoH_range is not None:
            self.set_zoH_range(zoH_range[0],zoH_range[1])

    def set_r_range(self,rmin,rmax):
        self.irmin = np.where(self.r >= rmin)[0][0]
        self.irmax = np.where(self.r >= rmax)[0][0]
        self.r_range = rmin,rmax

    def set_zoH_range(self,zoH_min,zoH_max):
        self.zoH_range = zoH_min,zoH_max
        self.zoH_uniform = np.linspace(zoH_min,zoH_max,
                                       self.hdr['N2'])
        return self.zoH_uniform

    def _interp_to_zoH(self,q,ir,iphi=0):
        if self.zoH_uniform is None:
            raise ValueError("You must set zoH_range "
                             +"before calling this function.")
        q_interp = interpolate.interp1d(self.zoH[ir,:,iphi],q[ir,:,iphi],
                                        fill_value='extrapolate')
        q_grid = q_interp(self.zoH_uniform)
        return q_grid

    def _integrate_radially(self,q):
        q_grid = np.empty((self.irmax-self.irmin,
                           len(self.zoH_uniform),
                           self.hdr['N3']))
        for ir in range(self.irmin,self.irmax):
            for iphi in range(self.hdr['N3']):
                q_grid[ir-self.irmin,:,iphi] = self._interp_to_zoH(q,ir,iphi)
        q_int = integrate.simps(q_grid,x=self.r[self.irmin:self.irmax],axis=0)
        norm = integrate.simps(np.ones_like(q_grid),
                               x=self.r[self.irmin:self.irmax],
                               axis=0)
        if self.hdr['N3'] > 1:
            q_int = integrate.simps(q_int,x=self.phi,axis=-1)
            norm = integrate.simps(norm,x=self.phi,axis=-1)
            q_avg = q_int/norm
        else:
            q_avg = q_int[...,0]/norm[...,0]
        return q_avg

    def __call__(self,q):
        return self._integrate_radially(q)
        
