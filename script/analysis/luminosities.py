################################################################################
#                                                                              #
#  Neutrino Luminosities                                                       #
#                                                                              #
################################################################################

import numpy as np
from scipy import integrate
import hdf5_to_dict as io
from units import cgs

def get_luminosity_from_dump(dump,hdr):
    lnu = np.linspace(hdr['lnumin'],
                      hdr['lnumax'],
                      dump['nuLnu'].shape[-1])
    nu = np.exp(lnu)
    nulnu_isotropic = np.sum(dump['nuLnu'],axis=(1,2))
    rate = nulnu_isotropic/cgs['HPL']/nu
    luminosity_tot = integrate.trapz(nulnu_isotropic,
                                     x=lnu,
                                     axis=1)
    rate_tot = integrate.trapz(rate,
                                     x=lnu,
                                     axis=1)
    return nu, rate_tot, luminosity_tot

def get_luminosity_from_fnam(fnam,hdr,geom):
    print(fnam)
    dump = io.load_dump(fnam,geom=geom)
    t = dump['t']
    nu,rate,lum = get_luminosity_from_dump(dump,hdr)
    return t,nu,rate,lum

class LuminosityGetter:
    def __init__(self,hdr,geom):
        self.hdr = hdr
        self.geom = geom

    def __call__(self,fnam):
        return get_luminosity_from_fnam(fnam,
                                        self.hdr,
                                        self.geom)
