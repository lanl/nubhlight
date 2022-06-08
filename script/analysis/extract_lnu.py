#!/usr/bin/env python

################################################################################
#                                                                              # 
#  Makes a text file from nuLnu from dumps                                     # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
import h5py, argparse
import numpy as np
from units import get_cgs
from scipy import integrate
from hdf5_to_dict import get_dumps_full
units = get_cgs()

def get_data_dump(dumpname):
    with h5py.File(dumpname,'r') as f:
        nuLnu = f['nuLnu'][()]
        nth = f['nth'][0]
        nphi = f['nphi'][0]
        nubins = f['nubins'][0]
        numin = f['numin'][0]
        numax = f['numax'][0]
        t_unit = f['T_unit'][0]
        t = f['t'][0]*t_unit
        R = f['Rout_rad'][0]*f['L_unit'][0]

    nu = np.zeros(nubins)
    lnumin = np.log(numin)
    lnumax = np.log(numax)
    dlnu = (lnumax - lnumin)/nubins
    for n in range(len(nu)):
      nu[n] = np.exp(lnumin + (0.5 + n)*dlnu)

    theta = np.linspace(0,np.pi,nth)
    phi = np.linspace(0,2*np.pi,nphi)

    eps = units['HPL']*nu/units['MEV']
    l10eps = np.log10(eps)
    Lnu = nuLnu/(nu[np.newaxis,np.newaxis,np.newaxis,:])
    integrand = R*R*Lnu*np.sin(theta[np.newaxis,:,np.newaxis,np.newaxis])
    Lnu_1 = integrate.simps(integrand,x=theta,axis=1)
    Lnu = integrate.simps(Lnu_1,x=phi,axis=1)

    return t,R,l10eps,Lnu

def get_time_series(folder,savepath='Lnu.h5'):
    dfnams = get_dumps_full(folder)
    t,R,l10eps,Lnu = get_data_dump(dfnams[0])

    times = np.empty(len(dfnams))
    Lnus  = np.empty((len(dfnams),*Lnu.shape))

    for i,fnam in enumerate(dfnams):
        print("...",fnam)
        times[i],R,l10eps,Lnus[i] = get_data_dump(fnam)

    with h5py.File(savepath,'w') as f:
        dset = f.create_dataset("logener",data=l10eps)
        dset.attrs['description'] = "log10 of neutrino energy in MeV"
        dset = f.create_dataset("t",data=times)
        dset.attrs['description'] = "time at each dump, in seconds"
        dset = f.create_dataset("Lnu",data=Lnus)
        dset.attrs['description'] = "Luminosity (in cgs)"
        dset.attrs['indices'] = "time; species: [electron, anti, heavy]; neutrino energy"
        dset = f.create_dataset('R',data=R)
        dset.attrs['description'] = "Radius at which Luminosity is measured"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract luminosity from a sequence of dumps.")
    parser.add_argument('directory',
                        type=str,
                        help='Directory to read from')
    parser.add_argument('-s','--save',
                        type=str,
                        default='Lnu.h5',
                        help='File to save to. Defaults to ./Lnu.h5')
    args = parser.parse_args()
    print("Extracting luminosities from...")
    get_time_series(args.directory,args.save)
    print("Done!")
