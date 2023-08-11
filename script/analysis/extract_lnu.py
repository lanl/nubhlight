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
from hdf5_to_dict import get_dumps_reduced
from multiprocessing import Pool
units = get_cgs()

def get_data_dump(dumpname):
    with h5py.File(dumpname,'r') as f:
        print("...",dumpname)
        nuLnu = f['nuLnu'][()]
        nth = f['nth'][0]
        nphi = f['nphi'][0]
        nubins = f['nubins'][0]
        numin = f['numin'][0]
        numax = f['numax'][0]
        t_unit = f['T_unit'][0]
        t = f['t'][0]*t_unit
        R = f['Rout_rad'][0]*f['L_unit'][0]

    return t,R,nth,nphi,nubins,numin,numax,nuLnu

def get_time_series(folder,twod=False,savepath='nuLnu.h5',nproc=None):
    dfnams = get_dumps_reduced(folder,twod)
    if len(dfnams) == 0:
        raise ValueError("No valid dumps in directory",folder)
    t,R,nth,nphi,nubins,numin,numax,nuLnu = get_data_dump(dfnams[0])

    nu = np.zeros(nubins)
    lnumin = np.log(numin)
    lnumax = np.log(numax)
    dlnu = (lnumax - lnumin)/nubins
    for n in range(len(nu)):
      nu[n] = np.exp(lnumin + (0.5 + n)*dlnu)
    lnu = np.log(nu)

    theta_f = np.linspace(0,np.pi,nth+1)
    phi_f = np.linspace(0,2*np.pi,nphi+1)
    
    theta_c = 0.5*(theta_f[1:] + theta_f[:-1])
    phi_c = 0.5*(phi_f[1:] + phi_f[:-1])

    eps = units['HPL']*nu/units['MEV']
    l10eps = np.log10(eps)

    times = np.empty(len(dfnams))
    nuLnus = np.empty((len(dfnams),*nuLnu.shape))

    p = Pool(processes=nproc)
    tuples = p.map(get_data_dump, dfnams)
    for i,t in enumerate(tuples):
        times[i] = t[0]
        nuLnus[i] = t[7]

    with h5py.File(savepath,'w') as f:
        dset = f.create_dataset("lnu", data = lnu)
        dset.attrs['description'] = "ln of neutrino frequency, in Hz"
        dset.attrs['units'] = 'Hz'
        dset = f.create_dataset("nu", data = nu)
        dset.attrs['description'] = "neutrino frequency, in Hz"
        dset.attrs['units'] = 'Hz'
        dset = f.create_dataset("logener",data=l10eps)
        dset.attrs['description'] = "log10 of neutrino energy in MeV"
        dset.attrs['units'] = 'MeV'
        dset = f.create_dataset("t",data=times)
        dset.attrs['description'] = "time at each dump, in seconds"
        dset.attrs['units'] = 's'
        dset = f.create_dataset("nuLnu",data=nuLnus)
        dset.attrs['description'] = "luminous Flux, nuLnu. To get bolumetric, sum over th,phi, integrate over lnu. To get luminosity vs frequency and angle, divide by nu"
        dset.attrs['units'] =  "Hz * ergs / solid angle / Hz."
        dset.attrs['indices'] = "time; species: [electron, anti, heavy]; theta; phi; neutrino energy"
        dset = f.create_dataset('R',data=R)
        dset.attrs['description'] = "radius at which Luminosity is measured, in cm"
        dset.attrs['units'] = 'cm'
        dset = f.create_dataset("theta_f", data=theta_f)
        dset.attrs['description'] = "edges of solid angle bins in theta"
        dset.attrs['units'] = "radians"
        dset = f.create_dataset("phi_f", data=phi_f)
        dset.attrs['description'] = "edges of solid angle bins in phi"
        dset.attrs['units'] = "radians"
        dset = f.create_dataset("theta_c", data=theta_c)
        dset.attrs['description'] = "centers of solid angle bins in theta"
        dset.attrs['units'] = "radians"
        dset = f.create_dataset("phi_c", data=phi_c)
        dset.attrs['description'] = "centers of solid angle bins in phi"
        dset.attrs['units'] = "radians"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Extract luminosity from a sequence of dumps.")
    parser.add_argument('directory',
                        type=str,
                        help='Directory to read from')
    parser.add_argument('-s','--save',
                        type=str,
                        default='nuLnu.h5',
                        help='File to save to. Defaults to ./nuLnu.h5')
    parser.add_argument("-n","--nprocs",type=int,default=None,
                        help="Number of processors to use. Defaults to all.")
    parser.add_argument("-t","--twod", action="store_true",
                        help="Use 2D dumps instead of 3D dumps.")
    args = parser.parse_args()
    print("Extracting luminosities from...")
    get_time_series(args.directory,args.twod,args.save,args.nprocs)
    print("Done!")
