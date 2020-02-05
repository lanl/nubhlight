#!/usr/bin/env python

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# A little example script showing how to plot some quantities in the EOS table
# Should work for python2 or python3

from __future__ import print_function,division
import numpy as np
import h5py
import matplotlib as mpl
from matplotlib import pyplot as plt

CL = 2.99792458e10
EOS_PATH = 'Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5'
mpl.rcParams.update({'font.size':18})

def load_eos(filename):
    """A little function that reads in the EOS
    and converts it into a dictionary.
    Also has a few names and aliases.
    """

    eos = {}
    # load file
    with h5py.File(filename,'r') as f:
        for k,v in f.items():
            eos[k] = v.value

    # derived quantities
    lrho = eos['logrho']
    lT = eos['logtemp']
    Ye = eos['ye']
    lP = eos['logpress']
    le = eos['logenergy']
    ent = eos['entropy']
    rho = 10.**lrho
    T = 10.**lT # temperature
    P = 10.**lP # pressure
    e = 10.**le - eos['energy_shift'] # specific internal energy
    w = rho*e + P # enthalpy by volume
    h = CL*CL + e + P/rho # specific enthalpy
    hgeom = h/(CL*CL) # specific enthalpy in natural units
    hm1 = hgeom -1
    lhm1 = np.log10(np.abs(hm1))

    # quick aliases
    eos['lrho'] = lrho
    eos['lT'] = lT
    eos['Ye'] = Ye
    eos['lP'] = lP
    eos['le'] = le
    eos['ent'] = ent
    
    eos['rho'] = rho
    eos['T'] = T
    eos['P'] = P
    eos['e'] = e
    eos['w'] = w
    eos['h'] = h
    eos['hgeom'] = hgeom
    eos['hm1'] = hm1
    eos['lhm1'] = lhm1    
    
    return eos

print("Loading eos name\n\t{}".format(EOS_PATH))
eos = load_eos(EOS_PATH)
print(("EOS tables are 3d and depend on:\n"
       +"\tlog_10(density) in g/cc\n"
       +"\tlog_10(temperature) in MeV\n"
       +"\tElectron fraction Ye\n"))
print("")
print("The tables are in column-major order")
print("")
print(("The independent variables are:\n"
       +"['logrho', 'logtemp', 'ye']\n"
       +"and have shapes:"))
for k in ['logrho','logtemp','ye']:
    print("\t{}: {}".format(k,len(eos[k])))
print("")
print("The have min and max values of:")
for k in ['logrho','logtemp','ye']:
    print("\t{}: [{}, {}]".format(k,eos[k].min(),eos[k].max()))
print("")      
print(("A dependent variable, "
       +"for example,\n"
       +"average mass of heavy ions, "
       +"has shape:\n"
       +"\t{}".format(eos['Abar'].shape)))
print("So you access Abar as:\n\tAbar(ye,logtemp,logrho)")
print("")
print("One subtlety in the tables is that\n"
      +"specific internal energy is stored\n"
      +"on a log scale. However, it can be negative\n"
      +"because binding energy can be larger than\n"
      +"thermal energy. Therefore you can get the\n"
      +"true value as:\n"
      +"\t10.**eos['logenergy'] - eos['energy_shift']\n"
      +"like so:\n"
      +"\tminimum energy is = {}".format(eos['e'].min()))
print("")
print("The relevant composition variables are:\n"
      +"\tAbar = Average heavy ion mass\n"
      +"\tZbar = Average heavy ion atomic number\n"
      +"\tXa   = Mass fraction of alpha particles\n"
      +"\tXh   = Mass fraction of heavy ions\n"
      +"\tXp   = Mass fraction of protons\n"
      +"\tXn   = Mass fraction of neutrons\n"
      +"\tXd   = Mass fraction of deuterium\n"
      +"\tX3he = Mass fraction of helium 3\n"
      +"\tX4li = Mass fraction of helium 4\n"
      +"\tXt   = Mass fraction of tau particles\n")
print("")
print("Let's try to make a 2d plot of neutron mass fraction\n"
      +"as a function of Ye and temperature,\n"
      +"for a density of about ~10^11 g/cm^3.")
irho = np.where(eos['logrho'] >= 11)[0][0]
print("The index where rho is this value is roughly {}".format(irho))
mesh = plt.pcolormesh(eos['logtemp'],eos['ye'],eos['Xn'][:,:,irho],
                      cmap='viridis',shading='gouraud')
plt.xlabel(r'$\log_{10}T$')
plt.ylabel(r'$Y_e$')
cbar = plt.colorbar(mesh,label=r'$X_n$')
plt.text(-1.9,0.5,
         r'$\log_{10}\rho \sim$ %.2f' % eos['logrho'][irho],
         color='r')
plt.savefig('Xn.png',bbox_inches='tight')
print("File saved as Xn.png.")
print("")
print("That about covers it!")
