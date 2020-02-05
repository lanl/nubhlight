#!/usr/bin/env python

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

# Takes spherically averaged data and generates a text file that
# summarizes it in the time-averaged sense.

from __future__ import print_function,division
import scipy as sp
import numpy as np

import sys,os,glob
import hdf5_to_dict as io

import h5py
from analyze_eos import load_eos
from scipy import interpolate,integrate

from argparse import ArgumentParser

SMALL = 1e-20

parser = ArgumentParser("Average spherically averaged data in time")
parser.add_argument('dumps',type=str,
                    help=('directory containing dumps to average. '
                          +'Assumed to contain 2d dump data.'))
parser.add_argument('eos',type=str,
                    help='Path to EOS file for evaluating pressure and sound speed')
parser.add_argument('--imin',type=int,
                    default=0,
                    help="Minimum index for dumps to average over")
parser.add_argument('--imax',type=int,
                    default=-1,
                    help="Maximum index for dumps to average over")
parser.add_argument('-s','--save',type=str,
                    default='averaged_quantities.txt',
                    help='Name to save file as')
parser.add_argument('--serial',action='store_true',
                    help='Run in serial')

args = parser.parse_args()

fnams = io.get_dumps_reduced(args.dumps,twod=True)
hdr = io.load_hdr(fnams[0])
eos = load_eos(args.eos)

cs_interp = interpolate.RegularGridInterpolator((eos['Ye'],eos['lT'],eos['lrho']),
                                                eos['cs'],
                                                method='linear',
                                                bounds_error=False,
                                                fill_value=None)
def get_cs(sadw):
    points = np.empty((len(sadw['r']),3))
    points[:,0] = sadw['Ye']
    points[:,1] = np.log10(sadw['TEMP'])
    points[:,2] = np.log10(sadw['RHO']*hdr['RHO_unit'])

    return  cs_interp(points)

P_interp = interpolate.RegularGridInterpolator((eos['Ye'],eos['lT'],eos['lrho']),
                                               eos['P'],
                                               method='linear',
                                               bounds_error=False,
                                               fill_value=None)
def get_PRESS(sadw):
    points = np.empty((len(sadw['r']),3))
    points[:,0] = sadw['Ye']
    points[:,1] = np.log10(sadw['TEMP'])
    points[:,2] = np.log10(sadw['RHO']*hdr['RHO_unit'])

    return  P_interp(points)

def get_quantities(fnam):
    print(fnam)
    t,sadw,sph,zoh = io.load_sadw(fnam)
    q = {}
    q['t'] = t*hdr['T_unit']*1e3
    q['r'] = sadw['r']*hdr['L_unit']
    q['rho']  = sadw['RHO']*hdr['RHO_unit']
    q['T'] = sadw['TEMP']
    q['Ye'] = sadw['Ye']
    q['u'] = sadw['UU']*hdr['U_unit']
    q['P'] = get_PRESS(sadw)
    q['thd'] = zoh['thd']
    q['H'] = zoh['H']*hdr['L_unit']
    q['cs'] = get_cs(sadw)
    return q

imin = args.imin
if imin < 0:
    imin = 0
imax = args.imax
if imax < 0:
    imax = len(fnams) + imax

if args.serial:
    qs = [get_quantities(f) for f in fnams[imin:imax+1]]
else:
    from concurrent import futures
    with futures.ProcessPoolExecutor() as p:
        qs = list(p.map(get_quantities,fnams[imin:imax+1]))

averages = np.zeros((len(qs[0]['r']),9))
averages[:,0] = qs[0]['r']
for j,name in enumerate(['rho','T','Ye','u','P','thd','H','cs']):
    print(j,name)
    for q in qs:
        averages[:,j+1] += q[name]
    averages[:,j+1] /= len(qs)

np.savetxt(args.save,
           averages,
           comments='# ',
           header='[0]:r (cm), [1]:rho (g/cc), [2]:T (MeV), [3]:Ye, [4]:u (erg/cc) [5]:P (erg/cc) [6]:scale angle (radians), [7]:scale height (cm), [8]:cs/c',
           fmt='%.14e')
