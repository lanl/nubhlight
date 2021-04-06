#!/usr/bin/env python

import sys, os
import re
from math import *
import numpy as np
sys.path.append("/home/korobkin/num/nubhlight/script/analysis")
import hdf5_to_dict as io
import check_transformation_matrices as chm

# --------------------------------------------------------------------------
# help message
helpstr= """
Calculate and output quantities averaged over extraction spheres at radii
{R1,R2,...Rn} for a sequence of dumps from nubhlight. For now, the set of
radii is fixed to approx 500, 1000, ... 3500 km (TODO: generalize)

For each radius, produces a separate ASCII file named integrated_ir{###}.dat

Usage:
  $ %s <dumpfile1> [<dumpfile2> [ ... ]]
""" % sys.argv[0]
if (len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] == '-h') ):
  print (helpstr)
  exit()
verbose = True

# --------------------------------------------------------------------------
# load grid data (geom)
dump_fname = sys.argv[1]
if not os.path.isfile(dump_fname):
  sys.stderr.write(f"ERROR: path {dump_fname} does not exist")
  sys.exit(1)

hdr = io.load_hdr(dump_fname)
geom = io.load_geom(hdr,recalc=True)

# extract coordinates
r   = geom['r'][:,:,0]
th  = geom['th'][:,:,0]
phi = geom['phi'][0,0,:]
dA  = (hdr['dx'][2]*hdr['dx'][3]*hdr['L_unit']**2  # area element
      *np.sqrt(geom['gcov'][:,:,2,2]*geom['gcov'][:,:,3,3]))

# compute hth and hphi, scaling factors for dth/dt-> vth and dphi/dt-> vphi
a = hdr['a']
a2= a*a
r2= r*r
sth2 = np.sin(th)**2
Sigma = r2 + a2*np.cos(th)**2
hth  = np.sqrt(Sigma)
hphi = np.sqrt(sth2 * (Sigma + a2*(1 + 2*r/Sigma)*sth2))

# TODO: automate
# list of radial indices -- roughly corresponds to
#      500, 1000, 1500, 2000, 2500, 3000, 3500 km
irs = [132,  152,  164,  173,  179,  184,  189]
output_fnames = {}
for ir in irs:
  output_dir = "." # TODO: add as an option
  output_fnames[ir] = f"{output_dir}/integrated_ir{ir}.dat"
  outfile = open(output_fnames[ir], 'w')
  Rext = r[ir,0]*hdr["L_unit"]
  outfile.write( f"""# Simulation: {hdr['PATH']}
# Extraction radius [cm]: R= {Rext}
# 1:time[s] 2:mdot[Msun/s] 3:<vr>/c 4:vth/c 5:vphi/c 6:P 7:T 8:Ye 9:beta
""")
  outfile.close()

# --------------------------------------------------------------------------
# Loop over the arguments and make sure file exists
for dump_fname in sys.argv[1:]:
  # ---
  if verbose: sys.stdout.write(f"Processing {dump_fname}...")
  if not os.path.isfile(dump_fname):
    sys.stderr.write(f"ERROR: path {dump_fname} does not exist")
    sys.exit(1)
  match = re.search("dump_([0-9]+).h5",dump_fname)
  if not match:
    sys.stderr.write(f"ERROR: {dump_fname} is not a valid dump file")
    sys.exit(1)
  dump_index = match.group(1)
  output_dir = "." # TODO: add as an option

  # load data
  hdr = io.load_hdr(dump_fname)
  data = io.load_dump(dump_fname, geom=geom)
  #chm.check_transformation_matrices(geom, hdr['a'], -1, 64)

  # shortcuts to various quantities
  rho = data['RHO']*hdr['RHO_unit']
  u   = data['UU']*hdr['U_unit']
  P   = data['PRESS']*hdr['U_unit']
  T   = data['TEMP']
  Ye  = data['Ye']
  bet = data['beta']
  tsec = data['t']*hdr['T_unit']

  # ------------------------------------------------------------------------
  # loop over extraction radii
  for ir in irs:
    # ---
    if verbose: sys.stdout.write(f" {ir}")

    # sample gcov and h2bl at the radial point `ir`
    gcov_R = geom['gcov'][ir,:]
    h2bl_con_R = geom['Lambda_h2bl_con'][ir,:]
    alp_R  = geom['alpha'][ir,:,:]
    ut_R   = data['ucon'][ir,:,:,0]

    # 4-velocity projections in code coordinates
    U1_R = data['U1'][ir,:,:]
    U2_R = data['U2'][ir,:,:]
    U3_R = data['U3'][ir,:,:]

    # Lorentz factor
    Gamma_R = alp_R[:,:]*ut_R[:,:]

    # 4-velocity projections in spherical KS coordinates
    u1pr_R = np.zeros_like(U1_R)
    u2pr_R = np.zeros_like(U1_R)
    u3pr_R = np.zeros_like(U1_R)
    for m in range(hdr['N3']):
      u1pr_R[:,m] = U1_R[:,m]*h2bl_con_R[:,1,1]
      u2pr_R[:,m] = U1_R[:,m]*h2bl_con_R[:,2,1] + U2_R[:,m]*h2bl_con_R[:,2,2]
    u3pr_R[:,:] = U3_R[:,:]

    # physical velocity in normal frame
    clight = 2.99792458e10 # [cm/s]
    Msun = 1.9891e33       # [g]
    vr   = u1pr_R/Gamma_R*clight
    vth  = u2pr_R/Gamma_R*clight
    vphi = u3pr_R/Gamma_R*clight
    for m in range(hdr['N3']):
      vth[:,m]  *=  hth[ir,:]
      vphi[:,m] *= hphi[ir,:]

    # mdot
    mdot = np.sum(rho[ir,:,:]*vr[:,:]*dA[ir,:,None], axis=(0,1))

    # areal mass
    A = np.sum(rho[ir,:,:]*dA[ir,:,None])

    # averages
    vthA = np.sum(rho[ir,:,:]*vth[:,:]*dA[ir,:,None], axis=(0,1))
    vphiA = np.sum(rho[ir,:,:]*vphi[:,:]*dA[ir,:,None], axis=(0,1))
    PA = np.sum(rho[ir,:,:]*P[ir,:,:]*dA[ir,:,None], axis=(0,1))
    TA = np.sum(rho[ir,:,:]*T[ir,:,:]*dA[ir,:,None], axis=(0,1))
    YeA = np.sum(rho[ir,:,:]*Ye[ir,:,:]*dA[ir,:,None], axis=(0,1))
    betA = np.sum(rho[ir,:,:]*bet[ir, :,:]*dA[ir,:,None], axis=(0,1))

    # output
    outfile = open(output_fnames[ir], 'a')
    outfile.write(
      "%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n"
      % (tsec, mdot/Msun, mdot/A/clight, vthA/A/clight, vphiA/A/clight,
         PA/A, TA/A, YeA/A, betA/A))
    outfile.close()

  # ---
  if verbose: sys.stdout.write("\n")


