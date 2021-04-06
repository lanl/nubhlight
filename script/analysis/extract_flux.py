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
Extracts physical variables (fluxes) from a SuperNu dump file at several
radii (for now, fixed at approx 500, 1000, ... 3500 km)
For each radius, produces a separate ASCII file "flux_{dump####}_ir{###}.dat
in a gnuplot format (single-line separated blocks)

Usage:
  $ %s <dumpfile> [<dumpfile2> [ ... ]]
""" % sys.argv[0]
if (len(sys.argv) < 2 or (len(sys.argv) == 2 and sys.argv[1] == '-h') ):
  print (helpstr)
  exit()
read_geom = True

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

# compute hth and hphi, scaling factors for dth/dt-> vth and dphi/dt-> vphi
a = hdr['a']
a2= a*a
r2= r*r
sth2 = np.sin(th)**2
Sigma = r2 + a2*np.cos(th)**2
hth  = np.sqrt(Sigma)
hphi = np.sqrt(sth2 * (Sigma + a2*(1 + 2*r/Sigma)*sth2))

# --------------------------------------------------------------------------
# Loop over the arguments and make sure file exists
for dump_fname in sys.argv[1:]:
  if not os.path.isfile(dump_fname):
    sys.stderr.write(f"ERROR: path {dump_fname} does not exist")
    sys.exit(1)
  match = re.search("dump_([0-9]+).h5",dump_fname)
  if not match:
    sys.stderr.write(f"ERROR: {dump_fname} is not a valid dump file")
    sys.exit(1)
  dump_index = match.group(1)
  output_dir = "." # TODO: add as an option

  # ------------------------------------------------------------------------
  # load data
  hdr = io.load_hdr(dump_fname)
  data = io.load_dump(dump_fname, geom=geom)
  #chm.check_transformation_matrices(geom, hdr['a'], -1, 64)


  # TODO: automate
  # list of radial indices -- roughly corresponds to
  #          500, 1000, 1500, 2000, 2500, 3000, 3500 km
  for ir in [132,  152,  164,  173,  179,  184,  189]:

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
    vr   = u1pr_R/Gamma_R
    vth  = u2pr_R/Gamma_R
    vphi = u3pr_R/Gamma_R
    for m in range(hdr['N3']):
      vth[:,m] *= hth[ir,:]
      vphi[:,m] *= hphi[ir,:]

    output_fname = f"{output_dir}/flux_{dump_index}_ir{ir}.dat"
    fluxfile = open(output_fname, 'w')

    # header
    tsec = data['t']*hdr['T_unit']
    Rext = r[ir,0]*hdr["L_unit"]
    fluxfile.write( f"""# Simulation: {hdr['PATH']}
# Produced from dump file: {dump_fname}
# Dump count: {dump_index}
# Time [s]: t= {tsec}
# Extraction radius [cm]: R= {Rext}
# 1:theta 2:phi      3:rho[g/cm3]   4:v_r[c]      5:v_theta[c]"""
"    6:v_phi[c]     7:u[erg/cm3]   8:P[dyne/cm2]  9:T[K]        10:Ye\n")

    # because in the data the conversion is into ergs somehow
    TEMP_unit = 1.1604525e+10 # MeV -> K

    # output loop
    for jth in range(0,hdr['N2']):
      for kph in range(0,hdr['N3']):
        rho = data['RHO'][ir,jth,kph]*hdr['RHO_unit']
        u   = data['UU'][ir,jth,kph]*hdr['U_unit']
        P   = data['PRESS'][ir,jth,kph]*hdr['U_unit']
        T   = data['TEMP'][ir,jth,kph]*TEMP_unit
        Ye  = data['Ye'][ir,jth,kph]
        bet = data['beta'][ir,jth,kph]
        fluxfile.write("%9.7f %9.7f "\
          "%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n"
          % (th[ir,jth], phi[kph],
            rho, vr[jth,kph], vth[jth,kph], vphi[jth,kph], u, P, T, Ye))
      fluxfile.write("\n") # empty line for gnuplot compatibility

    fluxfile.close()



