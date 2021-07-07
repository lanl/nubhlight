################################################################################
#                                                                              #
# BONDI INFLOW                                                                 #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
from shutil import copyfile
import glob
import numpy as np
#import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis/')
import util
import hdf5_to_dict as io
from bhlight import bcall
TMP_DIR = 'TMP'
TMP_BUILD = 'build_tmp.py'
util.safe_remove(TMP_DIR)

AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

RES = [32, 64, 128, 256]

util.make_dir(TMP_DIR)
os.chdir('../prob/bondi/')

copyfile('build.py', TMP_BUILD)
# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE

for n in range(len(RES)):
  util.change_cparm('N1TOT', RES[n], TMP_BUILD)
  util.change_cparm('N2TOT', RES[n], TMP_BUILD)
  call([sys.executable, TMP_BUILD, '-dir', TMP_DIR])
  call(['cp', os.path.join(os.getcwd(), TMP_DIR, 'bhlight'),
        '../../test/' + TMP_DIR + '/bhlight_' + str(RES[n])])
copyfile(os.path.join(os.getcwd(), TMP_DIR, 'param_template.dat'), '../../test/' +
         TMP_DIR + '/param.dat')
util.safe_remove(TMP_BUILD)
util.safe_remove(TMP_DIR)
os.chdir('../../test/')

NVAR = 8

L1 = np.zeros(len(RES))

os.chdir(TMP_DIR)

# RUN PROBLEM FOR EACH RESOLUTION AND ANALYZE RESULT
for m in range(len(RES)):
  call_string = ['./bhlight_' + str(RES[m]), '-p', 'param.dat']
  print(call_string)
  bcall(call_string)

  dfiles = np.sort(glob.glob('dumps/dump*.h5'))
  hdr = io.load_hdr(dfiles[0])
  geom = io.load_geom(hdr)
  dump0 = io.load_dump(dfiles[0], geom)
  dump1 = io.load_dump(dfiles[-1], geom)
  r = geom['r'][:,0,0]

  imin = 0
  while r[imin] < hdr['Reh']:
    imin += 1
  rho0 = np.mean(dump0['RHO'][imin:,:,0], axis=1)
  rho1 = np.mean(dump1['RHO'][imin:,:,0], axis=1)

  L1[m] = np.mean(np.fabs(rho1 - rho0))

  files = glob.glob('dumps/*')
  for f in files:
    os.remove(f)
  files = glob.glob('restarts/*')
  for f in files:
    os.remove(f)

# MEASURE CONVERGENCE
powerfit = np.polyfit(np.log(RES), np.log(L1), 1)[0]

os.chdir('../')

if not AUTO:
  # MAKE PLOTS
  fig = plt.figure(figsize=(16.18,10))

  ax = fig.add_subplot(1,1,1)
  ax.plot(RES, L1, marker='s', label='RHO')

  amp = 1.
  ax.plot([RES[0]/2., RES[-1]*2.],
    10.*amp*np.asarray([RES[0]/2., RES[-1]*2.])**-2.,
    color='k', linestyle='--', label='N^-2')
  plt.xscale('log', basex=2); plt.yscale('log')
  plt.xlim([RES[0]/np.sqrt(2.), RES[-1]*np.sqrt(2.)])
  plt.xlabel('N'); plt.ylabel('L1')
  #plt.title(NAMES[MODES[n]])
  plt.legend(loc=1)
  plt.savefig('bondi.png', bbox_inches='tight')

if AUTO:
  data = {}
  #data['SOL'] = -2.*np.zeros([len(MODES), NVAR])
  data['CODE'] = powerfit
  import pickle
  pickle.dump(data, open('data.p', 'wb'))

# CLEAN UP
util.safe_remove(TMP_DIR)

