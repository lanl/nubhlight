from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
from shutil import copyfile
import glob
import numpy as np
import multiprocessing
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

AUTO = '-auto' in sys.argv
TABLE = '-table' in sys.argv
MPI = "-mpi" in sys.argv

RES = [16, 32, 64, 128]

GAMMA = 4./3.
if TABLE:
    GAMMA_NAME='sc_eos_gamma_{}.h5'.format(str(GAMMA).replace('.','p'))

# Since this test is designed to run on a single machine (no batch scripts)
# set openmpi to only use a few threads. Let MPI handle the rest.
if MPI:
  nxcpu = 2
  num_mpi = 8
  num_cpus = multiprocessing.cpu_count()
  os.environ['OMP_NUM_THREADS'] = str(int(np.max([2,num_cpus/num_mpi])))
else:
  nxcpu = 1

util.make_dir(TMP_DIR)
os.chdir('../prob/mhdmodes3d/')

copyfile('build.py', TMP_BUILD)
# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE

for n in range(len(RES)):
  util.change_cparm('N1TOT', RES[n], TMP_BUILD)
  util.change_cparm('N2TOT', RES[n], TMP_BUILD)
  util.change_cparm('N3TOT', RES[n], TMP_BUILD)
  util.change_cparm('N1CPU', nxcpu,  TMP_BUILD)
  util.change_cparm('N2CPU', nxcpu,  TMP_BUILD)
  util.change_cparm('N3CPU', nxcpu,  TMP_BUILD)
  util.change_cparm('OPENMP', 'True', TMP_BUILD)
  args = [sys.executable, TMP_BUILD, '-dir', TMP_DIR]
  if TABLE:
    args.append('-table')
  call(args)
  call(['cp', os.path.join(os.getcwd(), TMP_DIR, 'bhlight'),
        '../../test/' + TMP_DIR + '/bhlight_' + str(RES[n])])
copyfile(os.path.join(os.getcwd(), TMP_DIR, 'param_template.dat'),
         '../../test/' + TMP_DIR + '/param.dat')
if TABLE:
    copyfile(os.path.join(os.getcwd(),GAMMA_NAME),
             os.path.join(os.getcwd(),'../../test/',TMP_DIR,GAMMA_NAME))
util.safe_remove(TMP_BUILD)
util.safe_remove(TMP_DIR)
os.chdir('../../test/')

# LOOP OVER EIGENMODES
MODES = [0, 1, 2, 3]
NAMES = ['ENTROPY', 'SLOW', 'ALFVEN', 'FAST']
MODES = [3]
NVAR = 8
VARS = ['rho', 'u', 'u1', 'u2', 'u3', 'B1', 'B2', 'B3']

amp = 1.e-4
k1 = 2.*np.pi
k2 = 2.*np.pi
k3 = 2.*np.pi
var0 = np.zeros(NVAR)
var0[0] = 1.
var0[1] = 1.
var0[5] = 1.
L1 = np.zeros([len(MODES), len(RES), NVAR])
powerfits = np.zeros([len(MODES), NVAR])

for n in range(len(MODES)):
  util.change_rparm('nmode', MODES[n], TMP_DIR + '/param.dat')
  os.chdir(TMP_DIR)

  # EIGENMODES
  dvar = np.zeros(NVAR)
  if MODES[n] == 0: # ENTROPY
    dvar[0] = 1.
  if MODES[n] == 1: # SLOW/SOUND
    dvar[0] = 0.556500332363
    dvar[1] = 0.742000443151
    dvar[2] = -0.282334999306
    dvar[3] = 0.0367010491491
    dvar[4] = 0.0367010491491
    dvar[5] = -0.195509141461
    dvar[6] = 0.0977545707307
    dvar[7] = 0.0977545707307
  if MODES[n] == 2: # ALFVEN
    dvar[3] = -0.339683110243
    dvar[4] = 0.339683110243
    dvar[6] = 0.620173672946
    dvar[7] = -0.620173672946
  if MODES[n] == 3: # FAST
    dvar[0] = 0.481846076323
    dvar[1] = 0.642461435098
    dvar[2] = -0.0832240462505
    dvar[3] = -0.224080007379
    dvar[4] = -0.224080007379
    dvar[5] = 0.406380545676
    dvar[6] = -0.203190272838
    dvar[7] = -0.203190272838
  dvar *= amp
  
  # RUN PROBLEM FOR EACH RESOLUTION AND ANALYZE RESULT
  for m in range(len(RES)):
    print("Res = {}".format(RES[m]))
    args = ['./bhlight_' + str(RES[m]), '-p', 'param.dat']
    if MPI:
      bcall(args,int(num_mpi))
    else:
      call(args)

    dfiles = np.sort(glob.glob('dumps/dump*.h5'))
    hdr = io.load_hdr(dfiles[-1])
    geom = io.load_geom(hdr, recalc=True)
    dump = io.load_dump(dfiles[-1], geom) 
    X1 = geom['x'][:,:,:]
    X2 = geom['y'][:,:,:]
    X3 = geom['z'][:,:,:]
    dvar_code = []
    dvar_code.append(dump['RHO'][:,:,:] - var0[0]) 
    dvar_code.append(dump['UU'][:,:,:]  - var0[1])
    dvar_code.append(dump['U1'][:,:,:]  - var0[2])
    dvar_code.append(dump['U2'][:,:,:]  - var0[3])
    dvar_code.append(dump['U3'][:,:,:]  - var0[4])
    dvar_code.append(dump['B1'][:,:,:]  - var0[5])
    dvar_code.append(dump['B2'][:,:,:]  - var0[6])
    dvar_code.append(dump['B3'][:,:,:]  - var0[7])

    #dvar_sol = []
    dvar_sol = np.zeros([RES[m], RES[m], RES[m]])
    for k in range(NVAR):
      #dvar_sol.append(np.real(dvar[k])*np.cos(k1*X1 + k2*X2))
      if abs(dvar[k]) != 0.:
        for i in range(RES[m]):
          for j in range(RES[m]):
            for kk in range(RES[m]):
              dvar_sol[i,j,kk] = np.real(dvar[k])*np.cos(k1*X1[i,j,kk] + 
                                                         k2*X2[i,j,kk] + 
                                                         k3*X3[i,j,kk])
              L1[n][m][k] = np.mean(np.fabs(dvar_code[k][i,j,kk] - 
                                            dvar_sol[i,j,kk]))

  # MEASURE CONVERGENCE
  for k in range(NVAR):
    if abs(dvar[k]) != 0.:
      powerfits[n,k] = np.polyfit(np.log(RES), np.log(L1[n,:,k]), 1)[0]
  
  os.chdir('../')

  if not AUTO:
    # MAKE PLOTS
    fig = plt.figure(figsize=(16.18,10))

    ax = fig.add_subplot(1,1,1)
    for k in range(NVAR):
      if abs(dvar[k]) != 0.:
        ax.plot(RES, L1[n,:,k], marker='s', label=VARS[k])
 
    ax.plot([RES[0]/2., RES[-1]*2.], 
      10.*amp*np.asarray([RES[0]/2., RES[-1]*2.])**-2.,
      color='k', linestyle='--', label='N^-2')
    plt.xscale('log', basex=2); plt.yscale('log')
    plt.xlim([RES[0]/np.sqrt(2.), RES[-1]*np.sqrt(2.)])
    plt.xlabel('N'); plt.ylabel('L1')
    plt.title(NAMES[MODES[n]])
    plt.legend(loc=1)
    plt.savefig('mhdmodes3d_' + NAMES[MODES[n]] + '.png', bbox_inches='tight')

if AUTO:
  data = {}
  data['SOL'] = -2.*np.zeros([len(MODES), NVAR])  
  data['CODE'] = powerfits
  import pickle
  pickle.dump(data, open('data.p', 'wb'))

# CLEAN UP
util.safe_remove(TMP_DIR)

