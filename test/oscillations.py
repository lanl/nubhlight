################################################################################
#                                                                              #
# UNIT TEST FOR Zaizen/Nakamura Oscillations                                   #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob, util, h5py
import numpy as np
import hdf5_to_dict as io
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'oscillations'
AUTO = '-auto' in sys.argv
MPI = '-mpi' in sys.argv
gam = 1.4

# devnull
try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
args = [sys.executable, 'build.py', '-dir', TMP_DIR]#, '-idim', IDIM]
if MPI:
  args += ['-mpi']
print(args)
call(args)
call(['mv', 'sc_eos_gamma_{}.h5'.format(str(gam).replace('.','p')),
      TMP_DIR])
os.chdir('../../test')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# Since this test is designed to run on a single machine (no batch scripts)
# set openmpi to only use a few threads. Let MPI handle the rest.
if MPI:
  import multiprocessing
  num_mpi = 8
  num_cpus = multiprocessing.cpu_count()
  os.environ['OMP_NUM_THREADS'] = str(int(np.max([2,num_cpus/num_mpi])))

# RUN EXECUTABLE
os.chdir(TMP_DIR)
args = ['./bhlight', '-p', 'param_template.dat']
if MPI:
  bcall(args,int(num_mpi),stdin=DEVNULL)
else:
  bcall(args)
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
geom = h5py.File(os.path.join(TMP_DIR,'')+'/dumps/grid.h5','r')
dumps = [h5py.File(f, 'r') for f in dfiles]

b = 0
mu = geom['local_angles_mu'][:]
angles = [d['local_angles'][b].sum(axis=(0,1)) for d in dumps]
Gnu = [(a[0] - a[1]) - (a[2] - a[3]) for a in angles]
diff = (Gnu[1].sum() - Gnu[0].sum())
reldiff = 2.0*diff/(np.abs(Gnu[1].sum()) + np.abs(Gnu[0].sum()) + 1e-20)
print("Relative diff in Gnu = {}".format(reldiff))

if AUTO:
    data = {}
    data['SOL'] = [0]
    data['CODE'] = diff
    data['THRESHOLD'] = 0.005
    import pickle
    pickle.dump(data, open('data.p', 'wb'))
    # clean up
    util.safe_remove(TMP_DIR)
    sys.exit()

# make figure
import matplotlib as mpl; mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':18})

fig,axarr = plt.subplots(1,3,sharex=True,figsize=(12,4))
axarr[0].plot(mu, Gnu[0]/1e53, label='initial')
axarr[0].plot(mu, Gnu[1]/1e53, label='final')
axarr[0].hlines(0,-1.1,1.1,color='k',linestyle='--')
axarr[0].legend()
axarr[0].set_ylabel(r'particle number$/10^{53}$')
axarr[0].set_xlim(-1.01, 1.01)

for i,ax in enumerate(axarr[1:]):
    ax.plot(mu, angles[i][0]/1e53,label=r'$\nu_e$')
    ax.plot(mu, angles[i][1]/1e53,linestyle='--',label=r'$\bar{\nu}_e$')
    ax.plot(mu, angles[i][2]/1e53,label=r'$\nu_X$')
    ax.plot(mu, angles[i][3]/1e53,linestyle='--',label=r'$\bar{\nu}_x$')
axarr[2].legend()

for ax in axarr:
    ax.set_xlabel(r'$\cos(\theta)$')

axarr[0].set_title(r'$G_\nu$')
axarr[1].set_title(r'$f(t = 0)$')
axarr[2].set_title(r'$f(t = t_{final})$')

plt.tight_layout()

plt.savefig('oscillations_1zone.png', bbox_inches='tight')
plt.savefig('oscillations_1zone.pdf', bbox_inches='tight')

# clean up
util.safe_remove(TMP_DIR)
