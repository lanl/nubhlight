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
restarts = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/restarts/restart*.h5'))
restart = h5py.File(restarts[-1], 'r')

b = 0
mu = geom['local_angles_mu'][:]
angles = [d['local_angles'][b].sum(axis=(0,1)) for d in dumps]
moments = dumps[0]['local_moments'][b].sum(axis=-1).sum(axis=-1)
ntot = [a.sum() for a in angles]
A = moments[0]
B = moments[1]

Gnu = [d['Gnu'][b].sum(axis=(0,1)) for d in dumps]
Gdiff = (Gnu[1].sum() - Gnu[0].sum())
Greldiff = 2.0*Gdiff/(np.abs(Gnu[1].sum()) + np.abs(Gnu[0].sum()) + 1e-20)

Nnu_tot = 1e55
dmu = mu[1] - mu[0]
fracs = [0.2, 0.4, 0.2, 0.2]
e0 = dmu*Nnu_tot*fracs[0]*(2./3.)*(1 - (3./4.)*mu*mu)
ebar0 = dmu*Nnu_tot*fracs[1]*((1./4.) + (3./4.)*mu*mu)
x0 = dmu*fracs[2]*Nnu_tot*np.ones_like(mu)/2
xbar0 = dmu*fracs[3]*Nnu_tot*np.ones_like(mu)/2

Amask = Gnu[0] < 0
Bmask = Gnu[0] > 0
e1 = e0.copy()
e1[Amask] = (1 - (2./3.)*(B/A))*e0[Amask] + (1./3.)*(B/A)*x0[Amask]
e1[Bmask] = (1./3.)*(e0[Bmask] + x0[Bmask])
x1 = x0.copy()
x1[Amask] = (2./3.)*(B/A)*e0[Amask] + (1 - (1./3.)*(B/A))*x0[Amask]
x1[Bmask] = (2./3.)*(e0[Bmask] + x0[Bmask])

# Check that osc_count is non-trivial
osc_count = np.array(restart['superphotons'].fields('osc_count'))
weights = np.array(restart['superphotons'].fields('w'))
print("Number of physical oscillations:",np.sum(osc_count * weights))

if AUTO:
    data = {}
    data['SOL'] = [0]
    data['CODE'] = Gdiff
    data['THRESHOLD'] = 0.005
    import pickle
    pickle.dump(data, open('data.p', 'wb'))
    # clean up
    util.safe_remove(TMP_DIR)
    sys.exit()

print("Ntot = ",ntot)
print("Ntot per species (initial) = ", angles[0].sum(axis=-1))
print("Ntot per species (final) = ", angles[1].sum(axis=-1))
print("Frac per species (initial) = ", angles[0].sum(axis=-1)/ntot[0])
print("Frac per species (final) = ", angles[1].sum(axis=-1)/ntot[1])
print("frac e+x, frac ebar+xbar = ",
      (angles[1].sum(axis=-1)/ntot[1])[0]+(angles[1].sum(axis=-1)/ntot[1])[2],
      (angles[1].sum(axis=-1)/ntot[1])[1]+(angles[1].sum(axis=-1)/ntot[1])[3])
print("Gnu start, Gnu end = ",Gnu[0].sum(), Gnu[1].sum())
print("Relative diff in Gnu = {}".format(Greldiff))

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

cmap = plt.get_cmap("tab10")
for i,ax in enumerate(axarr[1:]):
#     ax.plot(mu, angles[i][0]/1e53,
#             label=r'$\nu_e$')
#     ax.plot(mu, angles[i][1]/1e53,
#             #linestyle='--',
#             label=r'$\bar{\nu}_e$')
#     ax.plot(mu, angles[i][2]/1e53,
#             label=r'$\nu_X$')
#     ax.plot(mu, angles[i][3]/1e53,
#             #linestyle='--',
#             label=r'$\bar{\nu}_x$')
    ax.set_ylim(0,1.0)
axarr[1].plot(mu, angles[0][0]/1e53,
              label=r'$\nu_e$')
axarr[2].plot(mu, angles[1][0]/1e53)

axarr[1].plot(mu, angles[0][2]/1e53,
              label=r'$\nu_x$')
axarr[2].plot(mu, angles[1][2]/1e53)

axarr[2].plot(mu, e1/1e53, #color=cmap(1),
              linestyle=':', label=r'$\nu_e$ analytic')
axarr[2].plot(mu, x1/1e53, #color=cmap(1),
              linestyle=':', label=r'$\nu_x$ analytic')
axarr[1].legend(loc='lower center',
                ncol=1)
axarr[2].legend(loc='lower center',
                ncol=1)


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
