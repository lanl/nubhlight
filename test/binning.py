################################################################################
#                                                                              #
# UNIT TEST FOR SUPERPHOTON BINNING                                            #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
from scipy import optimize
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'binning'
AUTO = '-auto' in sys.argv
MPI = '-mpi' in sys.argv
if '-idim' in sys.argv:
  IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
  IDIM = 0

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
args = [sys.executable, 'build.py', '-dir', TMP_DIR, '-idim', str(IDIM)]
if MPI:
  args += ['-mpi']
call(args)
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
  bcall(args,str(num_mpi))
else:
  bcall(args)
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)

nuLnu = dump['nuLnu']
if np.any(np.abs(nuLnu[1:,:,:,:]) >= 1e-15):
  nuLnu = np.sum(nuLnu,axis=0)
else:
  nuLnu = nuLnu[0]

def make_1d(tensor,axis):
  assert axis in [0,1,2]
  assert len(tensor.shape) == 3
  if axis == 0:
    return np.sum(np.sum(tensor,axis=-1),axis=-1)
  if axis == 1:
    return np.sum(np.sum(tensor,axis=0),axis=-1)
  if axis == 2:
    return np.sum(np.sum(tensor,axis=0),axis=0)

# GET ANALYTIC SOLUTION
theta = np.linspace(0,np.pi,nuLnu.shape[0])
phi = np.linspace(0,2*np.pi,nuLnu.shape[1])
lnu = np.linspace(hdr['lnumin'],hdr['lnumax'],nuLnu.shape[2])
THETA,PHI = np.meshgrid(theta,phi,indexing='ij')

foo = lambda x,A,mu,sigma: A*np.exp(-(x-mu)**2/(2*(sigma**2)))
popt,pcov = optimize.curve_fit(foo,lnu,make_1d(nuLnu,2),p0=(1e56,30,3))
lnumu = 0.5*(hdr['lnumax']+hdr['lnumin'])
lnusigma = (hdr['lnumax']-hdr['lnumin'])/(2.*4)

if AUTO:
  data = {}
  data['SOL'] = [np.array([1.,2.]), np.array([lnumu,lnusigma])]
  data['CODE'] = [np.array([1.,2.]), np.array([popt[1],popt[2]])]
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit()

# MAKE FIGURE
errmu = 100.*np.abs(lnumu-popt[1])/popt[1]
errsigma = 100.*np.abs(lnusigma-popt[2])/popt[2]
print("mu error = {:.05}%\nsigma error = {:.03}%.".format(errmu,errsigma))

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':18})
dstring = ["isotropic","x","y","z"]

fig, (ax1,ax2) = plt.subplots(2,figsize=(12,16))
mesh = ax1.pcolormesh(PHI,THETA,nuLnu.sum(axis=-1)/1e54)
ax1.set_xlabel(r'$\phi$')
ax1.set_ylabel(r'$\theta$')
cbar = plt.colorbar(mesh,cmap='viridis',label='count'+r'$\times 10^{54}$')
ax2.bar(lnu,make_1d(nuLnu,2)/1e56)
ax2.set_xlabel(r'$\ln(\nu)$')
ax2.set_ylabel('count'+r'$\times 10^{56}$')
ax1.set_title("Bins In Direction and Frequency, Initial direction = {}".format(dstring[IDIM]))
plt.tight_layout()
plt.savefig('binning_idim_{}.png'.format(IDIM),bbox_inches='tight')
plt.savefig('binning_idim_{}.pdf'.format(IDIM),bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

