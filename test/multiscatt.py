################################################################################
#                                                                              #
# UNIT TEST FOR MULTIPLE SCATTERING BIASES                                     #
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
from scipy import integrate
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'multiscatt'
AUTO = '-auto' in sys.argv
MPI = '-mpi' in sys.argv
NOREBALANCE = '-norebalance' in sys.argv
gam = 1.4

# devnull
try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
args = [sys.executable, 'build.py', '-dir', TMP_DIR]
if MPI:
  args += ['-mpi']
if NOREBALANCE:
  args += ['-norebalance']
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
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)
nuLnu = dump['nuLnu']
onescatts = nuLnu[1].sum(axis = -1)
NTH,NPHI = onescatts.shape[0],onescatts.shape[1]
theta = np.linspace(0,np.pi,NTH+1)
theta = 0.5*(theta[1:] + theta[:-1])
onescatts_th = onescatts.sum(axis=-1)
norm_data = integrate.trapz(onescatts_th, theta)
onescatts_th /= norm_data

# GET ANALYTIC SOLUTION
imax = 3
ms_foo = lambda i,th: (2*2*i+1)*(1-np.cos(th)**(2*i+1))
ana = sum([ms_foo(i,theta) for i in range(imax)])*np.sin(theta)
norm_ana = integrate.trapz(ana,theta)
ana /= norm_ana

if AUTO:
  data = {}
  data['SOL'] = [theta, ana]
  data['CODE'] = [theta, onescatts_th]
  data['THRESHOLD'] = 0.1
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit()

# MAKE FIGURE
import matplotlib as mpl; mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':18})
plt.plot(theta,onescatts_th,'r-',lw=4,label='measured')
plt.plot(theta,ana,'k--',lw=2,label='analytic')
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$\int pdf(\theta,\phi) sin(\theta)d\phi$')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(8,5)
plt.savefig('multiscatt.png',bbox_inches='tight')
plt.savefig('multiscatt.pdf',bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)
