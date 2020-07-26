################################################################################
#                                                                              #
# ONE-ZONE OPTICALLY THIN BREMSSTRAHLUNG COOLING - MPI                         #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
from bhlight import bcall
import util

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM='brem'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call([sys.executable, 'build_mpi.py', '-dir', TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
bcall(['./bhlight','-p','param_template.dat'],8)
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
t_code = np.zeros(Nd)
Te_code = np.zeros(Nd)
for n in xrange(Nd):
  dump = io.load_dump(dfiles[n], geom)
  t_code[n] = dump['t']*hdr['T_unit']
  Te_avg = np.mean(dump['Thetae']*cgs['ME']*cgs['CL']**2/cgs['KBOL'])
  Te_code[n] = Te_avg

# GET ANALYTIC SOLUTION
tf = 1.e8
Te0 = 1.e8
dump = io.load_dump(dfiles[0], geom)
ne = dump['RHO'].mean()*hdr['Ne_unit']
gam = 5./3.
#N  = 5.4e-39 # cm^3 K^1/2 s^-1 Sr^-1 Hz^-1
t_sol = np.linspace(0, tf, 1024)
from scipy.integrate import odeint
def func(Te, t):
  gff = 1.2
  rel = 1. + 4.4e-10*Te
  J  = np.sqrt(2*np.pi*cgs['KBOL']*Te/(3*cgs['ME']))
  J *= 2**5*np.pi*cgs['QE']**6/(3*cgs['HPL']*cgs['ME']*cgs['CL']**3)
  J *= ne**2*gff*rel
  return -(gam-1.)*J/(2.*ne*cgs['KBOL'])
Te_sol = odeint(func, Te0, t_sol)

# MAKE FIGURE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,1,1)
ax.plot(t_code, Te_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(t_sol, Te_sol, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('t (s)'); plt.ylabel('Te (K)')
plt.xlim([0,256*hdr['T_unit']]); plt.ylim([0, 1.1e8])

plt.savefig('brem_mpi.png', bbox_inches='tight')

# CLEAN UP
call(['rm', '-rf', TMP_DIR])

