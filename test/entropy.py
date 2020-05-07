################################################################################
#                                                                              #
# ADVECTED ENTROPY MODE                                                        #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
import glob
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
sys.path.insert(0, '../script/analysis/')
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
sys.path.insert(0, '../script/')
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'entropy'

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call([sys.executable, 'build.py', '-dir', TMP_DIR])
os.chdir('../../test')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
executable = host.get('EXECUTABLE','')
call_string = ['./bhlight', '-p', 'param_template.dat']
bcall(call_string)
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)
x_code = geom['x'][:,0,0]
rho_code = dump['RHO'][:,0,0]
drho_code = rho_code - np.mean(rho_code)

# GET ANALYTIC SOLUTION
kwave = 2.*np.pi
omega = 0.
amp = 0.01
# rho, u, u1
mean = np.array([1., 0.01, 0.1])
delta = np.array([1., 0., 0.])

x_sol = np.linspace(0, 1, 1000)

mode = amp*np.cos(kwave*x_sol)
def mode(x, t):
  return np.real(amp*np.exp(1j*(kwave*(x-mean[2]*t) - omega*t)))
d_sol_init = []
d_sol = []
for n in xrange(len(mean)):
  d_sol_init.append(delta[n]*mode(x_sol, 0.))
  d_sol.append(delta[n]*mode(x_sol, dump['t']))

# MAKE FIGURE
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,1,1)
ax.plot(x_code, drho_code, color=code_col, linestyle=code_ls, marker=code_mrk, label='Code')
ax.plot(x_sol, d_sol_init[0], color=sol_col, linestyle=':', marker=sol_mrk, label='Initial')
ax.plot(x_sol, d_sol[0], color=sol_col, linestyle=sol_ls, marker=sol_mrk, label='Final')
plt.xlabel('x'); plt.ylabel('delta rho')
plt.xlim([0,hdr['startx1'] + hdr['dx1']*hdr['N1']]); 
plt.ylim([-1.5*amp*delta[0], 1.5*amp*delta[0]])
plt.legend(loc=1)

plt.savefig(PROBLEM + '.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

