################################################################################
#                                                                              #
# SOD SHOCKTUBE                                                                #
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
import hdf5_to_dict as io
import util
from bhlight import bcall
from scipy.optimize import root

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'sod'
AUTO = '-auto' in sys.argv
TABLE = '-table' in sys.argv
FORCE = '-force' in sys.argv
RELTABLE = '-reltable' in sys.argv
tscale = 1.e-2
gam = 1.4
if AUTO and not RELTABLE:
  import pickle  
else:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
args = [sys.executable, 'build.py', '-dir', TMP_DIR]
if TABLE:
  args.append('-table')
if FORCE:
  args.append('-force')
if RELTABLE:
  args.append('-reltable')
  tscale = 1.0
call(args)
if TABLE:
  call(['mv', 'sc_eos_gamma_{}.h5'.format(str(gam).replace('.','p')), TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
bcall(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
dump = io.load_dump(dfiles[-1], geom)

x_code   = geom['x'][:,0,0]
rho_code = dump['RHO'][:,0,0]
P_code   = dump['PRESS'][:,0,0]/(tscale*tscale)
u_code   = dump['U1'][:,0,0]/(tscale)
eps_code = (dump['UU'][:,0,0]/dump['RHO'][:,0,0])/(tscale*tscale)
if RELTABLE:
  ye_code = dump['Ye'][:,0,0]

if RELTABLE:
  import matplotlib as mpl; mpl.use('Agg')
  import matplotlib.pyplot as plt
  fig = plt.figure(figsize=(16.18,10))

  ax = fig.add_subplot(2,2,1)
  ax.plot(x_code,rho_code,label='density')
  ax.plot(x_code,ye_code,'r--',label=r'$Y_e$')
  plt.legend()
  plt.xlabel('x'); plt.ylabel('Density')

  ax = fig.add_subplot(2,2,2)
  ax.plot(x_code,P_code)
  plt.xlabel('x'); plt.ylabel('Pressure')

  ax = fig.add_subplot(2,2,3)
  ax.plot(x_code,u_code)
  plt.xlabel('x'); plt.ylabel('Velocity')

  ax = fig.add_subplot(2,2,4)
  ax.plot(x_code,eps_code)
  plt.xlabel('x'); plt.ylabel('Specific Internal Energy')

  plt.subplots_adjust(wspace=0.15)
  plt.savefig('sod_rel.png', bbox_inches='tight')
  util.safe_remove(TMP_DIR)
  sys.exit()

# GET ANALYTIC SOLUTION (ADAPTED FROM BRUCE FRYXELL'S exact_riemann.f)
x0 = 0.5
t = 0.25
rho1 = 1.; P1 = 1.; u1 = 0.
rho5 = 0.125; P5 = 0.1; u5 = 0.
cs1 = np.sqrt(gam*P1/rho1)
cs5 = np.sqrt(gam*P5/rho5)
Gam = (gam-1.)/(gam+1.)
beta = (gam-1.)/(2.*gam)
def func(x):
  P3 = x[0]
  u4 = (P3-P5)*np.sqrt((1.-Gam)/(rho5*(P3 + Gam*P5)))
  u2 = (P1**beta - P3**beta)*np.sqrt((1.-Gam**2.)*P1**(1./gam)/(Gam**2*rho1))
  return u2-u4
P3 = root(func, [(P1+P5)/2.]).x[0]
P4 = P3
rho3 = rho1*(P3/P1)**(1./gam)
rho4 = rho5*(P4 + Gam*P5)/(P5 + Gam*P4)
u4 = cs5*(P4/P5-1)/(gam*np.sqrt(1. + (gam+1.)/(2.*gam)*(P4/P5-1.)))
ushock = cs5*np.sqrt(1. + (gam+1.)/(2.*gam)*(P4/P5-1.))
u3 = u4
cs3 = np.sqrt(gam*P3/rho3)
xsh = x0 + ushock*t
xcd = x0 + u3*t
xft = 0.5 + (u3-cs3)*t
xhd = 0.5 - cs1*t
N = 1024
x = np.linspace(0, 1, 1024)
rho = np.zeros(N)
P   = np.zeros(N)
u   = np.zeros(N)
for n in range(N):
  if x[n] < xhd:
    rho[n] = rho1
    P[n]   = P1
    u[n]   = u1
  elif x[n] < xft:
    u[n]   = 2./(gam+1.)*(cs1 + (x[n] - x0)/t)
    fac    = 1. - 0.5*(gam-1.)*u[n]/cs1
    rho[n] = rho1*fac**(2./(gam-1.))
    P[n]   = P1*fac**(2.*gam/(gam-1.))
  elif x[n] < xcd:
    rho[n] = rho3
    P[n]   = P3
    u[n]   = u3
  elif x[n] < xsh:
    rho[n] = rho4
    P[n]   = P4
    u[n]   = u4
  else:
    rho[n] = rho5
    P[n]   = P5
    u[n]   = u5

if AUTO:
  data = {}
  data['SOL'] = [x, rho]
  data['CODE'] = [x_code, rho_code]
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit() 

# MAKE FIGURE
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(2,2,1)
ax.plot(x_code, rho_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(x, rho, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Density')
plt.ylim([0, 1.1])

ax = fig.add_subplot(2,2,2)
ax.plot(x_code, P_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(x, P, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Pressure')
plt.ylim([0, 1.1])

ax = fig.add_subplot(2,2,3)
ax.plot(x_code, u_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(x, u, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Velocity')
plt.ylim([0, 1.1])

ax = fig.add_subplot(2,2,4)
ax.plot(x_code, P_code/rho_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(x, P/rho, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('x'); plt.ylabel('Temperature')
plt.ylim([0.7, 1.2])

plt.subplots_adjust(wspace=0.15)

plt.savefig('sod.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

