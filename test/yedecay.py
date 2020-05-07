################################################################################
#                                                                              #
# ONE-ZONE OPTICALLY THIN NEUTRINO COOLING                                     #
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
from scipy import integrate,interpolate,optimize
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'yedecay'
AUTO = '-auto' in sys.argv
ANTINU = '-antinu' in sys.argv
MPI = '-mpi' in sys.argv
BOOST = '-boost' in sys.argv
LEPTON = '-lepton' in sys.argv
LABEL='electron_antinu' if ANTINU else 'electron_nu'
gam = 1.4

# devnull
try:
    from subprocess import DEVNULL # py3k
except ImportError:
    import os
    DEVNULL = open(os.devnull, 'wb')

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
compile_args = [sys.executable,'build.py','-dir',TMP_DIR]
if MPI:
    compile_args += ['-mpi']
if ANTINU:
    compile_args += ['-antinu']
if BOOST:
    compile_args += ['-boost']
if LEPTON:
    compile_args += ['-lepton']
call(compile_args)
call(['mv', 'sc_eos_gamma_{}.h5'.format(str(gam).replace('.','p')),
      TMP_DIR])
os.chdir('../../test')
call(['mv','../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
run_args = ['./bhlight', '-p', 'param_template.dat']
if MPI:
    NPROC = 8
    import psutil
    NCORE = psutil.cpu_count(logical=False)
    NTHREADS = (max(NCORE/NPROC, 1))
    os.environ['OMP_NUM_THREADS'] = '%d' % NTHREADS
    bcall(run_args,int(NPROC),stdin=DEVNULL)
else:
    bcall(run_args)
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = sorted(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
times = np.empty(Nd)
yes = np.empty(Nd)
us = np.empty(Nd)
rhos = np.empty(Nd)
press = np.empty(Nd)
for i,d in enumerate(dfiles):
    dump = io.load_dump(d)
    times[i] = dump['t']
    us[i] = dump['UU'].mean()
    yes[i] = dump['Ye'].mean()
    rhos[i] = dump['RHO'].mean()
    press[i] = dump['PRESS'].mean()
times *= hdr['T_unit']
rhos *= hdr['RHO_unit']
us *= hdr['U_unit']
press *= hdr['U_unit']
ucon0 = dump['ucon'][0,0,0,0]
ucov0 = dump['ucov'][0,0,0,0]

# GET ANALYTIC SOLUTION
#emissivity
cnu_flat = 1
j_unit = hdr['U_unit']
C = 4*np.pi*cnu_flat*j_unit/((hdr['numax']-hdr['numin'])*hdr['T_unit'])
# Ye
Cdecay = -2*C*(np.log(hdr['numax']/hdr['numin']))*((cgs['MP'])/cgs['HPL'])
Cs = (Cdecay/(rhos)).reshape(len(times))
if ANTINU:
    ana = 0.5*(1-np.exp(Cs*times))
else:
    ana = 0.5*np.exp(Cs*times)
# energy
us_ana = np.empty_like(us)
Cu = C*(hdr['numax']-hdr['numin'])
yefunc = interpolate.interp1d(times,yes,kind='cubic')
if ANTINU:
    urhs = lambda t: -Cu*(1-2.*yefunc(t))
else:
    urhs = lambda t: -2.*Cu*yefunc(t)
r = integrate.ode(urhs)
r.set_integrator('dopri5')
r.set_initial_value(us[0],times[0])
for i,t in enumerate(times):
    us_ana[i] = r.integrate(t)
#curve fit
if ANTINU:
    foo = lambda x,b: 0.5*(1-np.exp(b*x))
else:
    foo = lambda x,b: 0.5*np.exp(b*x)
popt,pcov = optimize.curve_fit(foo,times,yes)
slope_code = popt[0]
slope_sol = Cs.mean()

if AUTO:
    data = {}
    data['SOL'] = [np.array([0.0]), np.array([slope_sol])]
    data['CODE'] = [np.array([0.0]), np.array([slope_code])]
    data['THRESHOLD'] = 0.05
    import pickle
    pickle.dump(data, open('data.p','wb'))
    util.safe_remove(TMP_DIR)
    sys.exit()

import matplotlib as mpl; mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':22})
plt.plot(1e8*times,us/1e24,'r-',label = 'data')
plt.plot(1e8*times,us_ana/1e24,'b--',label='analytic')
plt.legend()
plt.xlabel('time (au)')
plt.ylabel(r'energy density (au)')
plt.savefig('yedecay_{}_energy.png'.format(LABEL),
            bbox_inches='tight',transparent=True)
plt.savefig('yedecay_{}_energy.pdf'.format(LABEL),
            bbox_inches='tight',transparent=True)
plt.cla()
plt.clf()

if ANTINU:
    plt.semilogy(1e8*times,0.5-yes,'r-',label='data')
    plt.semilogy(1e8*times,0.5-ana,'b--',label='analytic')
    plt.ylabel(r'$\frac{1}{2} - Y_e$')
else:
    plt.semilogy(1e8*times,yes,'r-',label='data')
    plt.semilogy(1e8*times,ana,'b--',label='analytic')
    plt.ylabel(r'$Y_e$')
plt.xlabel('time (s)')
plt.legend()
plt.xlabel('time (au)')
plt.savefig('yedecay_{}_ye.png'.format(LABEL),
            bbox_inches='tight')
plt.savefig('yedecay_{}_ye.pdf'.format(LABEL),
            bbox_inches='tight')
plt.cla()
plt.clf()

# CLEAN UP
util.safe_remove(TMP_DIR)
