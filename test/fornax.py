
#                                                                              #
# COMPARISON TO FORNAX                                                         #
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
PROBLEM = 'fornax_1zone'
MPI = '-mpi' in sys.argv
HDF = '-hdf' in sys.argv
AUTO = '-auto' in sys.argv
EQUIL = '-equil' in sys.argv
FAST = '-fast' in sys.argv
FORCE = '-force' in sys.argv

OUTNAME='equilibrium' if EQUIL else 'cooling'
INNAME ='fornax_equil.dat' if EQUIL else 'fornax_cooling.dat'

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
if EQUIL:
    compile_args += ['-equil']
if HDF:
    compile_args += ['-hdf']
if FAST:
    compile_args += ['-fast']
if FORCE:
    compile_args += ['-force']
call(compile_args)
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

# load fornax data
fornax = np.loadtxt('../data/{}'.format(INNAME))
ft, fT, fye = fornax.transpose()

# load simulation output
dumps = sorted(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr)
N1,N2,N3 = [hdr[i] for i in ['N1','N2','N3']]

times = np.empty(len(dumps))
yes = np.empty(len(dumps))
us = np.empty(len(dumps))
rhos = np.empty(len(dumps))
press = np.empty(len(dumps))
temps = np.empty(len(dumps))
for i,d in enumerate(dumps):
    dump = io.load_dump(d)
    times[i] = dump['t']
    us[i] = dump['UU'].mean()
    yes[i] = dump['Ye'].mean()
    rhos[i] = dump['RHO'].mean()
    press[i] = dump['PRESS'].mean()
    temps[i] = dump['TEMP'].mean()
times *= hdr['T_unit']
rhos *= hdr['RHO_unit']
us *= hdr['U_unit']
#temps /= hdr['TEMP_unit']
press *= hdr['U_unit']
ucon0 = dump['ucon'][0,0,0,0]
ucov0 = dump['ucov'][0,0,0,0]

if AUTO:
    data = {}
    data['SOL'] = [ft, fye]
    data['CODE'] = [times, yes]
    import pickle
    pickle.dump(data, open('data.p', 'wb'))
    util.safe_remove(TMP_DIR)
    sys.exit()

import matplotlib as mpl; mpl.use('Agg')
from matplotlib import pyplot as plt
mpl.rcParams.update({'font.size':18})


fT_int = interpolate.interp1d(ft,fT)
fye_int = interpolate.interp1d(ft,fye)
dTperc = 100.*np.abs(temps - fT_int(times))/temps
dYperc = 100.*np.abs(yes - fye_int(times))/yes

plt.plot(times,dTperc,label=r'$T$')
plt.plot(times,dYperc,label=r'$Y_e$')
plt.xlabel('time (s)')
plt.ylabel('% Difference')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(6,4,forward=True)
plt.savefig('fornax-{}-percent-difference.png'.format(OUTNAME),
            bbox_inches='tight')
plt.savefig('fornax-{}-percent-difference.pdf'.format(OUTNAME),
            bbox_inches='tight')

plt.clf()
plt.cla()

fig, (axT, axye) = plt.subplots(2,figsize=(6,8),sharex=True)
axT.plot(ft,fT,label='FORNAX')
axT.plot(times,temps,ls='--',label=r'$\nu {\tt bhlight}$')
axT.set_ylabel('T (MeV)')
axT.legend()
axye.plot(ft, fye,label='${\tt FORNAX}$')
axye.plot(times,yes,ls='--',label=r'$\nu {\tt bhlight}$')
axye.set_ylabel(r'$Y_e$')
axye.set_xlabel('time (s)')
plt.savefig('fornax-{}-comparison.png'.format(OUTNAME),
            bbox_inches='tight')
plt.savefig('fornax-{}-comparison.pdf'.format(OUTNAME),
            bbox_inches='tight')
