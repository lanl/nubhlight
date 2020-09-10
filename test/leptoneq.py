################################################################################
#                                                                              #
# Equilibrium Lepton Transport                                                 #
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
PROBLEM = 'leptoneq'

# MPI is default
NOMPI = '-nompi' in sys.argv
MPI = not NOMPI or '-mpi' in sys.argv 
AUTO = '-auto' in sys.argv

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
call(compile_args)
os.chdir('../../test')
call(['mv','../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
run_args = ['./bhlight', '-p', 'param_template.dat']
if MPI:
    NPROC = 4
    import psutil
    NCORE = psutil.cpu_count(logical=False)
    NTHREADS = (max(NCORE/NPROC, 1))
    os.environ['OMP_NUM_THREADS'] = '%d' % NTHREADS
    bcall(run_args,int(NPROC),stdin=DEVNULL)
else:
    bcall(run_args)
os.chdir('../')

# load simulation output
dumps = sorted(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr)
N1,N2,N3 = [hdr[i] for i in ['N1','N2','N3']]

times = np.empty(len(dumps))
Yes = np.empty((len(dumps),N1,N2))

dump = io.load_dump(dumps[0])
x,y = geom['xf'][:,:,0],geom['yf'][:,:,0]

#STOPPED EDITING HERE
for i,d in enumerate(dumps):
    dump = io.load_dump(d)
    times[i] = dump['t']
    Yes[i] = dump['Ye'].reshape(N1,N2)    
times *= hdr['T_unit']

Ymaxs = Yes.max(axis=(-2,-1))
Ymins = Yes.min(axis=(-2,-1))
Yavgs = Yes.mean(axis=(-2,-1))
Ydiffs = 100.*(Yes - 0.225)

dYmax = np.abs(Ymaxs - Yavgs)/Yavgs
dYmin = np.abs(Ymins - Yavgs)/Yavgs

if AUTO:
    data = {}
    data['SOL'] = [np.array([0.,1.]), np.array([0.,0.])]
    data['CODE'] = [np.array([0.,1.]), np.array([dYmin[-1],dYmax[-1]])]
    data['THRESHOLD'] = 0.05
    import pickle
    pickle.dump(data, open('data.p', 'wb'))
    util.safe_remove(TMP_DIR)
    sys.exit()

## STOPPED EDITING HERE

import matplotlib as mpl; mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import animation
mpl.rcParams.update({'font.size':18})

print("Plotting time evolution")
plt.plot(1e3*times,Ymaxs,lw=3,label='max')
plt.plot(1e3*times,Ymins,lw=3,label='min')
plt.plot(1e3*times,Yavgs,'k--',lw=3,label='avg')
plt.legend()
plt.xlabel(r'$t$ (ms)')
plt.ylabel(r'$Y_e$')
plt.savefig('leptoneq.pdf',bbox_inches='tight')
plt.savefig('leptoneq.png',bbox_inches='tight')
plt.cla()
plt.clf()

print("Plotting frames")
num_frames = 3
fig, axarr = plt.subplots(1,num_frames,
                          figsize=(4*num_frames+0.25*4,4),
                          sharex=True,
                          sharey=True)
pc = axarr[0].pcolormesh(x,y,Ydiffs[0],
                         cmap='viridis',
                         shading='gourad',
                         linewidth=0,
                         rasterized=True,
                         vmin=-15,vmax=15)
pc.set_edgecolor('face')
axarr[0].set_xlabel(r'$x$')
axarr[0].set_ylabel(r'$y$')
for i in range(1,num_frames-1):
    pc = axarr[i].pcolormesh(x,y,Ydiffs[i*len(times)//(num_frames-1)],
                             cmap='viridis',
                             shading='gourad',
                             linewidth=0,
                             rasterized=True,
                             vmin=-15,vmax=15)
    pc.set_edgecolor('face')
    axarr[i].set_xlabel(r'$x$')
mesh = axarr[-1].pcolormesh(x,y,Ydiffs[-1],cmap='viridis',
                            shading='gourad',
                            linewidth=0,
                            rasterized=True,
                            vmin=-15,vmax=15)
mesh.set_edgecolor('face')
axarr[-1].set_xlabel(r'$x$')

plt.tight_layout()

fig.colorbar(mesh,ax=axarr.ravel().tolist(),
             label=r'$100\times (Y_e - 0.225)$')
# plt.axes().set_aspect('equal','datalim')
plt.savefig('leptoneq_frames.png',bbox_inches='tight')
plt.savefig('leptoneq_frames.pdf',bbox_inches='tight',dpi=300)
plt.cla()
plt.clf()

print("Making movie")
fig,ax = plt.subplots()
ax.set_xlim((-1,1))
ax.set_ylim((-1,1))
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
pc = ax.pcolormesh(x,y,Ydiffs[0],
                   cmap='viridis',shading='gourad',
                   vmin=-15,vmax=15)
cbar = fig.colorbar(pc,label=r'$100\times (Y_e - 0.225)$')
def init():
    pc.set_array(Ydiffs[0].ravel())
    return pc
def animate(i):
    pc.set_array(Ydiffs[i].ravel())
    return pc
anim = animation.FuncAnimation(fig,animate,
                               init_func=init,
                               frames=len(Ydiffs),
                               interval=20,blit=False)
anim.save("leptoneq.mp4",writer='ffmpeg',
          extra_args=['-loglevel','verbose'])

print("Done.")
util.safe_remove(TMP_DIR)
