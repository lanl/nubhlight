################################################################################
#                                                                              #
# TRACER PARTICLES IN 1D                                                       #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
from shutil import copyfile
import numpy as np
import hdf5_to_dict as io
import util
from bhlight import bcall
import multiprocessing
import matplotlib as mpl; mpl.use('Agg')
from matplotlib import animation, rc
from matplotlib import pyplot as plt
rc('font',size=18)

num_mpi = 4
num_cpus = multiprocessing.cpu_count()
os.environ['OMP_NUM_THREADS'] = str(int(np.max([2,num_cpus/num_mpi])))

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)

PROBLEM = 'torus_cbc'
AUTO = '-auto' in sys.argv

if AUTO:
    raise ValueError("No auto for this test")

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call([sys.executable, 'build.py', '-dir', TMP_DIR,
      '-classic', '-nob', '-norenorm', '-tracertest'])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
bcall(['./bhlight', '-p', 'param_template.dat'],int(num_mpi))
os.chdir('../')

# READ DATA
tracers = io.TracerData.fromdir(os.path.join(TMP_DIR,'dumps','tracers'))

t0 = tracers.filter(tracers['time']==0.0)
rcyl0 = np.sqrt(t0['Xcart'][:,0]**2 + t0['Xcart'][:,1]**2)

print("Making movie")
fig,ax = plt.subplots()
ax.set_xlim(0,40)
ax.set_ylim(-17,17)
ax.set_xlabel(r'$r_{cyl}$')
ax.set_ylabel(r'$z$')
dots, = ax.plot(rcyl0,t0['Xcart'][:,2],'bo',ms=1)
plt.tight_layout()

def init():
    dots.set_data(rcyl0,t0['Xcart'][:,2])
    return dots,

def animate(i):
    trc = tracers.filter(tracers['time'] == tracers.times()[i])
    rcyl = np.sqrt(trc['Xcart'][:,0]**2 + trc['Xcart'][:,1]**2)
    dots.set_data(rcyl,trc['Xcart'][:,2])
    return dots,

anim = animation.FuncAnimation(fig,animate,init_func=init,
                               frames=len(tracers.times()),
                               interval=10,blit=True)

anim.save('tracers3d_relax.mp4')

plt.clf()
plt.cla()

print("Making orbit plot")
id_indcs = [200, 400, 600, 800]
ids = [tracers.ids()[i] for i in id_indcs]
traces = [tracers.get_trace(id) for id in ids]
for i,trace in zip(id_indcs,traces):
    plt.plot(trace['Xcart'][:,0],trace['Xcart'][:,1],label='id = {}'.format(i))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.axes().set_aspect('equal','datalim')
plt.savefig('tracers3d_orbit.png',bbox_inches='tight')
plt.savefig('tracers3d_orbit.pdf',bbox_inches='tight')
