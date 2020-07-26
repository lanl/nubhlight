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
import matplotlib as mpl; mpl.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plt
rc('font',size=18)

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)

PROBLEM = 'advection1d'
AUTO = '-auto' in sys.argv
if '-idim' in sys.argv:
  IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
  IDIM = 1
RES = 16

kwave = 2*np.pi
amp = 1.0
kcenter = 0.5
cadv = 0.5

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call([sys.executable, 'build.py', '-dir', TMP_DIR, '-idim',
      str(IDIM),'-ntot',str(RES),'-tracers'])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
bcall(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ DATA
tracers = io.TracerData.fromdir(os.path.join(TMP_DIR,'dumps','tracers'))
errtot = 0.0
for i in tracers.ids():
  single = tracers.get_trace(i)
  errlocal = (single['Xcart'][0,0] - single['Xcart'][-1,0])**2
  errtot += errlocal
errtot = np.sqrt(errtot)
print("Total error = %.2g %%" % (100*errtot))

if AUTO:
  data = {}
  data['SOL'] = [np.array([0.0]), np.array([errtot])]
  data['CODE'] = [np.array([0.0]), np.array([0.0])]
  data['THRESHOLD'] = 0.01
  print(data)
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  util.safe_remove(TMP_DIR)
  sys.exit()

print("Making sampling plot")
t0 = tracers.filter(tracers['time'] == 0.0)
plt.plot(t0['id'],t0['Xcart'][:,0],'bo',ms=2)
plt.xlabel('tracer id')
plt.ylabel(r'$x(t=0)$')
plt.savefig("tracers1d_sampled_position.pdf",
            bbox_inches='tight')
plt.savefig("tracers1d_sampled_position.png",
            bbox_inches='tight')
plt.clf()
plt.cla()

print("Making spacetime plot")
tid = len(tracers)//2
single = tracers.get_trace(tid)
plt.plot(single['Xcart'][:,0],single['time'],'bo',ms=2)
plt.xlabel(r'$x$')
plt.ylabel(r'$t$',rotation=0,labelpad=20)
plt.savefig('tracers1d_spacetime.pdf',bbox_inches='tight')
plt.savefig('tracers1d_spacetime.png',bbox_inches='tight')
plt.clf()
plt.cla()

# CLEAN UP
util.safe_remove(TMP_DIR)
