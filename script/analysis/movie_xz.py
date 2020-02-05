################################################################################
#                                                                              # 
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      # 
#                                                                              # 
################################################################################

import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import util
import glob
import os
import plot as bplt

FIGX = 13
FIGY = 10
SIZE = 40

if len(sys.argv) != 2:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit()

path = sys.argv[1]

#files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))
print 'Getting list of full dumps...'
#files = io.get_dumps_full(os.path.join(path, 'dumps/'))
files = io.get_dumps_reduced(os.path.join(path, 'dumps/'))

#FRAMEDIR = 'FRAMES'
FRAMEDIR = os.path.join(path, 'frames_xz/')
print FRAMEDIR
util.make_dir(FRAMEDIR)

hdr = io.load_hdr(files[0])
geom = io.load_geom(hdr)

print len(files)

def plot(args):
  n = args
  imname = 'frame_%08d.png' % n
  imname = os.path.join(FRAMEDIR, imname)
  print '%08d / ' % (n+1) + '%08d' % len(files) 
  print imname

  # Ignore if frame already exists
  if os.path.isfile(imname):
    return

  dump = io.load_dump(files[n], geom)
  fig = plt.figure(figsize=(FIGX, FIGY))
 
  fig.suptitle('t = %05.2g' % dump['t'])

  ax = plt.subplot(2,2,1)
  bplt.plot_xz(ax, geom, np.log10(dump['RHO']), dump, 
    vmin=-4, vmax = 0, label='RHO')
  bplt.overlay_field(ax, geom, dump, NLEV=10)
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,2,2)
  bplt.plot_xz(ax, geom, np.log10(dump['beta']), dump,
    vmin=-2, vmax=2, label='beta', cmap='RdBu_r')
  bplt.overlay_field(ax, geom, dump, NLEV=10)
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
 
  ax = plt.subplot(2,2,3)
  bplt.plot_xy(ax, geom, np.log10(dump['RHO']), dump,
     vmin=-4, vmax=0, label='RHO')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,2,4)
  bplt.plot_xy(ax, geom, np.log10(dump['beta']), dump,
     vmin=-2, vmax=2, label='beta', cmap='RdBu_r')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  #ax.pcolormesh(dump['X1'][:,:,0], dump['X2'][:,:,0], dump['RHO'][:,:,0])
  plt.savefig(imname, bbox_inches='tight', dpi=100)
  plt.close(fig)

import multiprocessing
import signal
import psutil

nthreads = psutil.cpu_count(logical=False)
nthreads = 4
print 'Number of CPUs: %i' % psutil.cpu_count(logical=False)

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  res = pool.map_async(plot, range(len(files)))
  res.get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
else:
  pool.close()
pool.join()

