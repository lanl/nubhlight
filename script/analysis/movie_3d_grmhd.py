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

FIGX = 14
FIGY = 10
SIZE = 40
#NTHREADS = 1

if len(sys.argv) != 2:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit()

path = sys.argv[1]

files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

hdr  = io.load_hdr(files[0]) 
geom = io.load_geom(hdr)

def plot(args):
  n = args
  print('%08d / ' % (n+1) + '%08d' % len(files))
  diag = io.load_diag(path)
  dump = io.load_dump(files[n], geom)
  #hdr = dump['hdr']
  fig = plt.figure(figsize=(FIGX, FIGY))

  # GET SHELL AVERAGES
  Thetae_sadw = np.zeros(hdr['N1'])
  sigma = np.zeros(hdr['N1'])
  mdot = np.zeros(hdr['N1'])

  for i in range(hdr['N1']):
    vol = 0.
    for j in range(hdr['N2']):
      for k in range(hdr['N3']):
        Thetae_sadw[i] += dump['Theta'][i,j,k]*dump['RHO'][i,j,k]*geom['gdet'][i,j]*hdr['dx1']*hdr['dx2']*hdr['dx3']
        sigma[i] += dump['RHO'][i,j,k]*hdr['dx2']*geom['gdet'][i,j]
        mdot[i] += dump['RHO'][i,j,k]*hdr['dx2']*hdr['dx3']*geom['gdet'][i,j]
        vol += dump['RHO'][i,j,k]*geom['gdet'][i,j]*hdr['dx1']*hdr['dx2']*hdr['dx3']
    sigma[i] /= (hdr['N2']*hdr['N3'])
    Thetae_sadw[i] /= vol
  
  ax = plt.subplot(3,3,1)
  bplt.plot_xz(ax, geom, np.log10(dump['RHO']), dump, 
    vmin=-4, vmax = 0, label='RHO', ticks=[-4,-3,-2,-1,0])
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  ax.set_xticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  ax.set_yticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  
  ax = plt.subplot(3,3,2)
  bplt.plot_xz(ax, geom, np.log10(dump['Theta']), dump,
    vmin=-2, vmax=2, label='Theta', cmap='RdBu_r', ticks=[-2,-1,0,1,2])
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  ax.set_xticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  ax.set_yticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])

  ax = plt.subplot(3,3,3)
  ax.plot(geom['r'][:,0,0], sigma, color='b', label='Sigma')
  ax.set_xscale('log'); ax.set_yscale('log')
  ax.set_xlim([1, SIZE]); ax.set_ylim([1.e-2, 1.e2])
  ax.set_aspect(np.log(SIZE/1)/np.log(1.e2/1.e-2))
  ax.plot(geom['r'][:,0,0], Thetae_sadw, color='r', label='SADW Thetae')
  ax.legend(loc=2, fontsize=10)
  ax.set_xlabel('r/M')

  ax = plt.subplot(3,3,4)
  bplt.plot_xy(ax, geom, np.log10(dump['RHO']), dump,
     vmin=-4, vmax=0, label='RHO', ticks=[-4,-3,-2,-1,0])
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  ax.set_xticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  ax.set_yticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  
  ax = plt.subplot(3,3,5)
  bplt.plot_xy(ax, geom, np.log10(dump['Theta']), dump,
     vmin=-2, vmax=2, label='Theta', cmap='RdBu_r',
     ticks=[-2,-1,0,1,2])
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  ax.set_xticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])
  ax.set_yticks([-SIZE,-SIZE/2,0,SIZE/2,SIZE])

  ax = plt.subplot(3,1,3)
  ax.plot(diag['t'], diag['mdot'], color='k')
  ax.axvline(dump['t'], color='k', linestyle='--')
  ax.set_xlim([0, dump['hdr']['tf']])
  ax.set_xlabel('t/M')
  ax.set_ylabel('mdot')
  ax2 = ax.twinx()
  ax2.set_ylabel('phi')

  plt.subplots_adjust(hspace=0.25)

  #ax.pcolormesh(dump['X1'][:,:,0], dump['X2'][:,:,0], dump['RHO'][:,:,0])
  plt.savefig(os.path.join(FRAMEDIR, 'frame_%08d.png' % n), 
      bbox_inches='tight', dpi=100)
  plt.close(fig)

import multiprocessing
import signal
import psutil

nthreads = psutil.cpu_count(logical=False)
print(psutil.cpu_count(logical=False))
NTHREADS = 10

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(NTHREADS)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  res = pool.map_async(plot, range(len(files)))
  res.get(720000)
except KeyboardInterrupt:
  pool.terminate()
else:
  pool.close()
pool.join()

