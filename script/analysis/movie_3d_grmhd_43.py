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
FIGY = 10.5
SIZE = 20
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

N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']

def plot(args):
  n = args
  print('%08d / ' % (n+1) + '%08d' % len(files))
  diag = io.load_diag(path)
  dump = io.load_dump(files[n], geom)
  #hdr = dump['hdr']
  fig = plt.figure(figsize=(FIGX, FIGY), frameon=False);
  ax = fig.gca()
  #ax = plt.subplot(1,1,1)

  x = geom['x']
  y = geom['y']
  z = geom['z']
  x = bplt.flatten_xz(x, hdr, flip=True)
  y = bplt.flatten_xz(y, hdr, flip=True)
  z = bplt.flatten_xz(z, hdr)

  # rcylindrical
  rcyl = np.sqrt(x**2 + y**2)

  lrho = bplt.flatten_xz(np.log10(dump['RHO']), hdr)
  mesh = ax.pcolormesh(rcyl[0:N1,0:N2], z[0:N1,0:N2], lrho[0:N1,0:N2], cmap='jet', vmin=-4, vmax=0, shading='gouraud')
  cax = fig.add_axes([0.135, 0.15, 0.025, 0.7])
  cb = fig.colorbar(mesh, cax=cax, ticks=[-4, -2, 0])
  cb.set_label(label='Density', fontsize=20)
  cax.tick_params(axis='both', which='major', labelsize=20)

  lbeta = bplt.flatten_xz(np.log10(dump['beta']), hdr)
  mesh = ax.pcolormesh(rcyl[N1:,0:N2], z[N1:,0:N2], lbeta[N1:,0:N2], cmap='inferno', vmin=-2, vmax=2, shading='gouraud')
  cax = fig.add_axes([0.865, 0.15, 0.025, 0.7])
  cb = fig.colorbar(mesh, cax=cax, ticks=[-2, 0, 2])
  cb.set_label(label='Plasma beta', fontsize=20)
  cb.ax.yaxis.set_label_position('left')
  cax.yaxis.set_ticks_position('left')
  cax.tick_params(axis='both', which='major', labelsize=20)

  ax.set_xlim([-4./3.*SIZE, 4./3.*SIZE]); ax.set_ylim([-SIZE, SIZE])
  ax.set_xticks([])
  ax.set_yticks([])
  circle1=plt.Circle((0,0),1. + np.sqrt(1.-hdr['a']**2),color='k',zorder=2)
  ax.add_artist(circle1)

  #ax.axvline(hdr['Risco'], color='k', linestyle='--', linewidth=2)
  #ax.axvline(-hdr['Risco'], color='k', linestyle='--', linewidth=2)

  bplt.overlay_field(ax, geom, dump, linestyle='-', linewidth=1, 
    linecolor='k', NLEV=10)

  plt.savefig(os.path.join(FRAMEDIR, 'frame_%08d.png' % n), 
      bbox_inches='tight', pad_inches=0, dpi=100)
  plt.close(fig)

import multiprocessing
import signal
import psutil

NTHREADS = psutil.cpu_count(logical=False)
print(psutil.cpu_count(logical=False))

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

