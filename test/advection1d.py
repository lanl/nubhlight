################################################################################
#                                                                              #
# Advection of Passive Scalars                                                 #
#                                                                              #
################################################################################

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
from shutil import copyfile
import glob
import numpy as np
import scipy as sp
from scipy import optimize
import hdf5_to_dict as io
import util
from bhlight import bcall
import matplotlib as mpl; mpl.use('Agg')
from matplotlib import animation, rc
from matplotlib import pyplot as plt
rc('font',size=18)

TMP_DIR = 'TMP'
TMP_BUILD = 'build_tmp.py'
util.safe_remove(TMP_DIR)

PROBLEM = 'advection1d'
AUTO = '-auto' in sys.argv
MOVIE = '-movie' in sys.argv
if '-idim' in sys.argv:
  IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
  IDIM = 1
  
RES = [128,256,512]

kwave = 2*np.pi
amp = 1.0
kcenter = 0.5
cadv = 0.5

util.make_dir(TMP_DIR)
os.chdir('../prob/' + PROBLEM)
copyfile('build.py', TMP_BUILD)

# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE
for n,res in enumerate(RES):
  util.change_cparm('N{}TOT'.format(IDIM), res, TMP_BUILD)
  call([sys.executable, TMP_BUILD, '-dir', TMP_DIR, '-idim', str(IDIM)])
  call(['cp', os.path.join(os.getcwd(), TMP_DIR, 'bhlight'),
        '../../test/' + TMP_DIR + '/bhlight_' + str(res)])
copyfile(os.path.join(os.getcwd(), TMP_DIR, 'param_template.dat'), '../../test/' + 
         TMP_DIR + '/param.dat')
util.safe_remove(TMP_BUILD)
util.safe_remove(TMP_DIR)
os.chdir('../../test/')
os.chdir(TMP_DIR)

def phi1_true(t,x):
  return amp*np.sin(2*np.pi*(x-kcenter - cadv*t))

# and convergence plot
print("Convergence test...")
errs = [None for res in RES]
for n,res in enumerate(RES):
  bcall(['./bhlight_' + str(res), '-p', 'param.dat'])
  dfiles = sorted(glob.glob('dumps/dump*.h5'))
  hdr = io.load_hdr(dfiles[-1])
  geom = io.load_geom(hdr, recalc=True)
  dump = io.load_dump(dfiles[-1],geom)
  N = hdr['N{}'.format(IDIM)]
  t = dump['t']
  x = geom[['x','y','z'][IDIM-1]].reshape(N)
  phi1 = dump['var1'].reshape(N)
  error = phi1 - phi1_true(t,x)
  max_error = np.max(np.abs(error))
  errs[n] = max_error

print("Richardson extrapolating...")
errf = lambda h, alpha, p: alpha*(h**p)
p0 = 2.0
h0 = 1.0/RES[0]
err0 = errs[0]
alpha0 = err0*h0*h0
hs = 1.0/np.array(RES)
(alpha,p),pcov = optimize.curve_fit(errf,hs,errs,p0=(alpha0,p0))
print("Convergence data:\nalpha = {}, p = {}\npcov = {}\n".format(
    alpha,p,pcov))

if AUTO:
  os.chdir("../")
  data = {}
  data['SOL'] = [np.array([0.0]), np.array([2.0])]
  data['CODE'] = [np.array([0.0]), np.array([p])]
  data['THRESHOLD'] = 0.03
  print(data)
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  util.safe_remove(TMP_DIR)
  sys.exit()

print("Plotting convergence...")
plt.loglog(RES,errf(hs,alpha,p),lw=2,ls='--',
           label=(r'$%.2f h^{%.2f}$' % (alpha,p)))
plt.loglog(RES,errs,'ro',ms=12,label='measured')
plt.xlabel('Resolution')
plt.ylabel(r'$\left| \phi_1^h - \phi_1 \right|_\infty$')
plt.legend()
plt.savefig('../advection1d.png',bbox_inches='tight')
plt.savefig('../advection1d.pdf',bbox_inches='tight')
plt.clf()

# Make one frame
print("Making frame")
dfiles = sorted(glob.glob('dumps/dump*.h5'))
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr, recalc=True)
N = hdr['N{}'.format(IDIM)]
x = geom[['x','y','z'][IDIM-1]].reshape(N)
dump = io.load_dump(dfiles[-1],geom)
phi0 = dump['var0'].reshape(N)
phi1 = dump['var1'].reshape(N)
fig, (ax0,ax1) = plt.subplots(2,sharex=True)
ax0.set_xlim((0,1))
ax0.set_ylim((-0.05,1.05))
ax1.set_ylim((-1.05,1.05))
ax1.set_xlabel('x')
ax0.set_ylabel(r'$\phi_0$')
ax1.set_ylabel(r'$\phi_1$')
ax0.plot(x,phi0)
ax1.plot(x,phi1)
plt.savefig('../advection1d-lastframe.png',bbox_inches='tight')
plt.savefig('../advection1d-lastframe.pdf',bbox_inches='tight')
plt.clf()

# Make movie
if MOVIE:
  print("Making movie")
  def get_phis(i):
    dump = io.load_dump(dfiles[i],geom) 
    phi0 = dump['var0'].reshape(N)
    phi1 = dump['var1'].reshape(N)
    return phi0,phi1
  fig, (ax0,ax1) = plt.subplots(2,sharex=True)
  ax0.set_xlim((0,1))
  ax0.set_ylim((-0.05,1.05))
  ax1.set_ylim((-1.05,1.05))
  ax1.set_xlabel('x')
  ax0.set_ylabel(r'$\phi_0$')
  ax1.set_ylabel(r'$\phi_1$')
  line0, = ax0.plot([],[],lw=2)
  line1, = ax1.plot([],[],lw=2)

  def init():
    line0.set_data([],[])
    line1.set_data([],[])
    return line0,line1
  def animate(i):
    phi0,phi1 = get_phis(i)
    line0.set_data(x,phi0)
    line1.set_data(x,phi1)
    return line0,line1

  anim = mpl.animation.FuncAnimation(fig,animate,
                                     init_func=init,
                                     frames=101,
                                     interval=20,blit=True)
  anim.save('../advection1d.mp4',
            writer='ffmpeg',
            extra_args=['-loglevel','verbose'])

# CLEAN UP
util.safe_remove(TMP_DIR)
