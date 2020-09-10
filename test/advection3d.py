################################################################################
#                                                                              #
# Advection in 3D of a Passive Scalar                                          #
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
import multiprocessing
from matplotlib import animation, rc
from matplotlib import pyplot as plt
rc('font',size=18)

TMP_DIR = 'TMP'
TMP_BUILD = 'build_tmp.py'
util.safe_remove(TMP_DIR)

PROBLEM = 'advection3d'
AUTO = '-auto' in sys.argv
BIGRUN = "-bigrun" in sys.argv
MPI = "-mpi" in sys.argv
MOVIE = "-movie" in sys.argv

if BIGRUN:
    RES = [64,128,256]
else:
    RES = [16,32,64]

kwave = 2*np.pi
amp = 1.0
nspace = 3.0
ddiag = np.sqrt(nspace)
cadv = 0.5*ddiag
csqr = cadv**2
gamma = np.sqrt(1./(1. - csqr))
qsqr = gamma**2 - 1.0
u1d = np.sqrt(qsqr/nspace)
c1d = cadv / ddiag

def phi_true(t,x,y,z):
    phi_x = np.cos(kwave*(x - c1d*t))
    phi_y = np.cos(kwave*(y - c1d*t))
    phi_z = np.cos(kwave*(z - c1d*t))
    return amp*phi_x*phi_y*phi_z

util.make_dir(TMP_DIR)
os.chdir('../prob/' + PROBLEM)
copyfile('build.py', TMP_BUILD)

# Since this test is designed to run on a single machine (no batch scripts)
# set openmpi to only use a few threads. Let MPI handle the rest.
if MPI:
    num_mpi = 8
    num_cpus = multiprocessing.cpu_count()
    os.environ['OMP_NUM_THREADS'] = str(int(np.max([2,num_cpus/num_mpi])))

# COMPILE CODE AT MULTIPLE RESOLUTIONS USING SEPARATE BUILD FILE
for n,res in enumerate(RES):
    for d in [1,2,3]:
        util.change_cparm('N{}TOT'.format(d), res, TMP_BUILD)
    if MPI:
        for d in [1,2,3]:
            util.change_cparm('N{}CPU'.format(d), 2, TMP_BUILD)
    call([sys.executable, TMP_BUILD, '-dir', TMP_DIR])
    parm_src = os.path.join(os.getcwd(), TMP_DIR, 'param_template.dat')
    parm_dest = '../../test/' +  TMP_DIR + '/param.dat'
    call(['cp', os.path.join(os.getcwd(), TMP_DIR, 'bhlight'),
          '../../test/' + TMP_DIR + '/bhlight_' + str(res)])
    copyfile(parm_src,parm_dest)
    if BIGRUN:
        util.change_rparm('tf',  0.1,  parm_dest)
        util.change_rparm('DTd', 0.05, parm_dest)
        util.change_rparm('DTl', 0.05, parm_dest)
util.safe_remove(TMP_BUILD)
util.safe_remove(TMP_DIR)
os.chdir('../../test/')
os.chdir(TMP_DIR)

# and convergence plot
print("Convergence test...")
errs = [None for res in RES]
for n,res in enumerate(RES):
    print("Res = {}".format(res))
    call_string = ['./bhlight_' + str(res), '-p', 'param.dat']
    if MPI:
        bcall(call_string,int(num_mpi))
    else:
        bcall(call_string)
    dfiles = sorted(glob.glob('dumps/dump*.h5'))
    hdr = io.load_hdr(dfiles[-1])
    geom = io.load_geom(hdr, recalc=True)
    dump = io.load_dump(dfiles[-1],geom)
    N1,N2,N3 = [hdr['N{}'.format(d)] for d in range(1,4)]
    mshape = (N1,N2,N3)
    t = dump['t']
    x,y,z = [geom[d].reshape(mshape) for d in ['x','y','z']]
    phi = dump['var0'].reshape(mshape)
    error = phi - phi_true(t,x,y,z)
    max_error = np.max(np.abs(error))
    errs[n] = max_error

print("Richardson extrapolating...")
errf = lambda h, alpha, p: alpha*(h**p)
p0 = 2.0
h0 = 1.0/RES[0]
err0 = errs[0]
alpha0 = 10.*err0*h0*h0
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
    import pickle
    pickle.dump(data, open('data.p', 'wb'))
    util.safe_remove(TMP_DIR)
    sys.exit()

print("Plotting pointwise Error")
errs = 100.0*np.array(errs)

fig, ax = plt.subplots()
pcm = ax.pcolormesh(x[:,:,0],y[:,:,0],
                    1000*(phi - phi_true(t,x,y,z))[:,:,0])
ax.set_xlabel('x')
ax.set_ylabel('y')
cbar = fig.colorbar(pcm)
plt.tight_layout()
cbar.set_label(r'$1000\times \phi$' + ' error')
plt.savefig("../{}-pointwise.png".format(PROBLEM),bbox_inches='tight')
plt.savefig("../{}-pointwise.pdf".format(PROBLEM),bbox_inches='tight')
plt.clf()

print("Plotting convergence...")
plt.loglog(RES,100.*errf(hs,alpha,p),lw=2,ls='--',
           label=(r'$\sim N^{- %.2f}$' % (p)))
plt.loglog(RES,errs,'ro',ms=12,label=r'$\nu {\tt bhlight}$')
plt.xlabel(r'Resolution $N$')
plt.ylabel(r'$100\times \left| \phi^h - \phi_{true} \right|_\infty$')
plt.legend()
plt.tight_layout()
ax = plt.gca()
ax.set_xticks(RES)
ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_yticks(errs)
ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
plt.savefig('../{}.png'.format(PROBLEM),bbox_inches='tight')
plt.savefig('../{}.pdf'.format(PROBLEM),bbox_inches='tight')
plt.clf()

if BIGRUN:
    print("Can't make a movie for big run. Quitting.")
    sys.exit()

if MOVIE:
    print("Making movie/2d slice plot")
    dfiles = sorted(glob.glob('dumps/dump*.h5'))
    hdr = io.load_hdr(dfiles[0])
    geom = io.load_geom(hdr)
    N1,N2,N3 = [hdr['N{}'.format(d)] for d in range(1,4)]
    mshape = (N1,N2,N3)
    x,y,z = [geom[d].reshape(mshape) for d in ['x','y','z']]
    def get_phi(i):
        dump = io.load_dump(dfiles[i],geom)
        phi = dump['var0'].reshape(mshape)
        return phi
    phi0 = get_phi(0)
    fig, ax = plt.subplots()
    ax.set_xlim((0,1))
    ax.set_ylim((0,1))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    pc = ax.pcolormesh(x[:,:,0],y[:,:,0],phi0[:,:,0],
                       cmap='viridis',shading='gouraud')
    cbar = fig.colorbar(pc)
    cbar.set_clim(-1.1,1.1)
    cbar.set_label(r'$\phi$')
    plt.tight_layout()
    plt.savefig('../{}-lastframe.png'.format(PROBLEM))
    plt.savefig('../{}-lastframe.pdf'.format(PROBLEM))

    def init():
        pc.set_array(phi0[:,:,0].ravel())
        return pc
    def animate(i):
        phi = get_phi(i)
        pc.set_array(phi[:,:,0].ravel())
        return pc

    anim = mpl.animation.FuncAnimation(fig,animate,
                                       init_func=init,
                                       frames=101,
                                       interval=20,blit=False)

    anim.save("../{}.mp4".format(PROBLEM),
              writer='ffmpeg',
              extra_args=['-loglevel','verbose'])

print("Done.")
util.safe_remove(TMP_DIR)
