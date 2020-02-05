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
import psutil
#import plot as bplt

FIGX = 14
FIGY = 10
#SIZE = 40
NTHREADS = psutil.cpu_count(logical=False)
#NTHREADS = 1
THMIN = np.pi/3.
THMAX = 2.*np.pi/3.

if len(sys.argv) != 3:
  util.warn('Format: python eht_analysis.py [dump path] [tavg]')
  sys.exit()

path = sys.argv[1]
tavg = float(sys.argv[2])

dumps = io.get_dumps_reduced(path)

hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr)

# Calculate jmin, jmax
ths = geom['th'][-1,:,0]
for n in xrange(len(ths)):
  if ths[n] > THMIN:
    jmin = n
    break
for n in xrange(len(ths)):
  if ths[n] > THMAX:
    jmax = n
    break

diag = io.load_diag(path)

dx1 = hdr['dx1']
dx2 = hdr['dx2']
dx3 = hdr['dx3']
N1 = hdr['N1']
N2 = hdr['N2']
N3 = hdr['N3']

r = geom['r'][:,N2/2,0]

ND = len(dumps)

t = np.zeros(ND)
rho_r = np.zeros([ND, N1])
Theta_r = np.zeros([ND, N1])
B_r = np.zeros([ND, N1])
Pg_r = np.zeros([ND, N1])
Ptot_r = np.zeros([ND, N1])
betainv_r = np.zeros([ND, N1])
uphi_r = np.zeros([ND, N1])
rho_SADW = np.zeros([ND, N1])

Mdot = np.zeros(ND)
Phi = np.zeros(ND)
Ldot = np.zeros(ND)
Edot = np.zeros(ND)
Lum = np.zeros(ND)

# EVALUATE DIAGNOSTICS
  
vol = (dx2*2.*np.pi*geom['gdet'][:,:]).sum(axis=-1)

def INT(var):
  return (dx2*dx3*geom['gdet'][:,:,None]*var[:,:,:]).sum(axis=-1).sum(axis=-1)

def WAVG(var, w):
  return INT(w*var)/INT(w)

def avg_dump(args):
  global t
  n = args
  print '%08d / ' % (n+1) + '%08d' % len(dumps) 
  dump = io.load_dump(dumps[n], geom)
  t[n] = dump['t']
  print dump['t']

  rho_SADW[n,:] = WAVG(dump['RHO'], dump['RHO'])

  # SHELL AVERAGES
  #vol = (dx2*dx3*geom['gdet'][:,:]).sum(axis=-1)
  rho_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*dump['RHO'][:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  Theta = (hdr['gam']-1.)*dump['UU'][:,:,:]/dump['RHO'][:,:,:]
  Theta_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Theta[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  B = np.sqrt(dump['bsq'][:,:,:])
  B_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*B[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Pg = (hdr['gam']-1.)*dump['UU'][:,:,:]
  Pg_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Pg[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Ptot = Pg + dump['bsq'][:,:,:]/2
  Ptot_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Ptot[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  betainv = (dump['bsq'][:,:,:]/2)/Pg
  betainv_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*betainv[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  uphi = (dump['ucon'][:,:,:,3])
  uphi_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*uphi[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1) 

  rho_r[n,:] /= vol
  Theta_r[n,:] /= vol
  B_r[n,:] /= vol
  Pg_r[n,:] /= vol
  Ptot_r[n,:] /= vol
  betainv_r[n,:] /= vol
  uphi_r[n,:] /= vol

  # FLUXES
  iF = 5
  Mdot[n] = (dump['RHO'][iF,:,:]*dump['ucon'][iF,:,:,1]*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  #Phi[n] = 0.5*(np.fabs(dump['bcon'][iF,:,:,1])*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  Phi[n] = 0.5*(np.fabs(dump['B1'][iF,:,:])*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  Trphi = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,3]
  Trphi -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,3]
  Trt = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,0]
  Trt -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,0]
  Ldot[n] = (Trphi[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  Edot[n] = -(Trt[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()

  rho = dump['RHO']
  P = (hdr['gam']-1.)*dump['UU']
  B = np.sqrt(dump['bsq'])
  C = 0.2
  j = rho**3*P**(-2)*np.exp(-C*(rho**2/(B*P**2))**(1./3.))
  Lum[n] = (j*geom['gdet'][:,:,None]*dx1*dx2*dx3).sum()

for n in xrange(len(dumps)):
  avg_dump(n)

# PARALLEL EXECUTION
#import multiprocessing
#import signal
#import psutil
#original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#pool = multiprocessing.Pool(NTHREADS)
#signal.signal(signal.SIGINT, original_sigint_handler)
#try:
#  res = pool.map_async(avg_dump, range(len(dumps)))
#  res.get(720000)
#except KeyboardInterrupt:
#  pool.terminate()
#else:
#  pool.close()
#pool.join()

# Time average
n = 0
for n in xrange(ND):
  if t[n] >= tavg:
    break

print 'nmin = %i' % n
avg = {}
rho = rho_r[n:,:].mean(axis=0) 
Theta = Theta_r[n:,:].mean(axis=0)
B = B_r[n:,:].mean(axis=0)
Pg = Pg_r[n:,:].mean(axis=0)
Ptot = Ptot_r[n:,:].mean(axis=0)
betainv = betainv_r[n:,:].mean(axis=0)
uphi = uphi_r[n:,:].mean(axis=0)

# OUTPUT
import pickle

out = {}
out['rho_SADW'] = rho_SADW[n:,:].mean(axis=0)
out['r'] = r
out['rho_r'] = rho
out['Theta_r'] = Theta
out['B_r'] = B
out['Pg_r'] = Pg
out['Ptot_r'] = Ptot
out['betainv_r'] = betainv
out['uphi_r'] = uphi
out['t'] = t
out['Mdot'] = Mdot
out['Phi'] = Phi
out['Ldot'] = Ldot
out['Edot'] = Edot
out['Lum'] = Lum

out['t_d'] = diag['t']
out['Mdot_d'] = diag['mdot']
out['Phi_d'] = diag['Phi']
out['Ldot_d'] = diag['ldot']
out['Edot_d'] = -diag['edot']
out['Lum_d'] = diag['lum_eht']
pickle.dump(out, open('eht_out.p', 'w'))




