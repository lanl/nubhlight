################################################################################
#                                                                              # 
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      # 
#                                                                              # 
################################################################################

import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import util
import os

if len(sys.argv) != 2:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit()

path = sys.argv[1]

files = io.get_dumps_full(path)

OUTDIR = os.path.join(path, '../analysis')
util.make_dir(OUTDIR)

hdr = io.load_hdr(files[0])
geom = io.load_geom(hdr)

NR = len(geom['r'][:,0,0])
NT = len(files)

t  = np.zeros(NT)
r  = geom['r'][:,hdr['N2']/2,0]
Phi = np.zeros([NT, NR])
mdot = np.zeros([NT, NR])
Sigma = np.zeros([NT, NR])
Thetae_sadw = np.zeros([NT, NR])
Thetap_sadw = np.zeros([NT, NR])
b2_sadw = np.zeros([NT, NR])
b2 = np.zeros([NT, NR])
betaz = np.zeros([NT, NR])

N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']
dx1 = hdr['dx1']; dx2 = hdr['dx2']; dx3 = hdr['dx3']
gdet = geom['gdet'][:,:,np.newaxis]

def get_ravg(n):
  print 'n = %i' % n

  dump = io.load_dump(files[n], geom)
  t[n] = dump['t']

  rho_av = (gdet*dx2*dx3*dump['RHO']).sum(axis=-1).sum(axis=-1)
  
  mdot[n,:] = -(gdet*dx2*dx3*dump['ucon'][:,:,:,1]*dump['RHO']).sum(axis=-1).sum(axis=-1)

  Phi[n,:] = 1./r*(gdet[:,N2/2,:]*dx1*dx3*np.fabs(dump['B2'][:,N2/2,:])).sum(axis=-1) 

  Sigma[n,:] = 1./(2.*np.pi*r[:])*(gdet*dx2*dx3*dump['RHO']).sum(axis=-1).sum(axis=-1)

  b2[n,:] = dump['bsq'][:,N2/2,:].mean(axis=-1)

  Thetae_sadw[n,:] = (gdet*dx2*dx3*dump['Thetae']*dump['RHO']).sum(axis=-1).sum(axis=-1)
  Thetap_sadw[n,:] = (gdet*dx2*dx3*dump['Thetap']*dump['RHO']).sum(axis=-1).sum(axis=-1)
  b2_sadw[n,:] = (gdet*dx2*dx3*dump['bsq']*dump['RHO']).sum(axis=-1).sum(axis=-1)

  ucon, ucov, bcon, bcov = io.get_state(dump, geom)
  bcon[:,:,:,0] = dump['B2'][:,:,:]*ucov[:,:,:,2]
  bcon[:,:,:,1] = bcon[:,:,:,0]*ucon[:,:,:,1]/ucon[:,:,:,0]
  bcon[:,:,:,2] = dump['B2'][:,:,:] + bcon[:,:,:,0]*ucon[:,:,:,2]/ucon[:,:,:,0]
  bcon[:,:,:,3] = bcon[:,:,:,0]*ucon[:,:,:,3]/ucon[:,:,:,0]
  for mu in xrange(4):
     bcov[:,:,:,mu] = (bcon[:,:,:,:]*geom['gcov'][:,:,None,mu,:]).sum(axis=-1)
  b2z = ((bcov*bcon).mean(axis=-1)).mean(axis=-1)
  P = ((hdr['gam']-1.)*dump['UU']).mean(axis=-1)
  betaz[n,:] = 2.*P[:,N2/2]/b2z[:,N2/2]

  Thetae_sadw[n,:] /= rho_av
  Thetap_sadw[n,:] /= rho_av
  b2_sadw[n,:] /= rho_av

for n in xrange(len(files)):
  get_ravg(n)

#import multiprocessing
#import signal
#import psutil

#nthreads = psutil.cpu_count(logical=False)
#print 'Number of CPUs: %i' % psutil.cpu_count(logical=False)

#original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
#pool = multiprocessing.Pool(nthreads)
#signal.signal(signal.SIGINT, original_sigint_handler)
#try:
  #res = pool.map_async(get_ravg, range(len(files)))
#  res = pool.map_async(get_ravg, range(10))
#  res.get(720000)
#except KeyboardInterrupt:
#  print 'Caught interrupt!'
#  pool.terminate()
#else:
#  pool.close()
#pool.join()


print t

out = {}
out['t'] = t
out['r'] = r
out['mdot'] = mdot.transpose()
out['Phi'] = Phi.transpose()
out['Sigma'] = Sigma.transpose()
out['b2'] = b2.transpose()
out['Thetae_sadw'] = Thetae_sadw.transpose()
out['Thetap_sadw'] = Thetap_sadw.transpose()
out['b2_sadw'] = b2_sadw.transpose()
out['betaz'] = betaz.transpose()
import pickle
pickle.dump(out, open(os.path.join(OUTDIR, 'ravgs.p'), 'wb'))

