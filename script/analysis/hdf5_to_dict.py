# ======================================================================
# copyright 2020. Triad National Security, LLC. All rights
# reserved. This program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
# is operated by Triad National Security, LLC for the U.S. Department
# of Energy/National Nuclear Security Administration. All rights in
# the program are reserved by Triad National Security, LLC, and the
# U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others
# acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
# license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display
# publicly, and to permit others to do so.
# ======================================================================

# Authors: Jonah Miller (jonahm@lanl.gov) and Ben Ryan (brryan@lanl.gov)
# PURPOSE:
# Provides interface to read tracer and dump files fron bhlight/nubhlight
# Mostly reads hdf5 files and converts to structured Python dictionary
# appropriate for use with associated plotting scripts.
# Performs some simple calculations like reading in coordinate
# system and using coordinate system to construct derived
# quantities like covariant and contravariant velocities

import sys; sys.dont_write_bytecode = True
import units
import os
import glob
units = units.get_cgs()
SMALL = 1.e-40
RADTYPE_LIGHT = 1
RADTYPE_NEUTRINOS = 2

def h5_to_string(char_array):
  """Python3 distinguishes between char strings and string objects.
  This causes incompatibilities with HDF5.
  """
  import numpy as np
  if type(char_array) in [bytes, np.bytes_]:
    return char_array.decode()
  if type(char_array) == str:
    return char_array
  raise TypeError("Char_array must be a string or byte array!\n"
                  +"Your type is: {}.\n".format(type(char_array)))

def get_dumps_reduced(folder,twod=False):
  if twod:
    return sorted(glob.glob(os.path.join(folder,'dump2d_*.h5')))
  return sorted(glob.glob(os.path.join(folder,'dump_*.h5')))

def get_dumps_full(folder,twod=False):
  import h5py
  alldumps = get_dumps_reduced(folder,twod)
  fulldumps = []

  for fname in alldumps:
    try:
      with h5py.File(fname, 'r') as dfile:
        if 'FULL_DUMP' not in dfile:
          fulldumps.append(fname)
        elif dfile['FULL_DUMP'][0]:
          fulldumps.append(fname)
    except:
      pass

  return sorted(fulldumps)

def get_tracer_fnams(tracer_dir):
  files1 = list(glob.glob(os.path.join(tracer_dir,'tracers_*.h5')))
  files2 = list(glob.glob(os.path.join(tracer_dir,'tracers_*.h5part')))
  files = sorted(files1 + files2)
  return files

def load_hdr(fname):
  import numpy as np
  import h5py

  dfile = h5py.File(fname, 'r')
  path = os.path.dirname(os.path.realpath(fname))  # more robust

  hdr = {}

  # closure
  def read_hdr_var(hdrname,varname,default = 0):
    do_string_conversion = type(default) is str
    if varname in dfile.keys():
      hdr[hdrname] = dfile[varname][0]
    else:
      hdr[hdrname] = default
    if do_string_conversion:
      hdr[hdrname] = h5_to_string(hdr[hdrname])

  hdr['PATH'] = path
  hdr['N1'] = dfile['N1tot'][0]
  hdr['N2'] = dfile['N2tot'][0]
  hdr['N3'] = dfile['N3tot'][0]
  hdr['METRIC'] = h5_to_string(dfile['metric'][0])
  hdr['ELECTRONS'] = dfile['electrons'][0]
  hdr['RADIATION'] = dfile['radiation'][0]
  hdr['OUTPUT_EOSVARS'] = dfile['output_eosvars'][0]
  hdr['NVAR'] = dfile['nvar'][0]
  hdr['startx'] = np.array([0, dfile['startx[1]'][0], dfile['startx[2]'][0],
    dfile['startx[3]'][0]])
  hdr['dx'] = np.array([0, dfile['dx[1]'][0], dfile['dx[2]'][0],
    dfile['dx[3]'][0]])
  hdr['stopx'] = hdr['startx'] + hdr['dx']*np.array([0, hdr['N1'],hdr['N2'],hdr['N3']])
  
  read_hdr_var('VERSION','version',None)
  read_hdr_var('NVAR_PASSIVE','nvar_passive',0)
  read_hdr_var('EOS','eos','GAMMA')
  read_hdr_var('TRACERS','tracers',0)
  read_hdr_var('nulnutype','nulnutype','camera')
  read_hdr_var('nth','nth',8)
  read_hdr_var('nphi','nphi',8)
  read_hdr_var('nubins_spec','nubins_spec',200)
  read_hdr_var('diagnostics_use_radtypes','diagnostics_use_radtypes',0)
  read_hdr_var('FULL_DUMP','FULL_DUMP',1)
  read_hdr_var("LOCAL_ANGULAR_DISTRIBUTIONS", "local_angular_distributions", 0)
  read_hdr_var("NEUTRINO_OSCILLATIONS", "neutrino_oscillations", 0)
  read_hdr_var("FORCE_EQUIPARTITION", "force_equipartition", 0)

  hdr['vnams'] = [h5_to_string(i) for i in dfile['P'].attrs['vnams']]

  keys = ['cour', 'DTd', 'DTl', 'DTp', 'DTr', 'tf']
  if hdr['ELECTRONS']:
    keys += ['game', 'gamp']
  hdr['NEED_UNITS'] = (hdr['EOS'] == 'TABLE' or hdr['RADIATION'])
  if hdr['NEED_UNITS']:
    keys += ['L_unit', 'T_unit', 'M_unit', 'RHO_unit', 'U_unit', 'B_unit']
  if hdr['RADIATION']:
    keys += ['tp_over_te', 'Ne_unit', 'Thetae_unit',
             'maxnscatt', 'nubins', 'numin', 'numax']
  if hdr['EOS'] == 'TABLE':
    keys += ['TEMP_unit']
  if hdr['METRIC'] == 'MKS':
    hdr['DEREFINE_POLES'] = dfile['derefine_poles'][0]
    keys += ['Rin', 'Rout','Reh', 'Risco', 'hslope', 'a', 'poly_xt',
             'poly_alpha', 'mks_smooth']
    if 'Rout_vis' in dfile.keys():
      keys += ['Rout_vis']
    if hdr['RADIATION']:
      keys += ['Mbh', 'mbh', 'Rout_rad']
  if hdr['EOS'] == 'GAMMA':
    keys += ['gam']
  if hdr['EOS'] == 'TABLE':
    keys += ['eospath']
  if hdr['EOS'] == 'POLYTROPE':
    keys += ['poly_K','poly_gam']

  for key in keys:
    hdr[key] = dfile[key][0]

  if hdr['METRIC'] == 'MKS' and hdr['RADIATION']:
    hdr['LEdd'] = 4.*np.pi*units['GNEWT']*hdr['Mbh']*units['MP']*units['CL']/units['THOMSON']
    hdr['nomEff'] = 0.1
    hdr['MdotEdd'] = hdr['LEdd']/(hdr['nomEff']*units['CL']**2)

  if hdr['RADIATION']:
    nu = np.zeros(hdr['nubins'])
    lnumin = np.log(hdr['numin'])
    lnumax = np.log(hdr['numax'])
    dlnu = (lnumax - lnumin)/hdr['nubins']
    for n in range(len(nu)):
      nu[n] = np.exp(lnumin + (0.5 + n)*dlnu)
    hdr['lnumin'] = lnumin
    hdr['lnumax'] = lnumax
    hdr['dlnu'] = dlnu
    hdr['nu'] = nu

  dfile.close()

  return hdr

def load_geom(hdr,
              recalc=False,
              use_3d_metrics=False):
  # IF RAD, CALCULATE DOMEGA!

  import numpy as np
  import h5py
  import os
  use_2d_metrics = not use_3d_metrics

  dfile = h5py.File(os.path.join(hdr['PATH'], 'grid.h5'),'r')

  geom = {}
  geom['gcov'] = np.array(dfile['gcov'])
  geom['gcon'] = np.array(dfile['gcon'])
  geom['gdet'] = np.array(dfile['gdet'])
  geom['alpha'] = np.array(dfile['alpha'])
  force_2d = hdr['N3'] == 1 and geom['gdet'].shape[2] > 1

  if use_2d_metrics or force_2d:
    geom['gcov'] = geom['gcov'][:,:,0]
    geom['gcon'] = geom['gcon'][:,:,0]
    geom['gdet'] = geom['gdet'][:,:,0]

  geom['X1'] = np.array(dfile['Xharm'][:,:,:,1])
  geom['X2'] = np.array(dfile['Xharm'][:,:,:,2])
  geom['X3'] = np.array(dfile['Xharm'][:,:,:,3])

  geom['x'] = np.array(dfile['Xcart'][:,:,:,1])
  geom['y'] = np.array(dfile['Xcart'][:,:,:,2])
  geom['z'] = np.array(dfile['Xcart'][:,:,:,3])

  geom['X1f'] = np.array(dfile['XFharm'][:,:,:,1])
  geom['X2f'] = np.array(dfile['XFharm'][:,:,:,2])
  geom['X3f'] = np.array(dfile['XFharm'][:,:,:,3])

  geom['xf'] = np.array(dfile['XFcart'][:,:,:,1])
  geom['yf'] = np.array(dfile['XFcart'][:,:,:,2])
  geom['zf'] = np.array(dfile['XFcart'][:,:,:,3])

  geom['Lambda_h2cart_con'] = np.array(dfile['Lambda_h2cart_con'])
  geom['Lambda_h2cart_cov'] = np.array(dfile['Lambda_h2cart_cov'])

  if force_2d:
    keys = ['alpha',
            'X1','X2','X3',
            'x','y','z',
            'X1f','X2f','X3f',
            'xf','yf','zf',
            'Lambda_h2cart_con',
            'Lambda_h2cart_cov']
    for k in keys:
      geom[k] = geom[k][:,:,0,np.newaxis]

  if hdr['METRIC'] == 'MKS':
    geom['r'] = np.array(dfile['Xbl'][:,:,:,1])
    geom['th'] = np.array(dfile['Xbl'][:,:,:,2])
    geom['phi'] = np.array(dfile['Xbl'][:,:,:,3])

    if hdr['N3'] == 1:
      geom['phi'][:,:,:] = 0. # TODO: is this actually right?

    geom['rcyl'] = geom['r']*np.sin(geom['th'])
    geom['rcyl'][:,0,:] = 0.
    geom['rcyl'][:,-1,:] = 0.

    if recalc:
      from jacobians import Jacobians
      jac = Jacobians.fromhdr(hdr)
      geom['Lambda_h2bl_con'] = np.empty((hdr['N1'],hdr['N2'],hdr['N3'],4,4))
      geom['Lambda_h2bl_cov'] = np.zeros_like(geom['Lambda_h2bl_con'])
      geom['Lambda_bl2cart_con'] = np.zeros_like(geom['Lambda_h2bl_con'])
      geom['Lambda_bl2cart_cov'] = np.zeros_like(geom['Lambda_h2bl_con'])
      for i in range(hdr['N1']):
        for j in range(hdr['N2']):
          for k in range(hdr['N3']):
            X1 = geom['X1'][i,j,k]
            X2 = geom['X2'][i,j,k]
            X3 = geom['X3'][i,j,k]
            geom['Lambda_h2bl_cov'][i,j,k],geom['Lambda_h2bl_con'][i,j,k] = jac.get_h2bl(X1,X2,X3)
            geom['Lambda_bl2cart_cov'][i,j,k],geom['Lambda_bl2cart_con'][i,j,k] = jac.get_bl2cart(X1,X2,X3)
    else:
      geom['Lambda_h2bl_con'] = np.array(dfile['Lambda_h2bl_con'])
      geom['Lambda_h2bl_cov'] = np.array(dfile['Lambda_h2bl_cov'])
      geom['Lambda_bl2cart_con'] = np.array(dfile['Lambda_bl2cart_con'])
      geom['Lambda_bl2cart_cov'] = np.array(dfile['Lambda_bl2cart_cov'])
    geom['Lambda_h2bl_con_3d'] = geom['Lambda_h2bl_con'].copy()
    geom['Lambda_h2bl_cov_3d'] = geom['Lambda_h2bl_cov'].copy()
    if force_2d:
      keys = ['r','th','phi','rcyl',
              'Lambda_h2bl_con',
              'Lambda_h2bl_cov',
              'Lambda_bl2cart_con',
              'Lambda_bl2cart_cov',
              'Lambda_h2bl_con_3d',
              'Lambda_h2bl_cov_3d']
      for k in keys:
        geom[k] = geom[k][:,:,0,np.newaxis]
    if use_2d_metrics or force_2d:
      geom['Lambda_h2bl_con'] = geom['Lambda_h2bl_con'][:,:,0]
      geom['Lambda_h2bl_cov'] = geom['Lambda_h2bl_cov'][:,:,0]
      
    if 'Gamma' in dfile.keys():
      geom['Gamma'] = np.array(dfile['Gamma'])
      if use_2d_metrics or force_2d:
        geom['Gamma'] = geom['Gamma'][:,:,0]

  # local angles stuff
  for key in ['local_angles_Xharm', 'local_angles_Xbl',
              'local_angles_Xcart', 'local_angles_mu']:
    if key in dfile.keys():
      geom[key] = np.array(dfile[key])

  dfile.close()

  return geom

def load_dump(fname, geom=None, nulegacy=False):
  if geom == None:
    hdr = load_hdr(fname)
    geom = load_geom(hdr)

  import h5py
  import numpy as np
  hdr = load_hdr(fname)

  dfile = h5py.File(fname, 'r')
  dump = {}
  dump['hdr'] = hdr
  dump['geom'] = geom
  dump['t'] = dfile['t'][0]
  if 'dump_cnt' in dfile.keys():
    dump['dump_cnt'] = dfile['dump_cnt'][0]
  else:
    dump['dump_cnt'] = None

  for n in range(hdr['NVAR']):
    #vnam = str(hdr['vnams'][n], STRTYPE)
    dump[hdr['vnams'][n]] = np.array(dfile['P'][:,:,:,n])

  dump['jcon'] = dfile['jcon']

  keys = []
  if hdr['FULL_DUMP']:
    keys += ['divb', 'fail_save']
    if hdr['ELECTRONS']:
      keys += ['Qvisc']
    if hdr['RADIATION']:
      keys += ['Rmunu', 'Nsph', 'nph', 'nuLnu', 'Jrad', 'Nem', 'Nabs', 'Nsc']
      if 'radG_int' in dfile.keys():
        keys += ['radG_int']
      if 'tau_cool' in dfile.keys():
        keys += ['tau_cool']
      if 'dtau_avg' in dfile.keys():
        keys += ['dtau_avg']
      if 'dtau_scatt' in dfile.keys():
        keys += ['dtau_scatt']
      if 'dtau_tot' in dfile.keys():
        keys += ['dtau_tot']
      if 'Nem_phys' in dfile.keys():
        keys += ['Nem_phys']
      if 'Nabs_phys' in dfile.keys():
        keys += ['Nabs_phys']
      if hdr['ELECTRONS'] and hdr['RADIATION']:
        keys += ['Qcoul']
    if hdr['OUTPUT_EOSVARS']:
      keys += ['PRESS', 'TEMP', 'ENT']
      if 'SC2' in dfile.keys():
        keys += ['SC2']
  for key in keys:
    dump[key] = np.array(dfile[key]) + SMALL

  # divb calculator is anomalous at polar axes
  if hdr['METRIC'] == 'MKS':
    dump['divb'][:,0,:] = SMALL
    dump['divb'][:,-1,:] = SMALL

  # microphysics
  if hdr['OUTPUT_EOSVARS']:
    dump['Theta'] = np.array(dfile['TEMP']) + SMALL
  else:
    dump['PRESS'] = (hdr['gam']-1.)*dump['UU']
    dump['Theta'] = dump['PRESS']/dump['RHO']
    dump['ENT'] = (hdr['gam']-1.)*dump['UU']*(dump['RHO']**(-1.0*hdr['gam']))

  dump['TEMP'] = dump['Theta']
  if hdr['EOS'] == 'TABLE': # dump[TEMP] is in MeV
    dump['THETA'] = (dump['TEMP']*hdr['TEMP_unit']
                     / (units['ME']*units['CL']*units['CL']))

  if hdr['ELECTRONS']:
    dump['Thetae'] = units['MP']/units['ME']*dump['KEL']*dump['RHO']**(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['up'] = dump['UU'] - dump['ue']
    dump['Thetap'] = (hdr['gamp']-1.)*dump['up']/dump['RHO']
    dump['TpTe'] = (hdr['gamp']-1.)*dump['up']/((hdr['game']-1.)*dump['ue'])
  elif hdr['RADIATION']:
    dump['Thetae'] = dump['UU']/dump['RHO']*hdr['Thetae_unit']

  if hdr['RADIATION']:
    if 'tau_cool' not in dump.keys():
      dump['tau_cool'] = np.abs((np.abs(dump['UU'])+SMALL)/(dump['Jrad'][0] + SMALL))
      dump['tau_cool'][dump['Jrad'][0] <= SMALL] = np.inf
    if 'tau_heat' not in dump.keys():
      dump['tau_heat'] = np.abs((np.abs(dump['UU'])+SMALL)/(dump['Jrad'][1] + SMALL))
    dump['Qrad'] = ((dump['Nem'] + dump['Nsc'])/hdr['DTd'])*dump['tau_cool']
    if 'Nem_phys' in dump.keys():
      if hdr['RADIATION'] == RADTYPE_NEUTRINOS and nulegacy == False:
        dump['Nem_e'] = dump['Nem_phys'][:,:,:,0]
        dump['Nem_anti'] = dump['Nem_phys'][:,:,:,1]
        dump['Nem_x'] = dump['Nem_phys'][:,:,:,2]
      else:
        dump['Nem_ph'] = dump['Nem_phys'][:,:,:,0]
    if 'Nabs_phys' in dump.keys():
      if hdr['RADIATION'] == RADTYPE_NEUTRINOS and nulegacy == False:
        dump['Nabs_e'] = dump['Nabs_phys'][:,:,:,0]
        dump['Nabs_anti'] = dump['Nabs_phys'][:,:,:,1]
        dump['Nabs_x'] = dump['Nabs_phys'][:,:,:,2]
      else:
        dump['Nabs_ph'] = dump['Nabs_phys'][:,:,:,0]
    if hdr['NEUTRINO_OSCILLATIONS'] or hdr['LOCAL_ANGULAR_DISTRIBUTIONS']:
      dump['local_angles'] = dfile['local_angles'][:]
      dump['Gnu'] = dfile['Gnu'][:]
      dump['local_Ns'] = dfile['local_Ns'][:]
      dump['local_wsqr'] = dfile['local_wsqr'][:]
      dump['local_moments'] = dfile['local_moments'][:]
      wmean = dump['local_angles'].sum(axis=3) / dump['local_Ns']
      Ns = dump['local_Ns'][:]
      wb2N = Ns*wmean*wmean
      w2 = dump['local_wsqr'][:]
      ones = np.ones_like(Ns)
      dump['local_stddev'] = np.sqrt(wb2N + Ns*(w2 - wb2N)/(np.minimum(ones, Ns - 1)))
      if hdr['NEUTRINO_OSCILLATIONS']:
        dump['local_osc_rate'] = dfile['local_osc_count'][:]

  ucon, ucov, bcon, bcov = get_state(dump, geom)

  dump['ucon'] = ucon
  dump['ucov'] = ucov
  dump['bcon'] = bcon
  dump['bcov'] = bcov

  dump['bsq'] = (bcon*bcov).sum(axis=-1)
  dump['beta'] = 2.*(dump['PRESS']/(dump['bsq'] + SMALL))

  dump['ucon'] = ucon
  dump['ucov'] = ucov
  dump['bcon'] = bcon
  dump['bcov'] = bcov

  # convenience
  dump['ut'] = ucon[:,:,:,0]
  dump['uX1'] = ucon[:,:,:,1]
  dump['uX2'] = ucon[:,:,:,2]
  dump['uX3'] = ucon[:,:,:,3]
  
  if 'jcon' in dfile.keys():
    jcov = np.zeros([hdr['N1'], hdr['N2'], hdr['N3'], 4])
    for mu in range(4):
      jcov[:,:,:,mu] = (dump['jcon'][:,:,:,:]*geom['gcov'][:,:,None,mu,:]).sum(axis=-1)
    dump['jcov'] = jcov
    dump['j2'] = (jcov[:,:,:,:]*dump['jcon'][:,:,:,:]).sum(axis=-1)

  if hdr['RADIATION']:
    dump['ur'] = (dump['ucon'][:,:,:,None,:]*dump['Rmunu'][:,:,:,:,:]).sum(axis=-1)
    dump['ur'] = (dump['ucov'][:,:,:,:]*dump['ur'][:,:,:,:]).sum(axis=-1)
    dump['ur'] = np.clip(dump['ur'], SMALL, None)
    dump['betar'] = dump['PRESS']/(1./3.*dump['ur'][:,:,:])
    dump['Jem'] = dump['Jrad'][0] # convenience
    dump['Jabs'] = dump['Jrad'][1] # convenience
    dump['Jsc'] = dump['Jrad'][2] + dump['Jrad'][3]
    if 'dtau_avg' in dump.keys():
      dump['dtau_abs'] = np.abs(dump['dtau_avg'][0])
      dump['dtau_dens'] = dump['dtau_tot']/(dump['RHO']*dump['ucon'][...,0])
    if 'Ye' in dump.keys():
      dump['dlepton_rad'] = dump['radG_int'][:,:,:,-2]
      dump['dyedt_rad'] = dump['dlepton_rad']/(dump['RHO']*dump['ucon'][:,:,:,0])

  if hdr['METRIC'] == 'MKS':
    dump['ucon_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_con_3d'],
                                           dump['ucon'])
    dump['ucov_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_cov_3d'],
                                           dump['ucov'])
    dump['bcon_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_con_3d'],
                                           dump['bcon'])
    dump['bcov_bl'] = grid_matrix_multiply(geom['Lambda_h2bl_cov_3d'],
                                           dump['bcov'])

  dump['ucon_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'],
                                           dump['ucon'])
  dump['ucov_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_cov'],
                                           dump['ucov'])
  dump['bcon_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'],
                                           dump['bcon'])
  dump['bcov_cart'] = grid_matrix_multiply(geom['Lambda_h2cart_con'],
                                           dump['bcov'])

  dfile.close()

  return dump

def load_hdr_geom(dumpdir):
  dumps = get_dumps_reduced(dumpdir)
  if dumps:
    hdr = load_hdr(dumps[0])
    geom = load_geom(hdr)
    return hdr,geom
  else:
    return None,None

def load_sadw(fname):
  import h5py
  sadw = {}
  sph = {}
  zoh = {}
  with h5py.File(fname,'r') as f:
    t = f['t'][()]
    for d,n in zip([sadw,sph,zoh],['sadw','sph_avg','zoh']):
      for k,v in f[n].items():
        d[k] = v[()]
  return t,sadw,sph,zoh

def load_diag(path, hdr = None, nulegacy = False, timers = True, twod = False):
  import numpy as np

  # load header
  if hdr is None:
    if twod:
      dfiles = sorted(glob.glob(os.path.join(path,'dump2d_*.h5')))
    else:
      dfiles = sorted(glob.glob(os.path.join(path,'dump_*.h5')))
    if len(dfiles) < 1:
      import util
      util.warn("Cannot read header. No dumps available.")
      sys.exit()
    hdr = load_hdr(dfiles[0])

  diag = {}
  dfile = np.loadtxt(os.path.join(path, 'diag.out')).transpose()
  diag['t']           = dfile[0]
  diag['rmed']        = dfile[1]
  diag['pp']          = dfile[2]
  diag['e']           = dfile[3]
  diag['adiabat rep'] = dfile[4]
  diag['u rep']       = dfile[5]
  diag['mdot']        = dfile[6]
  diag['edot']        = dfile[7]
  diag['ldot']        = dfile[8]
  diag['mass']        = dfile[9]
  diag['egas']        = dfile[10]
  diag['Phi']         = dfile[11]
  diag['phi']         = dfile[12]
  diag['jet_EM_flux'] = dfile[13]
  diag['divbmax']     = dfile[14]
  nbase = 14
  if hdr['RADIATION']:
    diag['step_made']      = dfile[nbase + 1]
    diag['step_abs']       = dfile[nbase + 2]
    diag['step_scatt']     = dfile[nbase + 3]
    diag['step_lost']      = dfile[nbase + 4]
    diag['step_rec']       = dfile[nbase + 5]
    diag['step_tot']       = dfile[nbase + 6]
    diag['step_sent']      = dfile[nbase + 7]
    diag['step_rcvd']      = dfile[nbase + 8]
    diag['step_made_all']  = dfile[nbase + 9]
    diag['step_abs_all']   = dfile[nbase + 10]
    diag['step_scatt_all'] = dfile[nbase + 11]
    diag['step_lost_all']  = dfile[nbase + 12]
    diag['step_rec_all']   = dfile[nbase + 13]
    diag['step_tot_all']   = dfile[nbase + 14]
    diag['step_sent_all']  = dfile[nbase + 15]
    diag['step_rcvd_all']  = dfile[nbase + 16]
    diag['load_imbalance'] = dfile[nbase + 17]
    diag['tune_emiss']     = dfile[nbase + 18]
    diag['tune_scatt']     = dfile[nbase + 19]
    diag['erad']           = dfile[nbase + 20]
    diag['lum']            = dfile[nbase + 21]
    diag['eff']            = dfile[nbase + 22]
    nbase += 22
    diag['Lum'] = diag['lum']*hdr['U_unit']*hdr['L_unit']**3/hdr['T_unit']
    diag['Mdot'] = diag['mdot']*hdr['M_unit']/hdr['T_unit']
    if hdr['ELECTRONS']:
      diag['num_super'] = dfile[nbase + 1]
      diag['lum_super'] = dfile[nbase + 2]
      nbase += 2
    if hdr['RADIATION'] == RADTYPE_NEUTRINOS and nulegacy == False:
      diag['lepton_tot']   = dfile[nbase + 1]
      diag['dlepton_tot']  = dfile[nbase + 2]
      diag['dlepton_perc'] = dfile[nbase + 3]
      nbase += 3
    if hdr['TRACERS'] and nulegacy == False:
      diag['tracer_tot']     = dfile[nbase + 1]
      diag['non_tracer_tot'] = dfile[nbase + 2]
      nbase += 2
  diag['lum_eht'] = dfile[nbase + 1]
  diag['mdot_eh'] = dfile[nbase + 2]
  diag['edot_eh'] = dfile[nbase + 3]
  diag['ldot_eh'] = dfile[nbase + 4]
  nbase += 4
  if timers:
    diag['TIMER_UPDATE']   = dfile[nbase + 1]
    diag['TIMER_FLUXCALC'] = dfile[nbase + 2]
    diag['TIMER_FIXUP']    = dfile[nbase + 3]
    diag['TIMER_BOUND']    = dfile[nbase + 4]
    diag['TIMER_DIAG']     = dfile[nbase + 5]
    diag['TIMER_OUT']      = dfile[nbase + 6]
    diag['TIMER_MAKE']     = dfile[nbase + 7]
    diag['TIMER_PUSH']     = dfile[nbase + 8]
    diag['TIMER_INTERACT'] = dfile[nbase + 9]
    diag['TIMER_MICRO']    = dfile[nbase + 10]
    diag['TIMER_ALL']      = dfile[nbase + 11]
    nbase += 11
    if hdr['ELECTRONS']:
      diag['TIMER_ELECTRON'] = dfile[nbase + 1]
      nbase += 1
    if hdr['NEUTRINO_OSCILLATIONS'] or hdr['LOCAL_ANGULAR_DISTRIBUTIONS']:
      diag['TIMER_OSCILLATIONS'] = dfile[nbase + 1]
      nbase += 1

  # Purge vanishing data due to restarts
  restart_mask = np.ones(len(diag['t']),dtype=bool)
  restart_mask[1:] = diag['t'][1:] - diag['t'][:-1] > 0
  for key in diag:
    diag[key] = diag[key][restart_mask]
  # not sure why I need to do this twice
  restart_mask = np.ones(len(diag['t']),dtype=bool)
  restart_mask[1:] = diag['t'][1:] - diag['t'][:-1] > 0
  for key in diag:
    diag[key] = diag[key][restart_mask]

  # Ignore old data due to restarts
  t = diag['t'][-1]
  ind = [len(diag['t'])-1]
  for i in reversed(range(ind[0])):
    if diag['t'][i] <= t:
      ind.append(i)
      t = diag['t'][i]
  ind = np.unique(np.array(ind))
  for key in diag:
    diag[key] = diag[key][ind]

  # remove points that suspiciously vanish
  mdot = np.abs(diag['mdot'])
  for i in range(1,len(mdot)):
    if mdot[i] < SMALL and mdot[i-1] > SMALL:
      for key in diag:
        diag[key][i] = diag[key][i-1]

  # header
  diag['hdr'] = hdr

  return diag

def delta(mu,nu):
  return 1. if mu == nu else 0.

def get_stresses(hdr,dump):
  "Gets Reynolds and Maxwell stress tensors and related"
  import numpy as np
  # Stresses defined one index up, one index down
  Reynolds = np.zeros((hdr['N1'],hdr['N2'],hdr['N3'],4,4))
  Maxwell = np.zeros_like(Reynolds)
  for mu in range(4):
    for nu in range(4):
      Reynolds[:,:,:,mu,nu] = ((dump['RHO']
                                + dump['UU']
                                + dump['PRESS'])
                               *dump['ucon'][:,:,:,mu]*dump['ucov'][:,:,:,nu]
                               + dump['PRESS']*delta(mu,nu))
      Maxwell[:,:,:,mu,nu] = ((dump['bsq']
                               *dump['ucon'][:,:,:,mu]
                               *dump['ucov'][:,:,:,nu])
                              + 0.5*dump['bsq']*delta(mu,nu)
                              - (dump['bcon'][:,:,:,mu]
                                 *dump['bcov'][:,:,:,nu]))

  TraceReynolds = np.zeros((hdr['N1'],hdr['N2'],hdr['N3']))
  TraceMaxwell = np.zeros_like(TraceReynolds)
  for mu in range(4):
    TraceReynolds += Reynolds[...,mu,mu]
    TraceMaxwell += Maxwell[...,mu,mu]
  
  ReynoldsTF = Reynolds.copy()
  MaxwellTF = Maxwell.copy()
  for mu in range(4):
    for nu in range(4):
      ReynoldsTF[...,mu,nu] -= (1./3.)*TraceReynolds*np.kron(mu,nu)
      MaxwellTF[...,mu,nu] -= (1./3.)*TraceMaxwell*np.kron(mu,nu)

  ReynoldsTFsq = np.zeros((hdr['N1'],hdr['N2'],hdr['N3']))
  MaxwellTFsq = np.zeros_like(ReynoldsTFsq)
  for mu in range(4):
    for nu in range(4):
      ReynoldsTFsq += ReynoldsTF[...,mu,nu]*ReynoldsTF[...,nu,mu]
      MaxwellTFsq += MaxwellTF[...,mu,nu]*MaxwellTF[...,nu,mu]

  ReynoldsTFsqrt = np.sqrt(ReynoldsTFsq)
  MaxwellTFsqrt = np.sqrt(MaxwellTFsq)

  stresses = {}
  stresses['Reynolds'] = Reynolds
  stresses['Maxwell'] = Maxwell
  stresses['T'] = Reynolds + Maxwell
  stresses['Reynolds/Trace'] = TraceReynolds
  stresses['Maxwell/Trace'] = TraceMaxwell
  stresses['Reynolds/TF'] = ReynoldsTF
  stresses['Maxwell/TF'] = MaxwellTF
  stresses['Reynolds/TF/sq'] = ReynoldsTFsq
  stresses['Maxwell/TF/sq'] = MaxwellTFsq
  stresses['Reynolds/TF/sqrt'] = ReynoldsTFsqrt
  stresses['Maxwell/TF/sqrt'] = MaxwellTFsqrt
  return stresses

def get_Lambda_bl2cyl(hdr,geom):
  """Transformation matrices for transforming
  spherical coordinates into cylindrical
  """
  import numpy as np
  r = geom['r']
  th = geom['th']
  cth = np.cos(th)
  sth = np.sin(th)
  Lambda_con = np.zeros((hdr['N1'],hdr['N2'],hdr['N3'],4,4))
  Lambda_cov = np.zeros_like(Lambda_con)

  Lambda_con[:,:,:,0,0] = 1.
  Lambda_con[:,:,:,1,1] = sth
  Lambda_con[:,:,:,1,2] = r*cth
  Lambda_con[:,:,:,2,3] = 1.
  Lambda_con[:,:,:,3,1] = cth
  Lambda_con[:,:,:,3,2] = -r*sth

  Lambda_cov[:,:,:,0,0] = 1.
  Lambda_cov[:,:,:,1,1] = sth
  Lambda_cov[:,:,:,1,3] = cth
  Lambda_cov[:,:,:,2,1] = cth/(r+SMALL)
  Lambda_cov[:,:,:,2,3] = -sth/(r+SMALL)
  Lambda_cov[:,:,:,3,2] = 1.

  return Lambda_cov,Lambda_con  

def get_state(dump, geom):
  import numpy as np
  hdr = dump['hdr']
  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']

  ucon = np.zeros([N1,N2,N3,4])
  ucov = np.zeros([N1,N2,N3,4])
  bcon = np.zeros([N1,N2,N3,4])
  bcov = np.zeros([N1,N2,N3,4])

  gcov = geom['gcov']
  gcon = geom['gcon']

  U1 = dump['U1']
  U2 = dump['U2']
  U3 = dump['U3']
  B1 = dump['B1']
  B2 = dump['B2']
  B3 = dump['B3']

  alpha = geom['alpha']
  qsq = (gcov[:,:,None,1,1]*U1**2 + gcov[:,:,None,2,2]*U2**2 +
         gcov[:,:,None,3,3]*U3**2 + 2.*(gcov[:,:,None,1,2]*U1*U2 +
                                        gcov[:,:,None,1,3]*U1*U3 +
                                        gcov[:,:,None,2,3]*U2*U3))
  gamma = np.sqrt(1. + qsq)
  ucon[:,:,:,0] = gamma/alpha
  ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,None,0,3]

  for mu in range(4):
    ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  for mu in range(4):
    bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  return ucon, ucov, bcon, bcov

def grid_matrix_multiply(M, v, transpose=False):
  import numpy as np
  # numpy crashes with a memory error unless
  # I use this horrible for loop. I don't know why.
  # TODO: fix this.
  out = np.empty_like(v)
  for i in range(M.shape[0]):
    for j in range(M.shape[1]):
      for k in range(M.shape[2]):
        if len(M.shape) > 4:
          if transpose:
            out[i,j,k] = np.dot(M[i,j,k].transpose(),v[i,j,k])
          else:
            out[i,j,k] = np.dot(M[i,j,k],v[i,j,k])
        else:
          if transpose:
            out[i,j,k] = np.dot(M[i,j].transpose(),v[i,j,k])
          else:
            out[i,j,k] = np.dot(M[i,j],v[i,j,k])
  return out


# ======================================================================
# Tracer stuff
# ======================================================================

# ----------------------------------------------------------------------
class TracerDataBase(object):
  def __init__(self,units,data):
    self.units = units
    self.data = data

  @classmethod
  def fromfile(cls, tracer_file, ids=None):
    import h5py
    import numpy as np
    units = {}
    data = {}
    with h5py.File(tracer_file,'r') as f:
      grp = f['units']
      for k in grp.keys():
        units[k] = grp[k][()]
      grp = f['data']
      for k in grp.keys():
        data[k] = grp[k][()]

    if ids is not None:
      mask = np.in1d(data['id'],ids)
      for k,v in data.items():
        data[k] = v[mask]
        
    return cls(units,data)

  @classmethod
  def fromtuple(cls,hdr,units,data):
    return cls(units,data)
  
  @classmethod
  def concatenate(cls,tdlist):
    import numpy as np
    # Empty list returns empty tracerdata
    if not tdlist:
      return cls({},{})
    # Concatenate only nonempty tracerdata
    tdlist = [td for td in tdlist if td]
    # Concatenate
    data = {}
    units = tdlist[0].units
    for k in tdlist[0].data.keys():
      data[k] = [None for td in tdlist]
    for i,td in enumerate(tdlist):
      for k,v in td.data.items():
        data[k][i] = v
    for k,v in data.items():
      data[k] = np.concatenate(v)
    return cls(units,data)

  @classmethod
  def sort_trace(cls,trace,key):
    sort_indcs = trace[key].argsort()
    out = {}
    for k,v in trace.items():
      out[k] = v[sort_indcs]
    return out

  def sort_data(self,key):
    sort_indcs = self.data[key].argsort()
    new_data = {}
    for k,v in self.data.items():
      new_data[k] = v[sort_indcs]
    return new_data

  def discrete_location(self,hdr):
    "Finds indexes i,j,k"
    import numpy as np
    i = ((self['Xharm'][:,0]-hdr['startx'][1])//hdr['dx'][1]).astype(int)
    j = ((self['Xharm'][:,1]-hdr['startx'][2])//hdr['dx'][2]).astype(int)
    k = ((self['Xharm'][:,2]-hdr['startx'][3])//hdr['dx'][3]).astype(int)

    i = np.maximum(0,np.minimum(i,hdr['N1']-1))
    j = np.maximum(0,np.minimum(j,hdr['N2']-1))
    k = np.maximum(0,np.minimum(k,hdr['N3']-1))

    return i,j,k
  
  def filter_data(self,mask):
    import numpy as np
    filtered_data = {}
    for k,v in self.data.items():
      filtered_data[k] = v[mask]
    return filtered_data

  def save(self,filename):
    import h5py
    with h5py.File(filename,'w') as f:
      grp = f.create_group('units')
      for k,v in self.units.items():
        grp.create_dataset(k,data=v)
      grp = f.create_group('data')
      for k,v in self.data.items():
        grp.create_dataset(k,data=v)

  def times(self):
    return sorted(set(self.data['time']))

  def keys(self):
    return self.data.keys()

  def values(self):
    return self.data.values()

  def items(self):
    return self.data.items()

  def __getitem__(self,key):
    return self.data[key]

  def __nonzero__(self):
    return bool(self.units)

  def __bool__(self):
    return bool(self.units)
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
class TracerData(TracerDataBase):
  varnames = {
    'hdr'   : ['dumptrace_id','nstep',
               'ntracers','Time','nstep'],
    'units' : ['B_unit','L_unit','M_unit',
               'RHO_unit','T_unit','U_unit'],
    'data'  : ['Press','T','rho','s',
               'Xcart','id','it','mass','time'],
    'net'   : ['Ye','rate_emitted','rate_absorbed'],
    'optional_data' : ['Ye','Ye_em','uu','Xharm',
                       'ucon','ucov',
                       'rate_emitted','rate_absorbed',
                       'bcon', 'bcov', 'tau_dyn']
    }

  @classmethod
  def fromdir(cls, tracer_dir,
                parallel=False,
                ids=None):
    units,data = load_tracers(tracer_dir,parallel,ids)
    return cls(units,data)

  @classmethod
  def frompath(cls, pth,
                 parallel=False,
                 ids=None):
    if os.path.isdir(pth):
      return cls.fromdir(pth, parallel, ids)
    elif os.path.isfile(pth):
      return cls.fromfile(pth, ids)
    else:
      raise IOError("Path is of unknown type")

  def ids(self):
    return sorted(set(self.data['id']))

  def sort(self, key='id'):
    sorted_data = self.sort_data(key)
    return TracerData(self.units, sorted_data)

  def filter(self,mask):
    filtered_data = self.filter_data(mask)
    return TracerData(self.units,filtered_data)

  def get_trace(self,t_id):
    import numpy as np
    # id
    mask = (self.data['id'] == t_id)
    trace_data = self.filter(mask)
    # uniqueness
    _, mask = np.unique(trace_data['time'], return_index=True)
    trace_data = trace_data.filter(mask)
    # sort
    trace_data = self.sort_trace(trace_data,'time')
    return Trace(self.units,trace_data)

  def remove_trace(self,t_id):
    import numpy as np
    mask = self.data['id'] != t_id
    return self.filter(mask)

  def __len__(self):
    return len(self.ids())
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
class Trace(TracerDataBase):
  def id(self):
    return self.data['id'][0]

  def sort(self):
    return super().sort_trace(self.data,'time')

  def filter(self,mask):
    filtered_data = self.filter_data(mask)
    return Trace(self.units,filtered_data)

  def save(self,filename=None):
    if filename is None:
      filename = 'trace_{:08d}.td'.format(self.id())
    super().save(filename)

  def proper_time(self):
    "Returns proper time tau(t)"
    import numpy as np
    from scipy.integrate import cumtrapz
    tau = np.zeros_like(self['time'])
    tau[1:] = cumtrapz(1/self['ucon'][...,0],self['time'])
    return tau

  def __len__(self):
    return len(self.data['time'])
    
# ----------------------------------------------------------------------


def load_tracers(tracer_dir,
                   parallel=False,
                   ids=None):
  import numpy as np
  import h5py, gc
  files = get_tracer_fnams(tracer_dir)
  tracers = {}
  if not files:
    raise ValueError("This directory is empty!")
  
  data_container = {}
  _,units,data = load_tracer(files[0])
  for d in data.keys():
    data_container[d] = [None for f in files]

  if parallel:
    from concurrent import futures
    if ids is not None:
      inp = [(f,ids) for f in files]
    else:
      inp = files
    with futures.ProcessPoolExecutor() as p:
      data_list = p.map(load_td_only,inp)
    for i,data in enumerate(data_list):
      for k,v in data.items():
        data_container[k][i] = v
  else:
    for i,f in enumerate(files):
      _,_,data = load_tracer(f,ids)
      for k,v in data.items():
        data_container[k][i] = v

  for k,v in data_container.items():
    tracers[k] = np.concatenate(v)

  del data_container
  gc.collect()

  return units,tracers

def load_td_only(inp):
  if type(inp) is str:
    tracer_file = inp
    ids = None
  else:
    tracer_file,ids = inp
  _,_,data = load_tracer(tracer_file,ids)
  return data

def load_tracer(tracer_file,ids=None):
  import h5py
  import numpy as np
  hdr = TracerData.varnames['hdr']
  units = TracerData.varnames['units']
  data = TracerData.varnames['data']
  optional_data = TracerData.varnames['optional_data']
  out_hdr,out_units,out_data = {},{},{}
  with h5py.File(tracer_file,'r') as f:
    for d in optional_data:
      if d in f.keys() or 'Step#0/{}'.format(d) in f.keys():
        data.append(d)
    for h in hdr:
      out_hdr[h] = f[h][0]
    for u in units:
      out_units[u] = f[u][0]
    for d in data:
      if 'Step#0' in f.keys():
        out_data[d] = f['Step#0/{}'.format(d)][()]
      else:
        out_data[d] = f[d][()]

  if ids is not None:
    mask = np.in1d(out_data['id'],ids)
    for k,v in out_data.items():
      out_data[k] = v[mask].copy()

  return out_hdr,out_units,out_data
