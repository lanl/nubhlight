#!/usr/bin/env python3

from __future__ import print_function, division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis/')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'thermalization'
AUTO = '-auto' in sys.argv
FAST = '-fast' in sys.argv
TF = 5.12 if FAST else 512.
if AUTO:
  import pickle
else:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  import pylab as pl

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call([sys.executable, 'build.py', '-dir', TMP_DIR])
os.chdir('../../test')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
if FAST:
  util.change_rparm('tf',str(TF),'param_template.dat')
bcall(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
geom = io.load_geom(hdr)
t_code = np.zeros(Nd)
Te_code = np.zeros(Nd)
Tr_code = np.zeros(Nd)
#Tr2_code = np.zeros(Nd)
Etot = np.zeros(Nd)
for n in range(Nd):
  dump = io.load_dump(dfiles[n], geom)
  t_code[n] = dump['t']*hdr['T_unit']
  Te_code[n] = dump['Thetae'][0][0][0]*cgs['ME']*cgs['CL']**2/cgs['KBOL']
  Te_code[n] = (hdr['gam']-1.)*dump['UU'][0,0,0]*hdr['U_unit']/(2.*dump['RHO']*hdr['Ne_unit']*cgs['KBOL'])
  Tr_code[n] = (-dump['Rmunu'][0,0,0,0,0]*hdr['U_unit']/cgs['AR'])**(1./4.)
  #Tr2_code[n] = (-dump['Rmunu'][0,0,0,0,0]*hdr['U_unit']/cgs['AR'])**(1./4.)
  ne = dump['RHO'][0,0,0]*hdr['Ne_unit']
  Etot[n] = 2.*ne*cgs['KBOL']*Te_code[n]/(hdr['gam']-1.) + cgs['AR']*Tr_code[n]**4.
print(ne)
print('gam = %e' % hdr['gam'])
diag = io.load_diag(os.path.join(TMP_DIR, 'dumps/'), hdr=hdr)

# SEMIANALYTIC SOLUTION
from scipy.integrate import odeint
tf = 1.e3
t = np.logspace(-4, np.log10(tf), 1000)
nubins = 100
numin = 1.e13
numax = 1.e18
dlnu = (np.log(numax) - np.log(numin))/nubins
nu = np.zeros(nubins)
for n in range(nubins):
  nu[n] = np.exp(np.log(numin) + (0.5 + n)*dlnu)

gamma_e = hdr['gam']
ne = hdr['Ne_unit']
T0 = 1.e6
kb = cgs['KBOL']
me = cgs['ME']
cl = cgs['CL']
ar = cgs['AR']
h = cgs['HPL']
qe = cgs['QE']
def jnu(nu, ne, thetae):
  Te = thetae*me*cl*cl/kb
  x = h*nu/(kb*Te)
  return 5.4e-39*ne**2.*Te**(-1./2.)*np.exp(-x)
def jnu(nu, ne, thetae):
  Te = thetae*me*cl*cl/kb
  x = h*nu/(kb*Te)
  rel = (1. + 4.4e-10*Te)
  gff = 1.2

  if x < 1.e-3:
    efac = (24 - 24*x + 12*x*x - 4.*x*x*x + x*x*x*x)/24.
  else:
    efac = np.exp(-x)

  jv = 1./(4.*np.pi)*pow(2.,5.)*np.pi*pow(qe,6)/(3.*me*pow(cl,3))
  jv *= pow(2.*np.pi/(3.*kb*me),1./2.)
  jv *= pow(Te,-1./2.)*ne*ne
  jv *= efac*rel*gff
  return jv

def Bnu_inv(nu, thetae):
  x = h*nu/(me*cl*cl*thetae)
  return (2.*h/(cl*cl))/(np.expm1(x))
def dydt(y, t0):
  source = np.zeros(len(y))

  Te = y[0]*(gamma_e - 1.)/(2.*ne*kb)
  thetae = kb*Te/(cgs['ME']*cgs['CL']**2)
  for n in range(nubins):
    Bnu_v = nu[n]**3.*Bnu_inv(nu[n], thetae)
    jnu_v = jnu(nu[n], ne, thetae)
    source[n+1] = cl*jnu_v*(4.*np.pi/cl - y[n+1]/Bnu_v)
    #source[n+1] = cl*jnu_v*(4.*np.pi/cl)
  source[0] = -sum(source[1:]*nu*dlnu)
  return source
y0 = np.zeros(nubins + 1)
y0[0] = 2.*ne*kb*T0/(gamma_e - 1.)
ans = odeint(dydt, y0, t)
ue_v = ans[:,0]
Te_v = ue_v*(gamma_e-1.)/(2.*ne*kb)
ur_v = np.zeros(len(t))
for n in range(len(t)):
  ur_v[n] = sum(ans[n,1:]*nu*dlnu)
Tr_v = (ur_v/ar)**(1./4.)

if AUTO:
  data = {}
  data['SOL'] = [t, Te_v]
  data['CODE'] = [t_code, Te_code]
  data['THRESHOLD'] = 0.02
  pickle.dump(data, open('data.p', 'wb'))
  # CLEAN UP
  util.safe_remove(TMP_DIR)
  sys.exit()

# MAKE FIGURE
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,1,1)
ax.plot(t_code, Te_code, color='r', linestyle=code_ls, marker=code_mrk,
  markersize=10)
ax.plot(t_code, Tr_code, color='b', linestyle=code_ls, marker=code_mrk,
  markersize=10)
ax.plot(t, Te_v, color='k', linestyle='--')
ax.plot(t, Tr_v, color='k', linestyle='--')
plt.yscale('log')
plt.ylim([0,1.e6])
plt.xlabel('t (s)'); plt.ylabel('Te (K)')
plt.xscale('log')
plt.xlim([1.e-3, 1.e0])
plt.savefig('thermalization.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

