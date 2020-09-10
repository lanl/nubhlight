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
PROBLEM = 'comptonization'
AUTO = '-auto' in sys.argv
FAST = '-fast' in sys.argv
TF = 100. if FAST else 10000.
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
call_string = ['./bhlight', '-p', 'param_template.dat']
bcall(call_string)
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
#ng = 2.38e18
ng = 2.536851e18
Tf = 4.895207e+06
for n in range(Nd):
  dump = io.load_dump(dfiles[n], geom)
  t_code[n] = dump['t']*hdr['T_unit']
  Te_code[n] = dump['Thetae'][0][0][0]*cgs['ME']*cgs['CL']**2/cgs['KBOL']
  #Te_code[n] = (hdr['gam']-1.)*dump['UU'][0,0,0]*hdr['U_unit']/(2.*dump['RHO']*hdr['Ne_unit']*cgs['KBOL'])
  #Tr_code[n] = (-dump['Rmunu'][0,0,0,0,0]*hdr['U_unit']/cgs['AR'])**(1./4.)
  Tr_code[n] = -dump['Rmunu'][0,0,0,0,0]*hdr['U_unit']/(3.*ng*cgs['KBOL'])
  #Tr2_code[n] = (-dump['Rmunu'][0,0,0,0,0]*hdr['U_unit']/cgs['AR'])**(1./4.)
  #ne = dump['RHO'][0,0,0]*hdr['Ne_unit']
  #Etot[n] = 2.*ne*cgs['KBOL']*Te_code[n]/(hdr['gam']-1.) + cgs['AR']*Tr_code[n]**4.
#print ne
#print 'gam = %e' % hdr['gam']
#diag = io.load_diag(os.path.join(TMP_DIR, 'dumps/'), hdr=hdr)

if AUTO:
  data = {}
  n0 = 3*len(t_code)//4
  data['SOL'] = [t_code[n0:], (t_code[n0:]*0. + 1.)*Tf]
  data['CODE'] = [t_code[n0:], Te_code[n0:]]
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
ax.axhline(Tf, color='k', linestyle='--')
#ax.plot(t, Te_v, color='k', linestyle='--')
#ax.plot(t, Tr_v, color='k', linestyle='--')
plt.yscale('log')
#plt.ylim([0,1.e6])
plt.xlabel('t (s)'); plt.ylabel('Te (K)')
plt.xscale('log')
plt.xlim([1.e-3, 1.e0])
plt.savefig(PROBLEM + '.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

