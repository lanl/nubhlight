#!/usr/bin/env python

################################################################################
#                                                                              # 
#  PLOT SEVERAL DIAG/IN-SITU QUANTITIES AT ONCE                                # 
#                                                                              # 
################################################################################

import matplotlib
import matplotlib.pyplot as plt
import plot as bplt
import sys, os
import hdf5_to_dict as io
import numpy as np
font = {'size' : 14}
matplotlib.rc('font', **font)

from argparse import ArgumentParser
parser = ArgumentParser(description='Summarize some ab-initio quantities.')
parser.add_argument('dumpfolder',type=str,
                    help='Folder containing dumps')
parser.add_argument('mode',type=str,
                    choices=['fluid','rad','perf'],
                    help=('Mode to use for output. '
                          +'Fluid, radiation, or performance.'))
parser.add_argument('-s','--save',type=str, metavar='plotname',
                    default=None,
                    help='Name to save your figure if you want to save.')
args = parser.parse_args()
folder = args.dumpfolder
if not os.path.exists(folder):
  print('ERROR File ' + folder + ' does not exist!')
  sys.exit()
mode = args.mode

diag = io.load_diag(folder)

if mode == 'perf':
  #raise ValueError("Performance not supported.")
  fig = plt.figure(figsize=(14,12))
  
  vals = ['TIMER_UPDATE', 'TIMER_FLUXCALC', 'TIMER_FIXUP', 'TIMER_BOUND', 
    'TIMER_DIAG', 'TIMER_OUT', 'TIMER_MAKE', 'TIMER_PUSH', 
    'TIMER_INTERACT', 'TIMER_ALL']
  if 'TIMER_MICRO' in diag.keys():
    vals.append('TIMER_MICRO')
  if diag['hdr']['ELECTRONS']:
    vals.append('TIMER_ELECTRON')
  for n, val in enumerate(vals):
    ax = fig.add_subplot(4, 3, n+1)
    ax.plot(diag['t'], diag[val])
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    if n < 9:
      ax.set_xticks([])
    else:
      ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.set_title(val)

elif mode == 'fluid':
  fig = plt.figure(figsize=(10,10))

  ax = fig.add_subplot(4,2,1)
  ax.plot(diag['t'], -np.array(diag['mdot']))
  #ax.set_yscale('log')
  ax.set_ylabel('mdot')

  ax = fig.add_subplot(4,2,2)
  ax.plot(diag['t'], np.array(diag['edot']))
  ax.set_yscale('log')
  ax.set_ylabel('edot')

  ax = fig.add_subplot(4,2,3)
  ax.plot(diag['t'], -np.array(diag['ldot']))
  ax.set_yscale('log')
  ax.set_ylabel('ldot')

  ax = fig.add_subplot(4,2,4)
  ax.plot(diag['t'], np.array(diag['mass']))
  ax.set_yscale('log')
  ax.set_ylabel('mass')

  ax = fig.add_subplot(4,2,5)
  ax.plot(diag['t'], np.array(diag['egas']))
  ax.set_yscale('log')
  ax.set_ylabel('egas')

  ax = fig.add_subplot(4,2,6)
  diag['phi'] = diag['Phi']/np.sqrt(-diag['mdot'])
  ax.plot(diag['t'], np.array(diag['phi']))
  #ax.set_yscale('log')
  ax.set_ylabel('phi')

  ax = fig.add_subplot(4,2,7)
  #ax.plot(diag['t'], np.array(diag['jet_EM_flux']))
  ax.plot(diag['t'], np.array(np.fabs(diag['jet_EM_flux']/diag['mdot'])))
  ax.set_yscale('log')
  ax.set_ylabel('jet')

  ax = fig.add_subplot(4,2,8)
  ax.plot(diag['t'], np.array(diag['divbmax']))
  ax.set_yscale('log')
  ax.set_ylabel('divbmax')

  for n in range(8):
    ax = fig.add_subplot(4,2,n+1)
    if n % 2 != 0:
      ax.yaxis.tick_right()
      ax.yaxis.set_label_position('right')
    if n < 6:
      ax.set_xticks([])
    else:
      ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

elif mode == 'rad':

  fig = plt.figure(figsize=(10,10))

  ax = fig.add_subplot(4,2,1)
  ax.plot(diag['t'], -np.array(diag['Mdot'])/diag['hdr']['MdotEdd'])
  ax.set_yscale('log')
  ax.set_ylabel('mdot')

  ax = fig.add_subplot(4,2,2)
  ax.plot(diag['t'], np.array(diag['Lum'])/diag['hdr']['LEdd'])
  ax.set_yscale('log')
  ax.set_ylabel('L/LEdd')

  ax = fig.add_subplot(4,2,3)
  ax.plot(diag['t'], diag['tune_emiss'])
  ax.set_yscale('log')
  ax.set_ylabel('tune_emiss')

  ax = fig.add_subplot(4,2,4)
  ax.plot(diag['t'], diag['step_made_all'])
  ax.set_yscale('log')
  ax.set_ylabel('made')

  ax = fig.add_subplot(4,2,5)
  ax.plot(diag['t'], diag['tune_scatt'])
  ax.set_yscale('log')
  ax.set_ylabel('tune_scatt')

  ax = fig.add_subplot(4,2,6)
  print(diag['step_scatt_all'].max())
  ax.plot(diag['t'], diag['step_scatt_all'])
  ax.set_yscale('log')
  ax.set_ylabel('scatt')

  ax = fig.add_subplot(4,2,7)
  ax.plot(diag['t'], diag['egas'])
  ax.plot(diag['t'], -diag['erad'])
  ax.set_yscale('log')
  ax.set_ylabel('egas erad')

  ax = fig.add_subplot(4,2,8)
  ax.plot(diag['t'], diag['step_tot_all'])
  ax.set_yscale('log')
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position('right')
  ax.set_ylabel('tot')
    
  for n in range(8):
    ax = fig.add_subplot(4,2,n+1)
    if n % 2 != 0:
      ax.yaxis.tick_right()
      ax.yaxis.set_label_position('right')
    if n < 6:
      ax.set_xticks([])
    else:
      ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

if args.save is not None:
  plt.savefig(args.save,bbox_inches='tight')
else:
  plt.show()
