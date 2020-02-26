#!/usr/bin/env python
import sys, os
from argparse import ArgumentParser
parser = ArgumentParser(
  description='Plot 2d slices of your simulation at a given time.')
parser.add_argument('dumpfile',type=str,
                    help='File name to plot')
parser.add_argument('variable',type=str,
                    help='Variable to plot')
parser.add_argument('--coords',type=str,
                    choices=['cart','mks'],default='cart',
                    help='Coordinate system. Cartesian or Modified Kerr-Schild')
parser.add_argument('-s','--size',
                    type=float,default=40,
                    help='Size of domain to plot')
parser.add_argument('--lin',
                    dest='log',
                    default=True,
                    action='store_false',
                    help='Sets scale to linear. Default is log.')
parser.add_argument('--vmin',
                    type=float,default=-4,
                    help='Colormap lower bound')
parser.add_argument('--vmax',
                    type=float,default=0,
                    help='Colormap upper bound')
parser.add_argument('-c','--cmap',
                    type=str,default='jet',
                    help='Colormap used')
parser.add_argument('--save',
                    type=str,default=None,
                    help='Figure filename if you want to save the figure')
parser.add_argument('--label',
                    type=str,default=None,
                    help='Label for colormap')
parser.add_argument('-i','--index',
                    type=int,default=None,
                    help='Index to plot in multi-d quantity')

def make_snap(dfnam,vnam,coords,size,cmap,logplot,
              savefig,label,
              vmin,vmax,
              index=None,
              geom=None):

  if not os.path.exists(dfnam):
    print('ERROR File ' + dfnam + ' does not exist!')
    sys.exit()

  import matplotlib
  if savefig is not None and savefig:
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  import plot as bplt
  import hdf5_to_dict as io
  import numpy as np
  font = {'size' : 16}
  matplotlib.rc('font', **font)

  hdr = io.load_hdr(dfnam)
  if geom is None:
    geom = io.load_geom(hdr)
  dump = io.load_dump(dfnam,geom=geom)

  if not vnam in dump.keys():
    print('ERROR Variable ' + vnam + ' is not contained in dump file!')
    print('Available variables are:')
    for key in dump.keys():
      if type(dump[key]) is np.ndarray:
        if (len(dump[key].shape) == 3 and dump[key].shape[0] == hdr['N1'] and
            dump[key].shape[1] == hdr['N2'] and dump[key].shape[2] == hdr['N3']):
          print(key, end=' ')
    print('')
    sys.exit()

  IS_3D = hdr['N3'] > 1

  var = dump[vnam]

  if index is not None:
    var = var[...,index]

  if logplot:
    var = np.log10(var)

  if label is None:
    label = vnam
    if index is not None:
      label += '[{}]'.format(index)

  if IS_3D:
    if coords == 'mks':
      fig, (a0, a1) = plt.subplots(1,2,gridspec_kw={'width_ratios':[1,1]}, figsize=(12,6))
    elif coords == 'cart':
      fig, (a0, a1) = plt.subplots(1,2,gridspec_kw={'width_ratios':[1,1]}, figsize=(13,7))
    ax = a0
    if coords == 'mks':
      bplt.plot_X1X2(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax, 
        cbar=False, label=label, ticks=None, shading='gouraud')
    elif coords == 'cart':
      bplt.plot_xz(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
        cbar=False, label=label, ticks=None, shading='gouraud')
      ax.set_xlim([-size,size]); ax.set_ylim([-size,size])
    ax = a1
    if coords == 'mks':
      bplt.plot_X1X3(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax, 
        cbar=True, label=label, ticks=None, shading='gouraud')
    elif coords == 'cart':
      bplt.plot_xy(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
        cbar=True, label=label, ticks=None, shading='gouraud')
      ax.set_xlim([-size,size]); ax.set_ylim([-size,size])
      
  else:
    if coords == 'mks':
      fig, ax = plt.subplots(1, 1, figsize=(10, 10))
      bplt.plot_X1X2(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
        cbar=True, label=label, ticks=None, shading='gouraud')
    elif coords == 'cart':
      fig, ax = plt.subplots(1, 1, figsize=(7, 10))
      bplt.plot_xz(ax, geom, var, dump, cmap=cmap, vmin=vmin, vmax=vmax,
        cbar=True, label=label, ticks=None, shading='gouraud')
      ax.set_xlim([0,size]); ax.set_ylim([-size,size])

  if savefig == False:
    plt.show()
  else:
    plt.savefig(savefig,bbox_inches='tight')
  plt.cla()
  plt.clf()
  plt.close()

if __name__ == "__main__":
  args = parser.parse_args()

  dfnam = args.dumpfile
  vnam = args.variable
  coords = args.coords
  size = args.size
  cmap = args.cmap
  logplot = args.log
  savefig = args.save if args.save is not None else False
  label = args.label
  vmin,vmax = args.vmin,args.vmax
  index = args.index

  make_snap(dfnam,vnam,coords,size,cmap,
            logplot,savefig,label,vmin,vmax,
            index=index)
