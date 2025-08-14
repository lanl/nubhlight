#!/usr/bin/env python

import sys, os
mlocation = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(mlocation,'..'))
import util
import hdf5_to_dict as io
from snap import make_snap
from multiprocessing import Pool
import gc

# Set your ffmpeg executable here. If not available, set codec to None
codec = 'ffmpeg'
mnam_base = 'anim'

from argparse import ArgumentParser
parser = ArgumentParser(
  description='Make a movie of your simulation.')
parser.add_argument('dumpfolder',type=str,
                    help='Folder containing data')
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
parser.add_argument('--label',
                    type=str,default=None,
                    help='Label for colormap')
parser.add_argument('-i','--index',
                    type=int,default=None,
                    help='Index to plot in multi-d quantity')
parser.add_argument('--cache',
                    default=None,
                    type=str,
                    help='Directory to cache frames in for restarts')
parser.add_argument('--serial',
                    dest='serial',
                    default=False,
                    action='store_true',
                    help='Run in serial')
parser.add_argument('--nproc',
                    dest='nproc',
                    default=None,
                    type=int,
                    help=('Number of parallel processe to use. '
                          +'If not set, defaults to all available cores.'))
parser.add_argument('--twod',
                    dest='twod',
                    default=False,
                    action='store_true',
                    help=('For 2d dump files'))


args = parser.parse_args()

dfold = util.sanitize_path(args.dumpfolder)
if not os.path.exists(dfold):
  print('ERROR Folder ' + dfnam + ' does not exist!')
  sys.exit()

if args.cache is not None:
  tmpdir = args.cache
  if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
else:
  tmpdir = 'FRAMES'
  util.safe_remove(tmpdir)
  os.mkdir(tmpdir)

dfnams = io.get_dumps_full(dfold, twod=args.twod)
hdr = io.load_hdr(dfnams[0])
geom = io.load_geom(hdr)
num_files = len(dfnams)

pairs = list(enumerate(dfnams))
if args.cache:
  fmap = {"frame_%08d.png" % i : (i, d) for (i,d) in pairs}
  all_frames = set(["frame_%08d.png" % i for i in range(num_files)])
  frames = set([f for f in os.listdir(tmpdir) if '.png' in f])
  frames_to_do = all_frames - frames
  pairs = [fmap[f] for f in frames_to_do]

def make_frame(pair):
  i,d = pair
  gc.collect()
  print("frame %d/%d" % (i,num_files))
  make_snap(d,args.variable,args.coords,
            args.size,args.cmap,args.log,
            os.path.join(tmpdir,'frame_%08d.png'  % i),
            args.label,
            args.vmin,args.vmax,
            index=args.index,
            geom=geom)

if args.serial:
  for pair in pairs:
    make_frame(pair)
else:
  p = Pool(processes=args.nproc)
  p.map(make_frame,pairs)

if codec is not None:
  from subprocess import call
  mnam = mnam_base + '_' + args.variable + '_' + args.coords + '.mp4'
  util.safe_remove(mnam)
  call([codec, '-i', os.path.join(tmpdir, 'frame_%08d.png'), mnam])
  util.safe_remove(tmpdir)
