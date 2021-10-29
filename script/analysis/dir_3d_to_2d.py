#!/usr/bin/env python

################################################################################
#                                                                              # 
#  Maps a directory of 3d dump files to 2d                                     # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
from argparse import ArgumentParser
import sys, os, h5py, re
import numpy as np
from multiprocessing import Pool
import hdf5_to_dict as io
from dump_3d_to_2d import avg_dump_file
SMALL = 1.e-40

parser = ArgumentParser('Map a directory of 3D files to 2D via averaging in phi.')
parser.add_argument('dir',type=str,
                    help='directory to map.')
parser.add_argument('--rmin',type=float,
                    default=5.5,
                    help='Min radius for Z/H averages')
parser.add_argument('--rmax',type=float,
                    default=25,
                    help='Max radius for Z/H averages')
parser.add_argument('--zmin',type=float,
                    default=-30,
                    help='Min value for z/H')
parser.add_argument('--zmax',type=float,
                    default=30,
                    help='Max value for z/H')
parser.add_argument('-f','--fixghosts',
                    action='store_true',
                    help='Set this flag to set boundary cells by interpolation.')
parser.add_argument('-N2',type=int,
                    default=2,
                    help='Num CPUs in X2. Only relevant if fixing ghosts.')
parser.add_argument('-N3',type=int,
                    default=11,
                    help='Num CPUS in X3. Only relevant if fixing ghosts.')
parser.add_argument('-n','--nproc',type=int,
                    default = None,
                    help = 'Number of parallel processes to use.')
parser.add_argument('-r','--restart', action='store_true',
                    help = 'Restart from where you left off, rather than overwriting old files')

def _avg_dump_file_worker(inp):
    infile = inp[0]
    rmin,rmax,zmin,zmax = inp[1],inp[2],inp[3],inp[4]
    fix_ghosts = inp[5]
    N2CPU,N3CPU = inp[6],inp[7]
    geom = inp[8]
    print("...{}".format(os.path.basename(infile)))
    avg_dump_file(infile,rmin,rmax,zmin,zmax,
                  fix_ghosts,N2CPU,N3CPU,
                  geom)

def avg_dump_dir(dirpath,rmin,rmax,zmin,zmax,
                 fix_ghosts = False,
                 N2CPU=2,N3CPU=11,
                 nproc = None,
                 restart = False):
    "Average dump data in directory. Save 2d files to same directory."
    print("Mapping dumps in {} from 3d to 2d...".format(dirpath))
    fnams = io.get_dumps_full(dirpath)
    if restart:
        fnams_finished = [name.replace('2d','') \
                          for name in io.get_dumps_reduced(dirpath,True)]
        fnams = list(set(fnams) - set(fnams_finished))
    hdr   = io.load_hdr(fnams[0])
    geom  = io.load_geom(hdr)
    inputs = [(fnam,rmin,rmax,zmin,zmax,fix_ghosts,N2CPU,N3CPU,geom) \
              for fnam in fnams]
    p = Pool(processes=nproc)
    p.map(_avg_dump_file_worker,inputs)

if __name__ == "__main__":
    args = parser.parse_args()
    avg_dump_dir(args.dir,
                 args.rmin,args.rmax,
                 args.zmin,args.zmax,
                 args.fixghosts,
                 args.N2,args.N3,
                 args.nproc,
                 args.restart)
