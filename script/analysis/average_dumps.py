#!/usr/bin/env python

################################################################################
#                                                                              # 
#  Average Dumps in Time                                                       # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
import sys, os, h5py, re
from argparse import ArgumentParser
import numpy as np
import hdf5_to_dict as io
from glob import glob
from dump_3d_to_2d import get_data, copy_hdr

SMALL = 1.e-40

parser = ArgumentParser('Average a directory of datafiles')
parser.add_argument('dir',type=str,
                    help='Directory containing dump files')
parser.add_argument('--imin',type=int,default=0,
                    help='Min index of dumps to start average')
parser.add_argument('--imax',type=int,default=-1,
                    help='Max index of dumps to start average')
parser.add_argument('-d','--twod',action='store_true',
                    help='Use two-dimensional dumps instead of 3d dumps')

def new_filename(imin,imax):
    return "tavg_dump_i_%08d_%08d.h5" % (imin, imax)

def avg_dump_files(fnams,imin=0,imax=-1,geom=None):
    if imax == -1:
        imax = len(fnams) - 1

    assert imin < imax
    assert len(fnams) >= imax+1
    
    out_fnam = new_filename(imin,imax)
    out_fnam = os.path.join(os.path.dirname(fnams[imin]),out_fnam)
    print("Saving file to ",out_fnam)

    avgs = {}

    print("Setup...")
    with h5py.File(out_fnam,'w') as dest:
        with h5py.File(fnams[imin],'r') as src:
            copy_hdr(src,dest,False)
            for name,var in src.items():
                if name in dest.keys():
                    continue
                var = get_data(var)
                if (var is not None
                    and type(var) is np.ndarray):
                    avgs[name] = np.zeros_like(var)
                        
        # Note: Assumes all dumps have same output cadence
        print("Averaging dumps...")
        for i in range(imin,imax):
            print("...",i)
            with h5py.File(fnams[i],'r') as src:
                for name in avgs.keys():
                    var = get_data(src[name])
                    if var is not None:
                        avgs[name] += var
        print("Saving to file...")
        for name in avgs.keys():
            avgs[name] /= float(imax - imin)
            dest.create_dataset(name,data=avgs[name])

        # Attributes need to be copied too!
        with h5py.File(fnams[imax],'r') as src:
            for n,v in src['P'].attrs.items():
                dest['P'].attrs.create(n,v)    

    print("Done.")
    return

def average_dir(dirname,imin=0,imax=-1,twod=False,geom=None):
    print("Averaging dump files in directory ",dirname)
    if twod:
        fnams = sorted(glob(os.path.join(dirname,'dump2d_*.h5')))
        print("Using 2d averaged dumps")
    else:
        fnams = sorted(glob(os.path.join(dirname,'dump_*.h5')))
        print("Using default dumps")
    avg_dump_files(fnams,imin,imax,geom)

if __name__ == "__main__":
    args = parser.parse_args()
    average_dir(args.dir,args.imin,args.imax,args.twod)
    
