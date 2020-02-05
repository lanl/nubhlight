#!/usr/bin/env python

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# A little script to accumulate skynet output

# Imports
# ----------------------------------------------------------------------
from __future__ import print_function, division

import numpy as np
import glob, os
from multiprocessing import Pool

try:
    from io import StringIO
except:
    from StringIO import StringIO
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
def get_header(filename):
    "Load header string for abundances file"
    from functools import reduce
    with open(filename,'r') as f:
        dstr = f.read().split('\n\n\n\n')[0]
        header = [l for l in dstr.split('\n') if '#' in l]
        header = reduce(lambda x,y: "{}\n{}".format(x,y),header)
    return header

def print_header(filename):
    "Print header string from abundances file"
    print(get_header(filename))

def get_data(filename):
    "Read data from file"
    with open(filename,'r') as f:
        data = np.loadtxt(StringIO(f.read().split('\n\n\n\n')[0]))
    return data

def get_metadata(root,filename):
    "Seeks out metadata file from root directory structure"
    dirname = "{:08d}".format(int(filename.split('/')[-1].split('.')[0]))
    metadata_path = os.path.join(root,dirname,'metadata.dat')
    metadata = np.loadtxt(metadata_path)
    return metadata

def get_mass_fractions(data):
    return data[-1,20:22]

def get_data_once(inp):
    fnam = inp[0]
    meta_root = inp[1]
    try:
        mass_fractions = get_mass_fractions(get_data(fnam))
        metadata = get_metadata(meta_root,fnam)
        idx = metadata[0]
        mass = metadata[1]
        out = np.array([idx,mass,mass_fractions[0],mass_fractions[1]])
        print("...{}".format(fnam))
    except:
        out = np.array([-1,-1,-1,-1],dtype=float)
        print("...{} failed".format(fnam))
    return out

def get_all_in_dir(skynet_root,meta_root,nproc=None):
    print("Accumulating in directory {}...".format(skynet_root))
    fnams = sorted(glob.glob(os.path.join(skynet_root,'*.h5.dat')))
    inp = [(f, meta_root) for f in fnams]
    p = Pool(processes=nproc)
    ma_array = p.map(get_data_once,inp)
    mass_fractions = np.vstack(ma_array).transpose()
    mass_fractions = mass_fractions[mass_fractions[:,0] > 0]
    return mass_fractions

def get_all_in_dirs(skynet_dirs,meta_dirs,nproc=None):
    arrays = [None for d in skynet_dirs]
    for i,(s,m) in enumerate(zip(skynet_dirs,meta_dirs)):
        arrays[i] = get_all_in_dir(s,m,nproc)
    mass_fractions = np.vstack(arrays)
    return mass_fractions

def make_file(skynet_dirs,meta_dirs,
              savename='mass_fractions.dat',
              nproc=None):
    mass_fractions = get_all_in_dirs(skynet_dirs,meta_dirs,nproc)
    np.savetxt(savename,mass_fractions,
               comments = '# ',
               header = '[1]:id [2]:mass [3]:X_lanthanides [4]:X_actinides')
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Get mass fractions from SkyNet data')
    parser.add_argument('-r','--root',type=str,nargs='+',
                        required=True,
                        help='List of directories containing SkyNet data')
    parser.add_argument('-m','--meta',type=str,nargs='+',
                        required=True,
                        help='List of directories containing metadata')
    parser.add_argument('-s','--save',type=str,
                        default='mass_fractions.dat',
                        help='Name for file to save')
    parser.add_argument('-n','--nproc',type=int,default=None,
                        help='Number of parallel processes to use')
    args = parser.parse_args()
    if len(args.root) != len(args.meta):
        raise ValueError("Must have equal number of root and meta directories")
    make_file(args.root,args.meta,args.save,args.nproc)
                        
