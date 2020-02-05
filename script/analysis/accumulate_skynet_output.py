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

# Functions
# ----------------------------------------------------------------------
def get_header(filename):
    "Load header string for abundances file"
    from functools import reduce
    with open(filename,'r') as f:
        dstr = f.read().split('\n\n\n\n')[-1]
        header = [l for l in dstr.split('\n') if '#' in l]
        header = reduce(lambda x,y: "{}\n{}".format(x,y),header)
    return header

def print_header(filename):
    "Print header string from abundances file"
    print(get_header(filename))

def get_data(filename):
    "Read data from file"
    with open(filename,'r') as f:
        data = np.loadtxt(StringIO(f.read().split('\n\n\n\n')[-1]))
    return data

def get_metadata(root,filename):
    "Seeks out metadata file from root directory structure"
    dirname = "{:08d}".format(int(filename.split('/')[-1].split('.')[0]))
    metadata_path = os.path.join(root,dirname,'metadata.dat')
    metadata = np.loadtxt(metadata_path)
    return metadata

def _accumulate_abundances_worker(inp):
    "Worker function"
    fnams = inp[0]
    meta_root = inp[1]
    cut = inp[2]
    data = get_data(fnams[0])
    A = data[:,0]
    Z = data[:,3]
    Abundances_A = np.zeros_like(A)
    Abundances_Z = np.zeros_like(Z)
    total_mass = 0.0
    for fnam in fnams:
        try:
            data = get_data(fnam)
            metadata = get_metadata(meta_root,fnam)
            idx = int(metadata[0])
            mass = metadata[1]
            if cut is None or idx in cut:
                Abundances_A += mass*data[:,1]
                Abundances_Z += mass*data[:,4]
                total_mass += mass
                print("...{} accumulated".format(fnam),idx)
            else:
                print("...{} unused".format(fnam),idx)
        except:
            print("...{} failed".format(fnam),idx)
            continue
    return Abundances_A,Abundances_Z,total_mass

def accumulate_abundances_parallel(skynet_root,meta_root,
                                   nproc=None,
                                   cut=None):
    "Accumulates abundances output from skynet root, using meta root"
    fnams = np.array(sorted(glob.glob(os.path.join(skynet_root,'*.h5.dat'))),
                     dtype=str)
    data = get_data(fnams[0])
    A = data[:,0]
    Z = data[:,3]
    Abundances_A = np.zeros_like(A)
    Abundances_Z = np.zeros_like(Z)
    total_mass = 0.0
    total_files = len(fnams)
    print("Files to accumulate = {}".format(total_files))
    p = Pool(processes=nproc)
    f_array = np.array_split(fnams,p._processes)
    inp = [(f,meta_root,cut) for f in f_array]
    work = p.map(_accumulate_abundances_worker,inp)
    for AA,AZ,m in work:
        Abundances_A += AA
        Abundances_Z += AZ
        total_mass += m
    Abundances_A /= total_mass
    Abundances_Z /= total_mass
    return A,Z,Abundances_A,Abundances_Z,total_mass

def gen_abundances_file(skynet_root,meta_root,savename,
                        nproc=None,
                        cut = None):
    A,Z,AA,AZ,m = accumulate_abundances_parallel(skynet_root,
                                                 meta_root,
                                                 nproc,
                                                 cut)
    data = np.vstack((A,Z,AA,AZ)).transpose()
    import h5py
    with h5py.File(savename + '.h5','w') as f:
        f.create_dataset('A',data=A)
        f.create_dataset('Z',data=Z)
        f.create_dataset('Abundance_A',data=AA)
        f.create_dataset('Abundance_Z',data=AZ)
        f.create_dataset('Mass',data=m)
    np.savetxt(savename + '.dat',
               data,
               comments = '# ',
               header = '[1]:A [2]:Z [3]:Abundance_A [4]:Abundance_Z')
    
# ----------------------------------------------------------------------

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description='Accumulate SkyNet output'
    )
    parser.add_argument('root',type=str,
                        help='Directory containing skynet output. Assumed flat.')
    parser.add_argument('meta',type=str,
                        help='Directory containing metadata. Assumed in PRISM format.')
    parser.add_argument('-s','--save',metavar='STR',
                        dest='savename',type=str,
                        default='abundances',
                        help='Saves data to STR.h5 and STR.dat')
    parser.add_argument('-n','--nproc',type=int,default=None,
                        help='Number of processes to use')
    parser.add_argument('-c','--cut',type=str,default=None,
                        help='File containing ids to use as cut')
    args = parser.parse_args()

    if args.cut is not None:
        cut = np.loadtxt(args.cut)
        cut = cut.astype(int)
    else:
        cut = args.cut

    gen_abundances_file(args.root,args.meta,args.savename,args.nproc,cut)

