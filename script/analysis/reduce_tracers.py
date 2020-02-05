#!/usr/bin/env python

import sys; sys.dont_write_bytecode = True
import hdf5_to_dict as io
import numpy as np
import os, glob, h5py

def reduce_directory(tracer_dir, out_dir = 'reduced'):
    print("Saving to directory: {}".format(out_dir))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    filenames = sorted(glob.glob(os.path.join(tracer_dir,'*tracers*h5part')))
    print("Iterating through tracer files in directory: {}".format(tracer_dir))
    for fnam in filenames:
        reduce_file(fnam,out_dir)

def reduce_file(infile_name,directory):
    outfile_name = os.path.join(directory,os.path.basename(infile_name))
    with h5py.File(outfile_name,'w') as dest:
        with h5py.File(infile_name,'r') as src:
            _copy_header(src,dest)
            _copy_data(src,dest)
    print("\t{} ... done".format(outfile_name))

def _copy_header(src,dest):
    variables = (io.TracerData.varnames['hdr']
                 +io.TracerData.varnames['units'])
    for v in variables:
        if v in src.keys():
            if v not in dest.keys():
                src.copy(v,dest)

def _copy_data(src,dest):
    variables = (io.TracerData.varnames['data']
                 + io.TracerData.varnames['net'])
    for v in variables:
        if _richin(v,src.keys()):
            vrich = _getrich(v,src.keys())
            if not _richin(v,dest.keys()):
                src.copy(vrich,dest)
                
def _richin(v,k):
    return v in k or 'Step#0/{}'.format(v) in k

def _getrich(v,k):
    if v in k:
        return v
    else:
        return 'Step#0/{}'.format(v)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description='Remove several extraneous datasets from tracer data')
    parser.add_argument('tracer_dir',type=str,
                        help='Directory containing tracer files you want to reduce')
    parser.add_argument('-d', '--outdir', dest='out_dir',
                        type=str,
                        default='reduced',
                        help="Directory to save files. Defaults to './reduced'")
    args = parser.parse_args()

    reduce_directory(args.tracer_dir,args.out_dir)
