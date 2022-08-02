#!/usr/bin/env python

import sys; sys.dont_write_bytecode = True
from hdf5_to_dict import TracerData,Trace
import numpy as np
import os, gc
from argparse import ArgumentParser

def get_trace_paths(trace_path_dir,partitioned=True):
  from functools import reduce
  tracer_files = []
  for r,d,f in os.walk(trace_path_dir):
    if f:
      tracer_files.append([os.path.join(r,fnam) for fnam in f])
  if partitioned:
    return tracer_files
  else:
    return reduce(lambda a,b: a+b, tracer_files)

def get_out_fnam(tid):
    return 'trace_{:08d}.td'.format(tid)

def get_subdirs(outdir,ndirs):
    if ndirs < 2:
        return [outdir]
    return [os.path.join(outdir,'{:04d}'.format(i)) \
            for i in range(ndirs)]

def make_if_needed(path):
    if not os.path.exists(path):
        os.makedirs(path)

def split_ids(ndirs,tids):
    if ndirs >= 2:
        np.random.shuffle(tids)
        tid_lists = np.array_split(tids,ndirs)
    else:
        tid_lists = [tids]
    return tid_lists

def save_traces_in_dir(i,nts,mdir,tidl,tracers):
    for tid in tidl:
        print("\t...{:08d}: {}/{}".format(tid,i,nts))
        savepath = os.path.join(mdir,get_out_fnam(tid))
        trace = tracers.get_trace(tid)
        trace.save(savepath)
        i += 1
    return i

if __name__ == "__main__":
    parser = ArgumentParser(
        description='Generate files containing individual traces from TracerData',
    )
    parser.add_argument('inpath',type=str,
                        help='Path to TracerData to split into traces'
    )
    parser.add_argument('-i','--ids',
                        type=str,
                        default=None,
                        help='Get ids from *.td file provided'
    )
    parser.add_argument('-d','--out_dir',type=str,
                        default='traces',
                        help='Directory to save new traces. Defaults to ./traces'
    )
    parser.add_argument('--serial',dest='parallel',
                        default=True,action='store_false',
                        help='Load data in serial'
    )
    parser.add_argument('-n','--ndirs',
                        type=int,
                        default=200,
                        help='Number of subdirectories to split traces into'
    )
    args = parser.parse_args()
    
    print("Saving traces to directory: {}".format(args.out_dir))
    print("Using {} subdirectories".format(args.ndirs))
    subdirs = get_subdirs(args.out_dir,args.ndirs)
    make_if_needed(args.out_dir)
    for dnam in subdirs:
        make_if_needed(dnam)

    print("Traces will be loaded from: {}".format(args.inpath))

    if args.ids is not None:
        tids = TracerData.frompath(args.ids).ids()
        nts = len(tids)
        print("There are {} traces to save".format(nts))
        
        print("Splitting tracers into subsets...")
        tid_lists = split_ids(args.ndirs,tids)
        
        print("Saving...")
        i = 0
        for mdir,tidl in zip(subdirs,tid_lists):
            print("...in directory {}...".format(mdir))
            print("...Loading traces...")
            tracers = TracerData.frompath(args.inpath,
                                              args.parallel,
                                              ids=tidl)
            gc.collect()
            i = save_traces_in_dir(i,nts,mdir,tidl,tracers)

    else:
        print("Loading traces...")
        tracers = TracerData.frompath(args.inpath,args.parallel)
        gc.collect()

        tids = tracers.ids()
        nts = len(tids)
        print("There are {} traces to save".format(nts))

        print("Splitting tracers into subsets...")
        tid_lists = split_ids(args.ndirs,tids)

        print("Saving...")
        i = 0
        for mdir,tidl in zip(subdirs,tid_lists):
            print("...in directory {}...".format(mdir))
            for tid in tidl:
                i = save_traces_in_dir(i,nts,mdir,tidl,tracers)
    
    print("Done!")    
