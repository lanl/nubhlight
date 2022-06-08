#!/usr/bin/env python

"""accumulate_tracers.py
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This little script creates a new tracer file containing tracers that
pass through some radius.
"""

from __future__ import print_function, division
import sys; sys.dont_write_bytecode = True
import hdf5_to_dict as io
import numpy as np
import os, glob, h5py
from multiprocessing import Pool
from units import cgs

def _get_r(data):
    X = data['Xcart']
    r = np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    return r

def _accumulate_files_worker(inpt):
    fnams,r_thresh = inpt
    ids = set([])
    tds = []
    for f in fnams:
        print("...",f,len(ids))
        try:
            hdr,units,data = io.load_tracer(f)
            td = io.TracerData.fromtuple(hdr,units,data)
            r = _get_r(data)
            to_add = td.filter(r >= r_thresh)
            if len(to_add) > 0:
                ids_to_add = set(to_add['id']) - ids
                if len(ids_to_add) > 0:
                    mask = np.in1d(to_add['id'],
                                   np.array(list(ids_to_add)))
                    to_add = to_add.filter(mask)
                    ids |= ids_to_add
                    tds.append(to_add)
        except:
            print("...error with",f,"skipping...")
            continue

    return io.TracerData.concatenate(tds)

def _accumulate_files_worker_t_cond(inpt):
    atol = 0.1
    fnams,t_cond = inpt
    ids = set([])
    tds = []
    for f in fnams:
        print("...",f,len(ids))
        try:
            hdr,units,data = io.load_tracer(f)
            td = io.TracerData.fromtuple(hdr,units,data)
            mask = np.isclose(td['T']*cgs['MEV']/cgs['GK'],
                              t_cond,atol=atol)
            to_add = td.filter(mask)
            if len(to_add) > 0:
                ids_to_add = set(to_add['id']) - ids
                if len(ids_to_add) > 0:
                    mask = np.in1d(to_add['id'],
                                   np.array(list(ids_to_add)))
                    to_add = to_add.filter(mask)
                    ids |= ids_to_add
                    tds.append(to_add)
        except:
            print("...error with",f,"skipping...")
            continue

    return io.TracerData.concatenate(tds)

def accumulate_files_parallel(fnams,
                              r_thresh=1000,
                              t_cond = None,
                              nproc=None):
    p = Pool(processes=nproc)
    fchunks = np.array_split(fnams,p._processes)
    if t_cond is not None:
        inp = [(f,t_cond) for f in fchunks]
        worker = _accumulate_files_worker_t_cond
    else:
        inp = [(f,r_thresh) for f in fchunks]
        worker = _accumulate_files_worker
    print("Accumulating for each file...")
    work = p.map(worker,inp)
    print("filtering...")
    work = [w for w in work if w]
    print("Sorting...")
    work.sort(key=lambda v: v.times()[-1])
    for i in range(len(work)-1):
        if work[i].times()[-1] > work[i+1].times()[0]:
            raise ValueError("Sorting failed!")
    print("Joining parallel work...")
    ids = set([])
    tds = []
    for td in work:
        try:
            ids_to_add = set(td['id']) - ids
            if len(ids_to_add) > 0:
                mask = np.in1d(td['id'],
                               np.array(list(ids_to_add)))
                to_add = td.filter(mask)
                ids |= ids_to_add
                tds.append(to_add)
        except:
            print("...There was an error with a worker. Skipping it.")
            continue
    if len(tds) > 0:
        tracers = io.TracerData.concatenate(tds)
    else:
        raise ValueError("There are no tracers to accumulate.")
    if len(tracers['id']) > len(set(tracers['id'])):
        raise ValueError("Some tracers are in multiple times.")
    return tracers

def accumulate_files(fnams,r_thresh=1000):
    ids = set([])
    data_container = {}
    data_cat = {}
    _,units,data = io.load_tracer(fnams[0])
    for k,v in data.items():
        data_container[k] = []

    for f in fnams:
        print("...",f,",",len(ids))
        try:
            hdr,units,data = io.load_tracer(f)
            td = io.TracerData.fromtuple(hdr,units,data)
            r = _get_r(data)
            to_add = td.filter(r >= r_thresh)
            if len(to_add) > 0:
                ids_to_add = set(to_add['id']) - ids
                if len(ids_to_add) > 0:
                    mask = np.in1d(to_add['id'],
                                   np.array(list(ids_to_add)))
                    to_add = to_add.filter(mask)
                    ids |= ids_to_add
                    for k,v in to_add.items():
                        data_container[k].append(v)
        except:
            print("...error with",f,"skipping...")
            continue

    for k,v in data_container.items():
        data_cat[k] = np.concatenate(v)

    return io.TracerData(units,data_cat)

def accumulate_directory(directory='./dumps/tracers',
                         r_thresh=1000,
                         t_cond=None,
                         savename=None,
                         nproc=None):
    if savename is None:
        savename = './tracers_accumulated_r_{}.td'.format(int(np.floor(r_thresh)))
    fnams = io.get_tracer_fnams(directory)
    print("Accumulating files in directory:"
          +"\n\t{}\n".format(directory)
          +"With radius R = {}\n".format(r_thresh)
          +"File,ntracers:")
    if nproc == 1:
        tracers = accumulate_files(fnams,r_thresh)
    else:
        tracers = accumulate_files_parallel(fnams,r_thresh,t_cond,nproc)
    print("Saving file to {}".format(savename))
    tracers.save(savename)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description=('Save each tracer particle '
                     +'as it passes througha sphere of radius R'))
    parser.add_argument('tracer_dir',type=str,
                        help='Directory containing tracer files')
    parser.add_argument('R',type=float,
                        help='radius at which to accumulate.')
    parser.add_argument('-T','--tthresh',type=float,default=None,
                        help=('If set, extract tracers when they drop below '
                              +'this temperature in GK.'))
    parser.add_argument('-s','--save',type=str,
                        default=None,
                        help=('Name to save file '
                              +'containing accumulated tracers. '
                              +'If not set, default name will be used.'))
    parser.add_argument('-n','--nproc',type=int,default=None,
                        help=("Number of parallel processes to use"))
    args = parser.parse_args()
    accumulate_directory(args.tracer_dir,
                         args.R,
                         args.tthresh,
                         args.save,
                         args.nproc)
