#!/usr/bin/env python

import sys; sys.dont_write_bytecode = True
import hdf5_to_dict as io
import numpy as np
import os, glob, h5py
from multiprocessing import Pool

def prune_directory(tracer_dir, base_file, conditions,
                    out_dir = 'prune',
                    serial = False,
                    nprocs = None):
    print("Saving to directory: {}".format(out_dir))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print("Using {} to get ids to keep.".format(base_file))
    ids = _get_ids_to_keep(base_file, conditions)
    filenames = sorted(glob.glob(os.path.join(tracer_dir,'*tracers*h5part')))
    print("Iterating through tracer files in directory: {}".format(tracer_dir))
    if serial:
        for fnam in filenames:
            prune_file(fnam,ids, out_dir)
    else:
        p = Pool(processes=nprocs)
        p.starmap(prune_file,zip(filenames,
                                 [ids for f in filenames],
                                 [out_dir for f in filenames]))

def prune_file(infile_name, ids_to_keep, directory,legacy_mode=False):
    outfile_name = os.path.join(directory,os.path.basename(infile_name))
    hdr,units,data = io.load_tracer(infile_name)
    mask = np.in1d(data['id'], ids_to_keep)
    ntracers = np.sum(mask)
    with h5py.File(outfile_name,'w') as dest:
        with h5py.File(infile_name,'r') as src:
            _copy_header(src,dest)
        for k, v in data.items():
            dsetname=k if legacy_mode else 'Step#0/{}'.format(k)
            dest.create_dataset(dsetname,data = v[mask])
        dest.create_dataset('ntracers', data = np.array([ntracers]))
    print("\t{} ... done".format(outfile_name))

def _get_ids_to_keep(filename, conditions):
    try:
        data = io.TracerData.fromfile(filename)
    except:
        hdr,units,data = io.load_tracer(filename)
    mask = conditions(data)
    return data['id'][mask]

def _get_r(data):
    X = data['Xcart']
    r = np.sqrt(X[:,0]**2 + X[:,1]**2 + X[:,2]**2)
    return r

def _get_theta(data):
    X = data['Xcart']
    rcyl = np.sqrt(X[:,0]**2 + X[:,1]**2)
    z = X[:,2]
    theta = np.arctan2(z,rcyl)
    return theta

def _cond_survived(data):
    return np.ones_like(data['id'],dtype=bool)

def _cond_rmin(data,rmin = 100):
    r = _get_r(data)
    return r >= rmin

def _cond_rmax(data,rmax):
    r = _get_r(data)
    return r <= rmax

def _cond_pole(data,tthresh):
    tthresh *= np.pi/180.
    theta = _get_theta(data)
    return np.abs(np.abs(theta) - np.pi) >= tthresh

def _cond_zmax(data,zmax):
    z = data['Xcart'][:,2]
    return np.abs(z) <= zmax

def _cond_theta(data,theta0,rng):
    theta = _get_theta(data)
    theta0 *= np.pi/180
    rng *= np.pi/180
    return np.logical_or(np.abs(theta-theta0) <= rng,
                         np.abs(theta+theta0) <= rng)
    
def _copy_header(src,dest):
    variables = (io.TracerData.varnames['hdr']
                 +io.TracerData.varnames['units'])
    variables.remove('ntracers')
    for v in variables:
        if v in src.keys():
            if v not in dest.keys():
                src.copy(v,dest)
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description='Remove tracers from files, based on conditions.')
    parser.add_argument('tracer_dir',type=str,
                        help='Directory containing tracer files you want to prune.')
    parser.add_argument('base_file',type=str,
                        help=('Tracer used to condition pruning.'
                              +'All tracers not present in this file will be pruned.'))
    parser.add_argument('-d', '--outdir', dest='out_dir',
                        type=str,
                        default='prune',
                        help="Directory to save pruned files. Defaults to './prune'")
    parser.add_argument('--rmin',dest='rmin',type=float,default=None,
                        help="If set, tracers with r <= RMIN are removed.")
    parser.add_argument('--rmax',dest='rmax',type=float,default=None,
                        help="If set, tracers with r >= RMAX are removed.")
    parser.add_argument('--pole', dest='tthresh',type=float,default=None,
                        help=("If set, tracers closer than "
                              +"TTHRESH degrees to the poles are removed."))
    parser.add_argument('--zmax',dest='zmax',type=float,default=None,
                        help=("If set, tracers with z > zmax are removed."))
    parser.add_argument('--theta',dest='theta',nargs=2,type=float,default=None,
                        action='store',metavar='TR',
                        help=("If set, the first argument is a value of theta."
                              +" The second is a range. Selects tracers "
                              +" within RANGE degrees of THETA."))
    parser.add_argument('-n','--nprocs',type=int,default=None,
                        help='Number of parallel procs to use')
    parser.add_argument('--serial',dest='serial',
                        default=False,action='store_true',
                        help='Run in serial')

    args = parser.parse_args()

    conds_list = []
    if args.rmin is not None:
        conds_list.append(lambda data: _cond_rmin(data,args.rmin))
    if args.rmax is not None:
        conds_list.append(lambda data: _cond_rmax(data,args.rmax))
    if args.tthresh is not None:
        conds_list.append(lambda data: _cond_pole(data,args.tthresh))
    if args.zmax is not None:
        conds_list.append(lambda data: _cond_zmax(data,args.zmax))
    if args.theta is not None:
        conds_list.append(lambda data: _cond_theta(data,
                                                   args.theta[0],
                                                   args.theta[1]))

    def conds(data):
        out = _cond_survived(data)
        for c in conds_list:
            out = np.logical_and(out,c(data))
        return out

    prune_directory(args.tracer_dir,args.base_file,
                    conds,args.out_dir,
                    args.serial,args.nprocs)
