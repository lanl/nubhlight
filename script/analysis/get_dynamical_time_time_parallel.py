#!/usr/bin/evn python

"""get_dynamical_time_time_parallel.py
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This little script computes the integrated characteristic dynamical
time of a tracer. This is performed time-paralle, using a set of
pruned tracer files.
"""

from __future__ import print_function, division
import sys; sys.dont_write_bytecode = True
import hdf5_to_dict as io
import numpy as np
import os, glob, h5py
from multiprocessing import Pool
from units import cgs
from plot_tracers import get_intersection

SAVE_TAU_IN_TRACERS=False
SMALL = 1e-20

def _get_id_worker(filename):
    with h5py.File(filename,'r') as f:
        dump_id = f['dumptrace_id'][0]
    return dump_id

def get_dump_ids(filenames,nproc=None):
    filenames = list(sorted(filenames))
    p = Pool(processes=nproc)
    ids = np.array(p.map(_get_id_worker, filenames))
    return ids

def generate_filename(parent_directory,nstep):
    template = 'tracers_{:08d}.h5part'
    fnam = template.format(nstep)
    return os.path.join(parent_directory, fnam)

def get_td_from_fnam(fnam):
    hdr,units,data = io.load_tracer(fnam)
    td = io.TracerData.fromtuple(hdr,units,data)
    return td

def set_tau_dyn(td_now, td_last):
    # now and last must share same ids, in same order
    td_now_ = get_intersection(td_now, td_last).sort()
    td_last_ = get_intersection(td_last, td_now).sort()
    # numerator
    rho_avg = 0.5*(td_now_['rho'] + td_last_['rho'])
    # denominator
    drho = td_now_['rho'] - td_last_['rho']
    dt = td_now_['time'] - td_last_['time']
    drhodt = np.divide(drho, dt)
    tau_dyn = np.divide(rho_avg, drhodt + SMALL)
    # add to td_now_
    td_now_.data['tau_dyn'] = tau_dyn
    return td_now_

def _tau_dyn_files_worker(inpt):
    parent_dir = inpt[0]
    trace_id_max = inpt[1]
    dump_ids = inpt[2]

    tau_dyns_int = np.zeros(trace_id_max+1)
    dt_int = np.zeros_like(tau_dyns_int)
    
    for dump_id in dump_ids:
        print("...",dump_id)
        fnam_now = generate_filename(parent_dir, dump_id)
        fnam_last = generate_filename(parent_dir, dump_id - 1)
        td_now = get_td_from_fnam(fnam_now)
        td_last = get_td_from_fnam(fnam_last)
        
        dt = td_now['time'][0] - td_last['time'][0]
        td_now_updated = set_tau_dyn(td_now, td_last)
        
        # save file for future reference
        if SAVE_TAU_IN_TRACERS:
            filename = os.path.join(parent_dir,
                                    'tracers_{:08d}.td'.format(dump_id))
            td_now_updated.save(filename)            # convert to cgs + GK
            
        time_unit = td_now_updated.units['T_unit']
        temp_unit = cgs['MEV']/cgs['GK']
        
        # generate arrays. Zero if no tracer is present. This
        # handles the fact that not all tracer files will have all
        # ids.  so we just make a flat index space for IDs.
        temps_now = np.zeros_like(tau_dyns_int)
        tau_dyns_now = np.zeros_like(tau_dyns_int)
        dt_contrib = np.zeros_like(tau_dyns_int)
    
        # Fill arrays
        temps_now[td_now_updated['id']] = td_now_updated['T']
        tau_dyns_now[td_now_updated['id']] = td_now_updated['tau_dyn']
        dt_contrib[td_now_updated['id']] = dt
        
        # unit conversion
        temps_now *= temp_unit
        tau_dyns_now *= time_unit
        
        # kill dynamical times if above 10 GK or below 1 GK
        mask = np.logical_or(temps_now > 10, temps_now < 1)
        tau_dyns_now[mask] = 0
        dt_contrib[mask] = 0
        
        # add to integral/sum (first-order)
        tau_dyns_int += tau_dyns_now*dt
        dt_int += dt_contrib
        # except:
        #     print("...dump id {} failed.".format(dump_id))
        #     continue
    return dt_int, tau_dyns_int

def compute_tau_dyn(directory,
                    accumulated,
                    nproc=None):
    # filenames and file dump ids
    name_template = 'tracers_*.h5part'
    filenames = glob.glob(os.path.join(directory, name_template))
    dump_ids = get_dump_ids(filenames, nproc)

    # Skips dump id 0, which we don't have enough info to process
    dids_hi, dids_lo = dump_ids[2::2], dump_ids[1::2]
    
    # tracer ids
    accumulated = accumulated.sort() # ensure data sorted by id
    trc_ids = accumulated['id']
    max_trc_id = trc_ids.max()

    print("Processing each time...")
    p = Pool(processes=nproc)
    dt_int = np.zeros(max_trc_id+1)
    tau_dyns_int = np.zeros(max_trc_id+1)
    for dids in dids_hi,dids_lo:
        chunks = np.array_split(dids, p._processes)
        inp = [(directory, max_trc_id, c) for c in chunks]
        result = p.map(_tau_dyn_files_worker, inp)
        dt_int += sum([r[0] for r in result])
        tau_dyns_int += sum([r[1] for r in result])

    # reduce to only ids that exist
    dt_int = dt_int[trc_ids]
    tau_dyns_int = tau_dyns_int[trc_ids]
    # take average
    tau_dyns_avg = np.abs(np.divide(tau_dyns_int, dt_int + SMALL))
    # extract masses. Assumes tracers are added to accumulated file
    # only once.
    mass = accumulated['mass']*accumulated.units['M_unit']/cgs['MSOLAR']
    return trc_ids, mass, tau_dyns_avg

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description=("Compute dynamical time scale "
                     +"while particles nucleosynthetically active."))
    parser.add_argument('tracer_dir',type=str,
                        help='Directory containing tracer files')
    parser.add_argument('accumulated',type=str,
                        help=('File containing 1 instance of '
                              +'each unbound tracer'))
    parser.add_argument('-s','--save',type=str,
                        default='dynamical_times.h5',
                    help='Filename to save output')
    parser.add_argument('-n','--nproc',type=int,default=None,
                        help='Number of parallel processes to use')
    args = parser.parse_args()
    print("Using {} for accumulated data.".format(args.accumulated))
    accumulated = io.TracerData.fromfile(args.accumulated)
    print("Walking trhough directory {} with {} procs.".format(args.tracer_dir,
                                                               args.nproc))
    trc_ids, mass, tau_dyn = compute_tau_dyn(args.tracer_dir,
                                             accumulated,
                                             args.nproc)
    print("tau_dyn computed. Saving to file {}".format(args.save))
    print("Units are: Msolar for mass, s for tau_dyn")
    with h5py.File(args.save,'w') as f:
        f.create_dataset('id',data=trc_ids)
        f.create_dataset('mass',data=mass)
        f.create_dataset('tau_dyn',data=tau_dyn)

