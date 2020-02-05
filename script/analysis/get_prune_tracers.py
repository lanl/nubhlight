import sys; sys.dont_write_bytecode = True
import hdf5_to_dict as io
import numpy as np
import h5py
import glob

def copy_header(src,dest):
    variables = ['B_unit','L_unit','M_unit',
                 'T_unit','U_unit',
                 'RHO_unit','TEMP_unit',
                 'dumptrace_id','nstep','t']
    for v in variables:
        if v in src.keys():
            if v not in dest.keys():
                src.copy(v,dest)

def prune_file(infile_name, ids_to_keep):
    outfile_name = "prune_" + infile_name
    hdr,units,data = io.load_tracer(infile_name)
    mask = np.in1d(data['id'], ids_to_keep)
    ntracers = np.sum(mask)
    with h5py.File(outfile_name,'w') as dest:
        with h5py.File(infile_name,'r') as src:
            copy_header(src,dest)
        for k, v in data.items():
            dest.create_dataset(k,data = v[mask])
        dest.create_dataset('ntracers', data = np.array([ntracers]))
    print("\t{} ... done".format(outfile_name))

def get_ids_to_keep(filename, conditions):
    hdr,units,data = io.load_tracer(filename)
    mask = conditions(data)
    return data['id'][mask]
    
