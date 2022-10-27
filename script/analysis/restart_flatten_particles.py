#!/usr/bin/env python
######################################################################
# Flattens particle dataset in restart file to new format
######################################################################

from __future__ import print_function,division
import h5py
import numpy as np
from argparse import ArgumentParser

def get_dsets_and_sizes(h5file):
    """gets photon datasets and sizes"""
    keys = []
    sizes = []
    for k,v in h5file.items():
        if 'photons_' in k:
            keys.append(k)
            sizes.append(v.shape[0])
    return keys,sizes

def update_file(h5file):
    "Updates hdf5 file object"
    keys,sizes = get_dsets_and_sizes(h5file)
    if len(keys) <= 0:
        raise ValueError("This file does not need conversion.")
    sizes = np.array(sizes)
    offsets = np.zeros(len(sizes)+1)
    offsets[1:] = np.cumsum(sizes)
    offsets = offsets[:-1]
    h5file.create_dataset('particle_counts',
                          data=sizes)
    h5file.create_dataset('particle_offsets',
                          data=offsets)
    data = list([None for k in keys])
    for i,k in enumerate(keys):
        print("\tReading dataset",
              k,"with size",
              h5file[k].shape)
        data[i] = h5file[k][:]
    data = np.hstack(data)
    print("\tWriting superphoton dataset.")
    h5file.create_dataset('superphotons',
                          data=data)

if __name__ == "__main__":
    parser = ArgumentParser("Flatten a particle dataset in a restart file to new restart format")
    parser.add_argument('restart', type=str,
                        help="Restart file to flatten")
    args = parser.parse_args()
    print("Updating file",args.restart)
    with h5py.File(args.restart, 'r+') as f:
        update_file(f)
    print("Done.")
