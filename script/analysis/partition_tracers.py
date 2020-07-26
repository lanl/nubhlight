#!/usr/bin/env python

################################################################################
# Author: Jonah Miller (jonahm@lanl.gov)
#
# PURPOSE:
# Partition a single tracer file into multiple subfiles to ease the process
# of converting tracers into traces.
################################################################################

from __future__ import print_function, division
import numpy as np
from numpy import random
from argparse import ArgumentParser
import hdf5_to_dict as io

def split_ids(ids,n):
    if n >= 2:
        np.random.shuffle(ids)
        id_lists = np.array_split(ids,n)
    else:
        id_lists = [ids]
    return id_lists

def partition(tracers,n):
    ids = tracers.ids()
    id_lists = split_ids(ids,n)
    partitions = [tracers.filter(np.in1d(tracers['id'],l))\
                  for l in id_lists]
    return partitions

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Partition tracer file into N new files"
    )
    parser.add_argument('infile', type=str,
                        help='File to partition')
    parser.add_argument('N',type=int,
                        help='Number of subfiles to partition into')

    args = parser.parse_args()
    tracers = io.TracerData.fromfile(args.infile)
    tpart = partition(tracers,args.N)
    for i,p in enumerate(tpart):
        outfile = args.infile + '{:03d}.td'.format(i)
        print("Saving file: {}".format(outfile))
        p.save(outfile)
