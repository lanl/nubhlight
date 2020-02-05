#!/usr/bin/env python

"""accumulate_and_prune_tracers.py
Author: Jonah Miller (jonahm@lanl.gov)

Combines accumulate_tracers and prune_tracers
to generate two tracer files and a pruned directory:

tracers_accumulated.td     -- Accumulated as they pass through some radius

prune                      -- Directory containing output
                              files for accumulated tracers
                              and thus full tracer histories

tracers_accumulated_nse.td -- Tracers from tracers_accumulated.td
                              but they are reccoreded when they
                              drop below a set temperature.
"""

from __future__ import print_function, division
import sys; sys.dont_write_bytecode = True
from accumulate_tracers import accumulate_directory
from prune_tracers import prune_directory,_cond_survived
from argparse import ArgumentParser

def conds(data):
    return _cond_survived(data)

parser = ArgumentParser(
    description='Accumulate and Prune Tracers'
)
parser.add_argument('base_dir',type=str,
                    help=("Directory containing original tracer files "
                          +"to accumulate/prune."))
parser.add_argument('accumulated',type=str,
                    help='Name of file to save accumulated tracers to')
parser.add_argument('R',type=float,
                    help='Radius at which to accumulate')
parser.add_argument('prune',type=str,
                    help='Name of new directory to save pruned tracers')
parser.add_argument('nse',type=str,
                    help='Name of file to save tracers out of NSE to')
parser.add_argument('T',type=float,
                    help=('Temperature to assume tracers '
                          +'drop out of NSE. In GK.'))
parser.add_argument('-n','--nprocs',type=int,default=None,
                    help='Number of parallel procs to use')

if __name__ == "__main__":
    args = parser.parse_args()

    print("First accumulate")
    accumulate_directory(args.base_dir,
                         args.R,
                         None,
                         args.accumulated,
                         args.nprocs)
    print("\n\n\n")
    print("Prune")
    prune_directory(args.base_dir,
                    args.accumulated,
                    conds,
                    args.prune,
                    False,
                    args.nprocs)
    print("\n\n\n")
    print("Second accumulate")
    accumulate_directory(args.prune,
                         args.R,
                         args.T,
                         args.nse,
                         args.nprocs)
