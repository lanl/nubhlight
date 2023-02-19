#!/usr/bin/env python

#"""
# Author: Kelsey Lund but based on Jonah Miller's analysis scripts
#
# Accumulate tracer data into a single file with values that are used as inputs for PRISM calculations.
# For now, this should only be used as a post-processing step in the sense that it assumes PRISM
# trajectory files already exist. 
# To do: Build this structure into the cleanup_trace function in tracers_to_PRISM so that one of the 
# outputs from that whole procedure is the accumulated file.
# """

import glob
import pickle
import traceback
import numpy as np
import hdf5_to_dict as io
from multiprocessing import Pool
from argparse import ArgumentParser
import units; units = units.get_cgs()

def find_t_closest(input,time):
	tlst = []
	for t in range(len(input)):
		tlst.append(abs(input[t] - time))
	return tlst.index(min(tlst))

def find_PRISM_start_t(trace_ID,dirt):
	# Pulls out starting time of PRISM trajectory file (in seconds, beware!)
	tstart = np.loadtxt(dirt + '/trace_%08d'%trace_ID + '/trajectory.dat')[0][0]
	return tstart 
    
def get_IDs(dirt):
	# dirt should be something like 'analysis/traces/' WITH trailing backslash(for now)
	# dirt should be traces directory that has (usually) 100 subdirectories
	IDs = []
	subdirs = sorted(glob.glob(dirt + '*/'))
	for subdir in subdirs:
		tr_ids = glob.glob(subdir + 'trace_*.td')
		for ID in tr_ids:
			IDs.append(ID)
	return IDs

def loop(inputs):#tr_id,dirt):
	try:
		tr_id,pr_dirt = inputs   # pr_dirt should be where PRISM trajectory files are- something like '/analysis/done_DND/'
		trace = io.TracerData.fromfile(tr_id)
		trace_times = trace['time'] * trace.units['T_unit']
		vals = {}
	
		start_t = find_PRISM_start_t(tr_id, prdirt)
		index = find_t_closest(trace_times,start_t)
	
		for key in trace.keys():
			vals[key] = [trace[key][index]]
		return vals
	except:
		traceback.print_exc()
#---------------------------------------------------------------------------------------#

parser = ArgumentParser()
parser.add_argument('trace_dir', type = str, help = ('Directory where subdirectories containing tracer files by IDs are located.'))
parser.add_argument('prism_dir', type=str, help = ('Directory where PRISM trajectory files are kept.'))
parser.add_argument('accumulated',type=str, help = ('Name of file to save accumulated tracers to.'))
parser.add_argument('-n','--nprocs',type=int, default = None, help=('Number of parallel procs to use.'))

if __name__ == "__main__":
	try:
		args = parser.parse_args()
		
		fnams = get_IDs(args.trace_dir) 
		p = Pool(args.nprocs)
		
		inputs = []
		vals = {}
		for f in fnams:
			inputs.append((f,args.prism_dir))
		vals_out = p.map(loop,inputs)
		for key in vals_out.keys():
			if key not in vals.keys():
				vals[key] = vals[key]
			else:
				vals[key] = np.concatenate((vals[key],vals_out[key]))
				
		with open(args.accumulated,'wb') as f:
			pickle.dump(vals, f)
	except:
		traceback.print_exc()

#trace_dir = '/users/klund/scratch4/nubhlight/torus_gw/analysis/traces/'
#prism_dir = '/users/klund/scratch4/nubhlight/torus_gw/analysis/done_DND/'
#nprocs = 10


