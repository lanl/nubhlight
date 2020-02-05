################################################################################
#                                                                              #
# BASE MODULE FOR PYTHON SCRIPTING                                             # 
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True

PATHS = {}
PATHS['BASE'] = os.path.join(os.path.abspath(__file__).rsplit('/', 2)[0], '')
PATHS['CORE'] = os.path.join(PATHS['BASE'], 'core')
PATHS['SCRIPT'] = os.path.join(PATHS['BASE'], 'script')
PATHS['ANALYSIS'] = os.path.join(PATHS['BASE'], 'script', 'analysis')
PATHS['MACHINE'] = os.path.join(PATHS['BASE'], 'script', 'machine')
PATHS['PROB'] = os.getcwd()
sys.path.insert(0, PATHS['SCRIPT'])
sys.path.insert(0, PATHS['ANALYSIS'])
sys.path.insert(0, PATHS['MACHINE'])
import units
cgs = units.get_cgs()
import util
import config
import hdf5_to_dict as io
from config import report_var

def build(PROBLEM):
  config.build(PROBLEM, PATHS)
  #config.build(PROBLEM, PATHS, COMPILE_PARAMS, RUNTIME_PARAMS)
  #config.copy_source_files()
  #config.collect_src()
  print(os.getcwd().split('/')[-1])

def get_host():
  return config.get_host(PATHS)

def bcall(call_string, nprocs:int = 1, **kwargs):
  if not call_string:
    raise ValueError("empty call")
  if type(nprocs) is not int or nprocs < 1:
    raise ValueError("positive integer processors required")

  from subprocess import call

  host = get_host()
  if nprocs > 1:
    call_string = (host.get('MPI_EXECUTABLE','mpirun').split(' ')
                   + ['-np',str(nprocs)]
                   + call_string)
  elif 'EXECUTABLE' in host:
      call_string = host['EXECUTABLE'].split(' ') + call_string

  util.gentle_warn(" ".join(call_string))
  call(call_string,**kwargs)

#check_params_set(COMPILE_PARAMS)
