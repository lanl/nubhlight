################################################################################
#                                                                              #
#  MACHINE-SPECIFIC FUNCTIONS                                                  #
#                                                                              #
#    OPTIONS:                                                                  #
#      COMPILER   : PATH TO COMPILER EXECUTABLE                                #
#      GSL_DIR    : PATH TO GSL INSTALLATION                                   #
#      MPI_DIR    : PATH TO MPI INSTALLATION                                   #
#      HDF5_DIR   : PATH TO HDF5 INSTALLATION                                  #
#      EXECUTABLE : BINARY WRAPPER USED TO LAUNCH BHLIGHT                      #
#                                                                              #
#    MPI_DIR AND HDF5_DIR ARE NOT REQUIRED IF COMPILER HANDLES HEADERS AND     #
#    LIBRARIES FOR THESE DEPENDENCIES                                          #
#                                                                              #
################################################################################

import util
import sys
import os
import re

# export PYTHONPATH=/users/jonahm/local-cray-mpich/lib/python3.7/site-packages:${PYTHONPATH}
# module load cce
# module load cray-mpich
# module load cray-hdf5-parallel
# module load cray-python/3.7.3.2

flags_base = '-fdiagnostics-color -fopenmp'
fflags_base = '-fdiagnostics-color -fopenmp -cpp'

GSL_NAME='local-cray-mpich'

def matches_host():
  host = os.uname()[1]
  frontend = 'cp-login'
  return frontend in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O2'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  # Change this to your locally installed GSL
  host['GSL_DIR']        = os.path.join(os.environ['HOME'],GSL_NAME)
  host['FORTRAN_COMP']   = 'h5pfc'
  host['FCFLAGS']        = fflags_base + ' ' + '-O2'
  host['FDEBUG_FLAGS']   = fflags_base + ' ' + '-g -O0'
  host['FORTLINK']       = '-lgfortran -lhdf5_fortran'
  host['FORTLIB']        = ''
  host['MEM_MODEL']      = 'small'
  host['EXECUTABLE']     = 'mpiexec -n 1'
  host['MPI_EXECUTABLE'] = 'mpiexec'

  return host

