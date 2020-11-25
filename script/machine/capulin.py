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

# export PYTHONPATH=${HOME}/local-cray-mpich/lib/python3.7/site-packages:${PYTHONPATH}
# module load cce
# module load cray-mpich
# module load cray-hdf5-parallel
# module load cray-python/3.7.3.2

# Also requires gsl somehow.
# Assumes gsl installed in ${HOME}/local-cray-mpich
# but gsl installed via spack also works

# DO NOT use fortran with this machine.
# It will fail.

# Some warnings:
# 1. Cray is a weird beast. Use cc, not the wrappers.
# 2. -O2 and higher cuase problems. Don't use them.
# 3. -fPIC is required to make the cray linker behave. Static linking only. No dynamic linking.
# 4. -fQunused-arguments is a clang nicety and makes some warnings go away.
# 5. Use system libraries when possible. I found it difficult to roll my own, e.g., HDF5.

flags_base = '-fdiagnostics-color -Qunused-arguments -fopenmp -fPIC'
fflags_base = ''

GSL_NAME='local-cray-mpich'

def matches_host():
  host = os.uname()[1]
  frontend = 'cp-login'
  return frontend in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'cc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O1'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  # Change this to your locally installed GSL
  host['GSL_DIR']        = os.path.join(os.environ['HOME'],GSL_NAME)
  host['FORTRAN_COMP']   = ''
  host['FCFLAGS']        = ''
  host['FDEBUG_FLAGS']   = ''
  host['FORTLINK']       = ''
  host['FORTLIB']        = ''
  host['MEM_MODEL']      = 'small'
  host['EXECUTABLE']     = 'mpiexec -n 1'
  host['MPI_EXECUTABLE'] = 'mpiexec'

  return host
