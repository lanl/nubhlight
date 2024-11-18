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

# module load gcc/4.9.3
# module load openmpi/1.10.0-gcc_4.9.3

# gsl must be installed by hand to
# ~/local/gcc/gsl

# hdf5-parallel must be installed by hand to
# ~/local/gcc/hdf5-parallel

flags_base = '-Wall -Werror -fdiagnostics-color -fopenmp'

def matches_host():
  host = os.uname()[1]
  return 'darwin-fe' in host# or 'cn' in host

def get_options():
  host = {}

  local_root = os.path.join(os.environ['HOME'],'local','skylake-gold')
  host['NAME']           = os.uname()[1]
  host['COMPILER']       = os.path.join(local_root,'bin','h5pcc')
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O3 -march=native'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  host['GSL_DIR']        = os.path.join(local_root,'lib')

  return host

