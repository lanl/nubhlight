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

flags_base = '-Wall -Werror=unknown-pragmas -Wstringop-overflow -fdiagnostics-color -fopenmp'

def matches_host():
  host = os.uname()[1]
  return 'sudi.local' in host # change to your username

def get_options():
  host = {}

  local_root = os.path.join(os.environ['HOME'],'local','gcc')
  host['NAME']           = os.uname()[1]
  host['COMPILER']       = os.path.join('/opt/homebrew/opt/gcc/','bin','gcc-13')
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O2'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  host['GSL_DIR']        = os.path.join(local_root,'gsl')
  host['GSL_DIR']        = os.path.join('/opt/homebrew/opt/','gsl')
  host['MPI_DIR']        = os.path.join('/opt/homebrew/opt/','open-mpi')
  host['HDF5_DIR']       = os.path.join('/opt/homebrew/opt/','hdf5-mpi')
  host['EXTRA_INCLUDES'] = "-I/opt/homebrew/Cellar/gcc/13.1.0/lib/gcc/current/gcc/aarch64-apple-darwin22/13/include/"
  host['MEM_MODEL']      = False
  host['USE_RPATH']      = False
  return host
 
