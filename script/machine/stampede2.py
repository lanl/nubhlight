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

import os

# module load phdf5
# module load gsl

def matches_host():
  host = os.uname()[1]
  return '.stampede2.tacc.utexas.edu' in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = '-O3 -fPIC -Wall -Werror -qopenmp'
  host['DEBUG_FLAGS']    = '-O0 -g -fPIC -Wall -Werror -qopenmp'
  host['HDF5_DIR']       = '/opt/apps/intel17/impi17_0/phdf5/1.8.16/x86_64/'
  host['GSL_DIR']        = '/opt/apps/intel17/gsl/2.3'

  return host

