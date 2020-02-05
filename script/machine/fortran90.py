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

#flags_base = '-Wall -Werror -fdiagnostics-color -fopenmp'
flags_base = '-Wall -fdiagnostics-color -fopenmp'

def matches_host():
  host = os.uname()[1]
  return host == 'fortran90'

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + '-O3'
  host['DEBUG_FLAGS']    = flags_base + ' ' + '-g -O0'
  host['GSL_DIR']        = '/home/brryan/software/gsl'

  return host

