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

def matches_host():                                                              
  host = os.uname()[1]                                                           
  return host == 'bh27'

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = '/home/bryan10/software/hdf5/bin/h5pcc'
  host['COMPILER_FLAGS'] = '-O3 -Werror -fopenmp'
  host['DEBUG_FLAGS']    = '-g -O0 -Wall -fopenmp'
  host['GSL_DIR']        = '/home/bryan10/software/gsl'

  return host

