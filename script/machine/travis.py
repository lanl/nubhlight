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

################################################################################
# This machine file is for continuous integration servers running on CI        #
################################################################################

import util
import sys
import os

flags_base = '-Wall -Werror -fopenmp'
fcflags = ''
fflags_base = '-fopenmp'

def matches_host():
  host = os.uname()[1]
  keywords = ['travis','job']
  kinh = [k in host for k in keywords]
  return any(kinh) and all(kinh)

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = '/usr/local/hdf5-parallel/bin/h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + fcflags + ' ' + '-O2 -march=native'
  host['DEBUG_FLAGS']    = flags_base + ' ' + fcflags + ' ' + '-g -O0'
  host['GSL_DIR']        = ''
  host['FORTRAN_COMP']   = '/usr/local/hdf5-parallel/bin/h5pfc'
  host['FCFLAGS']        = fflags_base + ' ' + '-cpp -O2'
  host['FORTLINK']       = '-lgfortran -lhdf5_fortran'
  host['FORTLIB']        = ''
  host['EXECUTABLE']     = 'mpirun'

  return host

