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
# This is an example machine file that shows how you might specify settings    #
# for a machine running Debian, assuming you have installed gsl, mpi, and hdf5 #
# parallel via the package manager.                                            #
################################################################################

import util
import sys
import os

hostnames = ['cnls-cf2', 'cnls-cf4', 'cnls-cf5']
flags_base = '-Wall -Werror -fopenmp'
fcflags = '-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh'
fflags_base = '-fopenmp'

local = '/netscratch/jonahm/local-gcc-8.2.0'

def matches_host():
  host = os.uname()[1]
  for h in hostnames:
    if h in host:
      return True
  return False

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = os.path.join(local,'bin','h5pcc')
  host['COMPILER_FLAGS'] = flags_base + ' ' + fcflags + ' ' + '-O2 -march=native'
  host['DEBUG_FLAGS']    = flags_base + ' ' + fcflags + ' ' + '-g -O0'
  host['GSL_DIR']        = local
  host['FORTRAN_COMP']   = os.path.join(local,'bin','h5pfc')
  host['FCFLAGS']        = fflags_base + ' ' + '-cpp -O2'
  host['FORTLINK']       = '-lgfortran -lhdf5_fortran'
  host['FORTLIB']        = ''
  host['EXECUTABLE']     = os.path.join(local,'bin','mpirun')

  return host

