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

# module purge
# module load gcc
# module load openmpi
# module load hdf5-parallel
# module load python/3.5-anaconda-4.1.1

flags_base = '-Wall -Werror -fdiagnostics-color -fopenmp'
fcflags = '-lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh'
fflags_base = '-fdiagnostics-color -fopenmp -cpp'

def matches_host():
  host = os.uname()[1]
  my_host = 'gr-fe'
  return my_host in host

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  #host['COMPILER']       = '/usr/local/hdf5-parallel/bin/h5pcc'
  host['COMPILER_FLAGS'] = flags_base + ' ' + fcflags + ' ' + '-O2 -march=native'
  host['DEBUG_FLAGS']    = flags_base + ' ' + fcflags + ' ' + '-g -O0'
  # or system GSL if it exists
  host['GSL_DIR']        = os.path.join(os.environ['HOME'],'local-gnu-openmpi')
  host['FORTRAN_COMP']   = 'h5pfc'
  #host['FORTRAN_COMP']   = '/usr/local/hdf5-parallel/bin/h5pfc'
  host['FCFLAGS']        = fflags_base + ' ' + '-O2'
  host['FDEBUG_FLAGS']   = fflags_base + ' ' + '-g -O0'
  host['FORTLINK']       = '-lgfortran -lhdf5_fortran'
  host['FORTLIB']        = ''
  host['EXECUTABLE']     = 'mpirun'

  return host

