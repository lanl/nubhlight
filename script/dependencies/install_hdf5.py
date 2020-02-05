################################################################################
#                                                                              #
# AUTOMATICALLY INSTALLS GSL TO bhlight/dependencies/gsl                       #
#                                                                              #
################################################################################

import os
import sys
import subprocess as sp

TMP_DIR = 'TMP/'
class color:
  BOLD = '\033[1m'
  WARNING = '\033[1;31m'
  NORMAL  = '\033[0m'

if not os.path.exists(TMP_DIR):
  os.makedirs(TMP_DIR)
else:
  print(color.WARNING + "  TMP_DIR " + TMP_DIR + " EXISTS." + color.NORMAL)
  print(color.WARNING + "  TO CONTINUE COULD CAUSE LOSS OF DATA!" + color.NORMAL)
  sys.exit()

os.chdir(TMP_DIR)

sp.call(['wget', 'https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.2/src/hdf5-1.10.2.tar.gz'])
sp.call(['tar', '-xvzf', 'hdf5-1.10.2.tar.gz'])
sp.call(['rm', 'hdf5-1.10.2.tar.gz'])

HDF5_DIR = os.listdir('.')[0]
os.chdir(HDF5_DIR)

HOME = os.path.expanduser("~")

os.environ['CC'] = HOME + '/software/openmpi/bin/mpicc'
sp.call(['./configure', '--enable-parallel', '--prefix=' + HOME + '/software/hdf5'])

#sp.call(['CC=' + HOME + '/software/openmpi/bin/mpicc', './configure', '--enable-parallel', 
#         '--prefix=' + HOME + '/software/hdf5/'])
sp.call(['make'])
sp.call(['make', 'install'])

os.chdir('..')
os.chdir('..')

sp.call(['rm', '-rf', TMP_DIR])

