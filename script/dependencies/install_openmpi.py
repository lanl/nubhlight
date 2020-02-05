# Convenience script for installing gsl to ~/software/gsl/

import os
import sys
import subprocess as sp

TMP_DIR = 'TMP/'
class color:
  BOLD = '\033[1m'
  WARNING = '\033[1;31m'
  NORMAL  = '\033[0m'

COMPILER = 'cc'
if len(sys.argv) == 2:
  COMPILER = sys.argv[1]
print COMPILER

if not os.path.exists(TMP_DIR):
  os.makedirs(TMP_DIR)
else:
  print(color.WARNING + "  TMP_DIR " + TMP_DIR + " EXISTS." + color.NORMAL)
  print(color.WARNING + "  TO CONTINUE COULD CAUSE LOSS OF DATA!" + color.NORMAL)
  sys.exit()

os.chdir(TMP_DIR)

sp.call(['wget', 'https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.gz'])
sp.call(['tar', '-xvf', 'openmpi-2.1.0.tar.gz'])
sp.call(['rm', 'openmpi-2.1.0.tar.gz'])

OPENMPI_DIR = os.listdir('.')[0]
os.chdir(OPENMPI_DIR)

HOME = os.path.expanduser("~")

sp.call(['./configure', 'CC=' + COMPILER, '--prefix=' + HOME + '/software/openmpi/'])
sp.call(['make', 'all', 'install'])

os.chdir('..')
os.chdir('..')

sp.call(['rm', '-rf', TMP_DIR])

