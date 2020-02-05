# Convenience script for installing gsl to ~/software/gsl/

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

sp.call(['wget', 'ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz'])
sp.call(['tar', '-xvzf', 'gsl-latest.tar.gz'])
os.remove('gsl-latest.tar.gz')

print os.walk('./')
print os.listdir('./')
GSL_DIR = os.listdir('.')[0]
#GSL_DIR = os.walk('./')
print GSL_DIR
os.chdir(GSL_DIR)

HOME = os.path.expanduser("~")
print HOME

sp.call(['./configure', '--prefix=' + HOME + '/software/gsl/'])
sp.call(['make'])
sp.call(['make', 'install'])

os.chdir('..')
os.chdir('..')

sp.call(['rm', '-rf', TMP_DIR])
