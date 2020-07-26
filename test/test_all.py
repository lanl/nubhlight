################################################################################ 
#                                                                              # 
# RUN ALL TESTS AND CHECK FOR ACCURACY                                         # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
import util
import subprocess as sp
import numpy as np
import glob
import pickle
from scipy.interpolate import interp1d as interp
import time

# SPECIFY EMAIL OPTIONS
TO = ['bryan10@illinois.edu']
CC = []
FROM = 'afd.illinois.testing@gmail.com'
PASSWORD = 'whatpasswordshouldichoose'
SUBJECT = 'BHLIGHT TESTING REPORT'
LOGNAME = 'test_all.txt'

SEND_REPORT = '-email' in sys.argv
TABLE = '-table' in sys.argv
FAST = '-fast' in sys.argv
VERBOSE = '-verbose' in sys.argv

LONG_TEST = '-long' in sys.argv
RAD_TEST = '-rad' in sys.argv
FLUID_TEST = '-fluid' in sys.argv
NU_TEST = '-nu' in sys.argv
LIGHT_TEST = '-light' in sys.argv

ERROR_THRESHOLD = 0.01

# REPORT TIME FOR TESTING

# GET ALL TEST SCRIPTS
M_FAST_TESTS = ['tracers1d.py', 'sod.py','table.py','advection1d.py','binning.py']
M_FLUID_TESTS = ['sod.py', 'table.py', 'advection2d.py', 'advection2d.py -mpi']
M_NU_TESTS = ['yedecay.py', 'yedecay.py -antinu',
              'yedecay.py -mpi', # 'yedecay.py -mpi -antinu',
              'multiscatt.py',
              'tracers1d.py']
M_LIGHT_TESTS = ['binning.py', 'brem.py',
                 'thermalization.py', 'thermalization_mpi.py',
                 'comptonization.py']
SKIP = ['generate_all_plots.py','test_all.py']
if LONG_TEST:
  TESTS = glob.glob('*.py')
  for test in SKIP:
    if test in TESTS:
      TESTS.remove(test)
elif FAST:
  TESTS = M_FAST_TESTS
elif LIGHT_TEST:
  TESTS = M_LIGHT_TESTS
elif RAD_TEST:
  TESTS = M_LIGHT_TESTS + M_NU_TESTS
elif FLUID_TEST:
  TESTS = M_FLUID_TESTS
elif NU_TEST:
  TESTS = M_NU_TESTS
else:
  TESTS = M_FLUID_TESTS + M_LIGHT_TESTS + M_NU_TESTS

print("")                                                                      
print("********************************************************************************")
print("")                                                                      
print("                                AUTOMATED TESTING")                    
print("")                                                                      
print("********************************************************************************")

util.log_output(sys, LOGNAME)

DATE = time.strftime('%Y/%m/%d')
TIME = time.strftime('%H:%M:%S')
MACHINE = os.uname()[1]
popen = sp.Popen(['git', 'show', '-s', '--format=%H'], stdout=sp.PIPE,
                   universal_newlines=True)
for line in iter(popen.stdout.readline, ""):
  HASH = line.lstrip().rstrip()
popen = sp.Popen(['git', 'branch'], stdout=sp.PIPE, universal_newlines=True)
for line in iter(popen.stdout.readline, ""):
  if line[0] == '*': 
    BRANCH = line[2:].rstrip()
print('\n  DATE:    ' + DATE)
print('  TIME:    ' + TIME)
print('  MACHINE: ' + MACHINE)
print('  BRANCH:  ' + BRANCH)
print('  COMMIT:  ' + HASH + '\n')

def name_to_args(namestring):
  """Takes a script name which may contain CLI args
  and splits out the args.

  Assumes string is generically of the form
  '<script name>.py -arg --arg'
  """
  namestring = namestring.rstrip().lstrip()
  if namestring[-3:] == '.py':
    return [namestring]
  args = namestring.split('.py ')
  args = ([args[0].lstrip().rstrip() + '.py']
          + [s.lstrip().rstrip() for s in args[1].split()])
  return args

# USE INTERPOLATION ON A (ANALYTIC SOLUTION) TO COMPARE TO B
def sanitize_array(a):
  a = np.array(a) # ensure a is a numpy array
  if len(a.shape) == 1:
    return a
  if np.prod(a.shape[1:]) > 1:
    raise ValueError(
      "Array should be 1d. Array shape = {}".format(a.shape)
    )
  return a.reshape(a.shape[0])

def L1_norm(xa, ya, xb, yb):
  # special case for 0d arrays
  if len(xa) == len(xb) == len(ya) == len(yb) == 1:
    if np.abs(yb[0]) <= 1e-12:
      return np.fabs(ya[0] - yb[0])
    return np.fabs((ya[0] - yb[0])/yb[0])
  xa,ya,xb,yb = [sanitize_array(a) for a in [xa,ya,xb,yb]]
  if xa[0] > xb[0]:
    xb = xb[1:]
    yb = yb[1:]
  fa = interp(xa, ya)
  norm = 0.
  nodenom = np.max(ya) <= 1e-12
  for n in range(len(xb)):
    num = np.fabs(yb[n] - fa(xb[n]))
    denom = np.fabs((yb[n] + fa(xb[n]))/2.)
    if nodenom:
      norm += num
    else:
      norm += num/denom

  return (norm/n)

FAIL = False
for TEST in TESTS:
  args = name_to_args(TEST)
  TESTNAME = args[0][:-3] if len(args) == 1 else TEST
  print('  ' + util.color.BOLD + TESTNAME + util.color.NORMAL)
  args = [sys.executable] + args + ['-auto']
  if TABLE:
    args += ['-table']
  if FAST:
    args += ['-fast']
  popen = sp.Popen(args,
                   stdout=sp.PIPE,
                   stderr=sp.PIPE,
                   universal_newlines=True)
  for line in iter(popen.stdout.readline, ""):
    if VERBOSE:
      print(line.rstrip())
    if line.lstrip().rstrip() == 'BUILD SUCCESSFUL':
      print('    BUILD SUCCESSFUL')
  print('    RUN FINISHED')
  popen.wait()

  if not os.path.isfile('data.p'):
    raise RuntimeError("Test did not succesfully complete.")
  
  with open('data.p', 'rb') as f:
    data = pickle.load(f)
    xa = data['SOL'][0]
    ya = data['SOL'][1]
    xb = data['CODE'][0]
    yb = data['CODE'][1]
    if 'THRESHOLD' in data.keys():
      error_threshold = data['THRESHOLD']
    else:
      error_threshold = ERROR_THRESHOLD

  norm = L1_norm(xa, ya, xb, yb)

  print('    ERROR: %.2g %%' % (100*norm))
  if norm < error_threshold:
    print(util.color.BOLD + '    PASS' + util.color.NORMAL + '\n')
  else:
    print(util.color.WARNING + '    FAIL' + util.color.NORMAL + '\n')
    FAIL = True

  sp.call(['rm', 'data.p'])
 
if not SEND_REPORT:
  sp.call(['rm', LOGNAME])
  if FAIL:
    raise RuntimeError("Tests failed!")
  else:
    print("All tests passed!")
  sys.exit()

import smtplib
if FAIL:
  SUBJECT += ' - FAIL'
else:
  SUBJECT += ' - PASS'
MESSAGE = ''
MFILE = open(LOGNAME, 'rb')
for line in MFILE:
  MESSAGE += line
MFILE.close()
EMAIL = ('From: %s\r\n' % FROM
         + 'To: %s\r\n' % ','.join(TO)
         + 'CC: %s\r\n' % ','.join(CC)
         + 'Subject: %s\r\n' % SUBJECT
         + '\r\n'
         + MESSAGE)
ADDRS = TO + CC
srvr = smtplib.SMTP('smtp.gmail.com', 587)
srvr.ehlo()
srvr.starttls()
srvr.ehlo()
srvr.login(FROM, PASSWORD)
srvr.sendmail(FROM, ADDRS, EMAIL)
srvr.close()
sp.call(['rm', LOGNAME])
