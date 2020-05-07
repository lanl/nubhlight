################################################################################ 
#                                                                              # 
# RUN ALL TESTS AND GENERATE RELEVANT PLOTS                                    # 
#                                                                              # 
################################################################################

import os, glob, time
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
import util
import subprocess as sp

LOGNAME = 'generate_all_plots.txt'
SKIP = ['generate_all_plots.py','test_all.py']
TESTS = glob.glob('*.py')
TESTS.remove('test_all.py')
for test in SKIP:
    if test in TESTS:
        TESTS.remove(test)

TABLE = '-table' in sys.argv

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
print '\n  DATE:    ' + DATE
print '  TIME:    ' + TIME
print '  MACHINE: ' + MACHINE
print '  BRANCH:  ' + BRANCH
print '  COMMIT:  ' + HASH + '\n'

for TEST in TESTS:
    print '  ' + util.color.BOLD + TEST + util.color.NORMAL
    args = [sys.executable, TEST]
    if TABLE:
        args += ['-table']
    popen = sp.Popen(args, stdout=sp.PIPE, stderr=sp.PIPE,
                     universal_newlines=True)
    for line in iter(popen.stdout.readline, ""):
        if line.lstrip().rstrip() == 'BUILD SUCCESSFUL':
            print '    BUILD SUCCESSFUL'  
    print '    RUN FINISHED' 
