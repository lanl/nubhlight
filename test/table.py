################################################################################
#                                                                              #
# TABULATED EOS READER                                                         #
#                                                                              #
################################################################################

from __future__ import print_function
import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import util
from bhlight import bcall

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'table'
AUTO = '-auto' in sys.argv
gam = 1.4

os.chdir('../prob/' + PROBLEM)
# COMPILE CODE
args = [sys.executable, 'build.py', '-dir', TMP_DIR]
call(args)
call(['mv', 'sc_eos_gamma_{}.h5'.format(str(gam).replace('.','p')), TMP_DIR])
os.chdir('../../test/')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
bcall(['./bhlight', '-p', 'param_template.dat'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfile = os.path.join(TMP_DIR,'dumps','table_errors.txt')
errs = np.loadtxt(dfile)
if AUTO:
    data = {}
    x = np.linspace(0,1,len(errs))
    data['CODE'] = [x, errs]
    data['SOL'] = [x, np.zeros_like(errs)]
    import pickle
    pickle.dump(data, open('data.p','wb'))
else:
    print("Max sum error = {:.5e}".format(errs.max()))

# CLEAN UP
util.safe_remove(TMP_DIR)
