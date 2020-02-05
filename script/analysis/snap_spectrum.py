#import matplotlib
#import matplotlib.pyplot as plt
from pylab import *
import plot as bplt
import sys
import hdf5_to_dict as io
#import numpy as np

if len(sys.argv) != 2:
  print 'ERROR Format is python snap_spectrum.py [filename]'
  sys.exit()

fnam = sys.argv[1]
dfolder = fnam.rsplit('/', 1)[0]
#dfolder = sys.argv[1]
#nd      = int(sys.argv[2])

#dumps = io.get_dumps_reduced(dfolder)
hdr = io.load_hdr(fnam)
dump = io.load_dump(fnam)

nu = hdr['nu']

ax = subplot(1,1,1)
nuLnu = dump['nuLnu'].sum(axis=0).sum(axis=0).sum(axis=0)
#nuLnu = dump['nuLnu'].sum(axis=1).sum(axis=1)
ax.step(nu, nuLnu, where='mid', color='k')
ax.set_xscale('log'); ax.set_yscale('log')
maxL = nuLnu.max()
ax.set_ylim([1.e-10*maxL, 1.e1*maxL])
ax.set_xlabel('nu (Hz)'); ax.set_ylabel('nuLnu (erg s^-1)')

show()

