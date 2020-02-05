from pylab import *
import plot as bplt
import sys
import hdf5_to_dict as io

if len(sys.argv) != 2:
  print 'ERROR format is'
  print '  python snap_diag.py [folder]'

folder = sys.argv[1]

diag = io.load_diag(folder, hdr = None)

ax = subplot(4,2,1)
ax.plot(diag['t'], -array(diag['Mdot'])/diag['hdr']['MdotEdd'])
ax.set_yscale('log')
ax.set_ylabel('mdot')
ax.set_xticks([])

ax = subplot(4,2,2)
ax.plot(diag['t'], array(diag['Lum'])/diag['hdr']['LEdd'])
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_ylabel('L/LEdd')
ax.set_xticks([])

ax = subplot(4,2,3)
ax.plot(diag['t'], diag['tune_emiss'])
ax.set_yscale('log')
ax.set_ylabel('tune_emiss')
ax.set_xticks([])

ax = subplot(4,2,4)
ax.plot(diag['t'], diag['step_made_all'])
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_ylabel('made')
ax.set_xticks([])

ax = subplot(4,2,5)
ax.plot(diag['t'], diag['tune_scatt'])
ax.set_yscale('log')
ax.set_ylabel('tune_scatt')
ax.set_xticks([])

ax = subplot(4,2,6)
ax.plot(diag['t'], diag['step_scatt_all'])
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_ylabel('scatt')
ax.set_xticks([])

ax = subplot(4,2,7)
ax.plot(diag['t'], diag['egas'])
ax.plot(diag['t'], -diag['erad'])
ax.set_yscale('log')
ax.set_ylabel('egas/erad')

ax = subplot(4,2,8)
ax.plot(diag['t'], diag['step_tot_all'])
ax.set_yscale('log')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.set_ylabel('tot')

show()

