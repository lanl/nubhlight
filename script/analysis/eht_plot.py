import sys; sys.dont_write_bytecode=True
sys.path.insert(0, '../')
import numpy as np
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import util
import pickle
import hdf5_to_dict as io

FIGX = 14
FIGY = 10
SIZE = 40

if len(sys.argv) != 3:
  util.warn('Format: python eht_plot.py [averages] [output]')
  sys.exit()

favg = sys.argv[1]
fout = sys.argv[2]

avg = pickle.load(open(favg, 'rb'))

fig = plt.figure(figsize=(FIGX, FIGY))

ax = plt.subplot(2,3,1)
ax.plot(avg['r'], avg['rho_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<rho>')
ax.set_ylim([1.e-2, 1.e0])

ax = plt.subplot(2,3,2)
ax.plot(avg['r'], avg['Pg_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<Pg>')
ax.set_ylim([1.e-6, 1.e-2])

ax = plt.subplot(2,3,3)
ax.plot(avg['r'], avg['B_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<|B|>')
ax.set_ylim([1.e-4, 1.e-1])

ax = plt.subplot(2,3,4)
ax.plot(avg['r'], avg['uphi_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<u^phi>')
ax.set_ylim([1.e-3, 1.e1])

ax = plt.subplot(2,3,5)
ax.plot(avg['r'], avg['Ptot_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<Ptot>')
ax.set_ylim([1.e-6, 1.e-2])

ax = plt.subplot(2,3,6)
ax.plot(avg['r'], avg['betainv_r'], color='k', linewidth=2)
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('<beta^-1>')
ax.set_ylim([1.e-2, 1.e1])

plt.savefig(fout + '_ravgs.png')
plt.close(fig)

# SADW
fig = plt.figure(figsize=(FIGX, FIGY))
ax = plt.subplot(2,3,1)
ax.plot(avg['r'], avg['rho_SADW'], color='k', linewidth=2)
ax.set_yscale('log')
plt.savefig(fout + '_sadw.png')
plt.close(fig)


fig = plt.figure(figsize=(FIGX, FIGY))

sx = ''
sx = '_d'

ax = plt.subplot(5,1,1)
ax.plot(avg['t'+sx], np.fabs(avg['Mdot'+sx]), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Mdot|')

ax = plt.subplot(5,1,2)
ax.plot(avg['t'+sx], avg['Phi'+sx], color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('Phi')

ax = plt.subplot(5,1,3)
ax.plot(avg['t'+sx], np.fabs(avg['Ldot'+sx]), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Ldot|')

ax = plt.subplot(5,1,4)
ax.plot(avg['t'+sx], np.fabs(avg['Edot'+sx] - avg['Mdot'+sx]), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Edot - Mdot|')

ax = plt.subplot(5,1,5)
ax.plot(avg['t'+sx], avg['Lum'+sx], color='k')
ax.set_xlim([0,1e4])
ax.set_xlabel('t/M')
ax.set_ylabel('Lum')

plt.savefig(fout + '_fluxes.png')
plt.close(fig)


fig = plt.figure(figsize=(FIGX, FIGY))

ax = plt.subplot(5,1,1)
ax.plot(avg['t'+sx], np.fabs(avg['Mdot'+sx]), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Mdot|')

print np.fabs(avg['Mdot'+sx]).max()

ax = plt.subplot(5,1,2)
ax.plot(avg['t'+sx], avg['Phi'+sx]/np.sqrt(np.fabs(avg['Mdot'+sx])), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('Phi/sqrt(|Mdot|)')

ax = plt.subplot(5,1,3)
ax.plot(avg['t'+sx], np.fabs(avg['Ldot'+sx])/(np.fabs(avg['Mdot'+sx])), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Ldot|/|Mdot|')

ax = plt.subplot(5,1,4)
ax.plot(avg['t'+sx], np.fabs(avg['Edot'+sx] - avg['Mdot'+sx])/(np.fabs(avg['Mdot'+sx])), color='k')
ax.set_xlim([0,1e4])
ax.set_xticklabels([])
ax.set_ylabel('|Edot - Mdot|/|Mdot|')

ax = plt.subplot(5,1,5)
ax.plot(avg['t'+sx], avg['Lum'+sx]/(np.fabs(avg['Mdot'+sx])), color='k')
ax.set_xlim([0,1e4])
ax.set_xlabel('t/M')
ax.set_ylabel('Lum/|Mdot|')

plt.savefig(fout + '_normfluxes.png')
plt.close(fig)

