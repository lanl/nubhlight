################################################################################
#                                                                              #
#  weighted radial profiles                                                    #
#                                                                              #
################################################################################

import numpy as np
from scipy import integrate
import matplotlib as mpl
#matplotlib.use('agg')
import matplotlib.pyplot as plt

def integrate_shell(hdr, geom, var):
    "Return (int var sqrt(g) dOmega)(r)"
    if hdr['N3'] > 1 and len(geom['gdet'].shape) < 3:
        integrand = var*geom['gdet'][:,:,np.newaxis]
    elif hdr['N3'] == 1:
        N1,N2,N3 = hdr['N1'],hdr['N2'],hdr['N3']
        integrand = var.reshape((N1,N2,N3))*geom['gdet'].reshape((N1,N2,N3))
    else:
        integrand = var*geom['gdet']
    if hdr['N3'] > 1: # assumes X1,X2,X3 are uniform. (They are!)
        intdphi = integrate.simps(integrand,dx = hdr['dx'][3], axis = 2)
    else:
        intdphi = 2*np.pi*integrand[:,:,0]
    intdomega = integrate.simps(intdphi,dx = hdr['dx'][2], axis = 1)
    return intdomega

def get_density_weighted_profile(dump, q):
    "Return ((int q rho sqrt(g) dOmega)/(int rho sqrt(g) dOmega))(r)"
    hdr = dump['hdr']
    geom = dump['geom']
    numerator = integrate_shell(hdr,geom,q*dump['RHO'])
    denomenator = integrate_shell(hdr,geom,dump['RHO'])
    return numerator/denomenator

def get_spherical_average(dump,q):
    "Return the spherical average, not weighted by density"
    hdr = dump['hdr']
    geom = dump['geom']
    numerator = integrate_shell(hdr,geom,q)
    denomenator = integrate_shell(hdr,geom,np.ones_like(q))
    return numerator/denomenator

def plot_density_weighted_profile(ax, geom, var, dump,
                                  lw=3, ylabel = None, label=None):
    r = geom['r'][:,0,0]
    profile = get_density_weighted_profile(dump, var)
    line = ax.plot(r,profile,lw=lw, label=label)
    ax.set_xlabel(r'$r/M$')
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    return line

if __name__ == "__main__":
    import argparse
    import hdf5_to_dict as io
    parser = argparse.ArgumentParser(
        description='Plot a radial average quantity in bhlight.')
    parser.add_argument('dump',
                        type=str,
                        help='Path to dump you want to use for the plot.')
    parser.add_argument('-q','--qname',
                        dest='qname',
                        type=str,
                        help='Name of qantity to plot.',
                        default='Theta')
    parser.add_argument('--ylabel',dest='ylabel',
                        type=str,
                        default=None,
                        help='Label for y axis')
    parser.add_argument('--rmin',
                        dest='rmin',
                        type=float,
                        help='Min radius',
                        default = None)
    parser.add_argument('--rmax'
                        ,dest='rmax',
                        type=float,
                        help='Max radius',
                        default=None)
    args = parser.parse_args()
        
    mpl.rcParams.update({'font.size':18})
    dump = io.load_dump(args.dump)
    geom = dump['geom']
    hdr = dump['hdr']
    varname = args.qname
    var = dump[varname]
    fig,ax = plt.subplots()
    if args.ylabel is not None:
        ylabel = args.ylabel
    else:
        ylabel = varname
    plot_density_weighted_profile(ax,geom,var,dump,ylabel=ylabel)
    if args.rmin is None:
        rmin = geom['r'].min()
    else:
        rmin = args.rmin
    if args.rmax is None:
        rmax = geom['r'].max()
    else:
        rmax = args.rmax
    ax.set_xlim(rmin,rmax)
    plt.savefig('{}.pdf'.format(varname),bbox_inches='tight')
    plt.savefig('{}.png'.format(varname),bbox_inches='tight')

