#!/usr/bin/env python

################################################################################
# Author: Jonah Miller (jonahm@lanl.gov)
#
# PURPOSE:
# To provide plotting and analysis routines for tracer particles for nubhlight.
# Performs simple binning. Saves plots, etc.
################################################################################

from __future__ import print_function, division
import numpy as np
import matplotlib as mpl
#matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from units import cgs

def get_r(tracers):
    "Returns spherical radius"
    return np.sqrt(tracers['Xcart'][:,0]**2
                   +tracers['Xcart'][:,1]**2
                   +tracers['Xcart'][:,2]**2)

def get_rcyl_z(tracers):
    "Returns cylindrical radius and Cartesian z"
    rcyl = np.sqrt(tracers['Xcart'][:,0]**2 + tracers['Xcart'][:,1]**2)
    z = tracers['Xcart'][:,2]
    return rcyl,z

def get_theta(tracers):
    "Returns theta_bl - pi/2"
    rcyl,z = get_rcyl_z(tracers)
    theta = np.arctan2(z,rcyl)
    return theta

def get_bernoulli(tracers):
    "Returns the Bernoulli parameter as defined in Novikov and Thorne"
    if 'bcon' in tracers.keys():
        bsq = (tracers['bcon']*tracers['bcov']).sum(axis=-1)
    else:
        bsq = np.zeros_like(tracers['rho'])
    h = (1
         + tracers['uu']/tracers['rho']
         + tracers['Press']/tracers['rho']
         + bsq/tracers['rho'])
    Be = -tracers['ucov'][:,0]*h - 1
    return Be

def get_vr(tracers,geom):
    """Assuming tracers all extracted at same extraction radius,
    get radial velocity. Requires metric information.
    """
    r_tcr = np.mean(np.sqrt(tracers['Xcart'][:,0]**2 
                        + tracers['Xcart'][:,1]**2 
                        + tracers['Xcart'][:,2]**2))
    iX1 = np.where(geom['r'][:,0,0] >= r_tcr)[0][0]
    Lambda_r = geom['Lambda_h2bl_con'][iX1,0,1,1]
    vr = tracers['ucon'][:,1]*Lambda_r/tracers['ucon'][:,0]
    return vr

def get_gamma(tracers,geom):
    """Assuming tracers all extracted at same extraction radius,
    get boost factor. Requires metric information.
    """
    r_tcr = np.mean(np.sqrt(tracers['Xcart'][:,0]**2 
                        + tracers['Xcart'][:,1]**2 
                        + tracers['Xcart'][:,2]**2))
    iX1 = np.where(geom['r'][:,0,0] >= r_tcr)[0][0]
    alpha = geom['alpha'][iX1,:,:].mean()
    return tracers['ucon'][:,0]*alpha

def get_vsq(tracers,geom):
    """Assuming tracers all extracted at same extraction radius,
    get square spatial velocity. Requires metric information.
    """
    gamma = get_gamma(tracers,geom)
    vsq = 1 - 1/(gamma*gamma)
    return vsq

def filter_angle(tracers,thmin,thmax):
    """Returns tracers with thmin <= |theta_bl - 90| <= thmax
    thmin and thmax assumed to be in degrees
    """
    theta = get_theta(tracers)
    mask = np.logical_and(180*np.abs(theta)/np.pi >= thmin,
                          180*np.abs(theta)/np.pi <= thmax)
    return tracers.filter(mask)

def get_three_angle_filters(tracers):
    "Returns three filters used for GW170817 letter"
    thsmall = filter_angle(tracers,0,15)
    thbig = filter_angle(tracers,50,90)
    thgw170817 = filter_angle(tracers,48,71)
    return thsmall,thbig,thgw170817

def get_intersection(tracers,tracers_ref):
    "Returns particles in tracers that also have ids in tracers_ref"
    mask = np.in1d(tracers['id'],tracers_ref['id'])
    return tracers.filter(mask)

def spacetime_bins(tracers,
                   nbins_theta=55, nbins_time=50,
                   do_nans=True):
    """Calculates the average electron fraction
    as a function of absolute value of angle off the midplane,
    | theta_bl - 90 |, and time.
    ----
    Inputs:
            tracers     -- tracer data set. Assumed to be accumulated.
            nbins_theta -- number of bins in angle
            nbins_time  -- number of bins in time

    Outputs:
            theta_range -- array of positions of theta bins
            t_range     -- array of position of time bins
            yebars      -- average electron fraction binned
            mtots       -- total mass in each bin
            mtot_v_t    -- total mass integrated over angle
    """
    theta = get_theta(tracers)
    theta_deg = 180*theta/np.pi
    theta_range = np.linspace(0,90,nbins_theta+1)
    t_range = np.linspace(0,tracers['time'].max(),nbins_time+1)
    TH,TI = np.meshgrid(theta_range,t_range,indexing='ij')
    yebars = np.zeros((nbins_theta,nbins_time))
    mtots = np.zeros_like(yebars)
    for i in range(nbins_theta):
        maski = np.logical_and(np.abs(theta_deg) >= theta_range[i],
                               np.abs(theta_deg) <= theta_range[i+1])
        ttempi = tracers.filter(maski)
        for j in range(nbins_time):
            maskj = np.logical_and(ttempi['time'] >= t_range[j],
                                   ttempi['time'] <= t_range[j+1])
            ttemp = ttempi.filter(maskj)
            if len(ttemp) > 0:
                num = np.sum(ttemp['Ye']*ttemp['mass'])
                den = np.sum(ttemp['mass'])
                yebars[i,j] = num/den
                mtots[i,j] = den
            else:
                if do_nans:
                    yebars[i,j] = np.nan
                    mtots[i,j] = np.nan
                else:
                    yebars[i,j] = 0.
                    mtots[i,j] = 0.
    mtot_v_t = np.empty(nbins_time)
    for j in range(nbins_time):
        maskj = tracers['time'] <= t_range[j+1]
        ttemp = tracers.filter(maskj)
        mtot_v_t[j] = np.sum(ttemp['mass'])

    return theta_range,t_range,yebars,mtots,mtot_v_t

def get_mass_mdot(tracers,
                      nbins_theta=5,
                      nbins_time=100):
    """Calculates total mass in wind and Mdot_wind
    as a function of time and angular bin.
    """
    theta_range,t_range,_,mtots,_ = spacetime_bins(tracers,
                                                   nbins_theta=nbins_theta,
                                                   nbins_time=nbins_time,
                                                   do_nans=False)
    mtots_sum = np.zeros_like(mtots)
    mdots = np.zeros_like(mtots)
    for j in range(1,len(t_range)-1):
        mtots_sum[:,j] = mtots_sum[:,j-1] + mtots[:,j-1]
    for j in range(len(t_range)-1):
        mdots[:,j] = mtots[:,j] / (t_range[j+1] - t_range[j])
    return theta_range,t_range,mdots,mtots_sum

def plot_yebar_spacetime(tracers,
                         nbins_theta=55, nbins_time=50,
                         vmin=0.15, vmax=0.45,
                         mmin=10.**(-7.5), mmax=10.**(-1.5),
                         cmap='inferno',
                         figsize=(12,6),
                         fontsize=24,
                         savename='ye-spacetime-horiz.pdf',
                         show=False):
    """Saves a contour plot of average electron fraction vs. angle and time.
    ----
    Inputs:
            tracers     -- tracer data set. Assumed to be accumulated.
            nbins_theta -- number of bins in angle
            nbins_time  -- number of bins in time
            figsize     -- figure size in inches
            savename    -- filename to save to

    Returns:
            mesh and contour plots
    """
    units = tracers.units
    theta_range,t_range,yebars,mtots,mtot_v_t = spacetime_bins(tracers,
                                                               nbins_theta,
                                                               nbins_time)
    theta_range_avg = 0.5*(theta_range[1:] + theta_range[:-1])
    t_range_avg = 0.5*(t_range[1:] + t_range[:-1])

    fig,axarr = plt.subplots(1,2,sharey=True,
                        figsize=figsize,
                        gridspec_kw = {'width_ratios':[1,3]})
    fig.subplots_adjust(wspace=0)

    mesh = axarr[1].contourf(theta_range_avg,
                             1e3*units['T_unit']*t_range_avg,
                             yebars.transpose(),
                             vmin=vmin,vmax=vmax,
                             cmap=cmap)
    cont = axarr[1].contour(theta_range_avg,
                            1e3*units['T_unit']*t_range_avg,
                            yebars.transpose(),
                            vmin=vmin,vmax=vmax,
                            colors='k')
    cbar = plt.colorbar(mesh,label=r'$\overline{Y}_e$')
    cbar.add_lines(cont)

    axarr[0].set_ylabel('time (ms)')

    axarr[1].set_xlabel(r'$|90^\circ - \theta_{bl}|$')
    axarr[1].set_xlim(1,85)
    axarr[1].set_ylim(0,0.99*tracers['time'].max()*units['T_unit']*1e3)

    axarr[0].fill_betweenx(1e3*units['T_unit']*t_range_avg,
                           units['M_unit']*mtot_v_t/cgs['MSOLAR'],
                           lw=3,color='b',alpha=0.5)
    axarr[0].plot(units['M_unit']*mtot_v_t/cgs['MSOLAR'],
                  1e3*units['T_unit']*t_range_avg,
                  lw=3,color='k')
    axarr[0].set_xscale('log')
    axarr[0].set_xlabel(r'Ejected Mass ($M_\odot$)')
    axarr[0].set_xlim(mmax,mmin)

    if fontsize is not None:
        for ax in fig.get_axes():
            for item in ([ax.title,ax.xaxis.label,ax.yaxis.label]
                         + ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(fontsize)
    
    if savename is not None:
        plt.savefig(savename,bbox_inches='tight')
    if show:
        plt.show()
    return fig,axarr,mesh,cont

def plot_minor_summary(tracers,geom,
                       vmin=10.**(-7.5),
                       vmax=10.**(-4),
                       clip=True,
                       cmap='viridis',
                       figsize=(12,8),
                       savename='accumualted-summary.pdf',
                       show=False):
    """Plots four interesting 2d histograms of tracer properties
    with log tracer mass as the color.
    -----
    Inputs:
            tracers  -- Tracer data set. Assumed to be accumulated.
            geom     -- Geometry object. Contains metric information.
            vmin     -- Bottom of colorbar
            vmax     -- Top of colorbar
            clip     -- Set values below vmin to NaN
            cmap     -- colormap
            figsize  -- figure size
            savename -- filename to save as
            show     -- display figure
    """
    units = tracers.units

    vr = get_vr(tracers,geom)
    theta = get_theta(tracers)

    norm = LogNorm(vmin=10.**(-7.5),vmax=10.**(-4),clip=True)
    weights = tracers['mass']*units['M_unit']/cgs['MSOLAR']

    fig,axarr = plt.subplots(2,2,figsize=(12,8))
    _,_,_,histY = axarr[0,0].hist2d(1e3*tracers.units['T_unit']*tracers['time'],
                                    tracers['Ye'],bins=50,
                                    weights=weights,
                                    cmap=cmap,
                                    norm=norm)
    axarr[0,0].set_xlabel('time (ms)')
    axarr[0,0].set_ylabel(r'$Y_e$')

    _,_,_,histTH = axarr[0,1].hist2d(1e3*tracers.units['T_unit']*tracers['time'],
                                     theta,bins=50,
                                     weights=weights,
                                     cmap=cmap,
                                     norm=norm)
    axarr[0,1].set_xlabel('time (ms)')
    axarr[0,1].set_ylabel(r'$\theta$')

    _,_,_,histS = axarr[1,0].hist2d(tracers['s'],tracers['Ye'],
                                    bins=100,range=[[10,40],[0.01,0.6]],
                                    weights=weights,
                                    cmap=cmap,
                                    norm=norm)
    axarr[1,0].set_xlabel(r'$s$ ($k_b/$b)')
    axarr[1,0].set_ylabel(r'$Y_e$')

    _,_,_,histV = axarr[1,1].hist2d(vr,tracers['Ye'],bins=50,
                                    weights=weights,
                                    cmap=cmap,
                                    norm=norm)
    axarr[1,1].set_xlabel(r'$v^r/c$')
    axarr[1,1].set_ylabel(r'$Y_e$')

    plt.tight_layout()

    cbar = fig.colorbar(histY,ax=axarr.ravel().tolist(),
                        label=r'Traced Mass ($M_{\odot}$)')

    if savename is not None:
        plt.savefig(savename,bbox_inches='tight')
    if show:
        plt.show()
    return histY,histTH,histS,histV

def plot_ye_th_hist(tracers,
                    tracers_nse,
                    figsize=(8,12),
                    fontsize=18,
                    vmin=10.**(-7.),
                    clip=True,
                    histmin=10.**(-7.),
                    histmax=10.**(-3.7),
                    xlabel=r'$Y_e|_{5GK}$',
                    savename='ye-vs-theta-folded.png',
                    show=False):
    "Plots 2d histogram of Ye vs theta, with 1d cuts"
    # fontsize assumed to be 18. Figsize assumed to be (8,12)
    from scipy import integrate

    solid1 = 2*np.pi*integrate.quad(np.cos,0,15*np.pi/180)[0]
    solid2 = 2*np.pi*integrate.quad(np.cos,50*np.pi/180,90*np.pi/180)[0]
    solid3 = 2*np.pi*integrate.quad(np.cos,48*np.pi/180,71*np.pi/180)[0]

    units = tracers.units
    theta = get_theta(tracers)

    thsmall,thbig,thgw170817 = get_three_angle_filters(tracers)
    thsmall_nse = get_intersection(tracers_nse,thsmall)
    thbig_nse = get_intersection(tracers_nse,thbig)
    thgw170817_nse = get_intersection(tracers_nse,thgw170817)

    M_UNIT = units['M_unit']/cgs['MSOLAR']
    norm = LogNorm(vmin=vmin,clip=clip)
    weights = tracers['mass']*M_UNIT

    fig,axarr = plt.subplots(2,sharex=True,
                             figsize=figsize,
                             gridspec_kw = {'height_ratios':[2,1]})
    fig.subplots_adjust(hspace=0)

    n,bins,patches,hist = axarr[0].hist2d(tracers['Ye'],
                                          180*np.abs(theta)/np.pi,bins=50,
                                          weights=weights,
                                          norm=norm,
                                          alpha=1.)
    axarr[0].set_ylabel(r'$|90^\circ - \theta_{bl}|$ (degrees)')
    axarr[0].set_ylim(0,90)

    ax1_divider = make_axes_locatable(axarr[0])
    cbar_ax= ax1_divider.append_axes("top",size="5%",pad="1%")
    cbar = fig.colorbar(hist,label=r'Traced Mass ($M_{\odot}$)',
                        orientation='horizontal',
                        cax=cbar_ax)
    cbar_ax.xaxis.set_ticks_position("top")
    cbar_ax.xaxis.set_label_position("top")

    axarr[0].hlines([48,71],axarr[0].get_xlim()[0],axarr[0].get_xlim()[-1],
                    colors='k',linestyles='--',linewidths=3)
    axarr[0].text(0.1575,65,'gw170817',fontsize=13)

    rect1 = mpl.patches.Rectangle((1.02*axarr[0].get_xlim()[0],
                                   axarr[0].get_ylim()[0]+0.5),
                                  0.15,15.,linewidth=3,
                                  edgecolor='r',facecolor='none')
    axarr[0].add_patch(rect1)
    rect2 = mpl.patches.Rectangle((0.21,50),
                                  0.99*axarr[0].get_xlim()[-1]-0.21,
                                  0.99*axarr[0].get_ylim()[-1]-50,
                                  linewidth=3,
                                  edgecolor='b',facecolor='none')
    axarr[0].add_patch(rect2)

    _,_,patches1 = axarr[1].hist(thsmall_nse['Ye'],60,
                                 weights=(1./solid1)*thsmall_nse['mass']*M_UNIT,
                                 alpha=1.,color='r',
                                 log=True, density=False,rasterized=True,
                                 label=r'$|90^\circ -\theta_{bl}| \leq 15^\circ$')
    _,_,patches2 = axarr[1].hist(thbig_nse['Ye'],60,
                                 weights=(1./solid2)*thbig_nse['mass']*M_UNIT,
                                 alpha=0.75,color='b',
                                 log=True, density=False,rasterized=True,
                                 label=r'$|90^\circ -\theta_{bl}| \geq 50^\circ$')
    _,_,patches3 = axarr[1].hist(thgw170817_nse['Ye'],60,
                                 weights=(1./solid3)*thgw170817_nse['mass']*M_UNIT,
                                 alpha=1,color='k',
                                 log=True, density=False,rasterized=True,
                                 histtype='step',lw=4,ls='-',
                                 label='gw170817')
    axarr[1].set_xlabel(xlabel)
    axarr[1].set_ylabel('Traced Mass per\nSolid Angle '+r'($M_{\odot}$sr$^{-1}$)')
    axarr[1].set_ylim(histmin,histmax)

    axarr[1].legend(bbox_to_anchor=(0.475,0.8),framealpha=1,fancybox=True)

    for ax in fig.get_axes():
        for item in ([ax.title,ax.xaxis.label,ax.yaxis.label]
                     + ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)

    if savename is not None:
        plt.savefig(savename,
                    bbox_inches='tight',dpi=200)
    if show:
        plt.show()

    return fig,axarr

def plot_mdot_mtot(tracers,
                   nbins_theta=5,
                   nbins_time=100,
                   figsize=(8,8),
                   savename='mass-mdot-v-t-th.pdf',
                   show=False):
    "Plots mass and mass outflow rate in the wind"
    units = tracers.units
    T_unit = units['T_unit']*1e3
    M_unit = units['M_unit'] / cgs['MSOLAR']
    MDOT_unit = M_unit / T_unit

    theta_range,t_range,mdots,mtots_sum = get_mass_mdot(tracers,
                                                        nbins_theta,
                                                        nbins_time)
    theta_range_avg = 0.5*(theta_range[1:] + theta_range[:-1])
    t_range_avg = 0.5*(t_range[1:] + t_range[:-1])

    fig,axarr = plt.subplots(2,sharex=True,figsize=figsize)

    lines = [None for i in range(len(theta_range)-1)]
    labels = [None for i in range(len(theta_range)-1)]
    for i in range(len(theta_range)-1):
        lines[i], = axarr[0].plot(T_unit*t_range_avg,mdots[i]*MDOT_unit)
        labels[i] = r'$%04.1f^\circ \leq \theta \leq %04.1f^\circ$' \
                    % (theta_range[i], theta_range[i+1])
        axarr[1].plot(T_unit*t_range_avg,mtots_sum[i]*M_unit)

    axarr[0].set_ylabel(r'$\dot{M}_w$ ($M_\odot /$ms)')
    axarr[1].set_xlabel(r'$t$ (ms)')
    axarr[1].set_ylabel(r'$M_w$ ($M_\odot$)')

    plt.legend(lines,labels,loc='lower left',
               bbox_to_anchor=(0.125,0.25),
               bbox_transform=plt.gcf().transFigure,
               fancybox=True,
               shadow=True)
    if savename is not None:
        plt.savefig(savename,bbox_inches='tight')
    if show:
        plt.show()
    return lines


if __name__ == "__main__":
    from argparse import ArgumentParser
    import hdf5_to_dict as io
    parser = ArgumentParser(
        description='Generate plots based on accumualted tracer data'
    )
    parser.add_argument('dumps',type=str,
                        help='directory containing simulation dumps')
    parser.add_argument('tracers',
                        type=str,
                        help=('File containing tracers '
                              +'extracted at some radius'))
    parser.add_argument('tracers_nse',
                        type=str,
                        help=('File containing same tracers, '
                              +'but extracted when they drop out of NSE.'))
    parser.add_argument('-s','--save',
                        dest='save',
                        action='store_true',
                        help='Save plots rather than display them')

    args = parser.parse_args()

    hdr,geom = io.load_hdr_geom(args.dumps)

    tracers = io.TracerData.fromfile(args.tracers)
    tracers_nse = io.TracerData.fromfile(args.tracers_nse)
    
    tracers = tracers.filter(get_bernoulli(tracers) > 0)
    tracers_nse = get_intersection(tracers_nse,tracers)

    if args.save:
        plot_yebar_spacetime(tracers)
        plot_minor_summary(tracers,geom)
        plot_ye_th_hist(tracers,tracers_nse)
        plot_mdot_mtot(tracers)
    else:
        plot_yebar_spacetime(tracers,fontsize=None,
                             savename=None,show=True)
        plot_minor_summary(tracers,geom,savename=None,show=True)
        plot_ye_th_hist(tracers,tracers_nse,savename=None,show=True)
        plot_mdot_mtot(tracers,savename=None,show=True)
