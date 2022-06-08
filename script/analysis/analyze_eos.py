#!/usr/bin/env python

################################################################################
#                                                                              # 
#  PLOTS SOME QUANTITIES FOR EOS                                               # 
#                                                                              # 
################################################################################

import numpy as np
import h5py
from scipy import interpolate,optimize
import matplotlib as mpl
from matplotlib import pyplot as plt
import sys,os
import argparse

CL = 2.99792458e10
SMALL = 1e-20
FONTSIZE = 18
LINEWIDTH=3

mpl.rcParams.update({'font.size':FONTSIZE})
# mpl.rcParams.update({'font.family':'serif'})
# mpl.rcParams.update({'mathtext.fontset':'stix'})

# ======================================================================
# Main function
# ======================================================================
# TODO: may need to be tuned per EOS
def analyze_eos(s, ye, infile_name, postfix='png'):
    print("Loading EOS")
    eos = load_eos(infile_name)
    outpath = get_outpath(infile_name)
    print("Will save to {}".format(outpath))
    print("Plotting entropy")
    plot_entropy(eos,ye,
                 figname = os.path.join(outpath,
                                        'entropy.'+postfix))
    plt.cla()
    plt.clf()
    print("Plotting enthalpy")
    plot_hm1(eos,ye,
             figname = os.path.join(outpath,
                                    'hm1.'+postfix))
    plt.cla()
    plt.clf()
    print("Making contour plots")
    contour_entropy(eos,ye,
                    figname = os.path.join(outpath,
                                           's_contours.'+postfix))
    plt.cla()
    plt.clf()
    contour_hm1_entropy(eos,ye,
                        figname = os.path.join(outpath,
                                               's_hm1_contours.'+postfix))
    plt.cla()
    plt.clf()

    print("Plotting adiabat")
    print("...summary")
    a = Adiabat(s,ye,eos)
    plot_adiabat(a,eos,
                 figname = os.path.join(outpath,'adiabat.'+postfix))
    plt.cla()
    plt.clf()

    print("...lP")
    plot_lP_adiabat(a,eos,
                    figname = os.path.join(outpath,'lP_adiabat.'+postfix))
    plt.cla()
    plt.clf()

    print("Plotting thermodynamics of disk edge")
    edge_map = map_disk_edge(eos,ye)
    plot_edge_map(*edge_map,
                  figname=os.path.join(outpath,'edge_map.'+postfix))
    plt.cla()
    plt.clf()

    return
    
# ======================================================================


# ======================================================================
# IO
# ======================================================================
def load_eos(filename):

    eos = {}
    # load file
    with h5py.File(filename,'r') as f:
        for k,v in f.items():
            eos[k] = v[()]

    # derived quantities
    lrho = eos['logrho']
    lT = eos['logtemp']
    Ye = eos['ye']
    lP = eos['logpress']
    le = eos['logenergy']
    ent = eos['entropy']
    rho = 10.**lrho
    T = 10.**lT
    P = 10.**lP
    e = 10.**le - eos['energy_shift']
    w = rho*e + P
    h = CL*CL + e + P/rho
    hgeom = h/(CL*CL)
    hm1 = hgeom -1
    lhm1 = np.log10(np.abs(hm1))

    # No maximum filter
    dpdrhoe = eos['dpdrhoe']
    dpderho = eos['dpderho']
    cs2 = (dpdrhoe + (P/(rho*rho))*dpderho)/h
    cs2 = np.minimum(cs2,CL*CL)
    cs2 = np.maximum(cs2,10*SMALL)
    cs = np.sqrt(cs2)

    # quick aliases
    eos['lrho'] = lrho
    eos['lT'] = lT
    eos['Ye'] = Ye
    eos['lP'] = lP
    eos['le'] = le
    eos['ent'] = ent
    
    eos['rho'] = rho
    eos['T'] = T
    eos['P'] = P
    eos['e'] = e
    eos['w'] = w
    eos['h'] = h
    eos['hgeom'] = hgeom
    eos['hm1'] = hm1
    eos['lhm1'] = lhm1

    eos['cs2'] = cs2
    eos['cs'] = cs
    
    return eos
# ======================================================================


# ======================================================================
# Adiabat workhorse
# ======================================================================
class Adiabat:
    def __init__(self, s, Ye, eos):
        self.s = s
        self.Ye = Ye
        self.eos = eos
        self.iYe = get_iYe(Ye,eos)
        self.lrho_bounds, self.ilrho_bounds = get_lrho_min_max(s,Ye,eos)
        self.lT_adiabat = get_adiabat(s,self.iYe, eos,
                                      self.lrho_bounds,
                                      self.ilrho_bounds)
        self.slc = np.s_[self.ilrho_bounds[0]:self.ilrho_bounds[1]]
        self.lrho = self.bnd1d(eos['lrho'])

    @classmethod
    def get_if_valid(cls, s, Ye, eos):
        try:
            adiabat = cls(s, Ye, eos)
        except:
            raise ValueError("Invalid Adiabat")
        if not adiabat.is_valid():
            raise ValueError("Invalid Adiabat")
        return adiabat

    def bnd1d(self,var):
        return var[self.slc]

    def bnd3d(self,var):
        lrho = self.eos['lrho']
        lT = self.eos['lT']
        var1d = np.empty_like(lrho)
        for i in range(self.ilrho_bounds[0],self.ilrho_bounds[1]):
            vartemp = var[self.iYe,:,i]
            varinterp = interpolate.interp1d(lT,vartemp)
            t = self.lT_adiabat[i]
            v = varinterp(t)
            var1d[i] = v
        return self.bnd1d(var1d)

    def lT(self):
        return self.bnd1d(self.lT_adiabat)

    def project(self,var):
        if len(var.shape) == 1:
            if len(var) == len(self.eos['lT']):
                return self.lT()
            return self.bnd1d(var)
        if len(var.shape) == 3:
            return self.bnd3d(var)
        raise ValueError("var must be 1d or 3d")

    def hm1bc(self):
        lP = self.project(self.eos['lP'])
        hm1 = self.project(self.eos['hm1'])
        lrho = self.project(self.eos['lrho'])
        dPdrhoe = self.project(self.eos['dpdrhoe'])
        rho = 10.**lrho
        P = 10.**lP
        rhointerp = interpolate.UnivariateSpline(rho,P)
        dPdrhos = rhointerp.derivative()
        hm1interp = interpolate.interp1d(lP,hm1)
        dPdrhoe_interp = interpolate.interp1d(lP,dPdrhoe)
        out = [lP.min(),
               float(dPdrhos(rho.min())),
               float(dPdrhoe_interp(lP.min())),
               float(hm1interp(lP.min()))]
        return out

    def is_valid(self):
        return np.all(np.gradient(self.project(self.eos['lT'])) >= 0)

    def get_rho_of_hm1_interp(self):
        rho_grid = self.project(self.eos['rho'])
        hm1_grid = self.project(self.eos['hm1'])
        fill_value = (0, self.project(self.eos['rho']).max())
        return interpolate.interp1d(hm1_grid,rho_grid,
                                    bounds_error = False,
                                    fill_value=fill_value)
    
    def __call__(self,var):
        return self.project(var)

def map_disk_edge(eos, Ye, entropies = None):
    if entropies is None:
        entropies = np.arange(2.,40.,0.5)

    lPmins = np.empty_like(entropies)
    hm1bcs = np.empty_like(entropies)
    dPdrhoss = np.empty_like(entropies)
    dPdrhoes = np.empty_like(entropies)

    for i,s in enumerate(entropies):
        try:
            a = Adiabat(s,Ye,eos)
            lPmins[i],dPdrhoss[i],dPdrhoes[i],hm1bcs[i] = a.hm1bc()
        except ValueError:
            lPmins[i],dPdrhoss[i],dPdrhoes[i],hm1bcs[i] = [np.NaN for i in range(4)]

    entropies = entropies[np.isfinite(hm1bcs)]
    lPmins = lPmins[np.isfinite(hm1bcs)]
    dPdrhoss = dPdrhoss[np.isfinite(hm1bcs)]
    dPdrhoes = dPdrhoes[np.isfinite(hm1bcs)]
    hm1bcs = hm1bcs[np.isfinite(hm1bcs)]

    return entropies,lPmins,dPdrhoss,dPdrhoes,hm1bcs
# ======================================================================


# ======================================================================
# Visualization
# ======================================================================
def plot_entropy(eos, ye, figname = None,
                 vmax = 4, levels = [1,3,7,10,20,100]):
    mpl.rcParams.update({'font.size':FONTSIZE})

    iYe = get_iYe(ye,eos)
    lT = eos['lT']
    lrho = eos['lrho']
    ent = eos['ent']

    mesh = plt.pcolormesh(lT,lrho,np.log10(ent[iYe,:,:]).transpose(),vmax=vmax)
    mesh.set_edgecolor('face')
    plt.colorbar(mesh,label=r'$\log_{10}$entropy ($k_b/$baryon)')
    CS = plt.contour(lT,lrho,ent[iYe,:,:].transpose(),
                     levels=levels,colors='k')
    plt.clabel(CS,inline=1)

    # fig = plt.gcf()
    # fig.set_size_inches(12,8)
    plt.ylabel(r'$\log_{10}\rho$ (cgs)')
    plt.xlabel(r'$\log_{10}T$ (MeV)')
    if figname is not None:
        plt.savefig(figname,bbox_inches='tight',rasterized=True)

def plot_hm1(eos, ye, figname = None,
             vmax = 0.025, levels = [0, 0.006, 0.009, 0.015, 0.02]):
    mpl.rcParams.update({'font.size':FONTSIZE})

    iYe = get_iYe(ye,eos)
    lT = eos['lT']
    lrho = eos['lrho']
    hm1 = eos['hm1']

    mesh = plt.pcolormesh(lT,lrho,hm1[iYe,:,:].transpose(),vmax=vmax)
    mesh.set_edgecolor('face')
    plt.colorbar(mesh,label=r'$\frac{h}{c^2}-1$')
    CS = plt.contour(lT,lrho,hm1[iYe,:,:].transpose(),colors='k',
                     levels=levels)
    plt.clabel(CS,inline=1)

    # fig = plt.gcf()
    # fig.set_size_inches(12,8)
    plt.ylabel(r'$\log_{10}\rho$ (cgs)')
    plt.xlabel(r'$\log_{10}T$ (MeV)')
    if figname is not None:
        plt.savefig(figname, bbox_inches='tight',rasterized=True)

def contour_entropy(eos, ye, figname = None,
                    vmax = 1, levels = [1,3,7,20,100]):
    mpl.rcParams.update({'font.size':FONTSIZE})

    iYe = get_iYe(ye,eos)
    lT = eos['lT']
    lrho = eos['lrho']
    ent = eos['ent']
    hm1 = eos['hm1']
    lhm1 = eos['lhm1']

    mesh = plt.pcolormesh(lT,lrho,lhm1[iYe,:,:].transpose(),vmax=vmax)
    mesh.set_edgecolor('face')
    plt.colorbar(mesh,label=r'$\log_{10}\left(\frac{h}{c^2}-1\right)$')
    CS = plt.contour(lT,lrho,ent[iYe,:,:].transpose(),
                     levels=levels,colors='r')
    plt.clabel(CS,inline=1)

    lines = [CS.collections[0]]
    labels = [r's ($k_b/$baryon)']
    plt.legend(lines,labels,loc = 'lower right')

    # fig = plt.gcf()
    # fig.set_size_inches(12,8)
    plt.ylabel(r'$\log_{10}\rho$ (cgs)')
    plt.xlabel(r'$\log_{10}T$ (MeV)')
    if figname is not None:
        plt.savefig(figname, bbox_inches='tight', rasterized=True)

def contour_hm1_entropy(eos, ye, figname = None,
                        levels_hm1 = [0.006,0.009, 0.015],
                        levels_ent = [1,4,8,12,16]):

    mpl.rcParams.update({'font.size':FONTSIZE})

    iYe = get_iYe(ye,eos)
    lT = eos['lT']
    lrho = eos['lrho']
    ent = eos['ent']
    hm1 = eos['hm1']
    
    CS1 = plt.contour(lT,lrho,hm1[iYe,:,:].transpose(),
                      levels=levels_hm1,colors='r')
    plt.clabel(CS1,inline=1)

    CS2 = plt.contour(lT,lrho,ent[iYe,:,:].transpose(),
                      levels=levels_ent,colors='k')
    plt.clabel(CS2,inline=1)

    lines = [CS1.collections[0], CS2.collections[0]]
    labels = [r'$\frac{h}{c^2}-1$', r'entropy ($k_b$/baryon)']
    plt.legend(lines,labels, loc = 'lower right')

    # fig = plt.gcf()
    # fig.set_size_inches(12,8)
    
    plt.xlabel(r'$\log_{10}T$ (MeV)')
    plt.ylabel(r'$\log_{10}\rho$ (cgs)')
    
    if figname is not None:
        plt.savefig(figname, bbox_inches='tight')

def plot_edge_map(entropies,lPmins,dPdrhoss,dPdrhoes,hm1bcs,
                  figname = None):
    mpl.rcParams.update({'font.size':FONTSIZE})

    fig, axarr = plt.subplots(2,2,sharex=True)
    lPax = axarr[1,0]
    hm1ax = axarr[0,0]
    dpdre_ax = axarr[1,1]
    dpdrs_ax = axarr[0,1]
    
    lPax.plot(entropies,lPmins)
    lPax.set_ylabel(r'$\min\left(\log_{10}P\right)$')
    
    hm1ax.plot(entropies,hm1bcs)
    hm1ax.set_ylabel(r'$\left[\frac{h}{c^2} - 1\right]_{\min\left(\log_{10}P\right)}$')
    
    dpdre_ax.plot(entropies,dPdrhoes/1e17)
    dpdre_ax.set_ylabel(r'$10^{17}\times \left[\frac{\partial P}{\partial\rho}\right]_\varepsilon^{\min P}$')
    
    dpdrs_ax.plot(entropies,dPdrhoss/1e17)
    dpdrs_ax.set_ylabel(r'$10^{17}\times \left[\frac{\partial P}{\partial\rho}\right]_s^{\min P}$')

    axarr[1,0].set_xlabel(r'entropy ($k_b$/baryon)')
    axarr[1,1].set_xlabel(r'entropy ($k_b$/baryon)')

    # fig.set_size_inches(12,8)
    # plt.tight_layout()

    if figname is not None:
        plt.savefig(figname,bbox_inches='tight')

def plot_adiabat(a, eos, figname=None):

    mpl.rcParams.update({'font.size':FONTSIZE})

    lP = eos['lP']
    hm1 = eos['hm1']
    lhm1 = eos['lhm1']
    lrho = eos['lrho']
    lT = eos['lT']
    le = eos['le']

    fig, axarr = plt.subplots(2,2)

    hm1_ax = axarr[0,0]
    hm1_ax.plot(a(lP),a(lhm1),lw=LINEWIDTH)
    hm1_ax.set_xlabel(r'$\log_{10}P$ (cgs)')
    hm1_ax.set_ylabel(r'$\log_{10}(h/c^2-1)$')

    lT_ax = axarr[0,1]
    lT_ax.plot(a(lrho),a(lT),lw=LINEWIDTH)
    lT_ax.set_xlabel(r'$\log_{10}\rho$ (cgs)')
    lT_ax.set_ylabel(r'$\log_{10}T$ (MeV)')

    P_ax = axarr[1,1]
    P_ax.plot(a(lrho),a(lP),lw=LINEWIDTH)
    P_ax.set_xlabel(r'$\log_{10}\rho$')
    P_ax.set_ylabel(r'$\log_{10}P$')

    eps_ax = axarr[1,0]
    eps_ax.plot(a(lrho),a(le),lw=LINEWIDTH)
    eps_ax.set_xlabel(r'$\log_{10}\rho$')
    eps_ax.set_ylabel(r'$\log_{10}\varepsilon$')

    plt.suptitle(r'$s = %s$ $k_b$/baryon, $Y_e = %s$' % (a.s, a.Ye),y=1.08)

    # fig.set_size_inches(12,8)
    # plt.tight_layout()

    if figname is not None:
        plt.savefig(figname,bbox_inches='tight')

def plot_lP_adiabat(a, eos, figname = None):

    mpl.rcParams.update({'font.size':FONTSIZE})

    lP = eos['lP']
    lhm1 = eos['lhm1']

    plt.plot(a(lP),a(lhm1), lw = LINEWIDTH)
    plt.xlabel(r'$\log_{10}P$ (cgs)')
    plt.ylabel(r'$\log_{10}(h/c^2-1)$')

    ax = plt.gca()
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(FONTSIZE)

    if figname is not None:
        plt.savefig(figname, bbox_inches = 'tight')

# ======================================================================


# ======================================================================
# Utility functions
# ======================================================================
def get_iYe(mYe,eos):
    iYe = np.where(eos['Ye'] >= mYe)[0][0]
    return iYe

def get_lrho_min_max(s, Ye, eos):
    iYe = get_iYe(Ye,eos)
    lrho = eos['lrho']
    ent = eos['ent']
    sinterp_minT = interpolate.interp1d(lrho,ent[iYe,0])
    sinterp_maxT = interpolate.interp1d(lrho,ent[iYe,-1])
    try:
        lrho_min = optimize.brentq(lambda r: sinterp_minT(r) - s,
                                   lrho.min(),lrho.max())
    except ValueError:
        lrho_min = lrho.min()
    try:
        lrho_max = optimize.brentq(lambda r: sinterp_maxT(r) - s,
                                   lrho.min(), lrho.max())
    except ValueError:
        lrho_max = lrho.max()

    ilrho_min = np.where(lrho >= lrho_min)[0][0]
    ilrho_max = np.where(lrho >= lrho_max)[0][0]

    return (lrho_min, lrho_max),(ilrho_min,ilrho_max)

def get_adiabat(s, iYe, eos, lrho_bounds, ilrho_bounds):
    lrho = eos['lrho']
    lT = eos['lT']
    ent = eos['ent']
    lrho_min,lrho_max = lrho_bounds
    ilrho_min,ilrho_max = ilrho_bounds

    lT_adiabat = -np.infty*np.ones_like(lrho)
    for i in range(ilrho_min,ilrho_max):
        s_interp = interpolate.interp1d(lT,ent[iYe,:,i])
        lTs = optimize.brentq(lambda t: s_interp(t) - s,lT.min(),lT.max())
        lT_adiabat[i] = lTs

    return lT_adiabat

def get_outpath(infile):
    outpath_name = infile.rstrip('.h5') + '_analysis'
    if not os.path.exists(outpath_name):
        os.makedirs(outpath_name)
    return outpath_name
# ======================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plot various properties of an EOS relevant for Fishbone-Moncrief disks')
    parser.add_argument('s',type=float,
                        help='Entropy for isocontours. In k_b/baryon')
    parser.add_argument('Ye',type=float,
                        help='Electron fraction. Assumed to be constant.')
    parser.add_argument('filename',type=str,
                        help='Name of EOS file to read. Assumed to be in Stellar Collapse format.')
    parser.add_argument('--pdf',dest='pdf',action='store_true',
                        help='Set to save plots to pdf format')
    args = parser.parse_args()
    print(("Analyzing EOS {}\n\t"
           +"Assuming:"
           +"\n\t\tEntropy s = {}\n\t\t"
           +"Electron fraction Ye = {}.").format(args.filename,
                                                 args.s,args.Ye))

    postfix = 'pdf' if args.pdf else 'png'
    analyze_eos(args.s,args.Ye,args.filename,postfix)
    print("Done!")
