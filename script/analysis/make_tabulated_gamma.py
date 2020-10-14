#!/usr/bin/env python

################################################################################
#                                                                              #
#  UTILITY---TABULATES AN IDEAL GAS LAW FOR STELLAR_CLLAPSE READER             #
#                                                                              #
################################################################################

from __future__ import print_function
import numpy as np
import h5py
from units import UnitSystem
import units as bunits

# constants, in CGS
GNEWT = bunits.cgs['GNEWT']  # Gravitational constant
ME    = bunits.cgs['ME']     # Electron mass
MP    = bunits.cgs['MP']     # Proton mass
CL    = bunits.cgs['CL']     # Speed of light
EV    = bunits.cgs['EV']     # Electron-volt
MEV   = bunits.cgs['MEV']    # Mega-Electron-Volt
MSUN  = bunits.cgs['MSOLAR'] # Solar mass
KBOL  = bunits.cgs['KBOL']   # boltzmann constant

def press(rho, eps, gam):
    return (gam - 1.)*rho*eps

def e2temp(eps, gam, M_unit = MP):
    #out = M_unit*(gam - 1.)*eps
    #out = MP*(gam - 1.)*eps
    out = (gam - 1.)*eps
    return out

def eps(u,rho):
    return u / rho

def t2eps(T, gam, M_unit = MP):
    #out = T / (M_unit*(gam - 1.))
    #out = T / ((gam - 1.))
    out = T*MEV / (MP*(gam - 1.))
    return out

def adiabatic_const(rho, eps, gam):
    u = rho*eps
    return u*rho**(-gam)

def entropy(rho, eps, gam):
    return (gam - 1.)*adiabatic_const(rho, eps, gam)

def sound_speed(rho, eps, gam):
    u = rho*eps
    ef = rho*CL*CL + u*gam
    press = (gam - 1.)*u
    cs2 = CL*CL*gam*press/ef
    # fac = (gam*(gam-1.)*eps)/(CL*CL+eps+(gam-1.)*eps)
    # cs2 = fac*CL*CL
    return np.sqrt(cs2)

def make_filename(gam):
    gamstring = str(gam)
    gamstring = gamstring.replace('.','p')
    name = "sc_eos_gamma_{}.h5".format(gamstring)
    return name

def make_table_u(rho_min, rho_max, n_rho,
                 u_min, u_max, n_u,
                 ye_min, ye_max, n_ye,
                 units, gam, filename = None,
                 crash_on_sound_speed = True):
    eps_min  = eps(u_min, rho_max)
    eps_max  = eps(u_max, rho_min)
    temp_min = e2temp(eps_min, gam, units.M)
    temp_max = e2temp(eps_max, gam, units.M)
    n_temp = n_u
    print("eps min, max = [%e, %e]"  % (eps_min, eps_max))
    print("Temp min, max = [%e, %e]" % (temp_min, temp_max))
    return make_table_temp(rho_min,rho_max,n_rho,
                           temp_min,temp_max,n_temp,
                           ye_min, ye_max, n_ye,
                           units, gam, filename,
                           crash_on_sound_speed)

def make_table_temp(rho_min, rho_max, n_rho,
                    temp_min, temp_max, n_temp,
                    ye_min, ye_max, n_ye,
                    units, gam, filename = None,
                    crash_on_sound_speed = True):
    assert gam > 1
    # from code units to CGS
    rho_min  *= units.RHO
    rho_max  *= units.RHO
    temp_min *= units.T # temperatures are in MeV
    temp_max *= units.T
    energy_shift = 0.0  # all energies positive for gamma law
    # push system into sub-luminal regime
    # 1D arrays
    lrho = np.linspace(np.log10(rho_min), np.log10(rho_max), n_rho)
    lt   = np.linspace(np.log10(temp_min), np.log10(temp_max), n_temp)
    ye   = np.linspace(ye_min, ye_max, n_ye)
    rho  = 10.**(lrho)
    temp = 10.**(lt)
    # 3D arrays
    YE,T,RHO = np.meshgrid(ye,temp,rho, indexing='ij')
    EPS = t2eps(T, gam, units.M)
    P   = press(RHO, EPS, gam)
    CS  = sound_speed(RHO, EPS, gam)
    if CS.max() >= CL and crash_on_sound_speed:
        raise RuntimeError("Sound speed exceeds speed of light!"
                           +" Max CS/CL = {}".format(CS.max()/CL))
    print("[cs_min/cl, cs_max/cl] = [{}, {}]".format(CS.min()/CL, CS.max()/CL))
    CS2 = CS*CS
    ENT = entropy(RHO, EPS, gam)
    #GAMMA = gam*np.ones_like(P)
    GAMMA = np.ones_like(P)
    DPDRHOE = (gam-1.)*EPS
    DPDERHO = (gam-1.)*RHO
    # logs
    LP = np.log10(P)
    LE = np.log10(EPS)
    # compositions and mass fractions
    # assume ionized hydrogen
    # Xa: mass fraction of alpha particles
    Xa = np.zeros_like(P)
    # Xh: mass fraction of heavy nuclei
    Xh = np.zeros_like(P)
    # Xp: mass fraction of protons
    Xp = YE
    # Xn: mass fraction of neutrons
    Xn = 1 - Xp
    # Abar: average atomic mass
    Abar = np.ones_like(P)
    # Zbar: average atomic number
    Zbar = np.ones_like(P)

    # make the hdf5 file
    
    if filename is None:
        filename = make_filename(gam)
    with h5py.File(filename,'w') as f:
        f.create_dataset('energy_shift', data = np.array([energy_shift]))
        f.create_dataset('pointsrho', data = np.array([n_rho]))
        f.create_dataset('pointstemp', data = np.array([n_temp]))
        f.create_dataset('pointsye', data = np.array([n_ye]))
        f.create_dataset('cs2', data = CS2)
        f.create_dataset('dpdrhoe', data=DPDRHOE)
        f.create_dataset('dpderho', data=DPDERHO)
        f.create_dataset('entropy', data = ENT)
        f.create_dataset('gamma', data = GAMMA)
        f.create_dataset('logenergy', data = LE)
        f.create_dataset('logpress', data = LP)
        f.create_dataset('logrho', data = lrho)
        f.create_dataset('logtemp', data = lt)
        f.create_dataset('ye', data = ye)
        f.create_dataset("Xa", data = Xa)
        f.create_dataset("Xh", data = Xh)
        f.create_dataset("Xp", data = Xp)
        f.create_dataset("Xn", data = Xn)
        f.create_dataset("Abar", data = Abar)
        f.create_dataset("Zbar", data = Zbar)
    print("Table created. Filename = {}.".format(filename))
    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description = ("Generate a table to be read by Stellar Collapse reader"
                       +" using tabulated equation of state."),
        epilog = "Defaults chosen to mimic Stellar Collapse tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--gamma', default = 1.4, type=float,
                        help = 'adiabatic index')
    parser.add_argument('-L', default=1e3, type=float,
                        help="Length unit in cgs. Can be used instead of --mbh")
    parser.add_argument('-M', default = 1e3, type=float,
                        help="Mass unit in cgs")
    parser.add_argument('--mbh', default=None, type=float,
                        help="Black hole mass. Can replace length unit.")
    parser.add_argument('--rhomin', default = 1e-4, type=float,
                        help="Minimum density. In code units.")
    parser.add_argument('--rhomax', default = 1e2, type=float,
                        help="Maximum density. In code units.")
    parser.add_argument('--nrho', default = 234, type=int,
                        help = "Number of points in rho")
    parser.add_argument('--umin', default = 1e-7, type=float,
                        help = "Minimum internal energy density. In code units.")
    parser.add_argument('--umax', default = 1e1, type=float,
                        help = "Maximum internal energy density. In code units.")
    parser.add_argument('--tempmin', default = None, type=float,
                        help = "Minimum temperature. Overwrites umin if set.")
    parser.add_argument('--tempmax', default = None, type=float,
                        help = "Maximum temperature. Overwrites umax if set.")
    parser.add_argument('--ntemp', default = 136, type=int,
                        help = "Number of points in temperature/internal energy")
    parser.add_argument('--yemin', default = 0.0, type=float,
                        help = "Minimum electron fraction.")
    parser.add_argument('--yemax', default = 0.55, type=float,
                        help = "Maximum electron fraction.")
    parser.add_argument('--nye', default = 50, type=int,
                        help = "Number of points in electron fraction")
    parser.add_argument('-o', '--output', default = None, type=str,
                        help = "Name of output file.")
    args = parser.parse_args()
    if args.L is None and args.mbh is None:
        raise ValueError("Either L or mbh must be set.")
    print("Welcome to the table generator.")
    if args.mbh is not None:
        print("\tUsing black hole mass = {}".format(args.mbh))
        units = UnitSystem(args.M, Mbh = args.mbh)
    else:
        units = UnitSystem(args.M, L_unit = args.L)
    print("Chosen units are:")
    print("\t[M, L, RHO, U] = [{}, {}, {}, {}]".format(units.M,
                                                       units.L,
                                                       units.RHO,
                                                       units.U))
    print("Making table!")
    print("\tgamma = {}".format(args.gamma))
    print("\tUsing [nrho, ntemp, nye] = [{}, {}, {}]".format(args.nrho,
                                                             args.ntemp,
                                                             args.nye))
    print("\tUsing [rhomin, rhomax]   = [{}, {}]".format(args.rhomin,
                                                         args.rhomax))
    print("\tUsing [yemin, yemax]     = [{}, {}]".format(args.yemin,
                                                         args.yemax))
    if args.tempmin is not None and args.tempmax is not None:
        print("\tUsing [tmin, tmax]       = [{}, {}]".format(args.tempmin,
                                                             args.tempmax))
        make_table_temp(args.rhomin,args.rhomax,args.nrho,
                        args.tempmin,args.tempmax,args.ntemp,
                        args.yemin,args.yemax,args.nye,
                        units, args.gamma, args.output)
    else:
        print("\tUsing [umin, umax]       = [{}, {}]".format(args.umin,
                                                             args.umax))
        make_table_u(args.rhomin,args.rhomax,args.nrho,
                     args.umin,args.umax,args.ntemp,
                     args.yemin,args.yemax,args.nye,
                     units, args.gamma, args.output)
    print("Done.")
