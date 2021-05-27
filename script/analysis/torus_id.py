#!/usr/bin/env python

# ======================================================================
# copyright 2020. Triad National Security, LLC. All rights
# reserved. This program was produced under U.S. Government contract
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
# operated by Triad National Security, LLC for the U.S. Department of
# Energy/National Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is granted
# for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable
# worldwide license in this material to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display publicly, and to
# permit others to do so.
# ======================================================================

"""Computes initial conditions for the Fishbone-Moncrief torus given a
a black hole mass and spin, an equation of state, and a desired disk
mass."""

from units import cgs
import numpy as np
import scipy as sp
from scipy import integrate, optimize

def r_isco(a):
    "Get the radius of the innermost stable circular orbit"
    # Assumes G = c = M = 1
    Z1 = 1 + (1-a*a)**(1./3.) * ((1 + a)**(1./3.) + (1 - a)**(1./3.))
    Z2 = np.sqrt(3*a*a + Z1*Z1)
    return 3 + Z2 - np.sign(a)*np.sqrt((3-Z1)*(3+Z1+2*Z2))

def lfish(a, r):
    """Get the specific angular momentum
    from Fishbone and Moncrief, 1976."""
    from numpy import sqrt
    try:
        return (
      ((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
          ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
                  sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
              ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) * (2. + r))) /
                  sqrt(1 + (2. * a) / pow(r, 1.5) - 3. / r))) /
      (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
          (pow(a, 2) + (-2. + r) * r)))
    except: 
        return 0

def get_lnh(a, rin, rmax, r, th):
    """Get the log of the specific enthalpy 
    from Fishbone, Moncrief, 1976."""
    # Assumes G = c = M = 1
    from numpy import log, sqrt
    
    l = lfish(a, rmax)
    sth = np.sin(th)
    cth = np.cos(th)
    
    DD = r * r - 2. * r + a * a
    AA = ((r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth)
    SS = r * r + a * a * cth * cth
    
    thin = np.pi/2.
    sthin = np.sin(thin)
    cthin = np.cos(thin)
    DDin  = rin * rin - 2. * rin + a * a;
    AAin  = ((rin * rin + a * a) * (rin * rin + a * a) 
             - DDin * a * a * sthin * sthin)
    SSin  = rin * rin + a * a * cthin * cthin
    
    try:
        return (0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
                                        (AA * sth * AA * sth))) /
                    (SS * DD / AA)) -
          0.5 * sqrt(1. + 4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth)) -
          2. * a * r * l / AA -
          (0.5 * log((1. + sqrt(1. + 4. * (l * l * SSin * SSin) * DDin /
                                         (AAin * AAin * sthin * sthin))) /
                     (SSin * DDin / AAin)) -
              0.5 * sqrt(1. + 4. * (l * l * SSin * SSin) * DDin /
                                  (AAin * AAin * sthin * sthin)) -
              2. * a * rin * l / AAin))
    except:
        return 0

def rho_from_hm1_ideal(hm1):
    "Dummy EOS call"
    kappa = 2e-7
    gam = 4./3.
    rho = pow(hm1 * (gam - 1.) / (kappa * gam), 1. / (gam - 1.))
    return rho

def compute_rho_from_coords(a, rin, rmax, eos, r, th):
    "Get the density from Fishbone, Moncrief, 1976."
    lnh = get_lnh(a, rin, rmax, r, th)
    if type(lnh) == np.ndarray:
        lnh[np.isnan(lnh)] = 0.
        lnh[r <= rin] = 0.
        lnh[lnh <= 0] = 0.
    elif np.isnan(lnh) or r <= rin or lnh <= 0:
        lnh = 0.
    hm1 = np.exp(lnh) - 1
    return eos(hm1)

def get_tot_mass(M, a, rin, rmax, eos, get_max = False):
    """Integrate the mass in a fishbone-moncrief disk.
    Returns mass in solar masses."""
    NX1 = 1024
    NX2 = 1024
    X1 = np.linspace(np.log(rin),np.log(max(1e3,2*rmax)), NX1)
    rgrid = np.exp(X1)
    thgrid = np.linspace(0, np.pi/2, NX2)
    RGRID,THGRID = np.meshgrid(rgrid,thgrid,indexing='ij')
    integrand = compute_rho_from_coords(a, rin, rmax, eos,
                                        RGRID, THGRID)

    # Normalization procedure
    rhomax = integrand.max()
    rhomin = integrand.min()
    RHO_unit = rhomax
    L_unit = M*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
    M_unit = RHO_unit*(L_unit**3)
    
    sigma = RGRID**2 + (a*np.cos(THGRID))**2
    delta = RGRID**2 - 2*RGRID + a**2
    detg = np.sqrt(sigma*((a**2 + RGRID**2)*sigma
                          + 2*a*a*RGRID*np.sin(THGRID)**2)/(delta + 1e-20))
    integrand[integrand <= rhomin] = 0
    integrand[integrand <= 1e-2*rhomax] = 0
    integrand *= detg
    integrand /= RHO_unit
    tot = 2*np.pi*2.*integrate.simpson(integrate.simpson(integrand,
                                               x=thgrid,axis=1),
                                       x = rgrid)
    tot *= M_unit/cgs['MSOLAR']
    if get_max:
        return RHO_unit, tot
    else:
        return tot

def solve_for_rmax_rho_unit(M, a, rin, eos, md_target, rmax_guess=None):
    """Given black hole mass, spin, inner disk radius,
    solve for rmax needed to produce a given disk mass.
    Also returns the peak density in the disk."""
    if rmax_guess is None:
        for g in np.linspace(rin,100,20):
            tot = get_tot_mass(M, a, rin, g, eos, False)
            if not np.isnan(tot):
                rmax_guess = g
                break
    if rmax_guess is None:
        raise ValueError("Could not find good initial guess")
    f = lambda rmax: get_tot_mass(M, a, rin, rmax, eos, False) - md_target
    sol = optimize.root_scalar(f, x0 = 1.5*rin, bracket = [rin, 100])
    if not sol.converged:
        print("Failure to converge! Solution object:")
        print(sol)
        raise ValueError("Convergence Failure")
    rmax = sol.root
    rho_max, md_measured = get_tot_mass(M, a, rin, rmax, eos, True)
    return md_measured, sol.root, rho_max

def solve_for_rmax_rho_unit_ideal_gas(M, a, md_target, rin_fac = 1.5,
                                      rmax_guess=None):
    """Sets rin = rin_fac*r_isco. Uses ideal gas."""    
    eos = rho_from_hm1_ideal
    rin = rin_fac*r_isco(a)
    md_measured, rmax, rho_unit = solve_for_rmax_rho_unit(M, a, rin, eos,
                                                          md_target,
                                                          rmax_guess)
    L_unit = M*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
    m_unit = rho_unit*(L*unit**3)
    return md_measured, rin, rmax, rho_unit, m_unit

def solve_for_rmax_rho_unit_tabulated(M, a, md_target, filepath, s, ye,
                                      rin_fac = 1.5,
                                      rmax_guess=None):
    """Sets rin = rin_fac*r_isco. Uses tabulated EOS
    Inputs:
    - M = mass of black hole, in solar masses
    - a = spin of the black hole
    - md_target = target mass of the disk, in solar masses
    - filepath = path to EOS hdf5 file in stelllar collapse format
    - s = entropy of the disk, in kb/baryon
    - ye = electron fraction of the disk
    - rin_fac = rin / r_isco. Recommended minimum value is 1.5

    Outputs:
    - md_measured = disk mass measured, in solar masses
    - rmax = The radius of maximum pressure in the disk, in rg
    - RHO_unit = The code unit for density required for
                 the maximum density in the disk to be unity
    - M_unit = The code unit for mass required for
               the maximum density in the disk to be unit

    Failure case:
    - If the equation of state has an invalid isoentropic curve at this radius,
      this function will raise an "Invalid Adiabat" error.
    """
    import analyze_eos
    eos_obj = analyze_eos.load_eos(filepath)
    adiabat = analyze_eos.Adiabat.get_if_valid(s, ye, eos_obj)
    eos = adiabat.get_rho_of_hm1_interp()
    rin = rin_fac*r_isco(a)
    md_measured, rmax, rho_unit =  solve_for_rmax_rho_unit(M, a, rin, eos,
                                                           md_target, rmax_guess)
    L_unit = M*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
    m_unit = rho_unit*(L_unit**3)
    return md_measured, rin, rmax, rho_unit, m_unit

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description=('Compute the values required for a Fishbone-Moncrief torus'
                     + ' around a black hole.'))
    parser.add_argument('M', type=float,
                        help='Mass of the black hole in solar masses')
    parser.add_argument('a', type=float,
                        help='Spin of the black hole in solar masses')
    parser.add_argument('md', type=float,
                        help='Target mass in the disk in solar masses')
    parser.add_argument('-f','--file',type=str,default=None,
                        help=('Path to stellar collapse file for EOS.'
                              + ' If none, defaults to ideal gas.'))
    parser.add_argument('-s','--entropy',type=float,default=4,
                        help='Entropy in the disk, in kb/baryon. Default is 4.')
    parser.add_argument('-Ye','--electronfraction',type=float,default=0.1,
                        help='Electron fraction in the disk. Default is 0.1.')
    parser.add_argument('-r','--rinfac',type=float,default=1.5,
                        help=('ratio ofinner radius of the disk'
                              + ' to radius of innermost stable circular orbit.'
                              + ' Defaults to 1.5.'))
    args = parser.parse_args()
    if args.file is None:
        print("Computing initial conditions for ideal gas...")
        md_measured, rin, rmax, rho_unit, m_unit \
            = solve_for_rmax_rho_unit_ideal_gas(args.M,
                                                args.a,
                                                args.md,
                                                args.rinfac)
    else:
        print("Computing initial conditions for tabulated eos...")
        md_measured, rin, rmax, rho_unit, m_unit \
            = solve_for_rmax_rho_unit_tabulated(args.M,
                                                args.a,
                                                args.md,
                                                args.file,
                                                args.entropy,
                                                args.electronfraction,
                                                args.rinfac)
    print(("\tmd_measured = {} solar masses\n"
           + "\terror       = {:e}\n"
           + "\trin         = {:f} rg\n"
           + "\trmax        = {:f} rg\n"
           + "\tRHO_unit    = {:e} g/cm^3\n"
           + "\tM_unit      = {:e} g\n").format(md_measured,
                                                md_measured-args.md,
                                                rin,
                                                rmax,
                                                rho_unit,
                                                m_unit))

