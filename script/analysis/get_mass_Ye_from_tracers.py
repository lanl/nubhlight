#!/usr/bin/env python

################################################################################
# Author: Jonah Miller (jonahm@lanl.gov)
#
# PURPOSE: Takes a *.td file containing tracer data and extracts the
# mass and Ye in the outflow in polar and equatorial bins and total
################################################################################

import numpy as np

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

def filter_angle(tracers,thmin,thmax):
    """Returns tracers with thmin <= |theta_bl - 90| <= thmax
    thmin and thmax assumed to be in degrees
    """
    theta = get_theta(tracers)
    mask = np.logical_and(180*np.abs(theta)/np.pi >= thmin,
                          180*np.abs(theta)/np.pi <= thmax)
    return tracers.filter(mask)

def get_bernoulli(tracers):
    "Returns the Bernoulli parameter as defined in Novikov and Thorne"
    if 'bcon' in tracers.keys():
        bsq = (tracers['bcon']*tracers['bcov']).sum(axis=-1)
    else:
        bsq = np.zeros_like(tracers['rho'])
    h = (1
         + tracers['uu']/tracers['rho']
         + tracers['Press']/tracers['rho']
         + bsq/tracers['rho']
         + bsq/tracers['rho'])
    Be = -tracers['ucov'][:,0]*h - 1
    return Be

def get_intersection(tracers,tracers_ref):
    "Returns particles in tracers that also have ids in tracers_ref"
    mask = np.in1d(tracers['id'],tracers_ref['id'])
    return tracers.filter(mask)

class OutflowProperties:
    "Container class containing outflow mass and Ye"
    def __init__(self, tracers, m_acc):
        from units import cgs
        # Mass in solar masses
        mass = tracers['mass'].sum()
        code2msolar = tracers.units['M_unit']/cgs['MSOLAR']
        self.mass = mass*code2msolar
        # mass divided by total mass accreted
        self.frac_mass = self.mass / m_acc
        # Mass-averaged Ye
        self.Ye = (tracers['mass']*tracers['Ye']).sum()/mass

    def __str__(self):
        return "{:.5e} {:.5e} {:.5e}".format(self.mass, self.frac_mass, self.Ye)

def make_table(MBH, a, Md, Ye, s, m_acc, sphere_files, nse_files):
    """Make formatted (ascii) table and return it as a string
    for each simulation with:
    - Black hole mass MBH in solar masses
    - Spin a (unitless)
    - Disk mass Md in solar masses
    - Electron fraction Ye (unitless)
    - Entropy s in K_b/baryon
    - m_acc is the amount of mass accreted by the end of the simulation
    - sphere_files are files of tracers extracted at a radius of 250 rg
    - nse_files are files of tracers extracted when they drop below 5GK
    Note some additional parameters that would be good to add:
    - Magnetic field strength
    - Magnetic field topology
    """
    from hdf5_to_dict import TracerData
    assert len(sphere_files) == len(nse_files) == len(m_acc)
    assert len(MBH) == len(a) == len(Md) == len(Ye) == len(s) == len(nse_files)

    out = """# Columns:
# [0]:  Mass of black hole (solar masses)
# [1]:  Spin of black hole (unitless)
# [2]:  Disk mass at the initial time (solar masses)
# [3]:  Disk Ye at the initial time (unitless)
# [4]:  Disk entropy at initial time (k_b/baryon)
# [5]:  Total mass in outflow (solar masses)
# [6]:  Ratio of total mass in outflow to mass accreted
# [7]:  Mass-averaged Ye in total outflow (unitless)
# [8]:  Extrapolated total mass in outflow to late times
# [9]:  Mass in polar region, >= 50 degress off equator (solar masses)
# [10]: Ratio of mass in polar region to mass accreted
# [11]: Mass-averaged Ye in polar region
# [12]: Extrapolated mass in polar region
# [13]: Mass in equatorial region, <= 15 degrees off equator (solar masses)
# [14]: Ratio of mass in equatorial region to mass accreted
# [15]: Mass-averaged Ye in equatorial region
# [16]: Extrapolated mass in equatorial region
# Note that polar + equatorial != total!
"""
    for i,(f_sph,f_nse) in enumerate(zip(sphere_files,nse_files)):
        tracers = TracerData.fromfile(f_sph) # load tracer data
        # ensure we only look at unbound material
        tracers = tracers.filter(get_bernoulli(tracers) > 0)

        # Get polar and equatorial outflows
        equatorial = filter_angle(tracers,0,15)
        polar = filter_angle(tracers,50,90)

        # get the nse versions
        tracers_nse = TracerData.fromfile(f_nse)
        tracers_nse = get_intersection(tracers_nse,tracers)
        equatorial_nse = get_intersection(tracers_nse,equatorial)
        polar_nse = get_intersection(tracers_nse, polar)

        # Get outflow properties
        prop_tot = OutflowProperties(tracers_nse,m_acc[i])
        prop_pol = OutflowProperties(polar_nse,m_acc[i])
        prop_equ = OutflowProperties(equatorial_nse,m_acc[i])
        
        out += "{:.3e} {:.3e} {:.3e} {:.2e} {:.2e} {} {:.3e} {} {:.3e} {} {:.3e}\n".format(
            MBH[i],a[i],
            Md[i],Ye[i],
            s[i],
            prop_tot,
            prop_tot.frac_mass*Md[i],
            prop_pol,
            prop_pol.frac_mass*Md[i],
            prop_equ,
            prop_equ.frac_mass*Md[i])

    return out
                                                                    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(
        description="Make table of outflow properties based on tracer data and print to screen"
    )
    parser.add_argument("--mbh",type=float,nargs="+",required=True,
                        help="black hole masses for each simulation (msolar)")
    parser.add_argument("-a","--spin",type=float,nargs="+",required=True,
                        help="black hole spins for each simulation (unitless")
    parser.add_argument("--mdisk",type=float,nargs="+",required=True,
                        help="disk masses for each simulation (msolar)")
    parser.add_argument("--ye",type=float,nargs='+',required=True,
                        help="initial Ye for each simulation (unitless)")
    parser.add_argument("-s","--entropy",type=float,nargs='+',required=True,
                        help="initial entropy for each simulation (k_b/baryon)")
    parser.add_argument("--macc",type=float,nargs='+',required=True,
                        help="mass accreted by end of simulation (msolar)")
    parser.add_argument("--tracers",type=str,nargs="+",required=True,
                        help="files with tracer data when tracers pass through extraction sphere")
    parser.add_argument("--nse",type=str,nargs="+",required=True,
                        help="files with tracer data when tracers fall below 5GK")
    args = parser.parse_args()
    table = make_table(args.mbh,args.spin,args.mdisk,args.ye,args.entropy,
                       args.macc,args.tracers,args.nse)
    print(table)
