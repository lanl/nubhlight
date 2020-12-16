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

class OutflowProperties:
    "Container class containing outflow mass and Ye"
    def __init__(self, tracers):
        from units import cgs
        # Mass in solar masses
        code2msolar = tracers.units['M_unit']/cgs['MSOLAR']
        self.mass = tracers['mass'].sum()*code2msolar
        # Mass-averaged Ye
        self.Ye = (tracers['mass']*tracers['Ye'])/self.mass

    def __str__(self):
        return "{:.5e} {.5e}".format(mass, Ye)

def make_table(MBH, a, Md, Ye, s, filenames):
    """Make formatted (ascii) table and return it as a string
    for each simulation with:
    - Black hole mass MBH in solar masses
    - Spin a (unitless)
    - Disk mass Md in solar masses
    - Electron fraction Ye (unitless)
    - Entropy s in K_b/baryon
    - a tracerdata file ending in *.td
    Note some additional parameters that would be good to add:
    - Magnetic field strength
    - Magnetic field topology
    """
    from hdf5_to_dict import TracerData
    assert len(MBH) == len(a) == len(Md) == len(Ye) == len(s) == len(filenames)

    out = """# Columns:
# [0]:  Mass of black hole (solar masses)
# [1]:  Spin of black hole (unitless)
# [2]:  Disk mass at the initial time (solar masses)
# [3]:  Disk Ye at the initial time (unitless)
# [4]:  Disk entropy at initial time (k_b/baryon)
# [5]:  Total mass in outflow (solar masses)
# [6]:  Mass-averaged Ye in total outflow (unitless)
# [7]:  Mass in polar region, >= 50 degress off equator (solar masses)
# [8]:  Mass-averaged Ye in polar region
# [9]:  Mass in equatorial region, <= 15 degrees off equator (solar masses)
# [10]: Mass-averaged Ye in equatorial region
# Note that polar + equatorial != total!
"""
    for i,f in enumerate(filenames):
        tracers = TracerData.fromfile(f) # load tracer data
        # ensure we only look at unbound material
        tracers = tracers.filter(get_bernoulli(tracers) > 0)

        # Get polar and equatorial outflows
        equatorial = filter_angle(tracers,0,15)
        polar = filter_angle(tracers,50,90)

        # Get outflow properties
        prop_tot = OutflowProperties(tracers)
        prop_pol = OutflowProperties(polar)
        prop_equ = OutflowProperties(equatorial)
        
        out += "{:.4e} {:.4e} {:.4e} {:.3e} {:.2e} {} {} {}".format(MBH[i],a[i],
                                                                    Md[i],Ye[i],
                                                                    s[i],
                                                                    prop_tot,
                                                                    prop_pol,
                                                                    prop_equ)

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
    parser.add_argument("-t","--tracers",type=str,nargs="+",required=True,
                        help="files with tracer data")
    args = parser.parse_args()
    table = make_table(args.mbh,args.spin,args.mdisk,args.ye,args.entropy,args.tracers)
    print(table)
