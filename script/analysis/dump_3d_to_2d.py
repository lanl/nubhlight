#!/usr/bin/env python

################################################################################
#                                                                              # 
#  Maps a 3d dump file to 2d                                                   # 
#                                                                              # 
################################################################################

from __future__ import print_function,division
import sys, os, h5py, re
from argparse import ArgumentParser
import numpy as np
from scipy import integrate,interpolate
import hdf5_to_dict as io
import sadw
from zoh import zoHProjector
SMALL = 1.e-40

parser = ArgumentParser('Map a 3D datafile to 2D via averaging in phi.')
parser.add_argument('dump',type=str,
                    help='Dump file to map.')
parser.add_argument('--rmin',type=float,
                    default=5.5,
                    help='Min radius for Z/H averages')
parser.add_argument('--rmax',type=float,
                    default=25,
                    help='Max radius for Z/H averages')
parser.add_argument('--zmin',type=float,
                    default=-30,
                    help='Min value for z/H')
parser.add_argument('--zmax',type=float,
                    default=30,
                    help='Max value for z/H')
parser.add_argument('-f','--fixghosts',
                    action='store_true',
                    help='Set this flag to set boundary cells by interpolation.')
parser.add_argument('-N2',type=int,
                    default=2,
                    help='Num CPUs in X2. Only relevant if fixing ghosts.')
parser.add_argument('-N3',type=int,
                    default=11,
                    help='Num CPUS in X3. Only relevant if fixing ghosts.')

def new_dumpname(count):
    if count is not None:
        return 'dump2d_%08d.h5' % count
    else:
        return 'dump2d.h5'

def get_data(h5val):
    """Returns data contained in hdf5 dataset.
    If h5val is not an hdf5 dataset, returns None.
    """
    try:
        data = h5val[()]
    except:
        data = None
    return data

def copy_hdr(src,dest,correct_meta=True):
    "Copies header data from h5py object src to h5py object dest."
    hdr_vars = ['N1tot','N2tot','N3tot',
                'metric','electrons','radiation','tracers',
                'output_eosvars',
                'nvar',
                'startx[1]','startx[2]','startx[3]',
                'dx[1]','dx[2]','dx[3]',
                'version',
                'nvar_passive',
                'eos',
                'nulnutype',
                'camera',
                'nth','nphi','nubins_spec',
                'vnams',
                'cour', 'DTd', 'DTf', 'DTl', 'DTp', 'DTr', 'DNr', 'tf',
                'game','gamp',
                'L_unit', 'T_unit', 'M_unit', 'RHO_unit', 'U_unit', 'B_unit',
                'tp_over_te', 'Ne_unit', 'Thetae_unit',
                'maxnscatt', 'nubins', 'numin', 'numax',
                'TEMP_unit',
                'derefine_poles',
                'Rin', 'Rout','Reh', 'Risco', 'hslope', 'a', 'poly_xt',
                'poly_alpha', 'mks_smooth',
                'Rout_vis',
                'Mbh', 'mbh', 'Rout_rad',
                'gam',
                'eospath',
                'poly_K','poly_gam']
    extra_vars = ['t','dt','dump_cnt',
                  'nuLnu','tp_over_te',
                  'failed','full_dump',
                  'gamma_fallback','nstep',
                  'passive_type']
    variables = list(set(hdr_vars + extra_vars))
    for v in variables:            
        if v in src.keys():
            if v not in dest.keys():
                src.copy(v,dest)
    if correct_meta:
        d3x_1d = src['dx[3]'][()]*src['N3tot'][()]
        dest['N3tot'][...] = 1
        dest['dx[3]'][...] = d3x_1d

def get_padded_metric(geom,q):
    gdet = geom['gdet'][:]
    istart = np.where(np.array(q.shape) == gdet.shape[0])[0][0]
    for i in range(istart):
        gdet = gdet[np.newaxis]
    iend = np.where(np.array(q.shape) == gdet.shape[-1])[0][-1]
    for i in range(len(q.shape) - iend -1):
        gdet = gdet[...,np.newaxis]
    return gdet

def put_along_axis(a,indices,values,axis):
    """copied from numpy docs
    https://docs.scipy.org/doc/numpy/reference/generated/numpy.put_along_axis.html
    """
    a = np.asarray(a)
    indices = np.asarray(indices)
    Ni, M, Nk = a.shape[:axis], a.shape[axis], a.shape[axis+1:]

    for ii in np.ndindex(Ni):
        for kk in np.ndindex(Nk):
            a_1d       = a      [ii + np.s_[:,] + kk]
            values_1d  = values [ii + np.s_[:,] + kk]
            a_1d[indices] = values_1d

def gen_idx(q,idx,axis):
    "gives 1d index along a specific axis"
    I = [slice(None)]*q.ndim
    I[axis] = idx
    return tuple(I)

def interp_ghosts(q,hdr,geom,NCPU,axis,direction):
    if direction == 0:
        NTOT = hdr['N1']
        x = geom['X1'][:,0,0]
    elif direction == 1:
        NTOT = hdr['N2']
        x = geom['X2'][0,:,0]
    elif direction == 2:
        NTOT = hdr['N3']
        x = geom['X3'][0,0,:]
    else:
        raise ValueError("invalid direction")
    NLOC = NTOT // NCPU
    qc = q.copy()
    qc[gen_idx(qc,0,axis)] = qc[gen_idx(qc,1,axis)]
    qc[gen_idx(qc,-1,axis)] = qc[gen_idx(qc,-2,axis)]
    for i in range(1,NCPU):
        idx = i*NLOC
        points = x.take(indices=[idx-2,idx+1])
        vals = q.take(indices=[idx-2,idx+1],axis=axis)
        q_interp = interpolate.interp1d(points,vals,axis=axis)
        put_along_axis(qc,range(idx-1,idx+1),q_interp(points),axis)
    return qc

def avg_phi(q,hdr,geom,axis=2):
    "Averages quantity q over azimuth."
    assert hdr['N3'] > 1
    dx = hdr['dx'][3]
    gdet = get_padded_metric(geom,q)
    num = integrate.simps(q*gdet,dx=dx,axis=axis)
    den = integrate.simps(gdet*np.ones_like(q),
                          dx=dx,axis=axis)
    avg = num/(den+SMALL)
    return avg

def avg_dump(src,dest,dump,
             rmin,rmax,zmin,zmax,
             fix_ghosts = False,
             N2CPU=2,N3CPU=11):
    "Copy dump data from src to dest"
    geom = dump['geom']
    hdr = dump['hdr']

    # prune aberrant dtau
    # WARNING: can be dangerous
    THRESH=1e2
    def isnan(a):
        return np.logical_or(np.isnan(a),np.abs(a) >= THRESH)
    if 'dtau_abs' in src.keys():
        dump['dtau_abs'][isnan(dump['dtau_abs'])] = 0.
    if 'dtau_scatt' in src.keys():
        dump['dtau_scatt'][isnan(dump['dtau_scatt'])] = 0.
    if 'dtau_tot' in src.keys():
        dump['dtau_tot'][isnan(dump['dtau_tot'])] = 0.
    if 'dtau_avg' in src.keys():
        dump['dtau_avg'][isnan(dump['dtau_avg'])] = 0.

    # avg in phi
    # ----------
    # Do the special cases first
    # Jrad
    if 'Jrad' in src.keys():
        Jrad = dump['Jrad']
        if fix_ghosts:
            Jrad = interp_ghosts(Jrad,hdr,geom,N2CPU,2,1)
            Jrad = interp_ghosts(Jrad,hdr,geom,N3CPU,3,2)
        Jrad_avg = avg_phi(Jrad,hdr,geom,-1)
        new_shape = list(Jrad.shape)
        new_shape[-1] = 1
        dest.create_dataset('Jrad',
                            data=Jrad_avg,
                            shape=new_shape)

    # Nabs/Nem phys
    for name in ['Nabs_phys','Nem_phys']:
        if name in src.keys():
            var = dump[name]
            if fix_ghosts:
                var = interp_ghosts(var,hdr,geom,N2CPU,1,1)
                var = interp_ghosts(var,hdr,geom,N3CPU,2,2)
            var_avg = avg_phi(var,hdr,geom)
            new_shape = list(var.shape)
            new_shape[2] = 1
            dest.create_dataset(name,
                                data=var_avg,
                                shape=new_shape)

    # dtau
    if 'dtau_avg' in src.keys():
        dtau_avg = dump['dtau_avg']
        dtau_avg_avg = avg_phi(dtau_avg,hdr,geom,-1)
        new_shape = list(dtau_avg.shape)
        new_shape[-1] = 1
        dest.create_dataset('dtau_avg',
                            data = dtau_avg_avg,
                            shape = new_shape)
        for name in ['dtau_tot','dtau_scatt']:
            var = dump[name]
            var_avg = avg_phi(var,hdr,geom)
            new_shape = list(var.shape)
            new_shape[2] = 1
            dest.create_dataset(name,
                                data=var_avg,
                                shape=new_shape)

    # Prim
    if 'P' in src.keys():
        P = src['P']
        P_avg = avg_phi(P,hdr,geom,2)
        new_shape = list(P.shape)
        new_shape[-2] = 1
        dest.create_dataset('P',
                            data=P_avg,
                            shape=new_shape)
        # Attributes need to be copied too!
        for n,v in src['P'].attrs.items():
            dest['P'].attrs.create(n,v)

    # All other cases
    skip = (['Jrad',
             'Nabs_phys','Nem_phys',
             'dtau_avg','dtau_tot','dtau_scat']
            + hdr['vnams'])
    for name,var in src.items():
        if name in dest.keys() or name in skip:
            continue
        var = get_data(var)
        if (var is not None
            and type(var) is np.ndarray
            and len(var.shape) >= 3
            and var.shape[0] == hdr['N1']
            and var.shape[1] == hdr['N2']
            and var.shape[2] == hdr['N3']):
            var_avg = avg_phi(var,hdr,geom)
            new_shape = list(var.shape)
            new_shape[2] = 1
            dest.create_dataset(name,
                                data=var_avg,
                                shape=new_shape)

    # Perform sadws
    variables = list(set(hdr['vnams']
                 + ['PRESS','TEMP','ENT',
                    'Thetae','Theta',
                    'Ye',
                    'bsq','beta','j2','betar']))
    sadw_grp = dest.create_group('sadw')
    sadw_grp.create_dataset('r',data=geom['r'][:,0,0])
    for v in variables:
        if v in dump.keys():
            prof = sadw.get_density_weighted_profile(dump,dump[v])
            sadw_grp.create_dataset(v,data=prof)

    # Perform spherical averages
    variables = list(set(['dtau_tot','dtau_abs','dtau_scatt']))
    sph_grp = dest.create_group('sph_avg')
    sph_grp.create_dataset('r',data=geom['r'][:,0,0])
    for v in variables:
        if v in dump.keys():
            prof = sadw.get_spherical_average(dump,dump[v])
            sph_grp.create_dataset(v,data=prof)
    # treat individual scattering cross-sections separately
    if 'dtau_avg' in dump.keys():
        dtau_avg_profs = [None for i in range(1,dump['dtau_avg'].shape[0])]
        for i in range(1,dump['dtau_avg'].shape[0]):
            dtau_avg_profs[i-1] = sadw.get_spherical_average(dump,
                                                             dump['dtau_avg'][i])
        dtau_avg = np.vstack(dtau_avg_profs)
        sph_grp.create_dataset('dtau_avg',data = dtau_avg)

    # Perform ZoHs
    variables = list(set(['Ye',
                 'dtau_tot','dtau_dens',
                 'dtau_abs','dtau_scatt',
                 'dlepton_rad','dyedt_rad']))
    pzoh = zoHProjector(hdr,geom,dump,(rmin,rmax),(zmin,zmax))
    zoh_grp = dest.create_group('zoh')
    zoh_grp.create_dataset('zoH_uniform',data=pzoh.zoH_uniform)
    zoh_grp.create_dataset('th0',data=pzoh.th0)
    zoh_grp.create_dataset('thd',data=pzoh.thd)
    zoh_grp.create_dataset('H',data=pzoh.H)
    zoh_grp.create_dataset('zoH',data=pzoh.zoH)
    for v in variables:
        if v in dump.keys():
            zoh_grp.create_dataset(v,data=pzoh(dump[v]))
    

def avg_dump_file(infile,rmin,rmax,zmin,zmax,
                  fix_ghosts = False,
                  N2CPU=2,N3CPU=11,
                  geom=None):
    "Average dump data in infile and save to new file."
    dump = io.load_dump(infile,geom=geom)
    if 'dump_id' in dump.keys():
        count = dump['dump_id']
    else:
        numbers = [int(n) for n in re.findall('\d+',infile)]
        if numbers:
            count = max(numbers)
        else:
            count = None
    outfile = new_dumpname(count)
    outfile = os.path.join(os.path.dirname(infile),outfile)
    with h5py.File(infile,'r') as src:
        with h5py.File(outfile,'w') as dest:
            copy_hdr(src,dest)
            avg_dump(src,dest,dump,rmin,rmax,zmin,zmax,
                     fix_ghosts,N2CPU,N3CPU)

if __name__ == "__main__":
    args = parser.parse_args()
    infile = args.dump
    rmin = args.rmin
    rmax = args.rmax
    zmin = args.zmin
    zmax = args.zmax
    N2CPU = args.N2
    N3CPU = args.N3
    fixgh = args.fixghosts
    if not os.path.exists(infile):
        print('ERROR File ' + infile + ' does not exist!')
        sys.exit()
    avg_dump_file(infile,rmin,rmax,zmin,zmax,
                  fix_ghosts = fixgh,
                  N2CPU = N2CPU,N3CPU = N3CPU)
