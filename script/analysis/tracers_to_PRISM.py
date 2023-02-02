#!/usr/bin/env python
# Little script to convert tracer data to a form suitable for PRISM
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Based on the work of Matt Mumpower

from __future__ import print_function, division
import sys; sys.dont_write_bytecode = True
import units; cgs = units.get_cgs()
import numpy as np
import os, shutil, stat
from itertools import chain
from scipy import integrate,interpolate
from subprocess import Popen
import hdf5_to_dict as io
SMALL = 1e-20

class HomologousExpansion:
    "Assumes everything in CGS!"
    power = -3.0
    gamma = 5./3.
    
    def __init__(self,tlast,rholast,kbTlast,Rlast):
        self.rho0 = self._get_rho0_homologous(rholast,tlast)
        self.C_adiabat = self._get_adiabatic_const(kbTlast,rholast)
        self.R0 = Rlast/tlast

    def rho(self,t):
        rho0 = self.rho0
        p = self.power
        return rho0*((t+SMALL)**(p))

    def R(self,t):
        return self.R0*t
    
    def kbT(self,t):
        rho = self.rho(t)
        return self.kbT_of_rho(rho)
    
    def kbT_of_rho(self,rho):
        C = self.C_adiabat
        gamma = self.gamma
        return C*(rho**(gamma-1.))

    @classmethod
    def _get_rho0_homologous(self,rholast,tlast):
        p = self.power
        return rholast*(tlast**(-p))
    
    @classmethod
    def _get_adiabatic_const(self,kbT,rho):
        gamma = self.gamma
        return kbT*(rho**(1.-gamma))

class TraceExtrapolator:

    t_transition = 0.1

    def __init__(self,trace):
        fv = 'extrapolate'
        kind = 'linear'
        be = False

        # internal arrays
        self._t = trace.units['T_unit']*trace['time']
        self._kbT = trace['T']*cgs['MEV']
        self._rho = trace.units['RHO_unit']*trace['rho']
        self._R = trace['R']*trace.units['L_unit']
        if 'Ye_avg' in trace.keys():
            self._Ye = trace['Ye_avg'][-1]
        else:
            self._Ye = trace['Ye'][-1]
        # interpolators
        self._kbT_int = interpolate.interp1d(self._t,
                                             self._kbT,
                                             kind=kind,
                                             bounds_error=be,
                                             fill_value=fv)
        self._rho_int = interpolate.interp1d(self._t,
                                             self._rho,
                                             kind=kind,
                                             bounds_error=be,
                                             fill_value=fv)
        self._R_int   = interpolate.interp1d(self._t,
                                             self._R,
                                             kind=kind,
                                             bounds_error=be,
                                             fill_value=fv)
        #params for Homologous Expansion
        self.t_last   = self._t[-1]
        self.rho_last = self._rho[-1]
        self.kbT_last = self._kbT[-1]
        self.R_last   = self._R[-1]
        # smoothing region
        self.t_first  = self._t[0]
        self.t_width  = self.t_last - self.t_first
        self.tt_start = (self.t_first
                         + (1. - self.t_transition)*self.t_width)
        self.tt_stop  = self.t_last
        # homologous expansion class
        self.HE = HomologousExpansion(self.t_last,  self.rho_last,
                                      self.kbT_last, self.R_last)

    def Ye(self,t):
        if type(t) is int or type(t) is float:
            return self._Ye
        else:
            return self._Ye*np.ones_like(t)
        
    def rho(self,t):
        return self._smoothed_q(t, self._rho_int, self.HE.rho)

    def kbT(self,t):
        return self._smoothed_q(t, self._kbT_int, self.HE.kbT)

    def R(self,t):
        return self._smoothed_q(t, self._R_int, self.HE.R)

    def _smoothed_q(self,t,q_int,q_HE):
        s = self._get_smoother(t)
        return (1.-s)*q_int(t) + s*q_HE(t)

    def _get_smoother(self,t):
        out = (t - self.tt_start) / (self.tt_stop - self.tt_start)
        out = np.minimum(1.,out)
        out = np.maximum(0.,out)
        return out

def read_files_in_dir(inpath):
    print("Reading directories")
    file_lists = []
    for root,dirs,files in os.walk(inpath):
        files = [os.path.join(root,f) for f in files \
                 if '.td' in f and 'trace_' in f]
        if len(files) > 0:
            file_lists.append(files)
    files = list(chain(*file_lists))
    return files

def copyComplete(source, target):
    # copy content, stat-info (mode too), timestamps...
    shutil.copy2(source, target)
    # copy owner and group
    st = os.stat(source)
    os.chown(target, st[stat.ST_UID], st[stat.ST_GID])

def get_dirname(basename):
    return basename[:-3]

def frdm_comp(Ye, T9, rho, filename, frdm_path):
    '''
    Construct an initial mass fraction file given:
    Ye  - Electron fraction
    T9  - Temperature in GK
    rho - Density in g/cm^3
    '''
    
    copyComplete(os.path.join(frdm_path,'sfho_frdm_comp'),
                 './sfho_frdm_comp')
    copyComplete(os.path.join(frdm_path,'sfho_frdm_comp_v1.03.bin'),
                 './sfho_frdm_comp_v1.03.bin')
    p = Popen(['./sfho_frdm_comp', '-ye', str(Ye), '-t', str(T9),
               '-rho', str(rho), '-f', filename])

    # Run SFHO_FRDM_COMP
    out, err = p.communicate()
    p.wait()
    
    if err != None:
        raise Exception("SFHO_FRDM_COMP: Failed with error:" + err)

    os.remove('sfho_frdm_comp')
    os.remove('sfho_frdm_comp_v1.03.bin')

def signal_smooth(x,window_len=11,window='hanning'):
    """From scipy cookbook:
    https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window;
                    should be an odd integer
        window: the type of window from
                 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be one of: 'flat', 'hanning', "
                         +"'hamming', 'bartlett', 'blackman'")

    if window_len % 2 == 0:
        window_len += 1

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    if len(x) < len(y):
        dlen = len(y) - len(x)
        if dlen % 2 == 0:
            y = y[dlen//2:-dlen//2]
        else:
            y = y[dlen//2:-dlen//2 + 1]
    return y

def log_smooth(x,window_len=11,window='hanning'):
    logx = np.log10(x)
    lxs = signal_smooth(logx,window_len,window)
    xsmooth = 10.**lxs
    if len(x) < len(xsmooth):
        dlen = len(xsmooth) - len(x)
        if dlen % 2 == 0:
            xsmooth = xsmooth[dlen//2:-dlen//2]
        else:
            xsmooth = xsmooth[dlen//2:-dlen//2 + 1]
    return xsmooth

def dydx(x,y):
    dx = x[1:] - x[:-1]
    xmid = 0.5*(x[1:] + x[:-1])
    dy = (y[1:] - y[:-1])/dx
    return xmid,dy

def get_cutoff(trace, p = 6, thresh = 1e-1):
    """Cutoff for where to end trace.
    p and thresh are set by eye.
    """
    tmid,dsdt = dydx(trace['time'],trace['s'])
    Q = signal_smooth(np.abs(dsdt**p))
    Qmid = len(Q) / 2.
    thresh *= 10**-p
    try:
        potential_cutoffs = np.where(Q > thresh)[0]
        cutoff = potential_cutoffs[potential_cutoffs > Qmid][0]
    except:
        cutoff = -1
    return cutoff

def value_crossing(a, v):
    residual = a - v
    first = residual[:-1]
    second = residual[1:]
    return np.logical_and((first*second)<=0,first > 0)
    

def cleanup_trace(trace,
                  T9 = 10, atol = 1,
                  p = 6, thresh = 1e-1):
  """Takes trace and cleans it up so it's suitable for PRISM"""
  # output trace
  trsm,trace_out = {},{}
  Tunit = cgs['MEV']/cgs['GK']
  
  trsm['T'],trsm['rho'] = log_smooth(trace['T']),log_smooth(trace['rho'])

  # starting temperature
  #Tgk   = trace['T']*cgs['MEV']/cgs['GK']
  Tgk   = trsm['T']*cgs['MEV']/cgs['GK']
  TgkT9 = value_crossing(Tgk, T9)
  sidx = np.where(TgkT9)[0][-1] + 1
  
  # Sometimes the temperature is not actually that close to T9=10,
  # so this will do a little interpolation to find a suitable
  # starting point. 
  
  if T9 > Tgk[sidx:][0] > T9-atol or Tgk[0] < T9:
      if Tgk[0] < T9:
          print("WARNING: Tracer {} starting at temperature {} < T9 = {}".format(trace['id'][0], Tgk[0], T9))
          sidx = 0
      for k in trace.keys() - ['T','rho']:
            trace_out[k] = trace[k][sidx:]
      trace_out['T'],trace_out['rho'] = trsm['T'][sidx:],trsm['rho'][sidx:]    
  else:
      interpsidx = interpolate.interp1d(trace['time'][sidx-1:sidx+1],Tgk[sidx-1:sidx+1])
      residual = lambda t: interpsidx(t) - T9
      root_results = optimize.root_scalar(residual, bracket=(trace['time'][sidx-1],trace['time'][sidx+1]))
      if not root_results.converged:
          raise ValueError("Root finding failed for tracer {}".format(trace['id'][0]))
      new_time = root_results.root
      new_val = interpsidx(new_time)
      trace_out['T'] = np.insert(trsm['T'][sidx:],0,new_val/Tunit)
      trace_out['time'] = np.insert(trace['time'][sidx:],0,new_time)
      trace_out['time'] = np.insert(trace['time'][sidx:],0,newx[ind])
      # Still need to get corresponding other key values
      for k in trace.keys() - ['time','T']:
          trace_in = trsm if k in trsm.keys() else trace
          interpsidx = interpolate.interp1d(trace['time'][sidx-1:sidx+1],trace_in[k][sidx-1:sidx+1],axis=0)
          trace_out[k] = np.insert(trace_in[k][sidx:],0,interpsidx(new_time),axis=0)

  # Cut off at the end
  cutoff = get_cutoff(trace_out,p,thresh)
  for k in trace.keys():
      trace_out[k] = trace_out[k][:cutoff]

  # smooth rho and T
  trace_out['Tgk'] = trace_out['T']*cgs['MEV']/cgs['GK']

  # heavies don't matter
  nu_ab = (np.abs(trace_out['rate_absorbed'][:,0])
           + np.abs(trace_out['rate_absorbed'][:,1]) + SMALL/2.)
  nu_em = (np.abs(trace_out['rate_emitted'][:,0])
           + np.abs(trace_out['rate_emitted'][:,1]) + SMALL/2.)
  nu_interact = nu_ab + nu_em

  # Ye
  try:
      ye_cutoff = np.where(nu_interact > np.sqrt(SMALL))[0][-1] + 1
  except:
      ye_cutoff = 0
  ye_avg = np.mean(trace_out['Ye'][ye_cutoff:])
  ye_std = np.std(trace_out['Ye'][ye_cutoff:])
  # trace_out['Ye_avg'] = trace_out['Ye'].copy()
  trace_out['Ye_avg'] = ye_avg*np.ones_like(trace_out['Ye'])
  # trace_out['Ye_avg'] = signal_smooth(trace_out['Ye_avg'],11)

  # R
  trace_out['R'] = np.sqrt(trace_out['Xcart'][:,0]**2
                           + trace_out['Xcart'][:,1]**2
                           + trace_out['Xcart'][:,2]**2)

  trace_out = io.TracerData.fromtuple(None,trace.units,trace_out)
    
  # Return
  return trace_out
    
def to_file(trace,
            filename = 'trajectory.dat'):
    "Assumes file has already been cleaned up"
    t   = trace['time']
    T   = trace['Tgk']
    rho = trace['rho']
    R   = trace['R']
    Ye  = trace['Ye_avg']
    dat = np.vstack((t,T,rho,R,Ye)).T
    tid = trace['id']
    header = ('gw170817 disk tracer {} from JMM: '.format(tid)
              + '1:time[s] 2:T[GK] 3:rho[g/cm3] 4:R[cm] 5:Ye')
    np.savetxt(filename,dat,
               fmt='%.6e',
               header=header,
               comments = '# ')

def meta_file(idx,mass,mass_unit,T_initx,rho_initx,total_mass = None):
    if total_mass is not None:
        mass_fraction = mass / total_mass
    mass *= mass_unit
    with open('metadata.dat','w') as f:
        if total_mass is not None:
            #f.write("# 1:id 2:mass[g] 3:mass fraction\n")
            #f.write("{} {:.6e} {:.6e}\n".format(idx,mass,mass_fraction))
            f.write("# 1:id 2:mass[g] 3:mass fraction 4:T_initx[GK] 5:rho_initx[g/cm3]\n")
            f.write("{} {:.6e} {:.6e} {:.6e} {:.6e}\n".format(idx,mass,mass_fraction,T_initx,rho_initx))
        else:
            f.write("# 1:id 2:mass[g] 3:T_smooth[GK] 4:rho_initx[g/cm3]\n")
            f.write("{} {:.6e} {:.6e} {:.6e}\n".format(idx,mass,T_initx,rho_initx))
            #f.write("# 1:id 2:mass[g]\n")
            #f.write("{} {:.6e}\n".format(idx,mass))            

def extrapolate_trace(trace,T9,atol,p,thresh,tfinal = 1e3):
    # extrapolation
    trace = cleanup_trace(trace,T9,atol,p,thresh)
    extrap = TraceExtrapolator(trace)
    # times
    t0 = trace['time'][0]*trace.units['T_unit']
    t1 = 1.1*trace['time'][-1]*trace.units['T_unit']
    tf = tfinal
    tl1,tlf = np.log10(t1),np.log10(tf)
    tlog = np.linspace(tl1,tlf,1000)
    t_extrap = 10.**tlog 
    tnew = np.concatenate((trace['time']*trace.units['T_unit'],t_extrap))
    # Interpolated output
    trace_out = {}
    trace_out['time'] = tnew
    trace_out['Tgk'] = extrap.kbT(tnew)/cgs['GK']
    trace_out['rho'] = extrap.rho(tnew)
    trace_out['R'] = extrap.R(tnew)/cgs['KM']
    trace_out['Ye_avg'] = extrap.Ye(tnew)
    trace_out['id'] = trace['id'][0]
    trace_out['mass'] = trace['mass'][0]
    trace_out['T_smooth'] = trace['Tgk'][0]
    trace_out['rho_smooth'] = trace['rho'][0]*trace.units['RHO_unit']
    return trace_out

def convert_file(tracer_file,outdir,
                 T9 = 10, atol = 0.3,
                 p = 6, thresh = 1e-1,
                 frdm_path = None,
                 tfinal = 1e3):
    currdir = os.getcwd()
    basename = os.path.basename(tracer_file)
    print("...",basename)
    directory = os.path.join(outdir,basename.rstrip('.td'))
    trace = io.Trace.fromfile(tracer_file)
    mass_unit = trace.units['M_unit']
    trace = extrapolate_trace(trace,T9,atol,p,thresh,
                              tfinal=tfinal)
    if not os.path.exists(directory):
        os.makedirs(directory,exist_ok=True)
    os.chdir(directory)
    to_file(trace)
    #meta_file(trace['id'],trace['mass'],mass_unit,trace['T_smooth'])
    meta_file(trace['id'],trace['mass'],mass_unit,trace['T_smooth'],trace['rho_smooth'])
    if frdm_path is not None:
        Ye = trace['Ye_avg'][0]
        #rho = trace['rho'][0]
        rhosmooth = trace['rho_smooth']
        Tsmooth = trace['T_smooth']
        #frdm_comp(Ye,T9,rho,'initx.dat',frdm_path)
        frdm_comp(Ye,Tsmooth,rhosmooth,'initx.dat',frdm_path)
        print('initx file written.')
    os.chdir(currdir)
    return

def convert_file_worker(inp):
    tracer_file = inp[0]
    outdir = inp[1]
    T9 = inp[2]
    atol = inp[3]
    p = inp[4]
    thresh = inp[5]
    frdm_path = inp[6]
    tfinal = inp[7]
    convert_file(tracer_file,outdir,
                 T9,atol,
                 p,thresh,
                 frdm_path,
                 tfinal)

def convert_ids(tracers,directory,
                ids = None,
                T9 = 10, atol = 0.3,
                p = 6, thresh = 1e-1,
                frdm_path=None,
                tfinal = 1e3):
    print("Generating PRISM files for i,id...")
    if ids is None:
        ids = tracers.ids()
    currdir = os.getcwd()
    initial = tracers.filter(tracers['time'] == tracers.times()[0])
    mass_unit = tracers.units['M_unit']
    total_mass = np.sum(initial['mass'])
    if not os.path.exists(directory):
        os.makedirs(directory,exist_ok=True)
    os.chdir(directory)
    for i,idx in enumerate(ids):
        print("\t...",i,idx)
        trace = extrapolate_trace(trace,T9,atol,p,thresh,
                                  tfinal)
        istr = "{:08d}".format(i)
        if not os.path.exists(istr):
            os.mkdir(istr)
        os.chdir(istr)
        to_file(trace)
        meta_file(idx,trace['mass'],mass_unit,total_mass)
        if frdm_path is not None:
            Ye = trace['Ye_avg'][0]
            rho = trace['rho'][0]
            frdm_comp(Ye,T9,rho,'initx.dat',frdm_path)
        os.chdir('..')
    os.chdir(currdir)

def convert_path_by_time(inpath,outdir,
                         T9 = 10, atol = 0.3,
                         p = 6, thresh = 1e-1,
                         frdm_path = None):
    tracers = io.TracerData.frompath(inpath)
    convert_ids(tracers,outdir,
                T9=T9,atol=atol,
                p=p,thresh=thresh,
                frdm_path=frdm_path)

def convert_path_by_files(inpath,outdir,
                          T9 = 10, atol = 0.3,
                          p = 6, thresh = 1e-1,
                          frdm_path = None,
                          tfinal = 1e3,
                          serial=False):
    print("Reading files in path ",inpath)
    files = read_files_in_dir(inpath)
    print("We have {} files to convert".format(len(files)))
    if not files:
        raise ValueError("No files to read!")
    if serial:
        print("...Running in serial")
        for f in files:
            convert_file(f,outdir,
                         T9,atol,
                         p,thresh,
                         frdm_path,
                         tfinal)
    else:
        print("...Running in parallel")
        from concurrent import futures
        inp = [(f,outdir,T9,atol,p,thresh,frdm_path,tfinal) \
               for f in files]
        with futures.ProcessPoolExecutor() as p:
            p.map(convert_file_worker,inp)
    print("Done!")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Convert tracer data to PRISM format',
                            epilog=('Tracer data saved '
                                    +'in one directory per trace'
                                    +'Method will fail if tracers '
                                    +'are not unbound.'))
    parser.add_argument('tracers',type=str,
                        help=('Tracer data to convert. '
                              +'Can be directory or *.td file.'))
    parser.add_argument('outdir',type=str,
                        help='Root directory to output PRISM data')
    parser.add_argument('-f','--frdm',type=str,
                        default=None,
                        help=('If you want an initx file with compositions, '
                              +'you must provide path to directory containing '
                              +'frdm code to extract it from EOS. '
                              +'Default is None, in which case no initx files '
                              +'are generated.'))
    parser.add_argument('-T9','--Tstart',type=float,
                        default=10,
                        help=('Temperature in GK at which traces start. '
                              +'Called T9 in PRISM. Default is 10.'))
    parser.add_argument('-a','--atol',type=float,
                        default=0.3,
                        help=('How close to Tstart we want to start. '
                              +'Default is 0.3'))
    parser.add_argument('-p','--power',type=float,
                        default=6.,
                        help=('Enhancement power factor '
                              +'for entropy trace stop criterion. '
                              +'Default is 6'))
    parser.add_argument('-t','--thresh',type=float,
                        default=1e-1,
                        help=('Threshold for enhanced entropy trace '
                              +'stop criterion. Default is 0.1.'))
    parser.add_argument('--tfinal',type=float,
                        default=1e3,
                        help='Final time to extrapolate to.')
    parser.add_argument('-s','--single',
                        dest='single',
                        action='store_true',
                        help='load tracers by id instad of by time')
    parser.add_argument('--serial',dest='serial',
                        action='store_true',
                        help=('Load data in serial. '
                              +'Only used if running in single id mode.'))
    parser.add_argument('--preserve',dest='preserve',
                        action='store_true',
                        help='Attempt to preserve input directory structure')
    args = parser.parse_args()

    print("PRISM converter")

    if args.frdm is not None:
        if not os.path.exists(args.frdm):
            raise IOError("frdm path does not exist")
        if not os.path.isdir(args.frdm):
            raise IOError("frdm path should be directory containing frdm")

    if not args.single:
        print("...Reading time dumps")
        print("...Reading tracer data")
        convert_path_by_time(args.tracers, args.outdir,
                             args.Tstart,  args.atol,
                             args.power,   args.thresh,
                             args.frdm,
                             tfinal=args.tfinal)
    else:
        print("...Reading single tracer files")
        if not os.path.isdir(args.tracers):
            raise ValueError("Tracers should be directory in this mode.")
        if args.preserve:
            subdirs = chain(*[d for r,d,f in os.walk(args.tracers)])
            for subdir in subdirs:
                inpath = os.path.join(args.tracers,subdir)
                outdir = os.path.join(args.outdir,subdir)
                if not os.path.exists(outdir):
                    os.makedirs(outdir,exist_ok=True)
                convert_path_by_files(inpath,outdir,
                                      args.Tstart,args.atol,
                                      args.power,args.thresh,
                                      args.frdm,
                                      tfinal=args.tfinal,
                                      serial=args.serial)
        else:
            convert_path_by_files(args.tracers, args.outdir,
                                  args.Tstart,  args.atol,
                                  args.power,   args.thresh,
                                  args.frdm,
                                  tfinal=args.tfinal,
                                  serial=args.serial)
                

