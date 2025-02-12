################################################################################
#                                                                              #
#  CONFIGURATION AND COMPILATION ROUTINE                                       #
#                                                                              #
################################################################################

import sys
import os
import util
from subprocess import call
import subprocess

PARAM_NAME   = 'param_template.dat'

HOST_OPTIONS = ['NAME',
                'COMPILER',
                'COMPILER_FLAGS',
                'DEBUG_FLAGS',
                'GSL_DIR',
                'DEBUG_FLAGS',
                'EXECUTABLE',
                'MPI_EXECUTABLE']

if sys.version_info <= (3,0):
  util.warn("ONLY TESTED FOR PYTHON 3.x")

CPARMS = {}
RPARMS = {}
REPVARS = {}

# CONTAINER CLASS FOR RUNTIME PARAMETER
class rparm:
  datatype = None
  value = None

# NICE FORMATTING
def print_config(key, var):
  print("    " + util.color.BOLD + "{:<15}".format(key) + util.color.NORMAL +
        str(var))

# SET COMPILE-TIME PARAMETER
def set_cparm(name, value):
  CPARMS[name] = value

def set_cparm_if_active(name, value):
  if util.parm_is_active(CPARMS, name):
    print_config(name + ' ', CPARMS[name])
  else:
    set_cparm(name, value)
  return

# SET RUNTIME PARAMETER. DO NOT OVERWRITE DEFAULT VALUES, AS THIS IS CALLED BY
# PROBLEM FILE BEFORE CORE ROUTINE
def set_rparm(name, datatype, default=None):
  if datatype == 'int':
    datatype = 'integer'
  
  if datatype != 'integer' and datatype != 'double' and datatype != 'string':
    util.warn('DATATYPE ' + datatype + ' NOT SUPPORTED')
    print('CHOICES: integer double string')
    sys.exit()
  
  obj = rparm()
  obj.datatype = datatype

  if name in RPARMS:
    obj.value = RPARMS[name].value
  else:
    obj.value = default
  RPARMS[name] = obj

def write_rparm(pf, name):
  if name not in RPARMS:
    util.warn('RUNTIME PARAMETER ' + name + ' NOT SET')
    sys.exit()

  datatype = RPARMS[name].datatype
  default = RPARMS[name].value
  if datatype == 'integer':
    if default == None:
      pf.write('[int] ' + name + ' = \n')
    else:
      pf.write('[int] ' + name + ' = %d\n' % default)
  elif datatype == 'double':
    if default == None:
      pf.write('[dbl] ' + name + ' = \n')
    else:
      pf.write('[dbl] ' + name + ' = %e\n' % default)
  elif datatype == 'string':
    if default == None:
      pf.write('[str] ' + name + ' = \n')
    else:
      pf.write('[str] ' + name + ' = %s\n' % default)
  else:
    print(name)
    util.warn("DATATYPE " + str(datatype) + " NOT RECOGNIZED")
    sys.exit()

  del RPARMS[name]

def report_var(name,value):
  REPVARS[name] = value

def is_user_debug():
  return '-debug' in sys.argv

def is_user_force():
  return '-force' in sys.argv

def is_user_noparam():
  return '-noparam' in sys.argv

def is_user_noclean():
  return '-noclean' in sys.argv

def is_user_help():
  return '-help' in sys.argv

def get_version():
    versionfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               '..','VERSION')
    with open(versionfile,'r') as f:
        version = f.read().rstrip().lstrip()
    return version

def set_dirs(PATHS):
  MOVEEXEC = True
  PATHS['BUILD'] = 'SRC'
  PATHS['PROB'] = os.getcwd()
  for n in range(len(sys.argv)):
    if sys.argv[n] == '-dir' and n < len(sys.argv) - 1:
      MOVEEXEC = False
      PATHS['BUILD'] = util.sanitize_path(sys.argv[n+1])
  PATHS['SRC'] = os.path.join(PATHS['BUILD'], 'source','')

  for key in list(PATHS):
    if PATHS[key][0] != '/':
      PATHS[key] = os.path.join(PATHS['PROB'], PATHS[key], '')
    util.make_dir(PATHS[key])

  return MOVEEXEC

def get_env_name(name:str):
  return 'BHLIGHT_{}'.format(name)

def get_host(PATHS):

  machines = util.get_files(PATHS['MACHINE'], '*')
  for n in range(len(machines)):
    machines[n] = machines[n].split('/')[-1].replace('.py', '')
    try:
      machine = __import__(machines[n])
    except ImportError:
      continue
    if machine.matches_host() == True:
      break
    del machine
    
  try: machine
  except NameError: util.warn("HOST " + os.uname()[1] + " UNKNOWN"); sys.exit()
  
  host = machine.get_options()

  # overwrite with environment variable if available
  for name in HOST_OPTIONS:
    env_name = get_env_name(name)
    if env_name in os.environ:
      util.gentle_warn(
        '{} environment variable overriding {} in machine file'.format(env_name,
                                                                       name)
      )
      host[name] = os.environ[env_name]

  return host

def build(PROBLEM, PATHS):
  print("")
  print("********************************************************************************")
  print("")
  print("                              BHLIGHT BUILD SCRIPT")
  print("")
  print("  OPTIONS:")
  print("    -help                (print this message and exit)")
  print("    -debug               (use debug compile params)")
  print("    -force               (do not abort upon compile-time warnings/errors)")
  print("    -noclean             (do not delete old source files)")
  print("    -noparam             (do not create new parameter file)")
  print("    -dir /path/to/target (create target dir and build there)")
  print("")
  print("********************************************************************************")
  print("")

  if is_user_help():
    sys.exit()

  # PROCESS USER INPUT
  DEBUG = is_user_debug()
  FORCE = is_user_force()
  MOVEEXEC = set_dirs(PATHS)
  NOPARAM = is_user_noparam()
  NOCLEAN = is_user_noclean()
  CLEAN   = not NOCLEAN
  WRITE_PARAM = not NOPARAM
  KEEP_SRC = '-src' in sys.argv
  NOKEEP_SRC = '-nosrc' in sys.argv
  REMOVE_SRC = NOKEEP_SRC

  # get version
  VERSION = get_version()

  # PRINT TO TERMINAL AND LOGFILE
  LOGFILE = os.path.join(PATHS['BUILD'], 'log_build')
  util.log_output(sys, LOGFILE)

  # SEARCH FOR MACHINE
  host = get_host(PATHS)

  if DEBUG and 'DEBUG_FLAGS' not in host.keys():
    util.warn("Debug compiler options not set! Using normal compiler flags.")
    host['DEBUG_FLAGS'] = host['COMPILER_FLAGS']

  C_FLAGS = '-std=c99 '
  if 'MEM_MODEL' in host:
    if host['MEM_MODEL']: # empty strings and False bools do not add
                          # this flag
      C_FLAGS += "-mcmodel=" + host['MEM_MODEL'] + ' '
  else:
    C_FLAGS += "-mcmodel=medium "
  if DEBUG:
    C_FLAGS += host['DEBUG_FLAGS']
  else:
    C_FLAGS += host['COMPILER_FLAGS']

  # MATH AND DYNAMIC LINKING
  LIB_FLAGS = '-lm -ldl'

  LIBRARIES = ''
  INCLUDES  = ''

  USE_RPATH = host.get('USE_RPATH', True)

  # GSL
  LIB_FLAGS += ' -lgsl -lgslcblas'
  if 'GSL_DIR' in host: 
    host['GSL_DIR'] = util.sanitize_path(host['GSL_DIR'])
    if USE_RPATH:
      LIB_FLAGS += (' -Wl,-rpath='
                    + host['GSL_DIR'] + 'lib/')
    LIBRARIES += '-L' + host['GSL_DIR'] + 'lib/'
    INCLUDES  += '-I' + host['GSL_DIR'] + 'include/'

  # MPI
  if 'MPI_DIR' in host:
    LIB_FLAGS += ' -lmpi'
    host['MPI_DIR'] = util.sanitize_path(host['MPI_DIR'])
    if USE_RPATH:
      LIB_FLAGS += (' -Wl,-rpath='
                    + host['MPI_DIR'] + 'lib/')
    LIBRARIES += ' -L' + host['MPI_DIR'] + 'lib/'
    INCLUDES  += ' -I' + host['MPI_DIR'] + 'include/'

  # HDF5
  if 'HDF5_DIR' in host:
    LIB_FLAGS += ' -lhdf5_hl -lhdf5'
    host['HDF5_DIR'] = util.sanitize_path(host['HDF5_DIR'])
    if USE_RPATH:
      LIB_FLAGS += (' -Wl,-rpath='
                    + host['HDF5_DIR'] + 'lib/')
    LIBRARIES += ' -L' + host['HDF5_DIR'] + 'lib/'
    INCLUDES  += ' -I' + host['HDF5_DIR'] + 'include/'

  if 'EXTRA_INCLUDES' in host:
    INCLUDES += (" " + host['EXTRA_INCLUDES'])
  if 'EXTRA_LIBRARIES' in host:
    LIBRARIES += (" " + host['EXTRA_LIBRARIES'])

  print("  CONFIGURATION\n")

  set_cparm("VERSION", '"{}"'.format(VERSION))

  print_config("VERSION",    VERSION)
  print_config("MACHINE",    host['NAME'])
  print_config("PROBLEM",    PROBLEM)
  print_config("BUILD DIR",  PATHS['BUILD'])
  print_config("COMPILER",   host['COMPILER'])
  print_config("GSL_DIR",    host['GSL_DIR'])
  if 'MPI_DIR' in host:
    print_config("MPI_DIR",  host['MPI_DIR'])
  if 'HDF5_DIR' in host:
    print_config("HDF5_DIR", host['HDF5_DIR'])
  if 'EXECUTABLE' in host:
    print_config("EXECUTABLE",     host['EXECUTABLE'])
  if 'MPI_EXECUTABLE' in host:
    print_config("MPI_EXECUTABLE", host['MPI_EXECUTABLE'])
  print_config("C_FLAGS",    C_FLAGS)
  print_config("LIB_FLAGS",  LIB_FLAGS)
  print_config("LIBRARIES",  LIBRARIES)
  print_config("INCLUDES",   INCLUDES)
  print_config("OPENMP", CPARMS['OPENMP'])

  print("\n  COMPILE-TIME PARAMETERS\n")
  print_config("N1TOT", CPARMS['N1TOT'])
  print_config("N2TOT", CPARMS['N2TOT'])
  print_config("N3TOT", CPARMS['N3TOT'])
  print_config("N1CPU", CPARMS['N1CPU'])
  print_config("N2CPU", CPARMS['N2CPU'])
  print_config("N3CPU", CPARMS['N3CPU'])
  print_config("METRIC", CPARMS['METRIC'])
  print_config("RECONSTRUCTION", CPARMS['RECONSTRUCTION'])
  if util.parm_is_active(CPARMS, 'EOS'):
    print_config("EOS", CPARMS['EOS'])
  else:
    set_cparm("EOS", 'EOS_TYPE_GAMMA')
  if util.parm_is_active(CPARMS, 'RADIATION'):
    print_config("RADIATION", CPARMS['RADIATION'])
    if util.parm_is_active(CPARMS, 'NU_BINS'):
      print_config('NU_BINS', CPARMS['NU_BINS'])
    else:
      set_cparm("NU_BINS", 200)
    if util.parm_is_active(CPARMS, 'NTH'):
      print_config("NTH", CPARMS["NTH"])
    else:
      set_cparm("NTH", 8)
    if util.parm_is_active(CPARMS, 'NPHI'):
      print_config("NPHI", CPARMS["NPHI"])
    else:
      set_cparm("NPHI", 8)
    if util.parm_is_active(CPARMS, "NU_BINS_SPEC"):
      print_config("NU_BINS_SPEC", CPARMS["NU_BINS_SPEC"])
    else:
      set_cparm("NU_BINS_SPEC", 200)
    if util.parm_is_active(CPARMS, "BURROWS_OPACITIES"):
      print_config("BURROWS_OPACITIES", CPARMS["BURROWS_OPACITIES"])
    else:
      set_cparm("BURROWS_OPACITIES", 0)
    if util.parm_is_active(CPARMS, "HDF5_OPACITIES"):
      print_config("HDF5_OPACITIES", CPARMS["HDF5_OPACITIES"])
    else:
      set_cparm("HDF5_OPACITIES", 0)
    if util.parm_is_active(CPARMS, "RAD_NUM_TYPES"):
      print_config("RAD_NUM_TYPES", CPARMS["RAD_NUM_TYPES"])
    else:
      if CPARMS['RADIATION'] == 1:
        set_cparm("RAD_NUM_TYPES", 1)
      else:
        set_cparm("RAD_NUM_TYPES", 3)
      print_config("RAD_NUM_TYPES", CPARMS["RAD_NUM_TYPES"])
    if util.parm_is_active(CPARMS, "NEUTRINO_OSCILLATIONS"):
      print_config("NEUTRINO_OSCILLATIONS ", CPARMS["NEUTRINO_OSCILLATIONS"])
      if util.parm_is_active(CPARMS, 'FORCE_EQUIPARTITION'):
        print_config("FORCE_EQUIPARTITION ", CPARMS["FORCE_EQUIPARTITION"])
      else:
        set_cparm('FORCE_EQUIPARTITION', 0)
    else:
      set_cparm('NEUTRINO_OSCILLATIONS', 0)
      set_cparm('FORCE_EQUIPARTITION', 0)
    if util.parm_is_active(CPARMS, "RZ_HISTOGRAMS"):
      print_config("RZ_HISTOGRAMS", CPARMS["RZ_HISTOGRAMS"])
      if not util.parm_is_active(CPARMS, "RZ_HISTOGRAMS_N"):
        set_cparm("RZ_HISTOGRAMS_N", 32)
      print_config("RZ_HISTOGRAMS_N", CPARMS["RZ_HISTOGRAMS_N"])
    else:
      set_cparm("RZ_HISTOGRAMS", 0)
  else:
    set_cparm("RADIATION", 0)
    set_cparm("RAD_NUM_TYPES", 0)
  if util.parm_is_active(CPARMS, 'ELECTRONS'):
    print_config("ELECTRONS", CPARMS['ELECTRONS'])
  else:
    set_cparm("ELECTRONS", 0)
  if util.parm_is_active(CPARMS,'NVAR_PASSIVE'):
    print_config("NVAR_PASSIVE", CPARMS["NVAR_PASSIVE"])
  else:
    set_cparm("NVAR_PASSIVE", 0)
  if util.parm_is_active(CPARMS, 'GAMMA_FALLBACK'):
    print_config('GAMMA_FALLBACK', CPARMS['GAMMA_FALLBACK'])
  else:
    set_cparm('GAMMA_FALLBACK', 0)
  if util.parm_is_active(CPARMS, 'EXIT_ON_INIT'):
    print_config('EXIT_ON_INIT', CPARMS['EXIT_ON_INIT'])
  else:
    set_cparm('EXIT_ON_INIT', 0)  
  if util.parm_is_active(CPARMS, 'OUTPUT_EOSVARS'):
    print_config("OUTPUT_EOSVARS", CPARMS["OUTPUT_EOSVARS"])
  else:
    set_cparm("OUTPUT_EOSVARS", 0)
  if CPARMS['RADIATION']:
    if 'EXPTAU_WEIGHTS' in CPARMS.keys():
      print_config('EXPTAU_WEIGHTS', CPARMS['EXPTAU_WEIGHTS'])
    else:
      set_cparm('EXPTAU_WEIGHTS', 1)

  if util.parm_is_active(CPARMS,'ELECTRONS')\
     and CPARMS['EOS'] != 'EOS_TYPE_GAMMA':
    raise ValueError("ELECTRONS only compatible with Gamma law EOS.\n"
                     +"Please set EOS = EOS_TYPE_GAMMA.\n")

  if CPARMS['EOS'] == 'EOS_TYPE_TABLE' and CPARMS['NVAR_PASSIVE'] < 2:
    raise ValueError("Tabulated EOS requires at least two passive scalars\n"
                     +"for the electron fraction Ye.\n"
                     +"Please set NVAR_PASSIVE >= 2\n"
                     +"and ensure your problem generator sets it appropriately.\n")
  if CPARMS['EOS'] == 'EOS_TYPE_TABLE' \
     and CPARMS['METRIC'] == 'MKS' \
     and CPARMS['NVAR_PASSIVE'] < 3:
    raise ValueError("Tabulated EOS and MKS metric requires at least three\n"
                     +"passive scalars, for Ye and atmosphere markers.\n"
                     +"Please set NVAR_PASSIVE >= 3\n"
                     +"and ensure your problem generator sets it appropriately.\n")
  if util.parm_is_active(CPARMS,'ESTIMATE_THETAE') \
     and (CPARMS['RADIATION'] == 'RADTYPE_NEUTRINOS'):
    raise ValueError("Neutrinos not compatible "
                     +"with estimating electron temperature.")
  if CPARMS['EOS'] == 'EOS_TYPE_POLYTROPE':
    util.gentle_warn("The polytropic EOS is totally untested. "
                     +"Use at your own risk!\n")

  if CPARMS['EOS'] != 'EOS_TYPE_GAMMA' and not CPARMS['OUTPUT_EOSVARS']:
    util.gentle_warn("Setting OUTPUT_EOSVARS = True.")
    set_cparm("OUTPUT_EOSVARS", 1)
    print_config("OUTPUT_EOSVARS", CPARMS["OUTPUT_EOSVARS"])

  if util.parm_is_active(CPARMS, 'RADIATION'):
    if util.parm_is_active(CPARMS, 'LOCAL_ANGULAR_DISTRIBUTIONS'):
      print_config('LOCAL_ANGULAR_DISTRIBUTIONS ',
                   CPARMS['LOCAL_ANGULAR_DISTRIBUTIONS'])
      # TODO(JMM): What should defaults be?
      set_cparm_if_active('LOCAL_ANGLES_NMU', 64)
      set_cparm_if_active('LOCAL_ANGLES_NX1', 64)
      set_cparm_if_active('LOCAL_ANGLES_NX2', 64)
    else:
      set_cparm('LOCAL_ANGULAR_DISTRIBUTIONS', 0)

  NEED_UNITS = (util.parm_is_active(CPARMS, 'RADIATION') or
                util.parm_is_active(CPARMS, 'COULOMB') or
                CPARMS['EOS'] == 'EOS_TYPE_TABLE')

  if util.parm_is_active(CPARMS, 'FLATEMISS'):
    if not util.parm_is_active(CPARMS, 'EMISSTYPE_FLAT'):
      util.gentle_warn("Flatemiss active, but not emission type.\n"
                       +"Setting EMISSTYPE_FLAT = ANTINU_ELECTRON.\n")
      set_cparm("EMISSTYPE_FLAT", "ANTINU_ELECTRON")

  if util.parm_is_active(CPARMS, 'RADIATION') \
     and 'X1R_RAD_BOUND' in CPARMS.keys():
    if CPARMS['X1R_RAD_BOUND'] == 'BC_CAMERA' \
       and CPARMS['METRIC'] == 'MINKOWSKI':
      util.warn("X1R_RAD_BOUND BC_CAMERA is "
                +"not supported for Minkowski metrics.")

  print("\n  EXTRA PARAMETERS\n")
  for k,v in REPVARS.items():
    print_config(k,v)

  # Set core runtime parameters
  set_rparm('tf', 'double')
  set_rparm('dt', 'double')
  if CPARMS['METRIC'] == 'MINKOWSKI':
    set_rparm('x1Min', 'double', default = 0.)
    set_rparm('x1Max', 'double', default = 1.)
    set_rparm('x2Min', 'double', default = 0.)
    set_rparm('x2Max', 'double', default = 1.)
    set_rparm('x3Min', 'double', default = 0.)
    set_rparm('x3Max', 'double', default = 1.)
  if CPARMS['METRIC'] == 'MKS':
    set_rparm('a', 'double', default = 0.5)
    set_rparm('hslope', 'double', default = 0.3)
    set_rparm('poly_xt', 'double', default = 0.82)
    set_rparm('poly_alpha', 'double', default = 14.)
    set_rparm('mks_smooth', 'double', default = 0.5)
    set_rparm('Rout', 'double', default = 40.)
    set_rparm('Rout_vis', 'double', default = 40.)
    #if util.parm_is_active(CPARMS, 'RADIATION'):
    #  set_rparm('Rout_rad', 'double')

  if NEED_UNITS:
    if CPARMS['METRIC'] == 'MINKOWSKI':
      set_rparm('L_unit', 'double')
      set_rparm('M_unit', 'double')
    if CPARMS['METRIC'] == 'MKS':
      set_rparm('M_unit', 'double')
      set_rparm('mbh', 'double', default = 1.989e34)
  
  if CPARMS['EOS'] == 'EOS_TYPE_GAMMA':
    set_rparm('gam', 'double', default = 5./3.)

  set_rparm('cour', 'double', default = 0.9)
  if util.parm_is_active(CPARMS, 'RADIATION'):
    set_rparm('cour_cool', 'double', default = 0.25)
  
  if util.parm_is_active(CPARMS, 'ELECTRONS'):
    set_rparm('game', 'double', default = 4./3.)
    set_rparm('gamp', 'double', default = 5./3.)
    set_rparm('fel0', 'double', default = 0.01)
    set_rparm('tptemin', 'double', default = 1.e-3)
    set_rparm('tptemax', 'double', default = 1.e3)
  
  if util.parm_is_active(CPARMS, 'RADIATION'):
    if not util.parm_is_active(CPARMS, 'ELECTRONS'):
      set_rparm('tp_over_te', 'double', default = 1.)
    set_rparm('nph_per_proc', 'double', default = 1.e5)
    set_rparm('numin', 'double', default = 1.e8)
    set_rparm('numax', 'double', default = 1.e20)
    set_rparm('tune_emiss', 'double', default = 1.)
    set_rparm('tune_scatt', 'double', default = 1.)
    set_rparm('t0_tune_emiss', 'double', default = -1.)
    set_rparm('t0_tune_scatt', 'double', default = -1.)
    set_rparm('thetae_max', 'double', default = 1.e3)
    set_rparm('sigma_max', 'double', default = 1.)
    set_rparm('kdotk_tol', 'double', default = 1.e-6)
    set_rparm('Nph_to_track', 'double', default = 0.);
    set_rparm('thbin', 'int', default = 8);
    set_rparm('phibin', 'int', default = 8);
    if util.parm_is_active(CPARMS, 'FLATEMISS'):
      set_rparm('cnu_flat', 'double', default = 1.)
    if util.parm_is_active(CPARMS, "RZ_HISTOGRAMS"):
      if CPARMS['METRIC'] == 'MINKOWSKI':
        set_rparm('rz_rmax', 'double', default = 1)
        set_rparm('rz_zmax', 'double', default = 1)
      if CPARMS['METRIC'] == 'MKS':
        set_rparm('rz_rmax', 'double', default = 40)
        set_rparm('rz_zmax', 'double', default = 40)
    
  set_rparm('init_from_grmhd', 'string', default = 'No')
  
  set_rparm('DTd', 'double', default = 0.5)
  set_rparm('DTl', 'double', default = 0.1)
  set_rparm('DTr', 'double', default = 1000.)
  set_rparm('DNr', 'integer', default = 1024)
  set_rparm('DTp', 'integer', default = 100)
  set_rparm('DTf', 'integer', default = 1)
  set_rparm('outputdir', 'string', default = './')

  # some fixes for fortran compilation
  USE_FORTRAN = util.parm_is_active(CPARMS,'BURROWS_OPACITIES')
  FORT_USE_MPI = (CPARMS['N1CPU'] > 1
                  or CPARMS['N2CPU'] > 1
                  or CPARMS['N3CPU'] > 1)
  FORT_USE_MPI_STR = 'TRUE' if FORT_USE_MPI else 'FALSE'
  if CPARMS['N1TOT'] > 1 and CPARMS['N2TOT'] > 1 and CPARMS['N2TOT'] > 1:
    FORT_NDIM=3
  elif ((CPARMS['N1TOT'] > 1 and CPARMS['N2TOT'] > 1)
        or (CPARMS['N1TOT'] > 1 and CPARMS['N3TOT'] > 1)
        or (CPARMS['N2TOT'] > 1 and CPARMS['N3TOT'] > 1)):
    FORT_NDIM=2
  else:
    FORT_NDIM=1
  if USE_FORTRAN:
    # -lgfortran for gcc -lifcore -limf for icc
    LIB_FLAGS += ' ' + host['FORTLINK']
    if len(host['FORTLIB']) > 0:
        LIBRARIES += ' -L' + host['FORTLIB']
    if DEBUG:
      if 'FDEBUG_FLAGS' not in host.keys():
        util.warn("Fortran debug options not set! Using normal fortran flags.")
        host['FDEBUG_FLAGS'] = host['FCFLAGS']
      FCFLAGS = host['FDEBUG_FLAGS']
    else:
      FCFLAGS = host['FCFLAGS']
    FCFLAGS += (' -DUSE_MPI=' + FORT_USE_MPI_STR
                + ' -DNDIM=' + str(FORT_NDIM))

  # GET ALL SOURCE FILES
  SRC_CORE = util.get_files(PATHS['CORE'], '*.c')
  INC_CORE = util.get_files(PATHS['CORE'], '*.h')
  if USE_FORTRAN:
    F90_CORE = util.get_files(PATHS['CORE'], '*.f90')
  else:
    F90_CORE = []
  SRC_PROB = util.get_files(PATHS['PROB'], '*.c')
  INC_PROB = util.get_files(PATHS['PROB'], '*.h')

  # Clean if necessary
  if CLEAN:
    util.make_clean(PATHS['SRC'])

  # COPY SOURCE FILES TO BUILD_DIR
  for src in SRC_CORE:
    call(['cp', src, PATHS['SRC'] + src.rsplit('/',1)[1]])
  for inc in INC_CORE:
    call(['cp', inc, PATHS['SRC'] + inc.rsplit('/',1)[1]])
  if USE_FORTRAN:
    for src in F90_CORE:
      call(['cp', src, PATHS['SRC'] + src.rsplit('/',1)[1]])
  for src in SRC_PROB:
    call(['cp', src, PATHS['SRC'] + src.rsplit('/',1)[1]])
  for inc in INC_PROB:
    call(['cp', inc, PATHS['SRC'] + inc.rsplit('/',1)[1]])

  # WRITE PARAMETERS FILE
  pf = open(PATHS['SRC'] + 'params.h', 'w')
  for KEY in CPARMS:
    if isinstance(CPARMS[KEY], str):
      pf.write("#define " + KEY + " (" + CPARMS[KEY] + ")\n")
    else: # True/False autocast to 1/0.
      pf.write("#define " + KEY + " (%g)\n" % CPARMS[KEY])
  pf.close()

  # GET SINGLE LISTS OF ALL SOURCE, OBJECT, AND HEADER FILES
  SRC_ALL = util.get_files(PATHS['SRC'], '*.c')
  INC_ALL = util.get_files(PATHS['SRC'], '*.h')
  SRC = ''
  OBJ = ''
  INC = ''
  for n in range(len(SRC_ALL)):
    SRC += '%s ' % os.path.basename(SRC_ALL[n])
    OBJ += '%s.o ' % os.path.basename(SRC_ALL[n])[:-2]
  for n in range(len(INC_ALL)):
    INC += '%s ' % os.path.basename(INC_ALL[n])
  if USE_FORTRAN:
    F90 = ''
    FOBJ = ''
    for src in F90_CORE:
      F90 += '%s ' % os.path.basename(src)
      FOBJ += '%s.o ' % os.path.basename(src)[:-4]

  # WRITE MAKEFILE
  os.chdir(PATHS['SRC'])
  mf = open('makefile', 'w')
  mf.write('CC = ' + host['COMPILER'] + '\n')
  if USE_FORTRAN:
    mf.write('F90 = ' + host['FORTRAN_COMP'] + '\n')
  mf.write('CCFLAGS = ' + C_FLAGS + ' ' + LIBRARIES + ' ' + INCLUDES + '\n')
  if USE_FORTRAN:
    mf.write('FCFLAGS = ' + FCFLAGS + '\n')
  mf.write('LIB_FLAGS = ' + LIB_FLAGS + '\n')
  mf.write('CC_COMPILE = $(CC) $(CCFLAGS) -c' + '\n')
  mf.write('CC_LOAD = $(CC) $(CCFLAGS)' + '\n')
  if USE_FORTRAN:
    mf.write('FSRC = ' + F90 + '\n')
    mf.write('FOBJ = ' + FOBJ + '\n')
  else:
    mf.write('FSRC = \n')
    mf.write('FOBJ = \n')
  mf.write('SRC = ' + SRC + '\n')
  mf.write('OBJ = ' + OBJ + '\n')
  mf.write('INC = ' + INC + '\n')
  mf.write('EXE = bhlight' + '\n')
  mf.write('.c.o:' + '\n')
  mf.write('\t$(CC_COMPILE) $*.c' + '\n')
  if USE_FORTRAN:
    mf.write('%.o: %.f90 makefile' + '\n')
    mf.write('\t$(F90) $(FCFLAGS) -c $<\n')
  mf.write('all: $(EXE)' + '\n')
  mf.write('$(OBJ): $(INC) makefile' + '\n')
  mf.write('$(EXE): $(OBJ) $(FOBJ) $(INC) makefile' + '\n')
  mf.write('\t$(CC_LOAD) $(OBJ) $(FOBJ) $(LIB_FLAGS) -o $(EXE)\n')
  mf.write('clean:\n')
  mf.write('\t$(RM) $(SRC) $(FSRC) $(OBJ) $(FOBJ) $(EXE) $(INC)\n')
  mf.close()

  print("\n  COMPILING SOURCE\n")

  ncomp = 0
  first_error = 1
  if DEBUG:
    popen = subprocess.Popen(['make'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
  else:
    popen = subprocess.Popen(['make','-j','10'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)

  for stdout_line in iter(popen.stdout.readline, ""):
    if stdout_line.rstrip()[-2:] == '.c' or stdout_line.rstrip()[-4:] == '.f90':
      print("    [" + util.color.BOLD + util.color.BLUE +
            "%2d%%" % (100.*float(ncomp)/len(SRC_ALL+F90_CORE)) + util.color.NORMAL +
            "] " + util.color.BOLD +
            stdout_line.rsplit(' -c ',1)[1].rstrip().lstrip().split('/')[-1] +
            util.color.NORMAL)
      ncomp += 1
  for stderr_line in iter(popen.stderr.readline, ""):
    # THIS ALSO FAILS FOR WARNINGS!!!
    if first_error == 1:
      util.warn("COMPILER ERROR")
      first_error = 0
    print(stderr_line.rstrip())

  if first_error != 1 and not FORCE:
    util.warn("COMPILATION FAILED")
    sys.exit()

  obj_files = util.get_files(PATHS['SRC'], '*.o')
  for f in obj_files:
    os.remove(f)
  os.rename(PATHS['SRC'] + 'bhlight', PATHS['BUILD'] + 'bhlight')

  if REMOVE_SRC:
    import shutil
    shutil.rmtree(PATHS['SRC'])

  print("\n  BUILD SUCCESSFUL")

  # CREATE RUNTIME PARAMETERS FILE
  PARAMFILE = PATHS['BUILD'] + PARAM_NAME
  if WRITE_PARAM:
    with open(PARAMFILE, 'w') as pf:
      pf.write("### RUNTIME PARAMETERS ###\n")
      pf.write("\n# COORDINATES\n")
      write_rparm(pf, 'tf')
      write_rparm(pf, 'dt')

      if CPARMS['METRIC'] == 'MINKOWSKI':
        write_rparm(pf, 'x1Min')
        write_rparm(pf, 'x1Max')
        write_rparm(pf, 'x2Min')
        write_rparm(pf, 'x2Max')
        write_rparm(pf, 'x3Min')
        write_rparm(pf, 'x3Max')
      if CPARMS['METRIC'] == 'MKS':
        write_rparm(pf, 'Rout')
        if util.parm_is_active(CPARMS, 'RADIATION'):
          write_rparm(pf, 'Rout_rad')

      if NEED_UNITS:
        pf.write("\n# UNITS\n")
        if CPARMS['METRIC'] == 'MINKOWSKI':
          write_rparm(pf, 'L_unit')
          write_rparm(pf, 'M_unit')
        if CPARMS['METRIC'] == 'MKS':
          write_rparm(pf, 'mbh')
          write_rparm(pf, 'M_unit')

      pf.write("\n# FLUID\n")
      write_rparm(pf, 'cour')
      if util.parm_is_active(CPARMS, 'RADIATION'):
        write_rparm(pf, 'cour_cool')
      if CPARMS['EOS'] == 'EOS_TYPE_GAMMA':
        write_rparm(pf, 'gam')
      if CPARMS['EOS'] == 'EOS_TYPE_TABLE':
        write_rparm(pf, 'eospath')

      if util.parm_is_active(CPARMS, 'ELECTRONS'):
        pf.write("\n# ELECTRONS\n")
        write_rparm(pf, 'game')
        write_rparm(pf, 'gamp')
        write_rparm(pf, 'fel0')
        write_rparm(pf, 'tptemin')
        write_rparm(pf, 'tptemax')

      if util.parm_is_active(CPARMS, 'RADIATION'):
        pf.write("\n# RADIATION\n")
        if not util.parm_is_active(CPARMS, 'ELECTRONS'):
          write_rparm(pf, 'tp_over_te')
        write_rparm(pf, 'nph_per_proc')
        write_rparm(pf, 'numin')
        write_rparm(pf, 'numax')
        write_rparm(pf, 'tune_emiss')
        write_rparm(pf, 'tune_scatt')
        write_rparm(pf, 't0_tune_emiss')
        write_rparm(pf, 't0_tune_scatt')
        write_rparm(pf, 'thetae_max')
        write_rparm(pf, 'sigma_max')
        write_rparm(pf, 'kdotk_tol')
        write_rparm(pf, 'Nph_to_track')
        write_rparm(pf, 'thbin')
        write_rparm(pf, 'phibin')
        if util.parm_is_active(CPARMS, 'BURROWS_OPACITIES'):
          write_rparm(pf, 'opac_param_file')
          write_rparm(pf, 'opac_file')
        if util.parm_is_active(CPARMS, 'HDF5_OPACITIES'):
          write_rparm(pf, 'opac_file')
        write_rparm(pf, 'init_from_grmhd')

      pf.write("\n# OUTPUT\n")
      write_rparm(pf, 'DTd')
      write_rparm(pf, 'DTl')
      write_rparm(pf, 'DTr')
      write_rparm(pf, 'DNr')
      write_rparm(pf, 'DTp')
      write_rparm(pf, 'DTf')
      write_rparm(pf, 'outputdir')
      if len(RPARMS.keys()) > 0:
        pf.write("\n# PROBLEM\n")
        prob_keys = RPARMS.keys()
        for key in list(prob_keys):
          write_rparm(pf, key)
    print("\n  RUNTIME PARAMETER FILE CREATED")

  if MOVEEXEC:
    os.rename(PATHS['BUILD'] + 'bhlight',   PATHS['BUILD'] + '../bhlight')
    if WRITE_PARAM:
      os.rename(PATHS['BUILD'] + PARAM_NAME, PATHS['BUILD'] + '../' + PARAM_NAME)

  print("")

  sys.exit()

