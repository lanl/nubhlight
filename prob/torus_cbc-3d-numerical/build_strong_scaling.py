################################################################################
#                                                                              #
#  TORUS, Tuned for Compact Binary Coalescence                                 #
#                                                                              #
################################################################################

import sys
sys.path.append('../../script')
sys.path.append('../../script/analysis')
sys.dont_write_bytecode = True
import bhlight as bhl
from bhlight import cgs
from units import UnitSystem
PROB = 'torus_cbc'

GAMTABLE   = False
RELTABLE   = True
INITIAL    = False # initial data only
NOB        = False # no B fields
NORENORM   = True
TWOD       = False
CLASSIC    = '-classic' in sys.argv # classic HARM disk vs. new data
DEBUG      = '-debug' in sys.argv
NONU       = '-nonu' in sys.argv
NOSCATT    = '-noscatt' in sys.argv
NOABS      = '-noabs' in sys.argv

if '-scale' not in sys.argv:
    raise IOError("Must set scale.")
SCALE      = int(sys.argv[sys.argv.index('-scale')+1])
if SCALE > 1 and SCALE % 2 != 0:
    raise ValueError("Scale must be a power of 2 in nodes.")

THREED     = not TWOD
NEUTRINOS  = not NONU
SCATT      = not NOSCATT
ABS        = not NOABS
FORTRAN    = NEUTRINOS
RADIATION  = 'RADTYPE_NEUTRINOS' if NEUTRINOS else False

USE_TABLE = GAMTABLE or RELTABLE
USE_GAMMA = GAMTABLE or not USE_TABLE
BFIELD = not NOB
RENORM = not NORENORM

TABLEPATH = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
TABLEPATH = "../../data/"+TABLEPATH

OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
OPACPARAM = "opacbin.LS220.evan.param"
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM

# Parameters from
# arXiv 1409.4426
# arXiv 1705.05473
MBH = 3.
ABH = 0.8
YE = 0.1
ENTROPY = 15
LRHO_GUESS = 5
LRHO_MIN = 0.0
LRHO_MAX = 14.0
if CLASSIC: # classic harm disk
    Rin = 6
    Rmax = 12
    BETA = 1e2
    RHO_unit = 1*10.**8
else:
    Rin = 3.7
    Rmax_km = 50.975
    Rmax_cm = 1e5*Rmax_km
    Rmax = Rmax_cm*(cgs['CL']**2)/(cgs['GNEWT']*cgs['MSOLAR']*MBH)
    print("Rmax = {}".format(Rmax))
    #Rmax = 9
    BETA = 1e2
    RHO_unit = 1*10.**8
# derived parameters
Rout = 1000. if THREED else 500.
GAMMA = 13./9.
KAPPA = 1.e-3
L_unit = cgs['GNEWT']*cgs['MSOLAR']*MBH/(cgs['CL']**2)
M_UNIT = RHO_unit*(L_unit**3)

# final time
TFINAL = 100.0

# output
#DTd = 0.1
#DTd = 1.e-6
DTd = 5.
DTl = 5.e-1
DTr = 25
DNr = 500

NPH_PER_CELL = 10.
N1TOT = 256
N2TOT = 128
N3TOT = 64
if SCALE == 1: # == number of nodes
    N1CPU = 1
    N2CPU = 1
    N3CPU = 1
if SCALE == 2:
    N1CPU = 1
    N2CPU = 1
    N3CPU = 2
if SCALE == 4:
    N1CPU = 1
    N2CPU = 1
    N3CPU = 4
if SCALE == 8:
    N1CPU = 1
    N2CPU = 2
    N3CPU = 4
if SCALE == 16:
    N1CPU = 1
    N2CPU = 2
    N3CPU = 8
if SCALE == 32:
    N1CPU = 1
    N2CPU = 4
    N3CPU = 8
if SCALE == 64:
    N1CPU = 1
    N2CPU = 8
    N3CPU = 8
if SCALE > 64:
    raise ValueError("Scale not supported.")

NCPU_TOT     = N1CPU*N2CPU*N3CPU
NCELL_TOT    = N1TOT*N2TOT*N3TOT
NPH_TOT      = NCELL_TOT*NPH_PER_CELL
NPH_PER_PROC = NPH_TOT / NCPU_TOT

EOS = "EOS_TYPE_TABLE" if USE_TABLE else "EOS_TYPE_GAMMA"
NVAR_PASSIVE = 3 if USE_TABLE else 0
GAMMA_FALLBACK = GAMTABLE
OPENMP = not DEBUG
#OPENMP = True

RHOMIN, RHOMAX, NRHO =  1e-16, 1e5, 400
UMIN, UMAX, NU = 1e-16, 1e5, 200
YEMIN, YEMAX, NYE = 0.0, 0.6, 50
CRASH_ON_SOUND_SPEED = False

if RELTABLE:
    tablepath = TABLEPATH
if GAMTABLE:
    import make_tabulated_gamma as tab
    tablepath = tab.make_filename(GAMMA)
    units = UnitSystem(M_UNIT, Mbh = MBH)
    tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                     UMIN,   UMAX,   NU,
                     YEMIN,  YEMAX,  NYE,
                     units,  GAMMA,  tablepath,
                     CRASH_ON_SOUND_SPEED)

# Radiation units
NUMIN_MEV = 0.1
NUMAX_MEV = 500
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
if RADIATION:
    print("Numin, Numax = [{}, {}]".format(NUMIN,NUMAX))

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', N1TOT)
bhl.config.set_cparm('N2TOT', N2TOT)
bhl.config.set_cparm('N3TOT', N3TOT)
bhl.config.set_cparm('N1CPU', N1CPU)
bhl.config.set_cparm('N2CPU', N2CPU)
bhl.config.set_cparm('N3CPU', N3CPU)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', OPENMP)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MKS')
bhl.config.set_cparm('DEREFINE_POLES', THREED)

# ELECTRONS. DO NOT USE THESE HERE!
bhl.config.set_cparm('ELECTRONS', False)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

# EOS
bhl.config.set_cparm("EOS", EOS)
bhl.config.set_cparm('NVAR_PASSIVE', NVAR_PASSIVE)
bhl.config.set_cparm('GAMMA_FALLBACK', GAMMA_FALLBACK)

# RADIATION
bhl.config.set_cparm('RADIATION', RADIATION)
bhl.config.set_cparm('EMISSION',   True)
bhl.config.set_cparm('ABSORPTION', ABS)
bhl.config.set_cparm('SCATTERING', SCATT)
bhl.config.set_cparm('BURROWS_OPACITIES', FORTRAN)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('GRAYABSORPTION',  False)
bhl.config.set_cparm('BREMSSTRAHLUNG',  False)
bhl.config.set_cparm('SYNCHROTRON',     False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

# generic
bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = Rout)
bhl.config.set_rparm('Rout_rad', 'double', default = 150.)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTr', 'double', default = DTr)
bhl.config.set_rparm('DNr', 'integer', default = DNr)
bhl.config.set_rparm('tune_emiss', 'double', 1e1)
bhl.config.set_rparm('tune_scatt', 'double', 1e1)
# bhl.config.set_rparm('t0_tune_emiss', 'double', 10.0) 

# problem
bhl.config.set_rparm('a', 'double', default = ABH)
bhl.config.set_rparm('mbh', 'double', default = MBH)
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm("bfield", 'int', default = int(BFIELD))
bhl.config.set_rparm("rin", "double", default = Rin)
bhl.config.set_rparm("rmax", "double", default = Rmax)
bhl.config.set_rparm("beta", "double", default = BETA)
bhl.config.set_rparm("renorm_dens", "int", default = int(RENORM))

if RADIATION:
    bhl.config.set_rparm('numin', 'double', default = NUMIN)
    bhl.config.set_rparm('numax', 'double', default = NUMAX)
    bhl.config.set_rparm('nph_per_proc', 'double', default = NPH_PER_PROC)

#EOS
if USE_TABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
    bhl.config.set_rparm('const_ye', 'double', YE)
    if RELTABLE:
        bhl.config.set_rparm('entropy','double',ENTROPY)
        bhl.config.set_rparm('lrho_guess','double',LRHO_GUESS)
        bhl.config.set_rparm('lrho_min', 'double', LRHO_MIN)
        bhl.config.set_rparm('lrho_max', 'double', LRHO_MAX)
if USE_GAMMA:
    bhl.config.set_rparm('gam', 'double', default = GAMMA)
    bhl.config.set_rparm("kappa", "double", default = KAPPA)

# Opacities
if FORTRAN:
    bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)
    bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
