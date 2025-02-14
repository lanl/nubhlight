################################################################################
#                                                                              #
#  TORUS, Tuned for Compact Binary Coalescence                                 #
#                                                                              #
################################################################################

# TODO: use argparse

from __future__ import print_function,division

import sys
sys.path.append('../../script')
sys.path.append('../../script/analysis')
sys.dont_write_bytecode = True
import bhlight as bhl
from bhlight import cgs
from units import UnitSystem
from math import log10,ceil
PROB = 'torus_cbc'

DO_GAMMA = '-gamma' in sys.argv # No table
GAMTABLE = '-gamtable' in sys.argv # fake table
RELTABLE = '-reltable' in sys.argv or (not DO_GAMMA or GAMTABLE)
INITIAL = '-initial' in sys.argv # initial data only
NOB = '-nob' in sys.argv # or INITIAL # no B fields
TOROIDALB = "-toroidalb" in sys.argv
RENORM = '-renorm' in sys.argv
NORENORM = '-norenorm' in sys.argv or (not RENORM)
CLASSIC = '-classic' in sys.argv
DEBUG = '-debug' in sys.argv
THREED = '-3d' in sys.argv
FAKETHREED = '-fake3d' in sys.argv
NEUTRINOS = '-nu' in sys.argv
NOEMISS = '-noemiss' in sys.argv
NOSCATT = '-noscatt' in sys.argv
NOABS = '-noabs' in sys.argv
DIAGNOSTIC = '-diag' in sys.argv
NOTRACE = '-notrace' in sys.argv
SMALL = '-small' in sys.argv or NOB
TRACERTEST = '-tracertest' in sys.argv
RESTARTTEST = '-restarttest' in sys.argv
HDF = '-hdf' in sys.argv
OSCILLATIONS = "-oscillations" in sys.argv
N1N2N3CPU_FROM_CLI = '-n1n2n3cpu' in sys.argv
N1N2N3TOT_FROM_CLI = '-n1n2n3tot' in sys.argv

LIMIT_RAD = '-limit' in sys.argv
MORE_RAD = '-morenu' in sys.argv

HPC = '-hpc' in sys.argv # used only for 2d runs
QUAD = '-quad' in sys.argv # quadrant symmetry
KILL = '-kill' in sys.argv # kill all packets

ENT_FROM_CLI  = '-ent' in sys.argv
RHO_FROM_CLI  = '-rho' in sys.argv
RIN_FROM_CLI  = '-rin' in sys.argv
RMAX_FROM_CLI = '-rmax' in sys.argv
MBH_FROM_CLI  = '-M' in sys.argv
ABH_FROM_CLI  = '-a' in sys.argv
YE_FROM_CLI   = '-Ye' in sys.argv
BETA_FROM_CLI = '-beta' in sys.argv
NPH_FROM_CLI  = '-nph' in sys.argv
TFINAL_FROM_CLI = '-tfinal' in sys.argv
N1CPU_FROM_CLI = '-n1cpu' in sys.argv
N2CPU_FROM_CLI = '-n2cpu' in sys.argv
N3CPU_FROM_CLI = '-n3cpu' in sys.argv
N1TOT_FROM_CLI = '-n1tot' in sys.argv
N2TOT_FROM_CLI = '-n2tot' in sys.argv
N3TOT_FROM_CLI = '-n3tot' in sys.argv

EMISS = not NOEMISS
SCATT = not (NOSCATT or KILL)
ABS = not (NOABS or KILL)
FORTRAN = NEUTRINOS and not HDF
TRACERS = not NOTRACE

USE_TABLE = GAMTABLE or RELTABLE
USE_GAMMA = GAMTABLE or not USE_TABLE
RENORM = not NORENORM
TWOD = not THREED

if any([N1CPU_FROM_CLI, N2CPU_FROM_CLI, N3CPU_FROM_CLI]) and not N1N2N3CPU_FROM_CLI:
    raise ValueError("Must turn on the -n1n2n3cpu flag if inputting -n1cpu, -n2cpu, -n3cpu")
if any([N1TOT_FROM_CLI, N2TOT_FROM_CLI, N3TOT_FROM_CLI]) and not N1N2N3TOT_FROM_CLI:
    raise ValueError("Must turn on the -n1n2n3tot flag if inputting -n1tot, -n2tot, -n3tot")

if NOB:
    BFIELD = "none"
elif TOROIDALB:
    BFIELD = "toroidal"
else:
    BFIELD = "classic"

if RELTABLE:
    bhl.report_var('EOS','RELTABLE')
elif GAMTABLE:
    bhl.report_var('EOS','GAMTABLE')
else:
    bhl.report_var('EOS','GAMMA')

if NEUTRINOS:
    RADIATION = 'RADTYPE_NEUTRINOS'
elif TRACERS and not NEUTRINOS:
    RADIATION = 'RADTYPE_LIGHT'
    EMISS = False
    SCATT = False
    ABS = False
else:
    RADIATION = False
    
TABLEPATH = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
#TABLEPATH = "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5"
#TABLEPATH = "HShenEOS_rho220_temp180_ye65_version_1.1_20120817.h5"
TABLEPATH = "../../data/"+TABLEPATH

if FORTRAN:
    OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
else:
    if OSCILLATIONS:
        OPACPATH = "NuLib_rho70_temp62_ye50_ng61_ns4_version1.0_20241120_bhlight.h5"
    else:
        OPACPATH = "opacity.SFHo.nohoro.juo.brem1.h5"

OPACPARAM = "opacbin.LS220.evan.param"
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM


# Parameters from
# arXiv 1409.4426
# arXiv 1705.05473
if MBH_FROM_CLI:
    MBH = float(sys.argv[sys.argv.index('-M') + 1])
else:
    MBH = 2.8
if ABH_FROM_CLI:
    ABH = float(sys.argv[sys.argv.index('-a') + 1])
else:
    ABH = 0.8
if YE_FROM_CLI:
    YE = float(sys.argv[sys.argv.index('-Ye') + 1])
else:
    YE = 0.1
if BETA_FROM_CLI:
    BETA = float(sys.argv[sys.argv.index('-beta') + 1])
else:
    BETA = 1e2

bhl.report_var('MBH',MBH)
bhl.report_var('ABH',ABH)
if USE_TABLE:
    bhl.report_var('Disk Ye',YE)

if ENT_FROM_CLI:
    ENTROPY = float(sys.argv[sys.argv.index('-ent') + 1])
else:
    ENTROPY = 4
bhl.report_var('ENTROPY',ENTROPY)

if CLASSIC: # classic harm disk
    # Rmax and rho fixed
    bhl.report_var('DISK_TYPE','CLASSIC')
    Rin = 6
    Rmax = 12
    #RHO_unit = 5.*10**6
    RHO_unit = 2*10.**12 # 1*10.**8
else:
    # Rin/Rmax chosen so disk mass is 5% Msun
    # RHO_unit chosen so peak desity is unity
    # Varies slightly between 2D and 3D.
    bhl.report_var('DISK_TYPE','CBC')
    if RIN_FROM_CLI:
        Rin = float(sys.argv[sys.argv.index('-rin') + 1])
    else:
        Rin = 3.7
    if RMAX_FROM_CLI:
        Rmax = float(sys.argv[sys.argv.index('-rmax') + 1])
    else:
        Rmax = 9.268 if THREED else 9.03
    if RHO_FROM_CLI:
        RHO_unit = float(sys.argv[sys.argv.index('-rho') + 1])
    else:
        RHO_unit = 6.9*10.**11 if THREED else 3.6*10.**11

# Note that if scattering is enabled,
# you'll get twice as many photons
# per cell as this number
# (unchanged for legacy reasons)
if NPH_FROM_CLI:
    NPH_PER_CELL = float(sys.argv[sys.argv.index('-nph') + 1])
else:
    NPH_PER_CELL = 10

bhl.report_var('Rin',Rin)
bhl.report_var('Rmax',Rmax)
bhl.report_var('RHO_unit',RHO_unit)

# IMPORTANT TO GET RIGHT!
# Chosen so that guess is right at rho = 1.0
# in code units
LRHO_GUESS = log10(RHO_unit)

# Example of how to chose Rmax by physical units:
# Rmax_km = 50.975
# Rmax_cm = 1e5*Rmax_km
# Rmax = Rmax_cm*(cgs['CL']**2)/(cgs['GNEWT']*cgs['MSOLAR']*MBH)

# derived parameters
if TRACERTEST or RESTARTTEST:
    Rout = 250.0
elif THREED or FAKETHREED:
    Rout = 1000.
else:
    Rout = 250.

# use selection instead
Rout_rad = ENTROPY*ceil(Rmax) # not safe to use 3x
Rout_vis = Rout
GAMMA = 13./9.
KAPPA = 1.e-3
L_unit = cgs['GNEWT']*cgs['MSOLAR']*MBH/(cgs['CL']**2)
M_UNIT = RHO_unit*(L_unit**3)

# final time
if TFINAL_FROM_CLI:
    TFINAL = float(sys.argv[sys.argv.index('-tfinal') + 1])
elif NOB and not INITIAL:
    TFINAL = 500.
elif INITIAL:
    TFINAL = 1.e-2
elif THREED or FAKETHREED:
    TFINAL = 10000.
else:
    TFINAL = 2000.
# Toroidal fields grow slowly
if TOROIDALB and not TFINAL_FROM_CLI:
    TFINAL *= 2

# time when radiation resolution controls turn on
t0_tune_emiss = Rout_rad if LIMIT_RAD else -1.
t0_tune_scatt = 2.*max(Rout_rad,t0_tune_emiss)
if LIMIT_RAD:
    tune_emiss = 1.
elif MORE_RAD:
    tune_emiss = 100.
else:
    tune_emiss = 1.

# output
DTd = 1.e-2 if INITIAL else 5.
DTl = 1.e-2 if INITIAL else 5.e-1
DTr = 100
DNr = 50 if RESTARTTEST else 1000

if N1N2N3TOT_FROM_CLI:
    if not all([N1TOT_FROM_CLI, N2TOT_FROM_CLI, N3TOT_FROM_CLI]):
        raise ValueError("If using -n1n2n3tot flag, should input -n1tot -n2tot -n3tot")
    N1TOT = int(sys.argv[sys.argv.index('-n1tot') + 1])
    N2TOT = int(sys.argv[sys.argv.index('-n2tot') + 1])
    N3TOT = int(sys.argv[sys.argv.index('-n3tot') + 1])
else:
    if TRACERTEST or RESTARTTEST:
        N1TOT = 36
        N2TOT = 36
        N3TOT = 24
    elif FAKETHREED:
        N1TOT = 96 if SMALL else 192 # 192
        N2TOT = 96 if SMALL else 128 # 192
        N3TOT = 1
    elif THREED: # 3D
        if QUAD:
            N1TOT = 96 if SMALL else 192 # 192
            N2TOT = 96 if SMALL else 128 # 192
            N3TOT = 16
        else:
            N1TOT = 96 if SMALL else 192 # 192
            N2TOT = 96 if SMALL else 128 # 192
            N3TOT = 64 if SMALL else 66  # 96
    else: # 2D
        N1TOT = 112 if SMALL else 256
        N2TOT = 112 if SMALL else 256
        N3TOT = 1

if N1N2N3CPU_FROM_CLI:
    if not all([N1CPU_FROM_CLI, N2CPU_FROM_CLI, N3CPU_FROM_CLI]):
        raise ValueError("If using -n1n2n3cpu flag, should input -n1cpu -n2cpu -n3cpu")
    N1CPU = int(sys.argv[sys.argv.index('-n1cpu') + 1])
    N2CPU = int(sys.argv[sys.argv.index('-n2cpu') + 1])
    N3CPU = int(sys.argv[sys.argv.index('-n3cpu') + 1])
else:
    if TRACERTEST:
        N1CPU = 1
        N2CPU = 1
        N3CPU = 4
    elif FAKETHREED:
        N1CPU = 1 if RADIATION else 2
        N2CPU = 1 if RADIATION and not HPC else 2
        N3CPU = 1
    elif THREED: # 3D
        if QUAD:
            N1CPU = 1
            N2CPU = 2
            N3CPU = 2
        else:
            N1CPU = 1
            N2CPU = 2  # 2
            N3CPU = 11 # 16
    else: # 2D
        N1CPU = 1 if RADIATION else 2
        N2CPU = 1 if RADIATION and not HPC else 2
        N3CPU = 1

NCPU_TOT     = N1CPU*N2CPU*N3CPU
NCELL_TOT    = N1TOT*N2TOT*N3TOT
NPH_TOT      = NCELL_TOT*NPH_PER_CELL
NPH_PER_PROC = NPH_TOT / NCPU_TOT

NTCR_TOT     = max(32000,NCELL_TOT)
#NTCR_TOT     = max(2000000,NCELL_TOT)
if TRACERS:
    bhl.report_var("NTCR_TOT",NTCR_TOT)
if not NEUTRINOS:
    NPH_PER_PROC = 0.0

EOS = "EOS_TYPE_TABLE" if USE_TABLE else "EOS_TYPE_GAMMA"
NVP0 = 3
NVAR_PASSIVE = NVP0 if USE_TABLE else 0
#GAMMA_FALLBACK = True
GAMMA_FALLBACK = GAMTABLE
#OPENMP = not DEBUG
OPENMP = True

RHOMIN, RHOMAX, NRHO =  1e-8, 1e8, 200
UMIN, UMAX, NU = 1e-8, 1e8, 200
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
print("Make sure to move the generated table to your run directory!")

# Radiation units
NUMIN_MEV = 1
NUMAX_MEV = 300
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
if RADIATION:
    bhl.report_var("[Numin, Numax] (meV)", [NUMIN_MEV,NUMAX_MEV])
    bhl.report_var("[Numin, Numax] (cgs)", [NUMIN,NUMAX])
    bhl.report_var('NPH_PER_CELL', NPH_PER_CELL)
    bhl.report_var('NPH_TOT',NPH_TOT)
    bhl.report_var('NPH_PER_PROC',NPH_PER_PROC)

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
bhl.config.set_cparm('DEREFINE_POLES', THREED or FAKETHREED)

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
bhl.config.set_cparm('QUADRANT_SYMMETRY', QUAD)

# EOS
bhl.config.set_cparm("EOS", EOS)
bhl.config.set_cparm('NVAR_PASSIVE', NVAR_PASSIVE)
bhl.config.set_cparm('GAMMA_FALLBACK', GAMMA_FALLBACK)

# RADIATION
bhl.config.set_cparm('RADIATION', RADIATION)
bhl.config.set_cparm('TRACERS', TRACERS)
bhl.config.set_cparm('EMISSION',   EMISS)
bhl.config.set_cparm('ABSORPTION', ABS)
bhl.config.set_cparm('SCATTERING', SCATT)
if KILL:
    bhl.config.set_cparm('KILL_ALL_PACKETS', True)
bhl.config.set_cparm('BURROWS_OPACITIES', FORTRAN)
bhl.config.set_cparm('HDF5_OPACITIES', HDF)
bhl.config.set_cparm('NU_BINS', 61)
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
bhl.config.set_cparm('DIAGNOSTICS_USE_RADTYPES', True)
bhl.config.set_cparm('RZ_HISTOGRAMS', True)
# bhl.config.set_cparm('RECORD_DT_MIN', True)

if OSCILLATIONS:
    bhl.config.set_cparm('LOCAL_ANGULAR_DISTRIBUTIONS', True)
    bhl.config.set_cparm('LOCAL_ANGLES_NMU', 32)
    bhl.config.set_cparm('LOCAL_ANGLES_NX1', 32)
    bhl.config.set_cparm('LOCAL_ANGLES_NX2', 32)
    bhl.config.set_cparm('RAD_NUM_TYPES', 4)
    bhl.config.set_cparm('NEUTRINO_OSCILLATIONS', True)

# Special. Don't turn this on if you don't need to
if DIAGNOSTIC:
    bhl.config.set_cparm("EXIT_ON_INIT", True)
    bhl.config.set_cparm("OUTPUT_BURROWS_OPACITIES", True)

                           ### RUNTIME PARAMETERS ###

# generic
bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = Rout)
# Maybe this should be something further out. Fine for test.
bhl.config.set_rparm('Rout_rad', 'double', default = Rout_rad)
bhl.config.set_rparm('Rout_vis', 'double', default = Rout_vis)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTr', 'double', default = DTr)
bhl.config.set_rparm('DNr', 'integer', default = DNr)
bhl.config.set_rparm('tune_emiss', 'double', default = tune_emiss)
bhl.config.set_rparm('t0_tune_emiss', 'double', default = t0_tune_emiss)
bhl.config.set_rparm('t0_tune_scatt', 'double', default = t0_tune_scatt) 

# problem
bhl.config.set_rparm('a', 'double', default = ABH)
bhl.config.set_rparm('mbh', 'double', default = MBH)
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm("bfield", 'string', default = BFIELD)
bhl.config.set_rparm("rin", "double", default = Rin)
bhl.config.set_rparm("rmax", "double", default = Rmax)
bhl.config.set_rparm("beta", "double", default = BETA)
bhl.config.set_rparm("renorm_dens", "int", default = int(RENORM))

if RADIATION:
    bhl.config.set_rparm('numin', 'double', default = NUMIN)
    bhl.config.set_rparm('numax', 'double', default = NUMAX)
    bhl.config.set_rparm('nph_per_proc', 'double', default = NPH_PER_PROC)
    bhl.config.set_rparm('rz_rmax', 'double', default = Rout_rad)
    bhl.config.set_rparm('rz_zmax', 'double', default = Rout_rad)

if TRACERS:
    bhl.config.set_rparm('ntracers', 'int', default = NTCR_TOT)

#EOS
if USE_TABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
    bhl.config.set_rparm('const_ye', 'double', YE)
    if RELTABLE:
        bhl.config.set_rparm('entropy','double',ENTROPY)
        bhl.config.set_rparm('lrho_guess','double',LRHO_GUESS)
if USE_GAMMA or GAMMA_FALLBACK:
    bhl.config.set_rparm('gam', 'double', default = GAMMA)
    bhl.config.set_rparm("kappa", "double", default = KAPPA)

# Opacities
if FORTRAN or HDF:
    bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)
if FORTRAN:
    bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)


                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
