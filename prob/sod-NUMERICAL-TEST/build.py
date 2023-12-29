################################################################################
#                                                                              #
#  SOD SHOCKTUBE                                                               #
#                                                                              #
################################################################################

import sys
sys.path.append('../../script/')
sys.dont_write_bytecode = True
import bhlight as bhl
PROB = 'sod-numerical'
MPI = '-mpi' in sys.argv

TABLE = '-table' in sys.argv
RELTABLE = '-reltable' in sys.argv

GAMMA = 1.4
TSCALE = 1.0 if RELTABLE else 1e-2
PSCALE = 0.25 if RELTABLE else 1.
TFINAL = 0.75 if RELTABLE else 0.25

CL = 2.99792458e10
RHO_UNIT = 6.17244*1e17   # c^6 /(G^3 * MSun^2) in cgs (2.99792458*10^10)^6/((6.6742*10^-8)^3*(1.989*10^33)^2)
T_UNIT = 1.e-3

#Grid extension
x1Min = -256.0
x1Max = 256.0
x2Min = -256.0
x2Max = 256.0
x3Min = -256.0
x3Max = 256.0

EOS = "EOS_TYPE_TABLE" if TABLE or RELTABLE else "EOS_TYPE_GAMMA"
NVAR_PASSIVE = 3 if TABLE or RELTABLE else 0
RELTABLENAME = "Hempel_DD2EOS_rho234_temp180_ye60_version_1.1_20120817.h5"
RELTABLEPATH = "../../data/"+RELTABLENAME

RHOMIN, RHOMAX, NRHO =  1e-4, 1e1, 234
UMIN, UMAX, NU = TSCALE/1e5, 1e4, 136
YEMIN, YEMAX, NYE = 0.0, 0.55, 50
CRASH_ON_SOUND_SPEED = False

# Tabulated EOS
if TABLE and RELTABLE:
    raise ValueError("Only -table XOR -reltable may be active at once.")
if TABLE:
    M_UNIT = 1.
    L_UNIT = 1./TSCALE
    sys.path.append('../../script/analysis')
    import make_tabulated_gamma as tab
    tablepath = tab.make_filename(GAMMA)
    units = tab.UnitSystem(M_UNIT, L_unit = L_UNIT)
    tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                     UMIN,   UMAX,   NU,
                     YEMIN,  YEMAX,  NYE,
                     units,  GAMMA,  tablepath,
                     CRASH_ON_SOUND_SPEED)
if RELTABLE:
    tablepath = RELTABLEPATH
    L_UNIT = T_UNIT*CL
    M_UNIT = RHO_UNIT*(L_UNIT**3)

# output and timestep parameters
DT0 = TFINAL/1.e6
DTd = TFINAL/10
DTl = TFINAL/100

FORCE = '-force' in sys.argv

args = [sys.executable, 'build.py', '-dir']
if FORCE:
  args.append('-force')

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 32)
bhl.config.set_cparm('N2TOT', 32)
bhl.config.set_cparm('N3TOT', 32)
bhl.config.set_cparm('N1CPU', 2 if MPI else 1)
bhl.config.set_cparm('N2CPU', 2 if MPI else 1)
bhl.config.set_cparm('N3CPU', 2 if MPI else 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)
#bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'NUMERICAL')
#bhl.config.set_cparm('METRIC', 'MINKOWSKI')

#EXIT AFTER GRID SETTING
bhl.config.set_cparm('EXIT_ON_INIT', True)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', False)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

# EOS
bhl.config.set_cparm("EOS", EOS)
bhl.config.set_cparm('NVAR_PASSIVE', NVAR_PASSIVE)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf',     'double', default = TFINAL)
bhl.config.set_rparm('dt',     'double', default = DT0)
bhl.config.set_rparm('DTd',    'double', default = DTd)
bhl.config.set_rparm('DTl',    'double', default = DTl)
bhl.config.set_rparm('DTr',    'double', default = 12837612)
bhl.config.set_rparm('tscale', 'double', default = TSCALE)
bhl.config.set_rparm('pscale', 'double', default = PSCALE)

bhl.config.set_rparm('x1Min', 'double', default = x1Min)
bhl.config.set_rparm('x1Max', 'double', default = x1Max)
bhl.config.set_rparm('x2Min', 'double', default = x2Min)
bhl.config.set_rparm('x2Max', 'double', default = x2Max)
bhl.config.set_rparm('x3Min', 'double', default = x3Min)
bhl.config.set_rparm('x3Max', 'double', default = x3Max)


#EOS
if TABLE or RELTABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
    bhl.config.set_rparm('RHO_unit',  'double', default = RHO_UNIT)
    
#    bhl.config.set_rparm('L_unit',  'double', default = L_UNIT)
#    bhl.config.set_rparm('M_unit',  'double', default = M_UNIT)
else:
    bhl.config.set_rparm('gam', 'double', default = GAMMA)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
