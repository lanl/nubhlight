################################################################################
#                                                                              #
#  3D LINEAR RMHD MODES                                                        #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True
import bhlight as bhl
PROB = 'mhdmodes3d'
MPI = '-mpi' in sys.argv

# Tabulated EOS stuff
GAMMA = 4./3.
TABLE = '-table' in sys.argv
EOS = "EOS_TYPE_TABLE" if TABLE else "EOS_TYPE_GAMMA"
NVAR_PASSIVE = 2 if TABLE else 0
M_UNIT = 1.
L_UNIT = 1.e-2
RHOMIN, RHOMAX, NRHO =  1e-7, 1e2, 234
UMIN, UMAX, NU = 1e-7, 1e2, 136
YEMIN, YEMAX, NYE = 0.0, 0.55, 50
CRASH_ON_SOUND_SPEED = False
if TABLE:
    sys.path.append('../../script/analysis')
    import make_tabulated_gamma as tab
    tablepath = tab.make_filename(GAMMA)
    units = tab.UnitSystem(M_UNIT, L_unit = L_UNIT)
    tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                     UMIN,   UMAX,   NU,
                     YEMIN,  YEMAX,  NYE,
                     units,  GAMMA,  tablepath,
                     CRASH_ON_SOUND_SPEED)

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 64)
bhl.config.set_cparm('N2TOT', 64)
bhl.config.set_cparm('N3TOT', 64)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_PERIODIC')
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

                           ### RUNTIME PARAMETERS ###

# TFINAL AND DTd ARE HARDCODED
bhl.config.set_rparm('tf', 'double', default = 5.)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('DTd', 'double', default = 5.e-1)
bhl.config.set_rparm('DTl', 'double', default = 5.e-1)
bhl.config.set_rparm('DTr', 'integer', default = 10000)
bhl.config.set_rparm('nmode', 'integer', default = 1)

#EOS
if TABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
    bhl.config.set_rparm('L_unit',  'double', default = L_UNIT)
    bhl.config.set_rparm('M_unit',  'double', default = M_UNIT)
else:
    bhl.config.set_rparm('gam', 'double', default = GAMMA)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

