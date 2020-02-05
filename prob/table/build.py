################################################################################
#                                                                              #
#  TABLE READER TEST                                                           #
#                                                                              #
################################################################################

import sys
sys.path.append('../../script/')
sys.path.append('../../script/analysis')
import make_tabulated_gamma as tab
sys.dont_write_bytecode = True
import bhlight as bhl
PROB = 'table'

GAMMA = 1.4
EOS = "EOS_TYPE_TABLE"
NVAR_PASSIVE = 16

# DON'T MAKE NRHO, NU, NYE too big
# or the the linker will complain.
# (or use MPI)
M_UNIT = 1.e2
L_UNIT = 1.e-2
LRHOMIN, LRHOMAX, NRHO =  -1, 2, 100
RHOMIN, RHOMAX = 10**LRHOMIN, 10**LRHOMAX
LUMIN, LUMAX, NU = -7, 0.5, 200
UMIN, UMAX = 10**LUMIN, 10**LUMAX
YEMIN, YEMAX, NYE = 0.0, 0.55, 5
CRASH_ON_SOUND_SPEED = False

# Tabulated EOS
tablepath = tab.make_filename(GAMMA)
units = tab.UnitSystem(M_UNIT, L_unit = L_UNIT)
tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                 UMIN,   UMAX,   NU,
                 YEMIN,  YEMAX,  NYE,
                 units,  GAMMA,  tablepath,
                 CRASH_ON_SOUND_SPEED)

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 3*NRHO)
bhl.config.set_cparm('N2TOT', 3*NU)
bhl.config.set_cparm('N3TOT', 3*NYE)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

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
bhl.config.set_cparm("GAMMA_FALLBACK", True)

# Exit on startup
bhl.config.set_cparm("EXIT_ON_INIT", True)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('x1Min', 'double', default = LRHOMIN)
bhl.config.set_rparm('x1Max', 'double', default = LRHOMAX)
bhl.config.set_rparm('x2Min', 'double', default = LUMIN)
bhl.config.set_rparm('x2Max', 'double', default = LUMAX)
bhl.config.set_rparm('x3Min', 'double', default = YEMIN)
bhl.config.set_rparm('x3Max', 'double', default = YEMAX)

bhl.config.set_rparm('tf',     'double', default = 1e-10)
bhl.config.set_rparm('dt',     'double', default = 0.25e-6)
bhl.config.set_rparm('DTd',    'double', default = 0.25e-1)
bhl.config.set_rparm('DTl',    'double', default = 0.25e-2)
bhl.config.set_rparm('DTr',    'double', default = 12837612)

#EOS
bhl.config.set_rparm('eospath', 'string', default = tablepath)
bhl.config.set_rparm('L_unit',  'double', default = L_UNIT)
bhl.config.set_rparm('M_unit',  'double', default = M_UNIT)
bhl.config.set_rparm('gam', 'double', default = GAMMA)

#PROBLEM
bhl.config.set_rparm("lrhomin", "double", default = LRHOMIN)
bhl.config.set_rparm("lrhomax", "double", default = LRHOMAX)
bhl.config.set_rparm("lumin", "double", default   = LUMIN)
bhl.config.set_rparm("lumax", "double", default   = LUMAX)
bhl.config.set_rparm("yemin", "double", default   = YEMIN)
bhl.config.set_rparm("yemax", "double", default   = YEMAX)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
