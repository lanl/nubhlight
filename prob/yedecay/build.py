################################################################################
#                                                                              #
#  NEUTRINO COOLING WITH FLAT EMISSIVITY                                       #
#                                                                              #
################################################################################

import sys
sys.dont_write_bytecode = True
sys.path.append('../../script/')
sys.path.append('../../script/analysis')
import bhlight as bhl
import make_tabulated_gamma as tab
PROB = 'yedecay'

ANTINU = '-antinu' in sys.argv
MPI = '-mpi' in sys.argv
BOOST = '-boost' in sys.argv
LEPTON = '-lepton' in sys.argv
EMISSTYPE='ANTINU_ELECTRON' if ANTINU else 'NU_ELECTRON'
YE0 = 0.0 if ANTINU else 0.5

# Use a fake table for this test
GAMMA = 1.4
TFINAL = 600

CL = 2.99792458e10
RHOMIN, RHOMAX, NRHO =  1e-4, 1e8, 234
UMIN, UMAX, NU = 1e-5, 1e8, 136
YEMIN, YEMAX, NYE = 0.0, 0.6, 50
CRASH_ON_SOUND_SPEED = False


M_UNIT = 1.0
L_UNIT = 1.0
tablepath = tab.make_filename(GAMMA)
units = tab.UnitSystem(M_UNIT, L_unit = L_UNIT)
tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                 UMIN,   UMAX,   NU,
                 YEMIN,  YEMAX,  NYE,
                 units,  GAMMA,  tablepath,
                 CRASH_ON_SOUND_SPEED)

# output and timestep parameters
DT0 = TFINAL/1.e6
DTd = TFINAL/100
DTl = TFINAL/100

# cells and CPUs
NTOT = 4 if MPI else 1
NCPU = 2 if MPI else 1
TUNE_EMISS = 7e-9 if MPI else 1e-9

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', NTOT)
bhl.config.set_cparm('N2TOT', NTOT)
bhl.config.set_cparm('N3TOT', NTOT)
bhl.config.set_cparm('N1CPU', NCPU)
bhl.config.set_cparm('N2CPU', NCPU)
bhl.config.set_cparm('N3CPU', NCPU)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC' if BOOST else 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC' if BOOST else 'BC_OUTFLOW')
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
bhl.config.set_cparm("EOS", "EOS_TYPE_TABLE")
bhl.config.set_cparm('NVAR_PASSIVE', 2)


# RADIATION
bhl.config.set_cparm('RADIATION', 'RADTYPE_NEUTRINOS')
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('FLATEMISS', True)
bhl.config.set_cparm('EXPTAU_WEIGHTS', False)
bhl.config.set_cparm('EMISSTYPE_FLAT', EMISSTYPE)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC' if BOOST else 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC' if BOOST else 'BC_ESCAPE')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_ESCAPE')
if LEPTON:
    bhl.config.set_cparm('COMPLAIN_ON_LEPTON_NONCON', True)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = DT0)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DNr', 'int', default = 1e7)
bhl.config.set_rparm('tune_emiss', 'double', default = TUNE_EMISS)
bhl.config.set_rparm('L_unit', 'double', default = L_UNIT)
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm('eospath', 'string', default = tablepath)
bhl.config.set_rparm('rho0', 'double', default = 1.e14)
bhl.config.set_rparm('uufactor','double', default = 8.e-11)
bhl.config.set_rparm('ye0', 'double', default = YE0)
bhl.config.set_rparm('cnu_flat','double',1)
bhl.config.set_rparm('numax', 'double', default = 7e13)
bhl.config.set_rparm('boost', 'int', 1 if BOOST else 0)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
