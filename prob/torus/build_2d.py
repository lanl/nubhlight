################################################################################
#                                                                              #
#  TORUS                                                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script');
sys.dont_write_bytecode = True; import bhlight as bhl
PROB = 'torus'

# HAVE build.py JUST BE COMPILE-TIME OPTIONS? CALL PYTHON BUILD.PY [prob] FROM
# BASE DIR AND HAVE BIN/ FOLDER BE CREATED? OR TWO FILES IN EACH FOLDER?

# SWITCH TO CPARMS AND RPARMS?

# set_cparm BECOMES cparm? config BECOMES cnfg (ELIMINATE config?)?

# DICTIONARIES cparm, rparm PASSED TO BUILD?

GAMTABLE = '-gamtable' in sys.argv
USE_TABLE = GAMTABLE
USE_GAMMA = GAMTABLE or not USE_TABLE

MBH = 1.e8
ABH = 0.9375
Rout = 1000.

GAMMA = 13./9.
TFINAL = 2000.
DTd = 5.
DTl = 5.e-1
DTr = 100
DNr = 1e10

M_UNIT = 8.e23
EOS = "EOS_TYPE_TABLE" if GAMTABLE else "EOS_TYPE_GAMMA"
NVAR_PASSIVE = 1 if GAMTABLE else 0
GAMMA_FALLBACK = GAMTABLE

RHOMIN, RHOMAX, NRHO =  1e-16, 1e5, 400
UMIN, UMAX, NU = 1e-16, 1e5, 200
YEMIN, YEMAX, NYE = 0.0, 0.6, 50
CRASH_ON_SOUND_SPEED = False

if GAMTABLE:
    sys.path.append('../../script/analysis')
    import make_tabulated_gamma as tab
    tablepath = tab.make_filename(GAMMA)
    units = tab.UnitSystem(M_UNIT, Mbh = MBH)
    tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                     UMIN,   UMAX,   NU,
                     YEMIN,  YEMAX,  NYE,
                     units,  GAMMA,  tablepath,
                     CRASH_ON_SOUND_SPEED)

                         ### COMPILE TIME PARAMETERS ###

NPH_TOT = 1.e5
N1CPU   = 1
N2CPU   = 4
N3CPU   = 1
NTOTCPU = N1CPU*N2CPU*N3CPU

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 96)
bhl.config.set_cparm('N2TOT', 96)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', N1CPU)
bhl.config.set_cparm('N2CPU', N2CPU)
bhl.config.set_cparm('N3CPU', N3CPU)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MKS')
bhl.config.set_cparm('DEREFINE_POLES', False)

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', False)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)
bhl.config.set_cparm('BETA_HEAT', True)
bhl.config.set_cparm('COULOMB', True)

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
bhl.config.set_cparm('OUTPUT_EOSVARS', USE_TABLE)

# RADIATION
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', True)
bhl.config.set_cparm('SYNCHROTRON', True)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('nph_per_proc', 'double', default = NPH_TOT/NTOTCPU)
bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = Rout)
bhl.config.set_rparm('Rout_rad', 'double', default = 40.)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTr', 'double', default = DTr)
bhl.config.set_rparm('DNr', 'integer', default = DNr)
bhl.config.set_rparm('a', 'double', default = ABH)
bhl.config.set_rparm('mbh', 'double', default = MBH)
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm('tune_emiss', 'double', 1.e0)
bhl.config.set_rparm('tune_scatt', 'double', 0.1)
bhl.config.set_rparm('t0_tune_emiss', 'double', 500)
bhl.config.set_rparm('t0_tune_scatt', 'double', 500)
bhl.config.set_rparm('MAD', 'int', default = 0)
bhl.config.set_rparm('BHflux', 'double', default = 0.)

#EOS
if USE_TABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
if USE_GAMMA:
    bhl.config.set_rparm('gam', 'double', default = GAMMA)
                         
                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

