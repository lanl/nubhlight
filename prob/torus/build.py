################################################################################
#                                                                              #
#  TORUS                                                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script');
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'torus'

# HAVE build.py JUST BE COMPILE-TIME OPTIONS? CALL PYTHON BUILD.PY [prob] FROM
# BASE DIR AND HAVE BIN/ FOLDER BE CREATED? OR TWO FILES IN EACH FOLDER?

# SWITCH TO CPARMS AND RPARMS?

# set_cparm BECOMES cparm? config BECOMES cnfg (ELIMINATE config?)?

# DICTIONARIES cparm, rparm PASSED TO BUILD?

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 96)
bhl.config.set_cparm('N2TOT', 96)
bhl.config.set_cparm('N3TOT', 96)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MKS')
bhl.config.set_cparm('DEREFINE_POLES', True)

# ELECTRONS
bhl.config.set_cparm('ELECTRONS', False)
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', False)
bhl.config.set_cparm('BETA_HEAT', True)
bhl.config.set_cparm('COULOMB', True)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
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
bhl.config.set_cparm('OUTPUT_EOSVARS', False)

# RADIATION
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', True)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 2000.)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = 40.)
bhl.config.set_rparm('Rout_rad', 'double', default = 40.)
bhl.config.set_rparm('gam', 'double', default = 13./9.)
bhl.config.set_rparm('DTd', 'double', default = 5.)
bhl.config.set_rparm('DTl', 'double', default = 5.e-1)
bhl.config.set_rparm('DTr', 'double', default = 1000.)
bhl.config.set_rparm('DNr', 'integer', default = 1000)
bhl.config.set_rparm('a', 'double', default = 0.9375)
bhl.config.set_rparm('mbh', 'double', default = 1.e8)
bhl.config.set_rparm('M_unit', 'double', default = 8.e23)
bhl.config.set_rparm('tune_emiss', 'double', 1.0)
bhl.config.set_rparm('tune_scatt', 'double', 0.1)
bhl.config.set_rparm('t0_tune_emiss', 'double', 0)
bhl.config.set_rparm('t0_tune_scatt', 'double', 0)
bhl.config.set_rparm('MAD', 'int', default = 0)
bhl.config.set_rparm('BHflux', 'double', default = 0.)
bhl.config.set_rparm('numin', 'double', default=1.e9) 
bhl.config.set_rparm('numax', 'double', default=1.e25) 

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

