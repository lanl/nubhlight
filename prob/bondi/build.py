################################################################################
#                                                                              #
#  BONDI                                                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script');
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'bondi'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 32)
bhl.config.set_cparm('N2TOT', 32)
bhl.config.set_cparm('N3TOT', 1)
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
bhl.config.set_cparm('SUPPRESS_HIGHB_HEAT', True)
bhl.config.set_cparm('BETA_HEAT', True)
bhl.config.set_cparm('COULOMB', False)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PROB')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_POLAR')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1L_INFLOW', False)
bhl.config.set_cparm('X1R_INFLOW', True)
bhl.config.set_cparm('X2L_INFLOW', False)
bhl.config.set_cparm('X2R_INFLOW', False)
bhl.config.set_cparm('X3L_INFLOW', False)
bhl.config.set_cparm('X3R_INFLOW', False)

# RADIATION
bhl.config.set_cparm('RADIATION', False)
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

bhl.config.set_rparm('tf', 'double', default = 100.)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = 20.)
bhl.config.set_rparm('gam', 'double', default = 5./3.)
bhl.config.set_rparm('DTd', 'double', default = 5.)
bhl.config.set_rparm('DTl', 'double', default = 5.e-1)
bhl.config.set_rparm('DTr', 'integer', default = 1000000)
bhl.config.set_rparm('a', 'double', default = 0.)
                         
                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

