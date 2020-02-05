################################################################################
#                                                                              #
#  OPTICALLY THIN BREMSSTRAHLUNG COOLING                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'thermalbrem'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 1)
bhl.config.set_cparm('N2TOT', 1)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_PERIODIC')

# RADIATION
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', True)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 512)
bhl.config.set_rparm('dt', 'double', default = 5120.e-6)
bhl.config.set_rparm('L_unit', 'double', default = 5.85532e7)
bhl.config.set_rparm('M_unit', 'double', default = 2.01466e15)
bhl.config.set_rparm('tune_emiss', 'double', default = 1.e-5)
bhl.config.set_rparm('DTd', 'double', default = 5.)
bhl.config.set_rparm('DTl', 'double', default = 500.)
bhl.config.set_rparm('DTr', 'integer', default = 500000.)
bhl.config.set_rparm('T0', 'double', default = 1.e6)
bhl.config.set_rparm('ne', 'double', default = 6.e15)
bhl.config.set_rparm('nph_per_proc', 'double', default=8.*4.e5)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

