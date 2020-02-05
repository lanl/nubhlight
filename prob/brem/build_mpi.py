################################################################################
#                                                                              #
#  OPTICALLY THIN BREMSSTRAHLUNG COOLING                                       #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl; del sys
PROB = 'brem_mpi'

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', 6)
bhl.config.set_cparm('N2TOT', 6)
bhl.config.set_cparm('N3TOT', 6)
bhl.config.set_cparm('N1CPU', 2)
bhl.config.set_cparm('N2CPU', 2)
bhl.config.set_cparm('N3CPU', 2)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', True)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
bhl.config.set_cparm('X1L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X1R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X2R_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3L_GAS_BOUND', 'BC_OUTFLOW')
bhl.config.set_cparm('X3R_GAS_BOUND', 'BC_OUTFLOW')

# RADIATION
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', True)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_ESCAPE')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = 256.)
bhl.config.set_rparm('dt', 'double', default = 256.e-6)
bhl.config.set_rparm('L_unit', 'double', default = 1.17106428906e16)
bhl.config.set_rparm('M_unit', 'double', default = 1.57378e32)
bhl.config.set_rparm('DTl', 'double', default = 10.)
bhl.config.set_rparm('T0', 'double', default = 1.e8)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

