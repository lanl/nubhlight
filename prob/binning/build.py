################################################################################
#                                                                              #
#  UNIT TEST FOR SUPERPHOTON BINNING                                           #
#                                                                              #
################################################################################

import sys; sys.path.append('../../script/'); 
sys.dont_write_bytecode = True; import bhlight as bhl;
PROB = 'binning'
MPI = '-mpi' in sys.argv
if '-idim' in sys.argv:
    IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
    IDIM = 0

NCPU = 2 if MPI else 1
NTOT = 12
TF = 1.
DTout = TF/2.
T0 = TF/1e6
NR = 3e6

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

# RADIATION
bhl.config.set_cparm('RADIATION', True)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('EMISSION', False)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('NU_BINS_SPEC', 200)
bhl.config.set_cparm("NTH", 100)
bhl.config.set_cparm("NPHI", 100)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = TF)
bhl.config.set_rparm('dt', 'double', default = T0)
bhl.config.set_rparm('tune_emiss', 'double', default = 100.)
bhl.config.set_rparm('L_unit', 'double', default = 1.17106428906e16)
bhl.config.set_rparm('M_unit', 'double', default = 1.57378e32)
bhl.config.set_rparm('DTl', 'double', default = DTout)
bhl.config.set_rparm('DTd', 'double', default = DTout)
bhl.config.set_rparm('DTr', 'double', default = 1e6)
bhl.config.set_rparm('direction','int', default = IDIM)
bhl.config.set_rparm('Nr0', 'int', default = NR)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)

