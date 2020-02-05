################################################################################
#                                                                              #
#  Advection of Passive Scalars in 1D                                          #
#                                                                              #
################################################################################

from __future__ import print_function
import sys
sys.path.append('../../script/')
sys.dont_write_bytecode = True; import bhlight as bhl

PROB = 'advection1d'
TRACERS = '-tracers' in sys.argv

if '-ntot' in sys.argv:
  NTOT = int(sys.argv[sys.argv.index('-ntot')+1])
else:
  NTOT = 512

if '-idim' in sys.argv:
  IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
  IDIM = 1
if IDIM == 1:
    N1 = NTOT
    N2 = 1
    N3 = 1
elif IDIM == 2:
    N1 = 1
    N2 = NTOT
    N3 = 1
elif IDIM == 3:
    N1 = 1
    N2 = 1
    N3 = NTOT
else:
    raise ValueError("Invalid IDIM. IDIM = {}".format(IDIM))

TF = 6 if TRACERS else 2

TRACERS_PER_CELL = 10
TRACERS_TOT = N1*N2*N3*TRACERS_PER_CELL

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', N1)
bhl.config.set_cparm('N2TOT', N2)
bhl.config.set_cparm('N3TOT', N3)

bhl.config.set_cparm('N1CPU', 1)
bhl.config.set_cparm('N2CPU', 1)
bhl.config.set_cparm('N3CPU', 1)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

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

# PASSIVE SCALARS
bhl.config.set_cparm('NVAR_PASSIVE', 2)

# Tracers
if TRACERS:
  bhl.config.set_cparm('RADIATION', True)
  bhl.config.set_cparm('ESTIMATE_THETAE', False)
  bhl.config.set_cparm('EMISSION', False)
  bhl.config.set_cparm('ABSORPTION', False)
  bhl.config.set_cparm('SCATTERING', False)
  bhl.config.set_cparm('TRACERS', True)
  bhl.config.set_cparm('NU_BINS', 200)
  bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
  bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
  bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
  bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
  bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
  bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = TF)
bhl.config.set_rparm('dt', 'double', default = 1e-6)
bhl.config.set_rparm('gam', 'double', default = 5./3.)
bhl.config.set_rparm('DTd', 'double', default = 0.02)
bhl.config.set_rparm('DTl', 'double', default = 0.02)
bhl.config.set_rparm('DTr', 'double', default = 10000)
bhl.config.set_rparm('cadv', 'double', default = 0.5)
bhl.config.set_rparm('idim', 'int', default = IDIM)

if TRACERS:
  bhl.config.set_rparm('ntracers', 'int', default=TRACERS_TOT)
  bhl.config.set_rparm('L_unit', 'double', default=1)
  bhl.config.set_rparm('M_unit', 'double', default=1)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
