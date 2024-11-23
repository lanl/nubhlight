################################################################################
#                                                                              #
#  UNIT TEST FOR MULTISCATT                                                    #
#                                                                              #
################################################################################

import sys
sys.dont_write_bytecode = True
sys.path.append('../../script/')
sys.path.append('../../script/analysis')
import bhlight as bhl
import make_tabulated_gamma as tab
from units import cgs
PROB = 'multiscatt'

MPI = '-mpi' in sys.argv
PHOTONS = '-photons' in sys.argv
NEUTRINOS = not PHOTONS
NOREBALANCE = '-norebalance' in sys.argv
if '-idim' in sys.argv:
  IDIM = int(sys.argv[sys.argv.index('-idim')+1])
else:
  IDIM = 3

print("MPI = {}".format(MPI))
print("PHOTONS = {}".format(PHOTONS))
print("NEUTRINOS = {}".format(NEUTRINOS))
print("IDIM = {}".format(IDIM))

NCPU = 2 if MPI else 1
NTOT = 12
TF = 1.
DTout = TF
T0 = TF/1e6
NR = 1e7 # if NEUTRINOS else 1e6
E_MEV = 25

# Use a fake table for this test
GAMMA = 1.4
TFINAL = 600
YE = 0.5

CL = cgs['CL']
RHOMIN, RHOMAX, NRHO =  1e-4, 1e20, 234
UMIN, UMAX, NU = 1e-8, 1e8, 136
YEMIN, YEMAX, NYE = 0.0, 0.6, 50
CRASH_ON_SOUND_SPEED = False

MS_THETA = (E_MEV*cgs['MEV']/(cgs['ME']*CL*CL))**2
MS_DELTA = -1.0

if PHOTONS:
  M_UNIT = 1.57378e32
  L_UNIT = 1.17106428906e16
else:
  dtau = 1.0
  rhol = dtau/((14./3.)*cgs['NUSIGMA0']*(MS_THETA**2)*(1.-YE)/cgs['MN'])
  RHO_UNIT = 2.8e14
  L_UNIT = rhol/RHO_UNIT
  M_UNIT = RHO_UNIT*(L_UNIT**3)

print("rhol = {:.5}".format(rhol))
print("RHO_UNIT = {:.5}".format(RHO_UNIT))
print("L_UNIT = {:.5}".format(L_UNIT))
print("M_UNIT = {:.5}".format(M_UNIT))

tablepath = tab.make_filename(GAMMA)
units = tab.UnitSystem(M_UNIT, L_unit = L_UNIT)
tab.make_table_u(RHOMIN, RHOMAX, NRHO,
                 UMIN,   UMAX,   NU,
                 YEMIN,  YEMAX,  NYE,
                 units,  GAMMA,  tablepath,
                 CRASH_ON_SOUND_SPEED)


RHO0 = 1.
#RHO0 = 1.e10 if NEUTRINOS else 1.e1
UU0 = 1.e-4*RHO0 if PHOTONS else 1e-6*RHO0
Theta0 = (GAMMA - 1.)*UU0/RHO0
print("rho_cgs = {}\nu_cgs = {}\n".format(RHO0,UU0))
print("Theta = {}".format(Theta0))
print("MS_THETA = {}".format(MS_THETA))
print("MS_DELTA = {}".format(MS_DELTA))


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

# EOS
bhl.config.set_cparm("EOS", "EOS_TYPE_TABLE")
bhl.config.set_cparm('NVAR_PASSIVE', 2)

# RADIATION
if NEUTRINOS:
    bhl.config.set_cparm('RADIATION', 'RADTYPE_NEUTRINOS')
else:
    bhl.config.set_cparm('RADIATION', True)
if NOREBALANCE:
    bhl.config.set_cparm('REBALANCE', False)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('EMISSION', False)
bhl.config.set_cparm('ABSORPTION', False)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('MULTISCATT_TEST', True)
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
if NEUTRINOS:
  if NOREBALANCE:
    bhl.config.set_rparm('tune_scatt', 'double', default = 20)
  else:
    bhl.config.set_rparm('tune_scatt', 'double', default = 20)#34)
bhl.config.set_rparm('t0_tune_emiss', 'double', default = 10.*TF)
bhl.config.set_rparm('t0_tune_scatt', 'double', default = 10.*TF)
bhl.config.set_rparm('L_unit', 'double', default = L_UNIT)
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm('DTl', 'double', default = DTout)
bhl.config.set_rparm('DTd', 'double', default = DTout)
bhl.config.set_rparm('DTr', 'double', default = 1e6)
bhl.config.set_rparm('Nr0', 'int', default = NR)
bhl.config.set_rparm('E_MEV', 'double', default = E_MEV)
bhl.config.set_rparm('rho0', 'double', default = RHO0)
bhl.config.set_rparm('uu0', 'double', default = UU0)
bhl.config.set_rparm('eospath', 'string', default = tablepath)
bhl.config.set_rparm('ms_theta_nu0', 'double', default = MS_THETA)
bhl.config.set_rparm('ms_delta0', 'double', default = MS_DELTA)
bhl.config.set_rparm('numin', 'double', default = 1.e8)
bhl.config.set_rparm('numax', 'double', default = 1.e30)
bhl.config.set_rparm('direction', 'int', default = IDIM)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
