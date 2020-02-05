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
from bhlight import cgs
PROB = 'f1z'

MPI   = '-mpi' in sys.argv
EQUIL = '-equil' in sys.argv
HDF   = '-hdf' in sys.argv
FAST  = '-fast' in sys.argv
YE0 = 0.1
T0 = 2.5 # mev
rho0 = 1.e9 # g/cm^3

BURROWS_OPACITIES = not HDF
KILL_ALL_PACKETS= not EQUIL

TABLEPATH = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
TABLEPATH = "../../data/"+TABLEPATH
OPACPARAM = "opacbin.LS220.evan.param"
if HDF:
    OPACPATH = "opacity.SFHo.nohoro.juo.brem1.h5"
else:
    OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM

GAS_BC = 'BC_PERIODIC' if EQUIL else 'BC_OUTFLOW'
RAD_BC = 'BC_PERIODIC' if EQUIL else 'BC_ESCAPE'
ABSORPTION = EQUIL

# radiation
NUMIN_MEV = 0.1
NUMAX_MEV = 500
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
print("Numin, Numax = [{}, {}]".format(NUMIN,NUMAX))
TUNE_EMISS = 1e-2 if EQUIL else 1.

# units
RHO_unit = rho0
L_unit = 10.*2.*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
M_unit = RHO_unit*(L_unit**3)
T_unit = L_unit / cgs['CL']

# output and timestep parameters
if FAST:
    TFINAL = 0.1 / T_unit # seconds
else:
    TFINAL = 0.5 / T_unit # seconds
DT0 = TFINAL/1.e6
DTd = TFINAL/100
DTl = TFINAL/100

# cells and CPUs
NTOT = 4 if MPI else 1
NCPU = 2 if MPI else 1

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
bhl.config.set_cparm('X1L_GAS_BOUND', GAS_BC)
bhl.config.set_cparm('X1R_GAS_BOUND', GAS_BC)
bhl.config.set_cparm('X2L_GAS_BOUND', GAS_BC)
bhl.config.set_cparm('X2R_GAS_BOUND', GAS_BC)
bhl.config.set_cparm('X3L_GAS_BOUND', GAS_BC)
bhl.config.set_cparm('X3R_GAS_BOUND', GAS_BC)
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
bhl.config.set_cparm('ABSORPTION', ABSORPTION)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('BURROWS_OPACITIES', BURROWS_OPACITIES)
bhl.config.set_cparm('HDF5_OPACITIES',    HDF)
bhl.config.set_cparm('GRAYABSORPTION', False)
bhl.config.set_cparm('BREMSSTRAHLUNG', False)
bhl.config.set_cparm('SYNCHROTRON', False)
bhl.config.set_cparm('EXPTAU_WEIGHTS', False)
bhl.config.set_cparm('X1L_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('X1R_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('X2L_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('X2R_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('X3L_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('X3R_RAD_BOUND', RAD_BC)
bhl.config.set_cparm('KILL_ALL_PACKETS', KILL_ALL_PACKETS)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = DT0)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DNr', 'int', default = 1e7)
bhl.config.set_rparm('tune_emiss', 'double', default = TUNE_EMISS)
bhl.config.set_rparm('L_unit', 'double', default = L_unit)
bhl.config.set_rparm('M_unit', 'double', default = M_unit)
bhl.config.set_rparm('eospath', 'string', default = TABLEPATH)
bhl.config.set_rparm('rho0', 'double', default = 1.0)
bhl.config.set_rparm('T0', 'double', default = 2.5)
bhl.config.set_rparm('ye0', 'double', default = YE0)
bhl.config.set_rparm('numin', 'double', default = NUMIN)
bhl.config.set_rparm('numax', 'double', default = NUMAX)
bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)
bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
