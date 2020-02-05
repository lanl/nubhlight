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
PROB = 'burrows'

TABLEPATH = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
TABLEPATH = "../../data/"+TABLEPATH
OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
OPACPARAM = "opacbin.LS220.evan.param"
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM
OPACOUT = "opacity.SFHo.nohoro.juo.brem1.h5"

GAS_BC = 'BC_PERIODIC' 
RAD_BC = 'BC_PERIODIC'

# radiation
NUMIN_MEV = 1
NUMAX_MEV = 320
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
print("Numin, Numax = [{}, {}]".format(NUMIN,NUMAX))

# units
RHO_unit = 1e12;
L_unit = 10.*2.*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
M_unit = RHO_unit*(L_unit**3)
T_unit = L_unit / cgs['CL']

# output and timestep parameters
TFINAL = 0.5 / T_unit # seconds
DT0 = TFINAL/1.e6
DTd = TFINAL/100
DTl = TFINAL/100

# cells and CPUs
NTOT = 1
NCPU = 1

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', NTOT)
bhl.config.set_cparm('N2TOT', NTOT)
bhl.config.set_cparm('N3TOT', NTOT)
bhl.config.set_cparm('N1CPU', NCPU)
bhl.config.set_cparm('N2CPU', NCPU)
bhl.config.set_cparm('N3CPU', NCPU)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', False)

# COORDINATES
bhl.config.set_cparm('METRIC', 'MINKOWSKI')

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'LINEAR')
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
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', True)
bhl.config.set_cparm('NU_BINS', 60)
bhl.config.set_cparm('BURROWS_OPACITIES', True)
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

# Exit on startup
bhl.config.set_cparm("EXIT_ON_INIT", True)

                           ### RUNTIME PARAMETERS ###

bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = DT0)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DNr', 'int', default = 1e7)
# bhl.config.set_rparm('tune_emiss', 'double', default = TUNE_EMISS)
bhl.config.set_rparm('L_unit', 'double', default = L_unit)
bhl.config.set_rparm('M_unit', 'double', default = M_unit)
bhl.config.set_rparm('eospath', 'string', default = TABLEPATH)
bhl.config.set_rparm('rho0', 'double', default = 1.0)
bhl.config.set_rparm('T0', 'double', default = 2.5)
bhl.config.set_rparm('ye0', 'double', default = 0.5)
bhl.config.set_rparm('numin', 'double', default = NUMIN)
bhl.config.set_rparm('numax', 'double', default = NUMAX)
bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)
bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)

# problem
bhl.config.set_rparm("outfile","string",default=OPACOUT)
bhl.config.set_rparm("lrho_min","double", default=5)
bhl.config.set_rparm("lrho_max","double", default=14)
bhl.config.set_rparm("numrho_out","int", default=48)
bhl.config.set_rparm("lT_min","double", default=-1.1)
bhl.config.set_rparm("lT_max","double", default=1.6)
bhl.config.set_rparm("numT_out","int",   default=52)
bhl.config.set_rparm("Ye_min","double", default=0.035)
bhl.config.set_rparm("Ye_max","double", default=0.56)
bhl.config.set_rparm("numYe_out","int",  default=30)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
