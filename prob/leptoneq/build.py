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
PROB = 'leptoneq'

MPI = '-mpi' in sys.argv
YE0 = 0.225
DYE = 0.125
T0 = 2.5 # MeV
rho0 = 1.e10
bhl.report_var('YE0',YE0)
bhl.report_var('DYE',DYE)
bhl.report_var('T0',T0)
bhl.report_var('rho0',rho0)

# Units
RHO_unit = rho0
L_unit = 10.*2.*cgs['GNEWT']*cgs['MSOLAR']/(cgs['CL']**2)
M_unit = RHO_unit*(L_unit**3)
T_unit = L_unit / cgs['CL']
bhl.report_var('RHO_unit', RHO_unit)
bhl.report_var('L_unit',   L_unit)
bhl.report_var('M_unit',   M_unit)
bhl.report_var('T_unit',   T_unit)

# output and timestep parameters
TFINAL = 1e-2 / T_unit # seconds
DT0    = TFINAL/1.e6
DTd    = TFINAL/100
DTl    = TFINAL/100
bhl.report_var('TFINAL',TFINAL)
bhl.report_var('DT0',   DT0)
bhl.report_var('DTd',   DTd)
bhl.report_var('DTl',   DTl)

# cells and CPUs
NTOT = 20
NCPU = 2 if MPI else 1

NPH_PER_CELL = 100
NCPU_TOT     = NCPU**2
NCELL_TOT    = NTOT**2
NPH_TOT      = NCELL_TOT*NPH_PER_CELL
NPH_PER_PROC = NPH_TOT / NCPU_TOT
bhl.report_var('NPH_PER_PROC',NPH_PER_PROC)

TUNE_EMISS   = 1.
bhl.report_var('TUNE_EMISS',TUNE_EMISS)

# Radiation units
NUMIN_MEV = 1e-1
NUMAX_MEV = 1e3
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
bhl.report_var("[Numin, Numax] (meV)", [NUMIN_MEV,NUMAX_MEV])
bhl.report_var("[Numin, Numax] (cgs)", [NUMIN,NUMAX])
bhl.report_var('NPH_TOT',NPH_TOT)
bhl.report_var('NPH_PER_PROC',NPH_PER_PROC)

TABLEPATH = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
OPACPARAM = "opacbin.LS220.evan.param"

TABLEPATH = "../../data/"+TABLEPATH
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM

bhl.report_var('TABLEPATH', TABLEPATH)
bhl.report_var('OPACPATH',  OPACPATH)
bhl.report_var('OPACPARAM', OPACPARAM)



                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', NTOT)
bhl.config.set_cparm('N2TOT', NTOT)
bhl.config.set_cparm('N3TOT', 1)
bhl.config.set_cparm('N1CPU', NCPU)
bhl.config.set_cparm('N2CPU', NCPU)
bhl.config.set_cparm('N3CPU', 1)

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
bhl.config.set_cparm('NO_GRMHD_UPDATE', True)

# EOS
bhl.config.set_cparm("EOS", "EOS_TYPE_TABLE")
bhl.config.set_cparm('NVAR_PASSIVE', 2)


# RADIATION
bhl.config.set_cparm('RADIATION', 'RADTYPE_NEUTRINOS')
bhl.config.set_cparm('EMISSION', True)
bhl.config.set_cparm('ABSORPTION', True)
bhl.config.set_cparm('SCATTERING', False)
bhl.config.set_cparm('BURROWS_OPACITIES', True)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('GRAYABSORPTION',  False)
bhl.config.set_cparm('BREMSSTRAHLUNG',  False)
bhl.config.set_cparm('SYNCHROTRON',     False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')

                           ### RUNTIME PARAMETERS ###

# Coordinates
# box [-1, 1]^2
bhl.config.set_rparm('x1Min', 'double', default=-1.)
bhl.config.set_rparm('x1Max', 'double', default=1.)
bhl.config.set_rparm('x2Min', 'double', default=-1.)
bhl.config.set_rparm('x2Max', 'double', default=1.)

# Generic
bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = DT0)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DNr', 'int', default = 1e7)
bhl.config.set_rparm('tune_emiss', 'double', default = TUNE_EMISS)
bhl.config.set_rparm('L_unit', 'double', default = L_unit)
bhl.config.set_rparm('M_unit', 'double', default = M_unit)

# Radiation
bhl.config.set_rparm('numin', 'double', default = NUMIN)
bhl.config.set_rparm('numax', 'double', default = NUMAX)
bhl.config.set_rparm('nph_per_proc', 'double', default = NPH_PER_PROC)

# EOS
bhl.config.set_rparm('eospath', 'string', default = TABLEPATH)

# opacities
bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)
bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)

# Problem
bhl.config.set_rparm('rho0', 'double', default = 1.0)
bhl.config.set_rparm('T0',   'double', default = T0)
bhl.config.set_rparm('ye0',  'double', default = YE0)
bhl.config.set_rparm('dye',  'double', default = DYE)
bhl.config.set_rparm('c1x',  'double', default = -0.5)
bhl.config.set_rparm('c1y',  'double', default = -0.5)
bhl.config.set_rparm('r1',   'double', default = 0.25)
bhl.config.set_rparm('c2x',  'double', default = 0.5)
bhl.config.set_rparm('c2y',  'double', default = 0.5)
bhl.config.set_rparm('r2',   'double', default = 0.25)

                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
