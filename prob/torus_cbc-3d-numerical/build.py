################################################################################
#                                                                              #
#  TORUS, Tuned for Compact Binary Coalescence                                 #
#                                                                              #
################################################################################

# TODO: use argparse

from __future__ import print_function,division

import sys
sys.path.append('../../script')
sys.path.append('../../script/analysis')
sys.dont_write_bytecode = True
import bhlight as bhl
from bhlight import cgs
from units import UnitSystem
from math import log10,ceil

PROB = 'torus_cbc-3dnumerical'

RELTABLE = '-reltable' in sys.argv
INITIAL = '-initial' in sys.argv # initial data only

NOB = '-nob' in sys.argv # or INITIAL # no B fields
TOROIDALB = '-toroidalb' in sys.argv
UNIFORMZB = '-uniformz' in sys.argv

RENORM = '-renorm' in sys.argv
NORENORM = '-norenorm' in sys.argv or (not RENORM)
DEBUG = '-debug' in sys.argv
THREED = '-3d' in sys.argv
FAKETHREED = '-fake3d' in sys.argv
NEUTRINOS = '-nu' in sys.argv
NOEMISS = '-noemiss' in sys.argv
NOSCATT = '-noscatt' in sys.argv
NOABS = '-noabs' in sys.argv
DIAGNOSTIC = '-diag' in sys.argv
RESTARTTEST = '-restarttest' in sys.argv
HDF = '-hdf' in sys.argv

LIMIT_RAD = '-limit' in sys.argv
MORE_RAD = '-morenu' in sys.argv

QUAD = '-quad' in sys.argv # quadrant symmetry
KILL = '-kill' in sys.argv # kill all packets

EMISS = not NOEMISS
SCATT = not (NOSCATT or KILL)
ABS = not (NOABS or KILL)
FORTRAN = NEUTRINOS and not HDF

RENORM = not NORENORM

if NOB:
    BFIELD = "none"
elif UNIFORMZB:
    BFIELD = "uniformz"
elif TOROIDALB:
    BFIELD = "toroidal"
else:
    BFIELD = "classic"

if RELTABLE:
    bhl.report_var('EOS','RELTABLE')

if NEUTRINOS:
    RADIATION = 'RADTYPE_NEUTRINOS'
else:
    RADIATION = False
    
EOSNAME = "Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
TABLEPATH = "../../data/"+EOSNAME
#TABLEPATH = EOSNAME

CARPETPROFNAME = "245414.h5"
CARPETPROFPATH = "../../data/"+CARPETPROFNAME

if FORTRAN:
    OPACPATH = "opacity.SFHo.nohoro.juo.brem1.bin"
else:
    OPACPATH = "opacity.SFHo.nohoro.juo.brem1.h5"
OPACPARAM = "opacbin.LS220.evan.param"
OPACPATH = "../../data/"+OPACPATH
OPACPARAM = "../../data/"+OPACPARAM

# Parameters from
# arXiv 1409.4426
# arXiv 1705.05473

BETA = 1e2
    
RHO_unit = 6.17244*1e17

# Note that if scattering is enabled,
# you'll get twice as many photons
# per cell as this number
# (unchanged for legacy reasons)
NPH_PER_CELL = 10

bhl.report_var('RHO_unit',RHO_unit)

# IMPORTANT TO GET RIGHT!
# Chosen so that guess is right at rho = 1.0
# in code units

L_UNIT = cgs['GNEWT']*cgs['MSOLAR']*1/(cgs['CL']**2) # M = 1*M_Sun
M_UNIT = RHO_unit*(L_UNIT**3)

# final time
if NOB and not INITIAL:
    TFINAL = 500.
elif INITIAL:
    TFINAL = 1.e-2
elif THREED or FAKETHREED:
    TFINAL = 10000.
else:
    TFINAL = 2000.
# Toroidal fields grow slowly
if TOROIDALB:
    TFINAL *= 2

Rout = 1000.
#Rout_rad = ENTROPY*ceil(Rmax) # not safe to use 3x
#Rout_vis = Rout
#
## time when radiation resolution controls turn on
#t0_tune_emiss = Rout_rad if LIMIT_RAD else -1.
#t0_tune_scatt = 2.*max(Rout_rad,t0_tune_emiss)

if LIMIT_RAD:
    tune_emiss = 1.
elif MORE_RAD:
    tune_emiss = 100.
else:
    tune_emiss = 1.

# output
DTd = 1.e-2 if INITIAL else 5.
DTl = 1.e-2 if INITIAL else 5.e-1
DTr = 100
DNr = 50 if RESTARTTEST else 1000

N1TOT = 32
N2TOT = 32
N3TOT = 1

N1CPU = 1
N2CPU = 1
N3CPU = 1
    
NCPU_TOT     = N1CPU*N2CPU*N3CPU
NCELL_TOT    = N1TOT*N2TOT*N3TOT
NPH_TOT      = NCELL_TOT*NPH_PER_CELL
NPH_PER_PROC = NPH_TOT / NCPU_TOT

if not NEUTRINOS:
    NPH_PER_PROC = 0.0

EOS = "EOS_TYPE_TABLE"
NVP0 = 3
NVAR_PASSIVE = NVP0 if RELTABLE else 0
#OPENMP = not DEBUG
OPENMP = True

tablepath = TABLEPATH

# Radiation units
NUMIN_MEV = 1.0
NUMAX_MEV = 100
NUMIN = NUMIN_MEV*cgs['MEV']/cgs['HPL']
NUMAX = NUMAX_MEV*cgs['MEV']/cgs['HPL']
if RADIATION:
    bhl.report_var("[Numin, Numax] (meV)", [NUMIN_MEV,NUMAX_MEV])
    bhl.report_var("[Numin, Numax] (cgs)", [NUMIN,NUMAX])
    bhl.report_var('NPH_PER_CELL', NPH_PER_CELL)
    bhl.report_var('NPH_TOT',NPH_TOT)
    bhl.report_var('NPH_PER_PROC',NPH_PER_PROC)

                         ### COMPILE TIME PARAMETERS ###

# SPATIAL RESOLUTION AND MPI DECOMPOSITION
bhl.config.set_cparm('N1TOT', N1TOT)
bhl.config.set_cparm('N2TOT', N2TOT)
bhl.config.set_cparm('N3TOT', N3TOT)
bhl.config.set_cparm('N1CPU', N1CPU)
bhl.config.set_cparm('N2CPU', N2CPU)
bhl.config.set_cparm('N3CPU', N3CPU)

# OPENMP PARALLELIZATION
bhl.config.set_cparm('OPENMP', OPENMP)

# COORDINATES
bhl.config.set_cparm('METRIC', 'NUMERICAL')
bhl.config.set_cparm('DEREFINE_POLES', THREED or FAKETHREED)

# ELECTRONS. DO NOT USE THESE HERE!
bhl.config.set_cparm('ELECTRONS', False)

# FLUID
bhl.config.set_cparm('RECONSTRUCTION', 'WENO')
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
bhl.config.set_cparm('QUADRANT_SYMMETRY', QUAD)

#EOS
bhl.config.set_cparm("EOS", EOS)
bhl.config.set_cparm('NVAR_PASSIVE', NVAR_PASSIVE)

# RADIATION
bhl.config.set_cparm('RADIATION', RADIATION)
bhl.config.set_cparm('EMISSION',   EMISS)
bhl.config.set_cparm('ABSORPTION', ABS)
bhl.config.set_cparm('SCATTERING', SCATT)
if KILL:
    bhl.config.set_cparm('KILL_ALL_PACKETS', True)
bhl.config.set_cparm('BURROWS_OPACITIES', FORTRAN)
bhl.config.set_cparm('HDF5_OPACITIES', HDF)
bhl.config.set_cparm('NU_BINS', 200)
bhl.config.set_cparm('ESTIMATE_THETAE', False)
bhl.config.set_cparm('GRAYABSORPTION',  False)
bhl.config.set_cparm('BREMSSTRAHLUNG',  False)
bhl.config.set_cparm('SYNCHROTRON',     False)
bhl.config.set_cparm('X1L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X1R_RAD_BOUND', 'BC_CAMERA')
bhl.config.set_cparm('X2L_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X2R_RAD_BOUND', 'BC_ESCAPE')
bhl.config.set_cparm('X3L_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('X3R_RAD_BOUND', 'BC_PERIODIC')
bhl.config.set_cparm('DIAGNOSTICS_USE_RADTYPES', True)

bhl.config.set_cparm("EXIT_ON_INIT", True)

# Special. Don't turn this on if you don't need to
if DIAGNOSTIC:
    bhl.config.set_cparm("EXIT_ON_INIT", True)
    bhl.config.set_cparm("OUTPUT_BURROWS_OPACITIES", True)

                           ### RUNTIME PARAMETERS ###
# Carpet profile path
bhl.config.set_rparm('carpetprofpath', 'string', default = CARPETPROFPATH)
# generic
bhl.config.set_rparm('tf', 'double', default = TFINAL)
bhl.config.set_rparm('dt', 'double', default = 1.e-6)
bhl.config.set_rparm('Rout', 'double', default = Rout)
# Maybe this should be something further out. Fine for test.
#bhl.config.set_rparm('Rout_rad', 'double', default = Rout_rad)
#bhl.config.set_rparm('Rout_vis', 'double', default = Rout_vis)
bhl.config.set_rparm('DTd', 'double', default = DTd)
bhl.config.set_rparm('DTl', 'double', default = DTl)
bhl.config.set_rparm('DTr', 'double', default = DTr)
bhl.config.set_rparm('DNr', 'integer', default = DNr)
bhl.config.set_rparm('tune_emiss', 'double', default = tune_emiss)
#bhl.config.set_rparm('t0_tune_emiss', 'double', default = t0_tune_emiss)
#bhl.config.set_rparm('t0_tune_scatt', 'double', default = t0_tune_scatt) 

# problem
bhl.config.set_rparm('M_unit', 'double', default = M_UNIT)
bhl.config.set_rparm('L_unit', 'double', default = L_UNIT)
bhl.config.set_rparm('RHO_unit', 'double', default = RHO_unit)
bhl.config.set_rparm("bfield", 'string', default = BFIELD)
bhl.config.set_rparm("beta", "double", default = BETA)
bhl.config.set_rparm("renorm_dens", "int", default = int(RENORM))

if RADIATION:
    bhl.config.set_rparm('numin', 'double', default = NUMIN)
    bhl.config.set_rparm('numax', 'double', default = NUMAX)
    bhl.config.set_rparm('nph_per_proc', 'double', default = NPH_PER_PROC)
    
#EOS
if RELTABLE:
    bhl.config.set_rparm('eospath', 'string', default = tablepath)
        
# Opacities
if FORTRAN or HDF:
    bhl.config.set_rparm('opac_file', 'string', default = OPACPATH)
if FORTRAN:
    bhl.config.set_rparm('opac_param_file', 'string', default = OPACPARAM)


                         ### CONFIGURE AND COMPILE  ###

bhl.build(PROB)
