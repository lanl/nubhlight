/******************************************************************************
 *                                                                            *
 * DECS.H                                                                     *
 *                                                                            *
 * GLOBAL MACROS, FUNCTION DEFINITIONS, INCLUDES, AND DECLARATIONS            *
 *                                                                            *
 ******************************************************************************/

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "constants.h"
#include "params.h"

/*******************************************************************************
      PRE-PROCESSOR MAGIC :
*******************************************************************************/
#define XSTR(x) STR(x)
#define STR(x) #x

/*******************************************************************************
      CONSTANTS :
*******************************************************************************/

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif

/*******************************************************************************
      COMPILE-TIME PARAMETERS :
*******************************************************************************/

// Number of active zones on each MPI process
#define N1 (N1TOT / N1CPU)
#define N2 (N2TOT / N2CPU)
#define N3 (N3TOT / N3CPU)

// Max size for 1D slice is NMAX
#define N12 (N1 > N2 ? N1 : N2)
#define NMAX (N12 > N3 ? N12 : N3)

#define NDIM (4) // Number of total dimensions
#define NPG (5)  // Number of positions on grid for grid functions
#define NG (3)   // Number of ghost zones
#define N123G                      \
  ((N1 + 2 * NG) * (N2 + 2 * NG) * \
      (N3 + 2 * NG)) // total number of cells on a cpu

// Fixup parameters
// rho_flor = RHOMIN / pow(r, -FLOOR_POWER)
#define FLR_POWER1 (2.)
#define FLR_POWER2 (2.)
#define FLR_R0 (50.) // TODO change this?
#define RHOMINLIMIT (1.e-17)
#define UUMINLIMIT (1.e-3 * RHOMINLIMIT)
#define RHOMIN (1.e-5)
#define UUMIN (1.e-3 * RHOMIN)
#define BSQORHOMAX (50.)
#define BSQOUMAX (2500.)
#define RHOEPS (0.9)
#define UORHOMAX (50.)
#define COLD_FLOORS (1)
#define ATM_THRESH (0.1)

// Numerical convenience to represent a small (<< 1) non-zero quantity
#define SMALL (1.e-20)
#define ABSURDLY_LARGE (1.e100)

// Maximum value of gamma, the Lorentz factor
#define GAMMAMAX (50.)

// Maximum fractional increase in timestep per timestep
#define SAFE (1.3)

// Superphoton diagnostics
#define MAXNSCATT (3)

// Chosen to give some wiggle room
// in scattering stability criterion
#define SCATT_BIAS_SAFETY (0.5) // (0.9/RAD_SCATT_TYPES)

//#define NUMIN     (1.e10)
//#define NUMAX     (1.e25)
//#define NU_BINS_SPEC    (200)
//#define NTH (8)
//#define NPHI (8)

// Whether to move polar axis slightly off of coordinate singularity
#define COORDSINGFIX 1
#define SINGSMALL (1.E-20)

// Whether to record a representative sample of superphoton positions
#define TRACK_PH (0)

// I/O format strings
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"
#define STRLEN (2048)

// Reconstruction algorithms
#define LINEAR (0)
#define PPM (1)
#define WENO (2)
#define MP5 (3)

// Primitive and conserved variables
#define RHO (0)
#define UU (1)
#define U1 (2)
#define U2 (3)
#define U3 (4)
#define B1 (5)
#define B2 (6)
#define B3 (7)
#define NVAR_BASE (B3 + 1)

// Passive variables (if present)
#define PASSIVE_START (NVAR_BASE)
#define PASSIVE_STOP (NVAR_BASE + NVAR_PASSIVE)
#define PASSTYPE_INTRINSIC (0)
#define PASSTYPE_NUMBER (1)
#if EOS == EOS_TYPE_TABLE
#define YE (PASSIVE_START)
#define YE_EM (YE + 1)
#if METRIC == MKS
#define ATM (YE + 2)
#define IS_ATM (0.0)
#define NOT_ATM (1.0)
#endif // METRIC
#endif // EOS_TYPE_TABLE

#if ELECTRONS
#define KEL (B3 + NVAR_PASSIVE + 1)
#define KTOT (B3 + NVAR_PASSIVE + 2)
#define NVAR (NVAR_BASE + NVAR_PASSIVE + 2)
#else
#define NVAR (NVAR_BASE + NVAR_PASSIVE)
#endif

// NETRINOS
#define RATYPE_NONE (0)
#define RADTYPE_LIGHT (1)
#define RADTYPE_NEUTRINOS (2)
#define RAD_TYPE_START (0)
#define TYPE_TRACER (-1)
#if RADIATION == RADTYPE_NEUTRINOS
#define NU_ELECTRON (0)
#define ANTINU_ELECTRON (1)
#define NU_HEAVY (2)
#define ANTINU_HEAVY (3)
#if MULTISCATT_TEST
#define RAD_SCATT_TYPES (3)
#else
#define RAD_SCATT_TYPES (4) // TODO: Should be 5, including electrons
#define RSCATT_TYPE_P (0)
#define RSCATT_TYPE_N (1)
#define RSCATT_TYPE_A (2)
#define RSCATT_TYPE_ALPHA (3)
#define RSCATT_TYPE_E (4) // TOOD: implement me
#endif
#define NRADCOMP (2)
#define RADG_YE (4)
#define RADG_YE_EM (5)
#elif RADIATION == RADTYPE_LIGHT
#define RAD_SCATT_TYPES (1)
#define NRADCOMP (0)
#define PHOTON (0)
#else
#define RAD_SCATT_TYPES (0)
#define NRADCOMP (0)
#endif
#define DLEPTON_THRESH (1e-10)

// EOS
#define EOS_TYPE_GAMMA (0)
#define EOS_TYPE_POLYTROPE (1)
#define EOS_TYPE_TABLE (2)
#if EOS == EOS_TYPE_GAMMA
#define EOS_NUM_EXTRA (0)
#define POLYTROPE_FALLBACK (0)
#elif EPS == EOS_TYPE_POLYTROPE
#define EOS_NUM_EXTRA (0)
#define POLYTROPE_FALLBACK (0)
#elif EOS == EOS_TYPE_TABLE
#if GAMMA_FALLBACK
#define POLYTROPE_FALLBACK (0)
// #define POLYTROPE_FALLBACK (1)
#else
#define POLYTROPE_FALLBACK (1)
#endif // GAMMA_FALLBACK
#define EOS_NUM_EXTRA (3)
#define EOS_LRHO (0)
#define EOS_LT (1)
#define EOS_YE (2)
// mass fractions
#define NUM_MASS_FRACTIONS (4)
#define MF_XA (0)
#define MF_XH (1)
#define MF_XN (2)
#define MF_XP (3)
#endif // EOS

// Centering of grid functions
#define FACE1 (0)
#define FACE2 (1)
#define FACE3 (2)
#define CORN (3)
#define CENT (4)
#define FACESTART (FACE1)
#define FACEEND (FACE3 + 1)
#define PGSTART (FACE1)
#define PGEND (NPG)

// Slope limiter
#define MC (0)
#define VANL (1)
#define MINM (2)

// Fluid and radiation boundaries
#define BC_OUTFLOW (0)
#define BC_PERIODIC (1)
#define BC_POLAR (2)
#define BC_PROB (3)
#define BC_ESCAPE (4)
#define BC_CAMERA (5)
#define BC_EQUILIB (6)

// Metric
#define MINKOWSKI (0)
#define SPHERICAL (1)
#define MKS (2)

// Diagnostic calls
#define DIAG_INIT (0)
#define DIAG_DUMP (1)
#define DIAG_LOG (2)
#define DIAG_FINAL (3)

// Types of restarts
#define RESTART_TEMP (0)
#define RESTART_PERM (1)

// Failure modes
#define FAIL_UTOPRIM (0)
#define FAIL_VCHAR_DISCR (1)
#define FAIL_COEFF_NEG (2)
#define FAIL_COEFF_SUP (3)
#define FAIL_GAMMA (4)
#define FAIL_METRIC (5)

// Geodesic integration and interpolation
#define PUSH_FAIL (0)
#define PUSH_SUCCESS (1)
#define SPH_INTERP_FAIL (0)
#define SPH_INTERP_SUCCESS (1)

// Root finding
#define ROOT_SUCCESS (1)
#define ROOT_FAIL (0)
#define FCOUNT_NBINS (6)
#define FCOUNT_MORE (FCOUNT_NBINS - 1)

// Timers
#define TIMER_UPDATE (0)
#define TIMER_FLUXCALC (1)
#define TIMER_FIXUP (2)
#define TIMER_BOUND (3)
#define TIMER_DIAG (4)
#define TIMER_OUT (5)
#define TIMER_ELECTRON (6)
#define TIMER_MAKE (7)
#define TIMER_PUSH (8)
#define TIMER_INTERACT (9)
#define TIMER_OSCILLATIONS (10)
#define TIMER_MICRO (11)
#define TIMER_ALL (12)
#define NUM_TIMERS (13)

// Units
#define NEED_UNITS (RADIATION || EOS == EOS_TYPE_TABLE)

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
typedef double grid_prim_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NVAR];
typedef double grid_double_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
typedef double grid_radtype_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG]
                                [RAD_NUM_TYPES];
typedef double grid_fourvector_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG]
                                   [NDIM];
typedef double grid_radg_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG]
                             [NDIM + NRADCOMP];
typedef double grid_tensor_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NDIM]
                               [NDIM];
typedef int    grid_int_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
typedef double grid_eosvar_type[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG]
                               [EOS_NUM_EXTRA];
#if RZ_HISTOGRAMS
typedef double rz_hist_type[RZ_HISTOGRAMS_N];
#endif // RZ_HISTOGRAMS
typedef void (*passive_init_ftype)(int, int, int, double *, double *);
typedef double (*hc_ftype)(double, double);

extern grid_prim_type       P;     // Primitive variables
extern grid_prim_type       F1;    // X1 fluxes
extern grid_prim_type       F2;    // X2 fluxes
extern grid_prim_type       F3;    // X3 fluxes
extern grid_prim_type       Ph;    // Half-step primitives
extern grid_prim_type       Psave; // Half-step primitives
extern grid_int_type        pflag; // Failure points
extern grid_int_type        fail_save;
extern grid_int_type        fixup_required;
extern grid_fourvector_type jcon;
extern grid_eosvar_type     extra; // extra variables needed by EOS
#if RADIATION
extern grid_radg_type   radG;     // Radiation four-force
extern grid_radg_type   radG_int; // ...integrated
extern grid_radg_type   radG_buf; // ...buffer for communication
extern grid_tensor_type Rmunu;    // Radiation stress-energy tensor
extern grid_int_type    Nsph;
extern grid_double_type nph;
#if RZ_HISTOGRAMS
extern rz_hist_type rz_r_orig_hist, rz_z_orig_hist;
#if NEUTRINO_OSCILLATIONS
extern rz_hist_type osc_rz_r_orig_hist, osc_rz_z_orig_hist;
#endif // NEUTRINO_OSCILLATIONS
#endif // RZ_HISTOGRAMS

extern struct of_photon **photon_lists;
extern struct of_photon **photon_mpi_lists;

#if DIAGNOSTICS_USE_RADTYPES
#define NULNU_IDX0 (RAD_NUM_TYPES)
#else
#define NULNU_IDX0 (MAXNSCATT + 1)
#endif
extern double nuLnu[NULNU_IDX0][NTH][NPHI][NU_BINS_SPEC];

#if RADIATION == RADTYPE_NEUTRINOS
extern double rad_type_counts[RAD_NUM_TYPES];
extern double lepton_tot, lepton_last, dlepton_tot, dlepton_perc;
extern double lepton_gas, lepton_rad;
extern double lepton_lost, lepton_lost_step, lepton_lost_local;
#endif

extern double Jrad[MAXNSCATT + 2][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
extern double Jrad_buf[MAXNSCATT + 2][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
extern double dtau_avg[RAD_SCATT_TYPES + 1][N1 + 2 * NG][N2 + 2 * NG]
                      [N3 + 2 * NG];
extern double en_int_avg[RAD_SCATT_TYPES + 1][N1 + 2 * NG][N2 + 2 * NG]
                        [N3 + 2 * NG];
extern grid_int_type     Nem, Nabs, Nsc;
extern grid_radtype_type Nem_phys, Nabs_phys, radtype_buf;

extern grid_int_type    Nsuper;
extern grid_double_type Esuper;
extern grid_prim_type   psupersave;

#if LOCAL_ANGULAR_DISTRIBUTIONS
#define LOCAL_NUM_BASES (2)
#define MOMENTS_A (0)
#define MOMENTS_B (1)
#define MOMENTS_DIFF (2)
#define LOCAL_NUM_MOMENTS (3)
typedef double grid_local_angles_type[LOCAL_NUM_BASES][LOCAL_ANGLES_NX1]
                                     [LOCAL_ANGLES_NX2][RAD_NUM_TYPES]
                                     [LOCAL_ANGLES_NMU];
extern grid_local_angles_type local_angles;
extern double                 local_dx1_rad, local_dx2_rad, local_dx_costh;

#if RAD_NUM_TYPES >= 4
typedef double grid_Gnu_type[LOCAL_NUM_BASES][LOCAL_ANGLES_NX1]
                            [LOCAL_ANGLES_NX2][LOCAL_ANGLES_NMU];
typedef double grid_local_moment_type[LOCAL_NUM_BASES][LOCAL_NUM_MOMENTS]
                                     [LOCAL_ANGLES_NX1][LOCAL_ANGLES_NX2];
typedef int grid_local_basis_idx_type[LOCAL_ANGLES_NX1][LOCAL_ANGLES_NX2];
typedef double grid_local_count_type[LOCAL_ANGLES_NX1][LOCAL_ANGLES_NX2];
extern grid_Gnu_type          Gnu, local_Ns, local_wsqr;
extern grid_local_moment_type local_moments;
extern grid_local_basis_idx_type local_b_osc;
extern grid_local_count_type local_osc_count;
#endif // #if RAD_NUM_TYPES >= 4
#endif // LOCAL_ANGULAR_DISTRIBUTIONS

#endif // RADIATION

// Default initialization is 0, which in this case is
// PASSTYPE_INTRINSIC, i.e., volume densities.
extern int                passive_type[NVAR_PASSIVE];
extern char               passive_name[NVAR_PASSIVE][STRLEN];
extern passive_init_ftype do_passive_fixup[NVAR_PASSIVE];

#if ELECTRONS
extern grid_double_type Qvisc, Qcoul;
#endif // ELECTRONS

/*******************************************************************************
    GLOBAL VARIABLES SECTION
*******************************************************************************/
// Command line arguments
extern char outputdir[STRLEN], dumpdir[STRLEN], restartdir[STRLEN];
extern char xmfdir[STRLEN];
#if RADIATION && TRACERS
extern char tracerdir[STRLEN];
#endif
extern char init_from_grmhd[STRLEN];
extern char metric[STRLEN], reconstruction[STRLEN];
extern char eos[STRLEN], nulnutype[STRLEN];
extern int  tracers;

// path for eos
#if EOS == EOS_TYPE_TABLE
extern char eospath[STRLEN];
#endif

// opacity table paths
#if RADIATION
#if (RADIATION == RADTYPE_NEUTRINOS)
#if BURROWS_OPACITIES
extern char opac_param_file[STRLEN];
extern char opac_file[STRLEN];
#endif // BURROWS
#if HDF5_OPACITIES
extern char opac_file[STRLEN];
#endif // hdf5
#endif // neutrinos
#endif // radiation

// Physics parameters
extern double a;
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
extern double gam;
#endif
#if EOS == EOS_TYPE_POLYTROPE
extern double poly_K, poly_gam;
#endif
#if POLYTROPE_FALLBACK
extern double rho_poly_thresh;
#endif
extern double M_unit;
extern double Reh;
extern double Risco;
#if NEED_UNITS
extern double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
#endif
#if EOS == EOS_TYPE_TABLE
extern double TEMP_unit;
#endif
#if RADIATION
extern double Ne_unit, Thetae_unit, kphys_to_num;
extern double tp_over_te, thetae_max, sigma_max, kdotk_tol;
#endif

// Numerical parameters
extern double Rin, Rout, Rout_vis, hslope;
extern double poly_norm, poly_xt, poly_alpha, mks_smooth;
#if RADIATION
extern double Rout_rad;
extern double nph_per_proc;
extern double tune_emiss, t_tune_emiss, dt_tune_emiss;
extern double tune_scatt, t_tune_scatt, dt_tune_scatt;
extern int    made_tune_proc, abs_tune_proc, scatt_tune_proc;
extern double numin, numax;
extern double kappa;
extern double startx_rad[NDIM], stopx_rad[NDIM];
extern double wgtC;
extern int    step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
extern int    tracer_tot, tracer_tot_all, ntcr_per_proc;
extern int    step_sent, step_rcvd, step_fail;
extern int    step_made_all, step_abs_all, step_scatt_all, step_lost_all;
extern int    step_rec_all, step_sent_all, step_rcvd_all, step_tot_all;
extern int    step_fail_all;
extern int    step_tot_max, step_tot_min;
extern double load_imbalance;
extern double Nph_to_track;
extern double sim_vol;
#if FLATEMISS
extern double cnu_flat;
#endif // FLATEMISS
#if MULTISCATT_TEST
extern double ms_theta_nu0, ms_delta0;
#endif // MULTISCATT_TEST
#if RZ_HISTOGRAMS
extern double rz_rmax, rz_zmax, delta_rcyl, delta_z;
#endif // RZ_HISTOGRAMS
#endif // RADIATION
#if ELECTRONS
extern double tptemin, tptemax;
#endif
extern double cour;
#if RADIATION
extern double cour_cool;
#endif
extern double dV, dx[NDIM], startx[NDIM], stopx[NDIM], startx_proc[NDIM],
    stopx_proc[NDIM];
extern double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
extern double dt, dtsave;
extern double t, tf;
extern int    nstep;
extern int    is_restart;

// Output parameters
extern double DTd;
extern double DTl;
extern double DTr;
extern int    DNr;
extern int    DTp;
extern int    DTf;
extern double DTw;
extern int    dump_cnt;
extern int    rdump_cnt;
extern double tdump, trestart, tlog;
extern double root_fcount[FCOUNT_NBINS];

// Global flags
extern int failed;
extern int lim;

// Diagnostics
extern double mdot, mdot_eh;
extern double edot, edot_eh;
extern double ldot, ldot_eh;
extern int    icurr, jcurr, kcurr;

// Parallelism
extern int nthreads;

// Electrons
#if ELECTRONS
extern double game, gamp;
extern double fel0;
#endif

// Set global variables that indicate current local metric, etc.
struct of_geom {
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double g;
  double alpha;
};

struct of_state {
  double ucon[NDIM];
  double ucov[NDIM];
  double bcon[NDIM];
  double bcov[NDIM];
};

#if RADIATION
// WARNING: if you change struct_of_photon, be sure to change
// the the photon hdf5 and MPI  types in io.c and mpi.c
#define NSUP 3
struct of_photon {
  // NSUP >=3 X^{\mu}, K^{\mu}, K_{\mu} so photon data always available anywhere
  // between n and n+1
  double X[NSUP][NDIM];
  double Kcov[NSUP][NDIM];
  double Kcon[NSUP][NDIM];
  double w;
  double KdotKprev;
  // radiation type. For neutrinos, flavor. Always active.
  // Not always important.
  // TODO: make sure to always set type when it's needed.
  int               type;
  int               nscatt;
  int               origin[NDIM];
  double            t0;
  int               is_tracked;
  // Only relevant for neutrino oscillations
  int               osc_count;
  struct of_photon *next;
};

#define PH_ELEM (12)
struct of_track_photon {
  double X1;
  double X2;
  double X3;
  int    nscatt;
};
struct of_microphysics {
#if RADIATION == RADTYPE_NEUTRINOS
  double rho;
  double T;
  double Ye;
  double Abar, Zbar;
  double Xi[NUM_MASS_FRACTIONS];
#if BURROWS_OPACITIES || HDF5_OPACITIES
  double jgrp[RAD_NUM_TYPES][NU_BINS + 1];
  double alphagrp[RAD_NUM_TYPES][NU_BINS + 1];
#endif
#else // RADTYPE_LIGHT
  double Thetae;
  double Ne;
#endif
  // needed for Bk-angle (which isn't really needed for neutrinos, but
  // this is the easy engineering solution)
  double B;
};
typedef double (*dsdom_ftype)(
    double, double, const struct of_microphysics *, int, int);
#endif // RADIATION

#if EOS == EOS_TYPE_TABLE
struct of_adiabat {
  double  s, ye;
  double  lrho_min, lrho_max;
  int     imin, imax;
  double  hm1_min, hm1_max;
  double *lT;
};
#endif

#if EOS == EOS_TYPE_TABLE || HDF_OPACITIES || BURROWS_OPACITIES
struct of_tablebounds {
  int    Nrho, NT, NYe;
  double lrho_min, lrho_max, dlrho;
  double lT_min, lT_max, dlT;
  double Ye_min, Ye_max, dYe;
};
#endif

// More grid functions. Axisymmetry assumed.
extern double         conn[N1 + 2 * NG][N2 + 2 * NG][NDIM][NDIM][NDIM];
extern struct of_geom ggeom[N1 + 2 * NG][N2 + 2 * NG][NPG];
#if RADIATION
// extern double dt_light, dt_light_min;
extern double                 dt_light[N1 + 2 * NG][N2 + 2 * NG], dt_light_min;
extern struct of_microphysics m_grd[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
extern grid_fourvector_type   Ucon_grd, Ucov_grd, Bcon_grd, Bcov_grd;
#endif

// MPI-specific stuff
extern int global_start[NDIM];
extern int global_stop[NDIM];

/*******************************************************************************
    MACROS
*******************************************************************************/
#define ILOOP for (int i = 0 + NG; i < N1 + NG; i++)
#define ILOOPALL for (int i = 0; i < N1 + 2 * NG; i++)
#define JLOOP for (int j = 0 + NG; j < N2 + NG; j++)
#define JLOOPALL for (int j = 0; j < N2 + 2 * NG; j++)
#define KLOOP for (int k = 0 + NG; k < N3 + NG; k++)
#define KLOOPALL for (int k = 0; k < N3 + 2 * NG; k++)
#define ZLOOP                              \
  for (int i = 0 + NG; i < N1 + NG; i++)   \
    for (int j = 0 + NG; j < N2 + NG; j++) \
      for (int k = 0 + NG; k < N3 + NG; k++)
#define ZLOOPALL ILOOPALL JLOOPALL KLOOPALL
#define ISLOOP(istart, istop) for (int i = istart + NG; i <= istop + NG; i++)
#define JSLOOP(jstart, jstop) for (int j = jstart + NG; j <= jstop + NG; j++)
#define KSLOOP(kstart, kstop) for (int k = kstart + NG; k <= kstop + NG; k++)
#define ZSLOOP(istart, istop, jstart, jstop, kstart, kstop) \
  for (int i = istart + NG; i <= istop + NG; i++)           \
    for (int j = jstart + NG; j <= jstop + NG; j++)         \
      for (int k = kstart + NG; k <= kstop + NG; k++)
// Loop over faces
#define FACELOOP for (int face = FACESTART; face < FACEEND; face++)
// Loop over all locations
#define LOCLOOP for (int loc = PGSTART; loc < PGEND; loc++)
// Loop over primitive variables
#define PLOOP for (int ip = 0; ip < NVAR; ip++)
#define BASELOOP for (int ip = 0; ip < NVAR_BASE; ip++)

// Loop over spacetime indices
#define DLOOP1 for (int mu = 0; mu < NDIM; mu++)
#define DLOOP2                      \
  for (int mu = 0; mu < NDIM; mu++) \
    for (int nu = 0; nu < NDIM; nu++)
#define DLOOP3                        \
  for (int mu = 0; mu < NDIM; mu++)   \
    for (int nu = 0; nu < NDIM; nu++) \
      for (int sigma = 0; sigma < NDIM; sigma++)
// spacelike indices
#define SDLOOP for (int mu = 1; mu < NDIM; mu++)

// Loop over extra variables
// TODO: Figure out how to make this conditionally defined. ~JMM
#define EOS_ELOOP for (int e = 0; e < EOS_NUM_EXTRA; e++)

// Loop over passive scalars
#define PASSLOOP for (int ipass = PASSIVE_START; ipass < PASSIVE_STOP; ipass++)
#define PASSELEM(ipass) (ipass - PASSIVE_START)
#define PASSTYPE(ipass) (passive_type[PASSELEM(ipass)])
#define PASSNAME(ipass) (passive_name[PASSELEM(ipass)])

#define TYPELOOP for (int itp = RAD_TYPE_START; itp < RAD_NUM_TYPES; itp++)
#define SCATTLOOP for (int iscatt = 0; iscatt < RAD_SCATT_TYPES; iscatt++)
#define JRADLOOP for (int n = 0; n < MAXNSCATT + 2; n++)
#define NULOOP for (int inu = 0; inu < NU_BINS + 1; inu++)

#define LOCALXLOOP                             \
  for (int i = 0; i < LOCAL_ANGLES_NX1; ++i)   \
    for (int j = 0; j < LOCAL_ANGLES_NX2; ++j)
#define LOCALMULOOP for (int imu = 0; imu < LOCAL_ANGLES_NMU; ++imu)
#define LOCALXMULOOP LOCALXLOOP LOCALMULOOP

#define MY_MIN(fval1, fval2) (((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1, fval2) (((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) (((fval) < 0.) ? -1. : 1.)

#define delta(i, j) ((i == j) ? 1. : 0.)
#define dot(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3])

/*******************************************************************************
    FUNCTION DECLARATIONS
*******************************************************************************/

// bounds.c
void bound_prim(grid_prim_type prim);
void fix_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
#if RADIATION
void bound_superphotons(grid_prim_type P, double t, double dt);
void polar_fix(double X[NDIM], struct of_photon *ph);
int  rad_error_check(struct of_photon **ph);
// int rad_boundary_transport(double X[NDIM], double Kcov[NDIM]);
void bound_rad_transport(
    double X[NDIM], struct of_photon *ph, int is_transported);
int bound_rad_isactive(double X[NDIM], struct of_photon *ph);
int rad_mpi_transport(struct of_photon **ph, struct of_photon **prev,
    struct of_photon **head, double X[NDIM], int active);
#endif

// coord.c
void   coord(int i, int j, int k, int loc, double X[NDIM]);
double r_of_X(const double X[NDIM]);
double th_of_X(const double X[NDIM]);
void   jac_harm_to_bl(
      const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void jac_bl_to_cart(
    const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void jac_harm_to_cart(
    const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]);
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]);
void bl_coord(const double X[NDIM], double *r, double *th);
int  bl_i_of_r(double r);
void cart_coord(const double X[NDIM], double Xcart[NDIM]);
void set_gcov(double X[NDIM], double gcov[NDIM][NDIM]);
void set_metric(double X[NDIM], struct of_geom *g);
void set_points();
void zero_arrays(void);
void set_grid(void);
void cart_to_sph(const double X[NDIM], double *r, double *th, double *ph);

// current.c
void current_calc();

// diag.c
void   reset_log_variables();
void   reset_dump_variables();
void   diag(int call_code);
void   fail(int fail_type);
void   area_map(int i, int j, int k, grid_prim_type prim);
void   diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
double flux_ct_divb(int i, int j, int k);
#if RADIATION
void record_superphoton(double X[NDIM], struct of_photon *ph);
void bin_all_superphotons();
void report_load_imbalance();
#if RADIATION == RADTYPE_NEUTRINOS
void print_rad_types();
void count_leptons(grid_prim_type P, double dt, int nstep);
#endif
#if RZ_HISTOGRAMS
void generate_rz_histograms();
#endif // RZ_HISTOGRAMS
#endif // RADIATION

// electrons.c
#if ELECTRONS
void init_electrons();
void heat_electrons(
    grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt);
double get_fel(int i, int j, int k, double p[NVAR]);
#if RADIATION
void coulomb(
    grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt);
void apply_rad_force_e(
    grid_prim_type Prh, grid_prim_type Pr, grid_radg_type radG, double Dt);
#endif // RADIATION
void fixup_electrons(grid_prim_type p);
#endif

// emissivity.c
#if RADIATION
double jnu(double nu, int type, const struct of_microphysics *m, double theta);
double Jnu(double nu, int type, const struct of_microphysics *m);
double get_J(struct of_microphysics *m);
void   init_emissivity();
#endif

// eos.c
void   init_EOS();
double EOS_bad_eos_error();
double EOS_get_gamma(const double *extra);
double EOS_pressure_rho0_u(double rho, double u, const double *extra);
double EOS_pressure_rho0_w(double rho, double w, double gamma,
    const struct of_geom *geom, double *extra);
double EOS_enthalpy_rho0_u(double rho, double u, const double *extra);
double EOS_sound_speed_rho0_u(double rho, double u, const double *extra);
void   EOS_set_floors(double scale, double rho, double u, double bsq,
      double *rhoflr, double *uflr, const double *extra);
double EOS_entropy_rho0_u(double rho, double u, const double *extra);
double EOS_adiabatic_constant(double rho, double u, const double *extra);
double EOS_temperature(double rho, double u, const double *extra);
double EOS_u_press(double press, double rho, double *extra);
#if RADIATION
double EOS_u_N_Theta(double rho, double N, double Theta, double *extra);
double EOS_pressure_N_Theta(double N, double Theta);
double EOS_Theta_unit();
#endif // RADIATION

// eos_gamma.c
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
double EOS_Gamma_pressure_rho0_u(double rho, double u);
double EOS_Gamma_pressure_rho0_w(double rho, double w);
double EOS_Gamma_entropy_rho0_u(double rho, double u);
double EOS_Gamma_enthalpy_rho0_u(double rho, double u);
double EOS_Gamma_adiabatic_constant_rho0_u(double rho, double u);
double EOS_Gamma_sound_speed_rho0_u(double rho, double u);
void   EOS_Gamma_set_floors(double scale, double rho, double u, double bsq,
      double *rhoflr, double *uflr);
double EOS_Gamma_rho_floor(double scale, double bsq);
double EOS_Gamma_u_floor(double scale, double bsq);
double EOS_Gamma_u_scale(double rho);
double EOS_Gamma_u_press(double press);
double EOS_Gamma_temp(double rho, double u);
#if RADIATION
double EOS_Gamma_Theta_unit();
#endif // RADIATION
#endif // EOS_TYPE_GAMMA

// eos_poly.c
#if EOS == EOS_TYPE_POLYTROPE || POLYTROPE_FALLBACK
double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam);
double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam);
double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam);
void   EOS_Poly_set_floors(double scale, double rho, double u, double bsq,
      double *rhoflr, double *uflr);
double EOS_Poly_rho_floor(double scale, double bsq);
double EOS_Poly_u_floor(double scale, double bsq);
double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam);
#endif // EOS_TYPE_POLY

// eos_stellar_collapse.c
#if EOS == EOS_TYPE_TABLE
void   EOS_SC_init(char *name);
void   EOS_SC_fill(double *restrict p, double *restrict eos);
double EOS_SC_pressure_rho0_u(double lrho, double lT, double ye);
double EOS_SC_pressure_rho0_w(double rho, double w, double ye, double *lTold);
double EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye);
double EOS_SC_sound_speed(double lrho, double lT, double ye);
double EOS_SC_entropy(double lrho, double lT, double ye);
double EOS_SC_gamma(double lrho, double lT, double ye);
double EOS_SC_temperature(double lT);
double EOS_SC_get_u_of_T(double rho, double T, double ye);
double EOS_SC_u_press(double press, double rho, double ye, double *lTold);
void EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double *extra);
void EOS_SC_avg_ions(double *Abar, double *Zbar, const double *extra);
void EOS_SC_set_floors(double scale, double rho, double u, double ye,
    double bsqr, double *rhoflr, double *uflr);
double EOS_SC_rho_floor(double scale, double bsq);
double EOS_SC_u_floor(double scale, double bsq, double ye);
double EOS_SC_get_min_lrho();
double EOS_SC_get_min_rho();
double EOS_SC_get_min_lT();
double EOS_SC_get_minu(double rho, double ye, double scale);
void   EOS_SC_get_polytrope(
      double lrho, double lT, double ye, double *poly_K, double *poly_gamma);
double EOS_SC_hm1_min_adiabat(const struct of_adiabat *a);
int    EOS_SC_find_adiabat_1d(double s, double ye, struct of_adiabat *a);
void   EOS_SC_print_adiabat(const struct of_adiabat *a);
void   EOS_SC_adiabat_free(struct of_adiabat *a);
void   EOS_SC_isoentropy_hm1(double hm1, const struct of_adiabat *a,
      double *lrho_guess, double *rho, double *u);
void   do_ye_fixup(
      int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]);
void do_ye_em_fixup(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]);
void do_atm_fixup(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]);
void EOS_SC_get_bounds(struct of_tablebounds *b);
void EOS_SC_print_table_mins();
#endif // EOS_TYPE_TABLE

// estimate_thetae.c
#if RADIATION
double get_Thetae_est(int i, int j, int k);
void   estimate_Thetae(
      grid_prim_type P, grid_eosvar_type extra, double t, double dt);
#endif

// fixup.c
void   fixup(grid_prim_type pv, grid_eosvar_type extra);
double get_scale(int i, int j, int k);
void   fixup1zone(
      int i, int j, int k, double pv[NVAR], double extra[EOS_NUM_EXTRA]);
void fixup_utoprim(grid_prim_type pv, grid_eosvar_type extra);

// interact.c
#if RADIATION
void interact(grid_prim_type P, grid_eosvar_type extra, double t, double dt);
#endif

// input.c
void init_params(char *pfname);
void set_param(char *key, void *data);

// io.c
// void set_core_params();
// void set_param(char *key, void *data);
// void read_params(char *pfname);
void init_io();
void init_fluid_restart();
#if RADIATION
void track_ph();
#if TRACERS
void dump_tracers();
#endif
#endif
void dump_grid();
void dump();
void restart_write(int restart_type);
void restart_read(char *fname);
int  restart_init();

// make_superphotons.c
#if RADIATION
void make_superphotons(
    grid_prim_type Prad, grid_eosvar_type extra, double t, double dt);
void set_weight(grid_prim_type Prad, grid_eosvar_type extra);
void get_dnz(grid_prim_type Prad, grid_eosvar_type extra);
#endif

// metric.c
double gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]);
void   conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM]);
void   lower(double ucon[NDIM], double gcov[NDIM][NDIM], double ucov[NDIM]);
void   raise(double ucov[NDIM], double gcon[NDIM][NDIM], double ucon[NDIM]);
struct of_geom *get_geometry(int ii, int jj, int kk, int loc);
void            blgset(int i, int j, struct of_geom *geom);
double          bl_gdet_func(double r, double th);
void            bl_gcov_func(double r, double th, double gcov[][NDIM]);
void            bl_gcon_func(double r, double th, double gcon[][NDIM]);
double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void   adjoint(double m[16], double adjOut[16]);
double determinant(double m[16]);
double invert(double *m, double *invOut);

// mpi.c
void init_mpi();
void sync_mpi_boundaries_X1L(grid_prim_type Pr);
void sync_mpi_boundaries_X1R(grid_prim_type Pr);
void sync_mpi_boundaries_X2L(grid_prim_type Pr);
void sync_mpi_boundaries_X2R(grid_prim_type Pr);
void sync_mpi_boundaries_X3L(grid_prim_type Pr);
void sync_mpi_boundaries_X3R(grid_prim_type Pr);
#if RADIATION
void sync_radG();
void sync_Jrad();
void sync_radtype_vec(grid_radtype_type v);
void sync_mpi_photons(
    struct of_photon **ph_mpi, grid_prim_type P, double t, double dt);
void mpi_reduce_nuLnu();
#endif
int    mpi_nprocs();
double mpi_max(double f);
double mpi_min(double f);
int    mpi_max_int(int f);
double mpi_reduce(double f);
int    mpi_reduce_int(int f);
void   mpi_dbl_allreduce_array(double *A, int size);
int    mpi_io_proc();
void   mpi_int_broadcast(int *val);
void   mpi_int_broadcast_array(int *val, int size);
void   mpi_int_broadcast_proc(int *val, int root);
void   mpi_dbl_broadcast(double *val);
void   mpi_allgather_int1(int *buffer, int element);
int    mpi_accumulate_int(int my_val);
double mpi_io_reduce(double val);
double mpi_io_max(double val);
int    mpi_myrank();
void   mpi_sync_output();
void   mpi_barrier();
int    mpi_is_periodic(int dir);

// opac_emis_neutrino.c
#if RADIATION
#if RADIATION == RADTYPE_NEUTRINOS && BURROWS_OPACITIES
void   init_opac_emis_burrows();
void   fill_opac_emis_burrows(struct of_microphysics *m);
double jnu_burrows(double nu, int type, const struct of_microphysics *m);
double Jnu_burrows(double nu, int type, const struct of_microphysics *m);
double int_jnudnudOmega_burrows(const struct of_microphysics *m);
double alpha_nu_burrows(double nu, int type, const struct of_microphysics *m);
void   test_opac_emis();
#if EOS == EOS_TYPE_TABLE
void opac_emis_to_hdf(const char *name, const struct of_tablebounds *b);
#endif // EOS
#endif // Burrows opacities
#if RADIATION == RADTYPE_NEUTRINOS && HDF5_OPACITIES
void   init_opac_emis_hdf(char *name);
void   fill_opac_emis_hdf(struct of_microphysics *m);
double jnu_hdf(double nu, int type, const struct of_microphysics *m);
double Jnu_hdf(double nu, int type, const struct of_microphysics *m);
double int_jnudnudOmega_hdf(const struct of_microphysics *m);
double alpha_nu_hdf(double nu, int type, const struct of_microphysics *m);
#endif // HDF opacities

// oscillations.c
#if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS
double get_dt_oscillations();
void get_local_angle_bins(
    struct of_photon *ph, int *pi, int *pj, int *pmu1, int *pmu2);
void accumulate_local_angles();
#if RAD_NUM_TYPES >= 4
void compute_local_gnu(grid_local_angles_type local_angles,
    grid_Gnu_type local_Ns, grid_Gnu_type local_wsqr, grid_Gnu_type gnu);
void compute_local_moments(grid_Gnu_type gnu, grid_local_moment_type moments);
void oscillate(grid_local_moment_type local_moments, grid_Gnu_type gnu);
#endif // RAD_NUM_TYPES >= 4
#endif // LOCAL_ANGULAR_DISTRIBUTIONS
#endif // RADIATION

// passive.c
//#if NVAR_PASSIVE > 0
void fixup_passive(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]);
void init_passives();
void name_passives();
//#endif

// phys.c
void   primtoflux(double *pr, struct of_state *q, int dir, int magnetic,
      struct of_geom *geom, double *flux);
void   bcon_calc(double *pr, double *ucon, double *ucov, double *bcon);
void   mhd_calc(double *pr, int dir, int magnetic, struct of_state *q,
      double *extra, double *mhd);
void   source(double *ph, struct of_geom *geom, int ii, int jj, double *dU,
      double Dt, double *extra);
double bsq_calc(double *pr, struct of_geom *geom);
void   get_state(double *pr, struct of_geom *geom, struct of_state *q);
void   ucon_calc(double *pr, struct of_geom *geom, double *ucon);
int    mhd_gamma_calc(double *pr, struct of_geom *geom, double *gamma);
void   mhd_vchar(double *pr, struct of_state *q, struct of_geom *geom, int js,
      double *vmax, double *vmin);

// problem.c
void set_problem_params();
void init_prob();
void bound_gas_prob_x1l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x1r(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x2l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x2r(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x3l(int i, int j, int k, grid_prim_type P);
void bound_gas_prob_x4r(int i, int j, int k, grid_prim_type P);

// push_superphotons.c
#if RADIATION
int push_X_K(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
    grid_prim_type P, grid_prim_type Prad, double KdotKprev, int type,
    double dtpush);
int push_superphoton(
    struct of_photon *ph, grid_prim_type P, grid_prim_type Prad, double dtpush);
void push_superphotons(grid_prim_type P, grid_prim_type Prad, double dt);
#endif

// rad_utils.c
#if RADIATION
unsigned long int count_particles_local();
void              init_rad(grid_prim_type Prad);
void              init_superphoton_resolution();
void update_superphoton_resolution(grid_prim_type Prad, grid_eosvar_type extra);
double linear_interp_log(double x, double *table, double lmin, double dl);
void   list_remove(struct of_photon **ph, struct of_photon **ph_head,
      struct of_photon **ph_prev);
double get_Thetae(double P[NVAR]);
double scatterer_dimensionless_temp(
    int radtype, int interaction, const struct of_microphysics *m);
double scatterer_number_density(
    int radtype, int interaction, const struct of_microphysics *m);
void precompute_microphysics();
void get_fluid_zone(int i, int j, int k, grid_prim_type Prad,
    grid_eosvar_type extra, struct of_microphysics *m, double Ucon[NDIM],
    double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM]);
int  is_null(double Kcov[NDIM], double Kcon[NDIM], double K0, double KdotKprev,
     double *KdotK);
void set_Rmunu();
// int is_null(struct of_photon *ph, double *KdotK);
void   Xtoijk(double X[NDIM], int *i, int *j, int *k);
void   copy_photon(struct of_photon *ph, struct of_photon *phc);
void   print_ph_diag(struct of_photon *ph);
int    get_X_K_interp(struct of_photon *ph, double t_interp, grid_prim_type P,
       double X[NDIM], double Kcov[NDIM], double Kcon[NDIM]);
int    to_be_pushed(double t, double dt, struct of_photon *ph);
double get_dtpush(struct of_photon *ph, double dt);
void   swap_ph(struct of_photon **donor, struct of_photon **recipient);
void   get_nuLnu_bin(double X[NDIM], int *thbin, int *phibin);
void   bin_superphoton_direction(const struct of_photon *ph);
double get_min_dt_cool(grid_prim_type P, grid_eosvar_type extra);
void   set_cooling_time(
      grid_double_type tau_cool, grid_prim_type P, grid_eosvar_type extra);
#if RADIATION == RADTYPE_NEUTRINOS
void record_lepton_flux(const struct of_photon *ph);
void check_nu_type(const char *location);
int  get_lepton_sign(const struct of_photon *ph);
int  nu_is_heavy(const int radtype);
#endif // NEUTRINOS
#endif // RADIATION

// radiation.c
#if RADIATION
double Bnu_inv(double nu, const struct of_microphysics *m); // TODO?
double jnu_inv(
    double nu, int type, const struct of_microphysics *m, double theta);
double alpha_inv_scatt(
    double nu, int type, int interaction, const struct of_microphysics *m);
double alpha_inv_abs(
    double nu, int type, const struct of_microphysics *m, double theta);
double get_fluid_nu(double X[NDIM], double Kcov[NDIM], double Ucon[NDIM]);
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
    double Bcov[NDIM], double B);
#endif

// random.c
void   init_random(int seed);
double get_rand();
double get_chisq(double nu);
double get_gaussian(double mu, double sigma);
void   get_ran_dir_3d(double *nx, double *ny, double *nz);

// reconstruction.c
void reconstruct(double ptmp[NMAX + 2 * NG][NVAR], int N,
    double p_l[NMAX + 2 * NG][NVAR], double p_r[NMAX + 2 * NG][NVAR]);

// root_finding.c
// counts root_finding iterations
void initialize_root_fcounts();
void print_root_fcounts();
// solves for f(x,params) - ytarget = 0
int find_root(double (*f)(const double, const void *), const void *params,
    const double ytarget, double xguess, const double xmin, const double xmax,
    const double xtol, const double ytol, double *xroot);

// scattering.c
#if RADIATION
int    scatt_temp_too_small(const struct of_microphysics *m);
int    scatter_superphoton(grid_prim_type P, grid_eosvar_type extra,
       struct of_photon *ph, double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
       int interaction);
void   init_all_hotcross();
double total_cross_lkup(
    double w, int type, int interaction, const struct of_microphysics *m);
#endif

// step.c
void step();

// tetrads.c
#if RADIATION
void coord_to_tetrad(
    double Ecov[NDIM][NDIM], double Kcoord[NDIM], double Ktetrad[NDIM]);
void tetrad_to_coord(
    double Econ[NDIM][NDIM], double Ktetrad[NDIM], double Kcoord[NDIM]);
void make_tetrad(int i, int j, int k, double Ucon[NDIM], double trial[NDIM],
    double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]);
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]);
void normalize_null_cov(double Gcon[NDIM][NDIM], double K[NDIM]);
#endif

// timing.c
void   time_set(int n, double val);
double time_read(int n);
void   time_init();
void   timer_start(int timerCode);
void   timer_stop(int timerCode);
double get_time_per_step(int timerCode);
void   report_performance();
void   timers_reset();

// tracers.c
#if RADIATION && TRACERS
double get_total_tracer_mass();
int    tracer_max_id(long unsigned int ntracers);
int    count_tracers_local();
void   lagrange_interp_grid_3d(
      double v_interp[], double X[NDIM], double *v, int size);
void lagrange_interp_prim_3d(
    double P_interp[NVAR], double X[NDIM], grid_prim_type P);
void   push_tracers(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
      grid_prim_type P, grid_prim_type Ph, double dt);
int    tracer_get_id(struct of_photon *ph);
double tracer_get_mass(struct of_photon *ph);
void   tracer_get_X_u(struct of_photon *ph, double t_interp, double X[NDIM],
      double ucov[NDIM], double ucon[NDIM]);
void   set_tracer(struct of_photon *ph, int id, int nstep, double t,
      double X[NDIM], double mass, grid_prim_type P);
void   make_tracer(struct of_photon **head, int id, int nstep, double t,
      double X[NDIM], double mass, grid_prim_type P);
void   sample_tracers_in_cell(struct of_photon **head, int nstep, int i, int j,
      int k, int *last_id, int ntracers_per_cell, double t, grid_prim_type P);
void   sample_all_tracers(long unsigned int ntracers_tot, int nstep, double t,
      grid_int_type tcrs_in_cell, grid_prim_type P);
void   prune_tracers();
#endif

// util.c
double find_min(const double *array, int size);
double find_max(const double *array, int size);
double find_median(double *array, int size);
double interp_1d(double x, const double xmin, const double xmax, const int imin,
    const int imax, const double *restrict tab_x, const double *restrict tab_y);
int    find_index(double value, const double *array, int size);
void * safe_malloc(size_t size);
void   safe_system(const char *command);
void   safe_fscanf(FILE *stream, const char *format, ...);
int    is_practically_nan(double v);
#if NEED_UNITS
void set_units();
#endif

// utop.c
int Utoprim(double U[NVAR], struct of_geom *geom, double prim[NVAR]);

// xdmf_output.c
void write_xml_file(int dump_id, double t, const char *vnams[NVAR]);
