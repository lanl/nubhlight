/******************************************************************************
 *                                                                            *
 * EOS_STELLAR_COLLAPSE.h                                                     *
 *                                                                            *
 * PROTOTYPES FOR READING EOS TABLES PROVIDED ON STELLARCOLLAPSE.ORG          *
 *                                                                            *
 ******************************************************************************/

/* The header for stellar collapse simply got too long, so I split it
 * into a separate file. These declarations belong here, not in decs.h
 * because they are intended to have private scope.
 * ~JMM
 */

#pragma once

// decs
#include "decs.h"

#if EOS == EOS_TYPE_TABLE

#define TABLE_TOL (1.e-10)
#define TABLE_FTOL (1.e-10)

// Debug settings
#define SC_DEBUG (0)
#define SC_MONOTONE_SAFE (1)

// Energy offset/norm settings
/*
 * If USE_ENERGY_SHIFT == 1,
 * then minimum specific internal energy in table can be negative
 * as defined by table. So can u = rho*epsilon

 * If OFFSET_U == 1
 * then shift total energy of system up by small amount
 * so that u cannot be negative.

 * If OFFSET_H_CS == 1
 * then use (possibly) negative specific internal energy
 * to calculate denominator for sound speed.
 * This results in positive u and epsilon, but C_s^2 is calculated
 * as if they could be negative.
 * Turn this on ONLY if USE_ENERGY_SHIFT == OFFSET_U == 0
 *
 * I don't know which of these settings is physically correct.
 * The consequences are non-negligible:
 *
 * - The sound speed can vary up to about 0.1%.
 *   but only when gas is cold.
 *   I don't know which sound speed is correct.
 *
 * - Negative u in the atmosphere, which you get with
 *   ENERGY_SHIFT == 1 and OFFSET_U == 0
 *   can cause (solveable) stability problems
 *   this can be mitigated by demanding u > 0
 *   no matter what the EOS says.
 *   to do this set OLD_U_FLORS == 1.
 *
 * - Setting ENERGY_SHIFT == OFFSET_U == 0
 *   changes the sound speed and amounts to shifting
 *   u by energy_shift*rho everywhere in the domain.
 *   For hot gas, this shift is negligible. For cold gas,
 *   it is of order 0.1%.
 *
 * - Enabling shifts in both u and epsilon is probably
 *   the most physically correct thing to do. This
 *   amounts to using negative minimum specific internal energy,
 *   as specified in the table, but shifting the total energy
 *   in the system, essentially resetting the zero point.
 *   Unfortunately, this produces
 *   floors for u which are perhaps unacceptably large,
 *   of order 10^{-7} throughout the entire domain.
 *
 * - Setting USE_ENERGY_SHIFT == 0 and OFFSET_U == 1
 *   changes the enthalpy, both specific and by-volume.
 *   A shift in the enthalpy produces a shift in
 *   functions used for root finding in utop.
 *   however this appears to be a wash.
 *
 * - These settings DO NOT change pressure or temperature in the gas.
 *   They also have no effect on the equilibrium torus initial data.
 *
 * For now, I've made what I believe is a conservative choice
 * which balances these considerations. The sound speed is
 * calculated as if the specific internal energy contains
 * a shift in it but otherwise there are no shifts
 * and specific internal energy and u have zero as their minima
 * to get these settings, use:
 * USE_ENERGY_SHIFT = 0
 * OFFSET_U         = 0
 * OFFSET_H_CS      = 1
 * OLD_U_FLOORS     = 0
 *
 * Another reasonable choice of settings is given below
 * This choice is safe if a sort of hot atmosphere
 * is acceptable. And/or if the temperature of the
 * wind is high.
 * USE_ENERGY_SHIFT = 1
 * OFFSET_U         = 1
 * OFFSET_H_CS      = 0
 * OLD_U_FLOORS     = 0
 *
 * A final reasonable choice is to use the energies
 * that the EOS suggests, but demand that u in the atmosphere
 * must be positive, even if the EOS doesn't like it.
 * To do this, set
 * USE_ENERGY_SHIFT = 1
 * OFFSET_U         = 0
 * OFFSET_H_CS      = 0
 * OLD_U_FLOORS     = 1
 *
 * ~JMM
 */
#define USE_ENERGY_SHIFT (0)
#define OFFSET_U (0)
#define OFFSET_H_CS (1)
// Demand u > 0.
// Good to set with USE_ENERGY_SHIFT=1 and
// OFFSET_U == OFFSET_H_CS == 0
#define OLD_U_FLOORS (0)

// Demand rho >= u/UORHOMAX
// honestly I have no idea
// if this should be on or off
#define BOUND_UORHO (0)

// Table cleanup settings
// Smooth out tables with median filter
#define DO_SMOOTH (1)
#define EPSSMOOTH (10)
#define MF_W (3)
#define MF_S ((2 * MF_W + 1) * (2 * MF_W + 1) * (2 * MF_W + 1))

#define EOS_ELEM(irho, iT, iY) (Nrho * ((iY)*NT + (iT)) + (irho))
#define MMA_ELEM(irho, iY) (Nrho * iY + irho)

static int     Nrho, NT, NYe;
static double *tab_lrho;
static double *tab_lT;
static double *tab_Ye;
static double *tab_tmp;
static double *tab_lP;
static double *tab_ent;
static double *tab_dpderho;
static double *tab_dpdrhoe;
static double *tab_cs2;
static double *tab_le;
static double *tab_Xa;
static double *tab_Xh;
static double *tab_Xn;
static double *tab_Xp;
static double *tab_Abar;
static double *tab_Zbar;
static double *tab_lwmrho;     // log enthalpy - rho, by volume
static double *tab_hm1;        // enthalpy - 1, by mass
static double *tab_lhm1;       // log enthalpy - 1, by mass
static double *tab_poly_gamma; // Polytrope gamma
static double *tab_poly_K;     // polytrope K

// min and max of wmrho given fixed ilrho and iY
static double *tab_le_min_2d;
static double *tab_le_max_2d;
static double *tab_lP_min_2d;
static double *tab_lP_max_2d;
static double *tab_lwmrho_min_2d;
static double *tab_lwmrho_max_2d;
// static double* tab_hm1_min_1d;

static double tab_lrho_min, tab_lrho_max;
static double tab_lT_min, tab_lT_max;
static double tab_Ye_min, tab_Ye_max;
static double tab_dlrho, tab_dlT, tab_dYe;

static double tab_lhm1_min, tab_lhm1_max;
static double tab_lP_min, tab_lP_max;
static double tab_ent_min, tab_ent_max;
static double tab_cs2_min, tab_cs2_max;
static double tab_le_min, tab_le_max;
static double tab_Xa_min, tab_Xa_max;
static double tab_Xh_min, tab_Xh_max;
static double tab_Xn_min, tab_Xn_max;
static double tab_Xp_min, tab_Xp_max;
static double tab_Abar_min, tab_Abar_max;
static double tab_Zbar_min, tab_Zbar_max;
static double tab_dpderho_min, tab_dpderho_max;
static double tab_dpdrhoe_min, tab_dpdrhoe_max;
static double tab_lwmrho_min, tab_lwmrho_max;

static double tab_rho_min, tab_rho_max;
static double tab_T_min, tab_T_max;
static double tab_e_min, tab_e_max;
static double tab_P_min, tab_P_max;
static double tab_wmrho_min, tab_wmrho_max;
static double tab_hm1_min, tab_hm1_max;
static double tab_u_min;

static double e_shift;
static double u_shift;
static double w_shift;
static double h_shift;
static double h_shift_cs;

// core
static double EOS_SC_interp(
    const double lrho, const double lT, const double Ye, double *restrict tab);

// max, min, etc.
// static void fill_min_1d(double* tab_min_1d, double* tab);
static void fill_max_min_2d(
    double *tab_max_2d, double *tab_min_2d, double *tab);
static double interp_2d(
    const double lrho, const double Ye, const double *tab_2d);
static void   temp_map(double lrho, double Ye, const double *tab);
static double catch_var_2d(const double lrho, const double Ye, const double var,
    const double *tab_min_2d, const double *tab_max_2d);
static double catch_var_2d_monotone(
    const double lrho, const double Ye, const double var, double *restrict tab);

// Root finding
static int  find_lT(const double lrho, double lTguess, const double ye,
     double *restrict tab, const double val, double *lT);
static void adiabat_lrho_min_max(
    double s, double ye, double *lrho_min, double *lrho_max);
static int find_adiabat_0d(
    double lrho, double lTguess, double ye, double s, double *lT);
static double lT_f(const double lT, const void *params);
struct of_lT_params {
  double  lrho, ye;
  double *restrict tab;
};
static double lT_f_adiabat(double lT, const void *params);
struct of_lT_adiabat_params {
  double *restrict         tab;
  const struct of_adiabat *a;
};
static double ent_f_adiabat(double lrho, const void *params);
struct of_ent_f_params {
  double lT, ye;
};

// fix borked tables
#if DO_SMOOTH
void fill_mdn_buffer(double buffer[], int width, int iYe, int iT, int irho,
    double *restrict tab);
static void median_filter(double *tab_in, double *tab_out);
#endif

// utilities
static double u2e(const double rho, const double u);
static double e2u(const double rho, const double e);
static double le2e(const double le);
static double e2le(const double e);
static double h2lh(const double h);
// static double lh2h(const double lh);
// static double lw2w(const double lw);
static double w2lw(const double w);
static double catch_rho(const double rho);
static double catch_lrho(const double lrho);
static double catch_e(const double e);
static double catch_press(const double press);
static double catch_w(const double w);
// static double catch_h(const double h);
static double catch_ye(const double ye);
// static double catch_temp(const double temp);
static double catch_lT(const double lT);
static double catch_s(const double s);
static double catch_hm1(const double hm1);

#endif // EOS_TYPE_TABLE
