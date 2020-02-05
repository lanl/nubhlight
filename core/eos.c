/******************************************************************************
 *                                                                            *
 * EOS.C                                                                      *
 *                                                                            *
 * DEFINES THE EQUATION OF STATE AND RELATED MICROPHYSICS                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

/*******************************************************************************
      EOS-Independent
*******************************************************************************/
#if RADIATION
// TODO: Is this right?
double EOS_pressure_N_Theta(double N, double Theta) {
  double press;
#if EOS == EOS_TYPE_GAMMA
  press = N * Theta * ME * CL * CL;
#elif EOS == EOS_TYPE_TABLE
  press = N * Theta * MP * CL * CL;
#else
  EOS_bad_eos_error();
#endif
  return press;
}

double EOS_u_N_Theta(double rho, double N, double Theta, double *extra) {
  double press = EOS_pressure_N_Theta(N, Theta);
  double u     = EOS_u_press(press, rho, extra);
  return u;
}
#endif

double EOS_bad_eos_error() {
  fprintf(stderr, "ERROR! UNKNOWN EOS TYPE! TYPE = %i\n", EOS);
  exit(1);
}
/******************************************************************************/

/*******************************************************************************
      Wrappers
*******************************************************************************/
void init_EOS() {
#if EOS == EOS_TYPE_GAMMA || EOS == EOS_TYPE_POLYTROPE
  return;
#elif EOS == EOS_TYPE_TABLE
  EOS_SC_init(eospath);
#else
  EOS_bad_eos_error();
#endif
}

double EOS_pressure_rho0_u(double rho, double u, const double *extra) {
  double press;
#if EOS == EOS_TYPE_GAMMA
  press = EOS_Gamma_pressure_rho0_u(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  press = EOS_Poly_pressure_rho0_u(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K, Gam;
#if COLD_FLOORS
    lrho = EOS_SC_get_min_lrho();
    lT   = EOS_SC_get_min_lT();
#endif
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    press = EOS_Poly_pressure_rho0_u(rho, u, K, Gam);
  } else {
    press = EOS_SC_pressure_rho0_u(lrho, lT, ye);
  }
#else
  press    = EOS_SC_pressure_rho0_u(lrho, lT, ye);
#endif // POLYTROPE_FALLBACK
#else
  EOS_bad_eos_error();
#endif
  return press;
}

/*
 * Here gamma is lorentz factor and geom is metric
 * added for enthalpy in particular because metric
 * factors may be needed to convert.
 * This one is special because it is used exclusively in Utoprim
 */
// WARNING: extra here is different than for other tabulated EOS funcs.
// Usually extra lrho, lT, ye.
// Here it contains NOTHING, lT, and
// the CONSERVED quantity u^0*Ye*rho
double EOS_pressure_rho0_w(double rho, double w, double gamma,
    const struct of_geom *geom, double *extra) {
  double press;
#if EOS == EOS_TYPE_GAMMA
  press = EOS_Gamma_pressure_rho0_w(rho, w);
#elif EOS == EOS_TYPE_POLYTROPE
  press = EOS_Poly_pressure_rho0_w(rho, w, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lapse  = geom->alpha;
  double detg   = geom->g;
  double ye_con = extra[EOS_YE];
  double yedens = lapse * ye_con / (gamma * detg);
  double ye     = fabs(yedens) / (fabs(rho) + SMALL);
  double lTold  = extra[EOS_LT];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    lTold       = EOS_SC_get_min_lT();
    double lrho = EOS_SC_get_min_lrho();
    double K, Gam;
    EOS_SC_get_polytrope(lrho, lTold, ye, &K, &Gam);
    press = EOS_Poly_pressure_rho0_w(rho, w, K, Gam);
  } else {
    press = EOS_SC_pressure_rho0_w(rho, w, ye, &lTold);
  }
#else
  press    = EOS_SC_pressure_rho0_w(rho, w, ye, &lTold);
#endif // POLYTROPE_FALLBACK
  extra[EOS_LT] = lTold;
#else
  EOS_bad_eos_error();
#endif
  return press;
}

double EOS_enthalpy_rho0_u(double rho, double u, const double *extra) {
  double enth;
#if EOS == EOS_TYPE_GAMMA
  enth = EOS_Gamma_enthalpy_rho0_u(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  enth  = EOS_Poly_enthalpy_rho0_u(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lrho   = extra[EOS_LRHO];
  double lT     = extra[EOS_LT];
  double ye     = extra[EOS_YE];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K, Gam;
    lrho = EOS_SC_get_min_lrho();
    lT   = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    enth = EOS_Poly_enthalpy_rho0_u(rho, MY_MAX(u, 0.0), K, Gam);
  } else {
    double h = EOS_SC_specific_enthalpy_rho0_u(lrho, lT, ye);
    enth     = h * rho;
  }
#else
  double h = EOS_SC_specific_enthalpy_rho0_u(lrho, lT, ye);
  enth     = h * rho;
#endif // POLYTROPE_FALLBACK
#else
  EOS_bad_eos_error();
#endif
  return enth;
}

double EOS_entropy_rho0_u(double rho, double u, const double *extra) {
  double ent;
#if EOS == EOS_TYPE_GAMMA
  ent = EOS_Gamma_entropy_rho0_u(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  ent   = EOS_Poly_entropy_rho0_u(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K, Gam;
    lrho = EOS_SC_get_min_lrho();
    lT   = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    ent = EOS_Poly_entropy_rho0_u(rho, u, K, Gam);
  } else {
    ent = EOS_SC_entropy(lrho, lT, ye);
  }
#else
  ent      = EOS_SC_entropy(lrho, lT, ye);
#endif // POLYTROPE_FALLBACK
#else
  EOS_bad_eos_error();
#endif
  return ent;
}

double EOS_sound_speed_rho0_u(double rho, double u, const double *extra) {
  double cs;
#if EOS == EOS_TYPE_GAMMA
  cs = EOS_Gamma_sound_speed_rho0_u(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  cs    = EOS_Poly_sound_speed_rho0_u(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    double K, Gam;
    lrho = EOS_SC_get_min_lrho();
    lT   = EOS_SC_get_min_lT();
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    cs = EOS_Poly_sound_speed_rho0_u(rho, u, K, Gam);
  } else {
    cs = EOS_SC_sound_speed(lrho, lT, ye);
  }
#else
  cs       = EOS_SC_sound_speed(lrho, lT, ye);
#endif // POLYTROPE_FALLBACK
#else
  EOS_bad_eos_error();
#endif
  return cs;
}

void EOS_set_floors(double scale, double rho, double u, double bsq,
    double *rhoflr, double *uflr, const double *extra) {
#if EOS == EOS_TYPE_GAMMA
  EOS_Gamma_set_floors(scale, rho, u, bsq, rhoflr, uflr);
#elif EOS == EOS_TYPE_POLYTROPE
  EOS_Poly_set_floors(scale, rho, u, bsq, rhoflr, uflr);
#elif EOS == EOS_TYPE_TABLE
  double ye = extra[EOS_YE];
  EOS_SC_set_floors(scale, rho, u, ye, bsq, rhoflr, uflr);
#else
  EOS_bad_eos_error();
#endif
}

double EOS_adiabatic_constant(double rho, double u, const double *extra) {
  double cad;
#if EOS == EOS_TYPE_GAMMA
  cad = EOS_Gamma_adiabatic_constant_rho0_u(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  EOS_Poly_adiabatic_constant(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double gam  = EOS_get_gamma(extra);
  cad         = u * pow(rho, -gam);
#else
  EOS_bad_eos_error();
#endif
  return cad;
}

double EOS_get_gamma(const double *extra) {
#if EOS == EOS_TYPE_GAMMA
  return gam;
#elif EOS == EOS_TYPE_POLYTROPE
  return poly_gam;
#elif EOS == EOS_TYPE_TABLE
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double ye   = extra[EOS_YE];
  double gam  = EOS_SC_gamma(lrho, lT, ye);
  return gam;
#else
  EOS_bad_eos_error();
#endif
}

double EOS_temperature(double rho, double u, const double *extra) {
#if EOS == EOS_TYPE_GAMMA
  return EOS_Gamma_temp(rho, u);
#elif EOS == EOS_TYPE_POLYTROPE
  return EOS_Poly_temperature(rho, u, poly_K, poly_gam);
#elif EOS == EOS_TYPE_TABLE
  double lT = extra[EOS_LT];
  return EOS_SC_temperature(lT);
#else
  EOS_bad_eos_error();
#endif
}

double EOS_u_press(double press, double rho, double *extra) {
  double u;
#if EOS == EOS_TYPE_GAMMA
  u = EOS_Gamma_u_press(press);
#elif EOS == EOS_TYPE_PLYTROPE
  u = EOS_Poly_u_press(press, rho, poly_K, poly_Gam);
#elif EOS == EOS_TYPE_TABLE
  double ye = extra[EOS_YE];
  // double yedens = extra[EOS_YE];
  // double ye     = fabs(yedens) / (fabs(rho) + SMALL);
  double lTold = extra[EOS_LT];
#if POLYTROPE_FALLBACK
  if (rho < rho_poly_thresh) {
    u     = EOS_SC_get_minu(rho, ye, 0.0);
    lTold = EOS_SC_get_min_lT();
  } else {
    u = EOS_SC_u_press(press, rho, ye, &lTold);
  }
#else
  u        = EOS_SC_u_press(press, rho, ye, &lTold);
#endif // POLYTROPE_FALLBACK
  extra[EOS_LT] = lTold;
#else
  EOS_bad_eos_error();
#endif // EOS
  return u;
}

#if RADIATION
double EOS_Theta_unit() {
#if EOS == EOS_TYPE_GAMMA
  return EOS_Gamma_Theta_unit();
#elif EOS == EOS_TYPE_TABLE
  return MEV / (MP * CL * CL); // TODO: is this right?
#else
  EOS_bad_eos_error();
#endif // EOS
}
#endif // RADIATION

/******************************************************************************/
