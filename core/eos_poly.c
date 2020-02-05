/******************************************************************************
 *                                                                            *
 * EOS_POLY.C                                                                 *
 *                                                                            *
 * DEFINES THE EQUATION OF STATE AND RELATED MICROPHYSICS FOR POLYTROPE EOS   *
 *                                                                            *
 ******************************************************************************/

/* TODO: For now this is a pure polytrope. But it could be mixed
 *       with an ideal gas (say a photon gas) if this helps.
 */

#include "decs.h"

#if EOS == EOS_TYPE_POLYTROPE || POLYTROPE_FALLBACK

double EOS_Poly_pressure_rho0_u(double rho, double u, double K, double Gam) {
  rho = fabs(rho + SMALL);
  return K * pow(rho, Gam);
}

double EOS_Poly_pressure_rho0_w(double rho, double w, double K, double Gam) {
  rho = fabs(rho + SMALL);
  return K * pow(rho, Gam);
}

double EOS_Poly_enthalpy_rho0_u(double rho, double u, double K, double Gam) {
  double P = EOS_Poly_pressure_rho0_u(rho, u, K, Gam);
  return rho + u + P;
}

double EOS_Poly_entropy_rho0_u(double rho, double u, double K, double Gam) {
  // TODO: units?
  return K;
}

double EOS_Poly_sound_speed_rho0_u(double rho, double u, double K, double Gam) {
  rho            = MY_MAX(0.0, rho);
  u              = MY_MAX(0.0, u);
  double rhogam  = K * pow(rho, Gam - 1);
  double cs2_num = Gam * rhogam;
  double cs2_den = 1 + u + rhogam;
  double cs2     = cs2_num / cs2_den;
  return sqrt(cs2);
}

void EOS_Poly_set_floors(double scale, double rho, double u, double bsq,
    double *rhoflr, double *uflr) {
  *rhoflr = EOS_Poly_rho_floor(scale, bsq);
  *uflr   = EOS_Poly_u_floor(scale, bsq);
}

double EOS_Poly_rho_floor(double scale, double bsq) {
  double rhoflr = RHOMIN * scale;
  rhoflr        = MY_MAX(rhoflr, RHOMINLIMIT);
  rhoflr        = MY_MAX(rhoflr, bsq / BSQORHOMAX);

  return rhoflr;
}

double EOS_Poly_u_floor(double scale, double bsq) {
  double uflr = UUMIN * scale;
  uflr        = MY_MAX(uflr, UUMINLIMIT);
  uflr        = MY_MAX(uflr, bsq / BSQOUMAX);

  return uflr;
}

double EOS_Poly_adiabatic_constant(double rho, double u, double K, double Gam) {
  return K / (Gam - 1.);
}

double EOS_Poly_temperature(double rho, double u, double K, double Gam) {
  // TODO: is this right?
  return 0.0; // Polytrope is at zero internal energy/temperature
}

double EOS_Poly_u_press(double press, double rho, double K, double Gam) {
  return 0.0; // pressure and internal energy independent
}

double EOS_Poly_Theta_unit() { return MP / ME; }

#endif
