/******************************************************************************
 *                                                                            *
 * EOS_GAMMA.C                                                                *
 *                                                                            *
 * DEFINES THE EQUATION OF STATE FOR IDEAL GAS                                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
double EOS_Gamma_pressure_rho0_u(double rho, double u) {
  return ((gam - 1.) * u);
}

double EOS_Gamma_pressure_rho0_w(double rho, double w) {
  return ((w - rho) * (gam - 1.) / gam);
}

double EOS_Gamma_Gaspress_entropy_rho0_u(double rho, double u) {  
  return (gam - 1.0) * EOS_Gamma_adiabatic_constant_rho0_u(rho, u);
}

double EOS_Gamma_Radpress_entropy_rho0_u(double rho, double u) {
  //double a = AR * pow(T_unit,2) * L_unit * pow(TEMP_unit, 4) * pow(M_unit, -2);
  double u_cgs = u * U_unit;
  double rho_cgs = rho * RHO_unit;
  double entropy_cgs = pow((64./3) * AR * (pow(u_cgs, 3)*pow(gam, 3)/pow(rho_cgs, 4)) * pow((gam -1)/gam, 3), 1./4);
  //return pow((64./3) * a * (pow(u, 3)*pow(gam, 3)/pow(rho, 4)) * pow((gam -1)/gam, 3), 1./4);
  return entropy_cgs / (KBOL / MP);
}

double EOS_Gamma_enthalpy_rho0_u(double rho, double u) { return rho + u * gam; }

double EOS_Gamma_adiabatic_constant_rho0_u(double rho, double u) {
  return u * pow(rho, -gam);
}

double EOS_Gamma_u_scale(double rho) { return pow(rho, gam); }

void EOS_Gamma_set_floors(double scale, double rho, double u, double bsq,
    double *rhoflr, double *uflr) {
  *rhoflr = EOS_Gamma_rho_floor(scale, bsq);
  *uflr   = EOS_Gamma_u_floor(scale, bsq);
  *rhoflr = MY_MAX(*rhoflr, u / UORHOMAX);
}

double EOS_Gamma_rho_floor(double scale, double bsq) {
  double rhoflr = RHOMIN * scale;
  rhoflr        = MY_MAX(rhoflr, RHOMINLIMIT);
  rhoflr        = MY_MAX(rhoflr, bsq / BSQORHOMAX);

  return rhoflr;
}

double EOS_Gamma_u_floor(double scale, double bsq) {
  double uscal = EOS_Gamma_u_scale(scale);
  double uflr  = UUMIN * uscal;
  uflr         = MY_MAX(uflr, UUMINLIMIT);
  uflr         = MY_MAX(uflr, bsq / BSQOUMAX);

  return uflr;
}

double EOS_Gamma_sound_speed_rho0_u(double rho, double u) {
  double ef    = EOS_Gamma_enthalpy_rho0_u(rho, u);
  double press = EOS_Gamma_pressure_rho0_u(rho, u);
  return sqrt(gam * press / ef);
}

double EOS_Gamma_u_press(double press) { return press / (gam - 1.); }

double EOS_Gamma_temp(double rho, double u) {
  #if EOS_GAMMA == RADPRESS
    double press_cgs = (gam - 1.) * u * U_unit;
    double temp_cgs = pow((3 * press_cgs)/AR, 1./4);
    double temp = temp_cgs * KBOL / MEV;
    //double temp = pow((3 * press), 1./4);
  #else
    // return TEMP_unit*(gam - 1.)*u/rho;
    double temp = (gam - 1.) * u / rho;
  #endif
  return temp;
}

#if RADIATION
double EOS_Gamma_Theta_unit() {
  return (gam - 1.) * MP / ME / (1. + tp_over_te);
}
#endif // RADIATION
#endif // EOS_TYPE_GAMMA
