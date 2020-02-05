/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR OPTICALLY THIN BREMSSTRALUNG COOLING                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Problem-specific variables to set at runtime
static double T0;
void set_problem_params()
{
  set_param("T0", &T0);
}

void init_prob()
{
  double N  = 5.4e-39; // cm^3 K^1/2 s^-1 Sr^-1 Hz^-1
  double n0 = HPL*pow(T0, 1./2.)/((gam-1.)*M_PI*N*tf*T_unit);
  double u0 = n0*KBOL*T0/(gam-1.);
  
  ZLOOP {
    P[i][j][k][RHO] = 1.e0;
    P[i][j][k][UU]  = 2.*u0/U_unit;
    P[i][j][k][U1]  = 0.;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
  }
}

