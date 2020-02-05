/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR OPTICALLY THIN BREMSSTRALUNG COOLING                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Problem-specific variables to set at runtime
static double T0, ne;
void set_problem_params()
{
  set_param("T0", &T0);
  set_param("ne", &ne);
}

// Initialize dynamical variables 
void init_prob()
{
  double u0 = ne*KBOL*T0/(gam-1.);
  
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

