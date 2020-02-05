/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NEUTRINO COOLING WITH FLAT EMISSIVITY               *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define NONE    (0)
#define CENTER1 (1)
#define CENTER2 (2)

// Problem-specific variables to set at runtime
static double rho0,T0;
static double ye0;
static double dye;
static double c1x, c1y, r1, c2x, c2y, r2;

int center(double X[NDIM])
{
  if (X[1] >= c1x - r1 && X[1] <= c1x + r1
      && X[2] >= c1y - r1 && X[2] <= c1y + r1) {
    return CENTER1;
  }
  if (X[1] >= c2x - r2 && X[1] <= c2x + r2
      && X[2] >= c2y - r2 && X[2] <= c2y + r2) {
    return CENTER2;
  }
  // if (pow(X[1]-c1x,2.) + pow(X[2]-c1y,2.) <= r1*r1) return CENTER1;
  // if (pow(X[1]-c2x,2.) + pow(X[2]-c2y,2.) <= r2*r2) return CENTER2;
  return NONE;
}

void set_problem_params()
{
  set_param("rho0", &rho0);
  set_param("T0",&T0);
  set_param("ye0", &ye0);
  set_param("dye", &dye);
  set_param("c1x",&c1x);
  set_param("c1y",&c1y);
  set_param("r1",&r1);
  set_param("c2x",&c2x);
  set_param("c2y",&c2y);
  set_param("r2",&r2);
}

void init_prob()
{
  if (NVAR_PASSIVE <= 0) {
    fprintf(stderr,"Problem must be run with electron fraction!\n");
    exit(1);
  }
  
  #if NVAR_PASSIVE > 0
  PASSTYPE(YE) = PASSTYPE_NUMBER;
  #endif
  double X[NDIM];
  double u0;

  ZLOOP {
    coord(i,j,k,CENT,X);
    int loc = center(X);
    double ye = ye0;
    if (loc == CENTER1) ye -= dye;
    if (loc == CENTER2) ye += dye;
    u0 = EOS_SC_get_u_of_T(rho0*RHO_unit,T0,ye);

    
    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU]  = u0;
    P[i][j][k][U1]  = 0.;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
    #if EOS == EOS_TYPE_TABLE
    P[i][j][k][YE]  = ye;
    P[i][j][k][YE_EM]  = ye;
    #endif
  }
}
