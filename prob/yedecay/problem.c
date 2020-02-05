/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NEUTRINO COOLING WITH FLAT EMISSIVITY               *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Problem-specific variables to set at runtime
static double rho0,uufactor;
static int boost;
#if NVAR_PASSIVE > 0
static double ye0;
#endif
void set_problem_params()
{
  set_param("rho0", &rho0);
  set_param("uufactor",&uufactor);
  set_param("boost", &boost);
  #if NVAR_PASSIVE > 0
  set_param("ye0", &ye0);
  #endif
}

void init_prob()
{
  double u0 = boost ? 0.5 : 0.0;

  ZLOOP {
    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU]  = uufactor*rho0;
    P[i][j][k][U1]  = u0;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
    #if EOS == EOS_TYPE_TABLE
    P[i][j][k][YE]     = ye0;
    P[i][j][k][YE_EM]  = ye0;
    #if METRIC == MKS
    P[i][j][k][ATM]    = NOT_ATM;
    #endif
    #endif
  }
}

