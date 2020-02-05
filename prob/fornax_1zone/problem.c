/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NEUTRINO COOLING WITH FLAT EMISSIVITY               *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Problem-specific variables to set at runtime
static double rho0,T0;
#if NVAR_PASSIVE > 0
static double ye0;
#endif
void set_problem_params()
{
  set_param("rho0", &rho0);
  set_param("T0",&T0);
  #if NVAR_PASSIVE > 0
  set_param("ye0", &ye0);
  #endif
}

void init_prob()
{
  #if EOS != EOS_TYPE_TABLE
  fprintf(stderr, "Must use tabulated EOS!\n");
  exit(1);
  #endif

  #if NVAR_PASSIVE > 0
  PASSTYPE(YE) = PASSTYPE_NUMBER;
  #endif

  double u0 = EOS_SC_get_u_of_T(rho0*RHO_unit,
				T0,ye0);

  ZLOOP {
    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU]  = u0;
    P[i][j][k][U1]  = 0.;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
    #if EOS == EOS_TYPE_TABLE
    P[i][j][k][YE]  = ye0;
    P[i][j][k][YE_EM]  = ye0;
    #endif
  }
}

