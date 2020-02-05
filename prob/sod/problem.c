/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOD SHOCKTUBE                                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double tscale;
static double pscale;
void set_problem_params()
{
  set_param("tscale", &tscale);
  set_param("pscale", &pscale);
}

void init_prob()
{
  double X[NDIM];
  double extra[EOS_NUM_EXTRA];

  // Make problem nonrelativistic
  // double tscale = 1.e-2;
  tf /= tscale;
  dt /= tscale;
  DTd /= tscale;
  DTl /= tscale;

  ZLOOP {
    coord(i, j, k, CENT, X);

    double rhogas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;

    // rescale pressure in relativistic case since
    // P must produce subluminal speeds
    double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
    pgas *= pscale;

    #if NVAR_PASSIVE > 0 || EOS == EOS_TYPE_TABLE
    double ye = (X[1] < 0.5 || X[1] > 1.5) ? 0.45 : 0.25;
    double yedens = ye*rhogas;
    #endif

    #if EOS == EOS_TYPE_TABLE
    {
      P[i][j][k][YE] = ye;
      PASSTYPE(YE) = PASSTYPE_NUMBER;
      P[i][j][k][YE_EM] = ye;
    
      extra[EOS_YE] = ye;

      #if NVAR_PASSIVE >= 4
      P[i][j][k][PASSIVE_START+3] = yedens;
      PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
      #endif
    }
    #endif

    P[i][j][k][RHO] = rhogas;
    P[i][j][k][UU] = EOS_u_press(pgas, rhogas, extra);
    P[i][j][k][U1] = 0.;
    P[i][j][k][U2] = 0.;
    P[i][j][k][U3] = 0.;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;

  } // ZLOOP

  // Rescale to make problem nonrelativistic
  ZLOOP {
    P[i][j][k][UU] *= tscale*tscale;
    P[i][j][k][U1] *= tscale;
    P[i][j][k][U2] *= tscale;
    P[i][j][k][U3] *= tscale;
    P[i][j][k][B1] *= tscale;
    P[i][j][k][B2] *= tscale;
    P[i][j][k][B3] *= tscale;
  }
}

