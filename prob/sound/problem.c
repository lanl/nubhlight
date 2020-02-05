/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOUND WAVE                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

//static double tscale;
void set_problem_params()
{
  //set_param("tscale", &tscale);
}

void init_prob()
{
	double X[NDIM];

  // Make problem nonrelativistic
  //tf /= tscale;
  //dt /= tscale;
  //DTd /= tscale;
  //DTl /= tscale;

  // omega^2 - cs^2 k^2 = 0, cs^2 = gam P / rho

  double rho0 = 1., u0 = 0.01, u10 = 0.;
  double delta_rho = 0.994443291303;
  double delta_u = 0.016574054855;
  double delta_u1 = -0.103960767065;
  //double delta_rho = 0.324033853851;
  //double delta_u = 0.928101794093;
  //double delta_u1 = -0.183382445616;
  
  double kwave = 2.*M_PI;
  double amp = 0.01;

  ZLOOP {
    coord(i, j, k, CENT, X);
    double mode = amp*cos(kwave*X[1]);

    P[i][j][k][RHO] = rho0 + delta_rho*mode;
    P[i][j][k][UU] = u0 + delta_u*mode;
    P[i][j][k][U1] = u10 + delta_u1*mode;
    P[i][j][k][U2] = 0.;
    P[i][j][k][U3] = 0.;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
  } // ZLOOP
}

