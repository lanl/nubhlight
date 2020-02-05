/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ENTROPY WAVE                                        *
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

  double rho0 = 1., u0 = 0.01, u10 = 0.1;
  double delta_rho = 1.;//0.324033853851;
  double delta_u = 0.;//0.928101794093;
  double delta_u1 = 0.;//-0.183382445616;
  
  double kwave = 2.*M_PI;
  double amp = 0.01;

  // Rotate solution if possible
  double theta = 0., phi = 0.;
  if (N2 > 1 && N3 == 1) {
    theta = M_PI/2.;
    phi = M_PI/4.;
    kwave *= sqrt(2.);
  }
    //kwave *= 2.*cos(theta);
  else if (N3 > 1) {
    theta = atan(sqrt(2.));
    phi = M_PI/4.;
    kwave *= sqrt(3.);
    tf = 17.3210000e+00;
    DTd = tf/200.; 
  }

    //kwave *= 2.*cos(phi);
    /*if (theta <= M_PI/4.)
      kwave /= (1./(cos(theta)));
    else 
      kwave /= (1./(cos(M_PI/2. - theta)));*/
  double k1 = kwave*sin(theta)*cos(phi);
  double k2 = kwave*sin(theta)*sin(phi);
  double k3 = kwave*cos(theta);

  if (N2 == 1 && N3 == 1) {
    k1 = 2.*M_PI;
    k2 = 0.;
    k3 = 0.;
  }

  ZLOOP {
    coord(i, j, k, CENT, X);
    double mode = amp*cos(k1*X[1] + k2*X[2] + k3*X[3]);

    P[i][j][k][RHO] = rho0 + delta_rho*mode;
    P[i][j][k][UU] = u0 + delta_u*mode;
    P[i][j][k][U1] = (u10 + delta_u1*mode)*sin(theta)*cos(phi);
    P[i][j][k][U2] = (u10 + delta_u1*mode)*sin(theta)*sin(phi);
    P[i][j][k][U3] = (u10 + delta_u1*mode)*cos(theta);
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
    if (N2 == 1 && N3 == 1) {
      P[i][j][k][U1] = u10;
      P[i][j][k][U2] = 0.;
      P[i][j][k][U3] = 0.;
    }
  } // ZLOOP
}

