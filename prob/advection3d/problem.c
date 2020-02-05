/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ADVECTION OF A PASSIVE SCALAR in 3D                 *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void set_problem_params()
{
  // stub
}

void init_prob()
{
  #if NVAR_PASSIVE < 1
  fprintf(stderr,
    "ERROR: This problem must be run with at least one passive scalar!\n");
  exit(1);
  #endif // NVAR_PASSIVE < 1

  #if METRIC != MINKOWKI
  fprintf(stderr,
    "ERROR: This problem must be run with the Minkowski metric!\n");
  exit(1);
  #endif // METRIC != MINKOWSKI

  double X[NDIM];
  double nspace = 3.0;

  // diagonal distance on the cube
  double ddiag = sqrt(nspace);
  // advection speed
  double cadv = 0.5*ddiag;
  double csqr = cadv*cadv;
  if (csqr >= 1) {
    fprintf(stderr,
      "ERROR: Advection speed must be less than speed of light!\n");
    exit(1);
  }

  // lorentz factor
  double gamma = sqrt(1.0 / (1 - csqr));
  // four-speed
  double qsqr = gamma*gamma - 1.0;
  // speed in x/y/z
  double ux = sqrt(qsqr/nspace);
  double uy = ux;
  double uz = ux;
  double threevsqr = ux*ux + uy*uy + uz*uz;
  double threev = sqrt(threevsqr);
  double u0sqr = -(threevsqr + 1.);
  if ( (threevsqr + u0sqr + 1) > SMALL ) {
    fprintf(stderr,
      "ERROR: Fluid speed must be less than speed of light!\n");
    fprintf(stderr,"\tSpeeds are: [%f, %f, %f].\n",ux,uy,uz);
    fprintf(stderr,"\tNorm speed is: %f.\n",threev);
    fprintf(stderr,"\tcadv = %f.\n\tGamma = %f.\n\tq^2 = %f.\n\t(u^0)^2 = %f.\n",
      cadv,gamma,qsqr,u0sqr);
    fprintf(stderr,"\tu^mu u_mu = %f.\n",u0sqr + threevsqr);
    exit(1);
  }

  // no fluid gradients
  // rest mass density
  double rho0 = 1.0;
  // energy. 1/10th rest energy
  double uu = 0.1;

  // wave data
  double amp = 1.0;
  double kwave = 2*M_PI;
  double wave;

  ZLOOP {
    coord(i, j, k, CENT, X);
    wave = amp*cos(kwave*X[1])*cos(kwave*X[2])*cos(kwave*X[3]);

    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU] = uu;
    P[i][j][k][U1] = ux;
    P[i][j][k][U2] = uy;
    P[i][j][k][U3] = uz;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
    P[i][j][k][PASSIVE_START] = wave;
  } // ZLOOP

  PASSLOOP { PASSTYPE(ipass) = PASSTYPE_INTRINSIC; }
}
