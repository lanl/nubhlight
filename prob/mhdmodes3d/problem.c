/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ENTROPY WAVE                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <complex.h>

static int nmode;
void set_problem_params()
{
  set_param("nmode", &nmode);
}

void init_prob()
{
	double X[NDIM];

  // Mean state
  double rho0 = 1.;
  double u0 = 1.;
  double B10 = 1.;
  double B20 = 0.;
  double B30 = 0.;

  // Wavevector
  double k1 = 2.*M_PI;
  double k2 = 2.*M_PI;
  double k3 = 2.*M_PI;
  double amp = 1.e-4;

  complex omega, drho, du, du1, du2, du3, dB1, dB2, dB3;

  // Eigenmode
  if (nmode == 0) { // Entropy
    omega = 2.*M_PI/5.*I;
    drho = 1.;
    du   = 0.;
    du1  = 0.;
    du2  = 0.;
    du3  = 0.;
    dB1  = 0.;
    dB2  = 0.;
    dB3  = 0.;
  } else if (nmode == 1) { // Slow
    omega = 2.35896379113*I;
    drho = 0.556500332363;
    du   = 0.742000443151;
    du1  = -0.282334999306;
    du2  = 0.0367010491491;
    du3  = 0.0367010491491;
    dB1  = -0.195509141461;
    dB2  = 0.0977545707307;
    dB3  = 0.0977545707307;
  } else if (nmode == 2) { // Alfven
    omega = 3.44144232573*I;
    drho = 0.;
    du   = 0.;
    du1  = 0.;
    du2  = -0.339683110243;
    du3  = 0.339683110243;
    dB1  = 0.;
    dB2  = 0.620173672946;
    dB3  = -0.620173672946;
  } else { // Fast
    omega = 6.92915162882*I;
    drho = 0.481846076323;
    du   = 0.642461435098;
    du1  = -0.0832240462505;
    du2  = -0.224080007379;
    du3  = -0.224080007379;
    dB1  = 0.406380545676;
    dB2  = -0.203190272838;
    dB3  = -0.203190272838;
  }
    
    /*omega = 3.87806616532*I;
    drho = 0.580429492464;
    du   = 0.773905989952;
    du1  = -0.179124430208;
    du2  = -0.179124430208;
    du3  = 0.;
    dB1  = 0.;
    dB2  = 0.;
    dB3  = 0.;*/

  tf = 2.*M_PI/cimag(omega);
  DTd = tf/10.;

  ZLOOP {
    coord(i, j, k, CENT, X);
    double mode = amp*cos(k1*X[1] + k2*X[2] + k3*X[3]);

    P[i][j][k][RHO] = rho0 + creal(drho*mode);
    P[i][j][k][UU] = u0 + creal(du*mode);
    P[i][j][k][U1] = creal(du1*mode);
    P[i][j][k][U2] = creal(du2*mode);
    P[i][j][k][U3] = creal(du3*mode);
    P[i][j][k][B1] = B10 + creal(dB1*mode);
    P[i][j][k][B2] = B20 + creal(dB2*mode);
    P[i][j][k][B3] = B30 + creal(dB3*mode);

    // Allow for passive scalars if needed
    #if EOS == EOS_TYPE_TABLE
    double ye = 0.25;
    P[i][j][k][YE] = ye;
    P[i][j][k][YE_EM] = ye;
    #endif
  } // ZLOOP
}

