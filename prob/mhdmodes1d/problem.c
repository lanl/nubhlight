
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
static int idim;
void set_problem_params()
{
  set_param("nmode", &nmode);
  set_param("idim",&idim);
}

void init_prob()
{
  double X[NDIM];

  if ( idim != 1 && idim != 2 && idim != 3) {
    fprintf(stderr,"Primary axis must be 1, 2, or 3!\n");
    exit(1);
  }
  if ( nmode != 0 && nmode != 1 && nmode != 2 && nmode != 3 ) {
    fprintf(stderr,
      "ERROR. Not a valid mode. Modes are:\n"
      "\t0: ENTROPY\n"
      "\t1: SLOW\n"
      "\t2: ALFVEN\n"
      "\t3: FAST\n"
      "\tYou chose: %d.\n",nmode);
    exit(1);
  }

  // Mean state
  double rho0 = 1.;
  double u0 = 1.;
  double B10, B20,  B30;

  if (idim == 1) {
    B10 = 1.;
    B20 = 0.;
    B30 = 0.;
  } else if (idim == 2) {
    B10 = 0.;
    B20 = 1.;
    B30 = 0.;
  } else {
    B10 = 0.;
    B20 = 0.;
    B30 = 1.;
  }
  
  // Wavevector
  double k1 = 2.*M_PI;
  double amp = 1.e-4;

  complex omega, drho, du, du1, du2, du3, dB1, dB2, dB3, dutemp, dBtemp;

  // Eigenmode
  if (nmode == 0) { // Entropy
    omega = 2.*M_PI/5.*I; // To get tf
    drho  = 1.;
    du    = 0.;
    du1   = 0.;
    du2   = 0.;
    du3   = 0.;
    dB1   = 0.;
    dB2   = 0.;
    dB3   = 0.;
  } else if (nmode == 1) { // Slow
    omega = 2.74220688339*I;
    drho  = 0.580429492464;
    du    = 0.773905989952;
    du1   = -0.253320198552;
    du2   = 0.;
    du3   = 0.;
    dB1   = 0.;
    dB2   = 0.;
    dB3   = 0.;
  } else if (nmode == 2) { // Alfven
    omega = 3.44144232573*I;
    drho  = 0.;
    du    = 0.;
    du1   = 0.;
    du2   = 0.480384461415;
    du3   = 0.;
    dB1   = 0.;
    dB2   = 0.877058019307;
    dB3   = 0.;
  } else { // Fast
    omega = 3.44144232573*I;
    drho  = 0.;
    du    = 0.;
    du1   = 0.;
    du2   = 0.;
    du3   = 0.480384461415;
    dB1   = 0.;
    dB2   = 0.;
    dB3   = 0.877058019307;
  }

  // if idim =! 1, permute
  if ( idim == 2 ) {
    dutemp = du3;
    du3 = du2;
    du2 = du1;
    du1 = dutemp;
    dBtemp = dB3;
    dB3 = dB2;
    dB2 = dB1;
    dB1 = dBtemp;
  }
  if ( idim == 3 ) {
    dutemp = du1;
    du1 = du2;
    du2 = du3;
    du3 = dutemp;
    dBtemp = dB1;
    dB1 = dB2;
    dB2 = dB3;
    dB3 = dBtemp;
  }
    
  tf = 2.*M_PI/cimag(omega);
  DTd = tf/10.;

  ZLOOP {
    coord(i, j, k, CENT, X);
    double mode = amp*cos(k1*X[idim]);

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

