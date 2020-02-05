/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ADVECTION OF A PASSIVE SCALAR                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double cadv;
static int idim;
#if RADIATION && TRACERS
static int ntracers;
grid_int_type tracers_in_cell;
#endif
void set_problem_params()
{
  // advection speed = fluid 3-speed
  set_param("cadv", &cadv);
  set_param("idim", &idim);
  #if RADIATION && TRACERS
  set_param("ntracers", &ntracers);
  #endif
}

void init_prob()
{
  if (cadv >= 1) {
    fprintf(stderr,
      "ERROR: Advection speed must be less than speed of light.\n");
    exit(1);
  }
  if (idim != 1 && idim != 2 && idim !=3) {
    fprintf(stderr,
      "ERROR: idim must be 1, 2, or 3!. idim = %d. cadv = %f.\n",idim,cadv);
    exit(1);
  }

  #if NVAR_PASSIVE < 2
  fprintf(stderr,
    "ERROR: This problem must be run with at least two passive scalars!\n");
  exit(1);
  #endif // NVAR_PASSIVE < 2

  #if METRIC != MINKOWKI
  fprintf(stderr,
    "ERROR: This problem must be run with the Minkowski metric!\n");
  exit(1);
  #endif // METRIC != MINKOWSKI

  double X[NDIM];

  // lorentz factor
  double gamma = sqrt(1.0 / (1 - cadv*cadv));
  // four-velocity
  double qsqr = gamma*gamma - 1;
  // x-component
  double u1 = sqrt(qsqr);
  if ( u1 >= 1 ) {
    fprintf(stderr,
      "ERROR: Fluid speed must be less than speed of light!\n");
    exit(1);
  }

  // no fluid gradients
  // rest mass density
  double rho0 = 1.0;
  // energy. 1/10th rest energy
  double uu = 0.1;
  
  double Xmid = 0.5;
  double amp = 1.0;
  double kwave = 2.*M_PI;
  double kcenter = 0.5;

  double ux,uy,uz;
  if (idim == 1) {
    ux = u1; uy = 0.; uz = 0.;
  } else if (idim == 2) {
    ux = 0.; uy = u1; uz = 0.;
  } else if (idim == 3) {
    ux = 0.; uy = 0.; uz = u1;
  } else {
    fprintf(stderr,
      "ERROR: idim must be 1, 2, or 3!. idim = %d. cadv = %f.\n",idim,cadv);
    exit(1);
  }

  if (idim != 1 && idim != 2 && idim !=3) {
    fprintf(stderr,
      "ERROR: idim must be 1, 2, or 3!. idim = %d. cadv = %f.\n",idim,cadv);
    exit(1);
  }

  printf("Beginning passive setup.\n");

  ZLOOP {
    coord(i, j, k, CENT, X);

    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU] = uu;
    P[i][j][k][U1] = ux;
    P[i][j][k][U2] = uy;
    P[i][j][k][U3] = uz;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
    P[i][j][k][PASSIVE_START] = X[idim] <= Xmid ? 4*X[idim]*X[idim] : 0.0;
    P[i][j][k][PASSIVE_START + 1] = amp*sin(kwave*(X[idim] - kcenter));    

  } // ZLOOP

  // Ensure passives are intrinsic
  PASSLOOP { PASSTYPE(ipass) = PASSTYPE_INTRINSIC; }

  // Tracers
  #if RADIATION && TRACERS
  printf("Beginning tracer setup.\n");

  ZLOOP tracers_in_cell[i][j][k] = 1;

  sample_all_tracers(ntracers, 0, 0.0, tracers_in_cell, P);

  printf("Tracer count = %d\n",count_tracers_local());
  #endif // TRACERS

  printf("Problem setup complete.\n");
}
