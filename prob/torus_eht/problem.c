/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Local functions
void coord_transform(double *Pr, int i, int j);
double lfish_calc(double rmax) ;

static int MAD;
static double BHflux;
void set_problem_params()
{
  set_param("MAD", &MAD);
  set_param("BHflux", &BHflux);
}

void init_prob()
{
  double r, th, sth, cth, ur, uh, up, u, rho, X[NDIM];
  struct of_geom *geom;

  // Disk interior
  double l, rin, lnh, expm2chi, up1, DD, AA, SS, thin, sthin, cthin, DDin, AAin;
  double SSin, kappa, hm1;

  // Magnetic field
  static double A[N1+2*NG][N2+2*NG];
  double rho_av, rhomax, umax, beta, bsq_ij, bsq_max, q, rmax;

  rin = 6.;
  rmax = 12.;

  l = lfish_calc(rmax);
  kappa = 1.e-3;

  // Plasma beta for initial magnetic field
  beta = 1.e2;

  rhomax = 0.;
  umax = 0.;
  int passed = 0;
  ZSLOOP(-1, N1, -1, N2, -1, N3) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
    if (j == NG && k == NG && passed == 0 && r > 50.) {
      printf("r = %e, i = %i\n", r, i-NG);
      passed = 1;
    }

    sth = sin(th);
    cth = cos(th);

    // Calculate lnh
    DD = r * r - 2. * r + a * a;
    AA = (r * r + a * a) * (r * r + a * a) -
        DD * a * a * sth * sth;
    SS = r * r + a * a * cth * cth;

    thin = M_PI / 2.;

    sthin = sin(thin);
    cthin = cos(thin);
    DDin = rin * rin - 2. * rin + a * a;
    AAin = (rin * rin + a * a) * (rin * rin + a * a)
        - DDin * a * a * sthin * sthin;
    SSin = rin * rin + a * a * cthin * cthin;

    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA *
                     sth *
                     AA *
                     sth)))
        / (SS * DD / AA))
          - 0.5 * sqrt(1. +
           4. * (l * l * SS * SS) * DD /
           (AA * AA * sth * sth))
          - 2. * a * r * l / AA -
          (0.5 *
           log((1. +
          sqrt(1. +
               4. * (l * l * SSin * SSin) * DDin /
               (AAin * AAin * sthin * sthin))) /
         (SSin * DDin / AAin))
           - 0.5 * sqrt(1. +
            4. * (l * l * SSin * SSin) *
            DDin / (AAin * AAin * sthin *
              sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // regions outside torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by fixup
      rho = 1.e-7 * RHOMIN;
      u = 1.e-7 * UUMIN;

      ur = 0.;
      uh = 0.;
      up = 0.;

      P[i][j][k][RHO] = rho;
      P[i][j][k][UU] = u;
      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      hm1 = exp(lnh) - 1.;
      rho = pow(hm1 * (gam - 1.) / (kappa * gam),
          1. / (gam - 1.));
      u = kappa * pow(rho, gam) / (gam - 1.);
      ur = 0.;
      uh = 0.;

      // Calculate u^phi
      expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;


      P[i][j][k][RHO] = rho;
      if (rho > rhomax) rhomax = rho;
      P[i][j][k][UU] = u* (1. + 4.e-2 * (get_rand() - 0.5));
      //P[i][j][k][UU] = 1.05*u;
      if (u > umax && r > rin) umax = u;
      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(P[i][j][k], i, j);
    }

    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
  } // ZSLOOP

  // Normalize the densities so that max(rho) = 1
  umax = mpi_max(umax);
  rhomax = mpi_max(rhomax);
  /*ZSLOOP(-1, N1, -1, N2, -1, N3) {
    P[i][j][k][RHO] /= rhomax;
    P[i][j][k][UU] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;
  fixup(P);
  bound_prim(P);*/


    // Normalize densities
    ZSLOOP(-1, N1, -1, N2, -1, N3) {
      P[i][j][k][RHO] /= rhomax;
      P[i][j][k][UU] /= rhomax;
    }
    umax /= rhomax;
    rhomax = 1.;
    fixup(P);
    bound_prim(P);
    //return;

    // Find vector potential at corners
    ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
    ZSLOOP(0, N1, 0, N2, 0, 0) {
      rho_av = 0.25*(P[i][j  ][0][RHO] + P[i-1][j  ][0][RHO] +
                     P[i][j-1][0][RHO] + P[i-1][j-1][0][RHO]);

      coord(i, j, k, CORN, X);
      bl_coord(X, &r, &th);

      q = rho_av/rhomax - 0.2;
      if (q > 0.) A[i][j] = q;
    }

    // Differentiate to find cell-centered B, and begin normalization
    bsq_max = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT) ;

      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1])/
                        (2.*dx[2]*geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/
                       (2.*dx[1]*geom->g);
      P[i][j][k][B3] = 0.;

      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    }
    bsq_max = mpi_max(bsq_max);

    // Normalize to set field strength
    double beta_act = (gam-1.)*umax/(0.5*bsq_max);
    double norm = sqrt(beta_act/beta);
    ZLOOP {
      P[i][j][k][B1] *= norm;
      P[i][j][k][B2] *= norm;
    }

    // Check normalization
    bsq_max = 0.;
    umax = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT);
      bsq_ij = bsq_calc(P[i][j][k], geom);
      
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;
      
      if (P[i][j][k][UU] > umax) umax = P[i][j][k][UU];
    }
    printf("MAX bsq = %e Pmax = %e\n", bsq_max, (gam-1.)*umax);
    printf("FINAL BETA: %e\n", (gam-1.)*umax/(0.5*bsq_max));
    
    
    /*ZLOOP {
      P[i][j][k][B1] = 0.;
      P[i][j][k][B2] = 0.;
      P[i][j][k][B3] = 0.;
    }*/
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(double *Pr, int ii, int jj)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct of_geom *geom, blgeom;

  coord(ii, jj, 0, CENT, X);
  bl_coord(X, &r, &th);
  blgset(ii, jj, &blgeom);

  ucon[1] = Pr[U1];
  ucon[2] = Pr[U2];
  ucon[3] = Pr[U3];

  AA = blgeom.gcov[0][0];
  BB = 2. * (blgeom.gcov[0][1] * ucon[1] +
       blgeom.gcov[0][2] * ucon[2] +
       blgeom.gcov[0][3] * ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1] * ucon[1] * ucon[1] +
      blgeom.gcov[2][2] * ucon[2] * ucon[2] +
      blgeom.gcov[3][3] * ucon[3] * ucon[3] +
      2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
      blgeom.gcov[1][3] * ucon[1] * ucon[3] +
      blgeom.gcov[2][3] * ucon[2] * ucon[3]);

  discr = BB * BB - 4. * AA * CC;
  ucon[0] = (-BB - sqrt(discr)) / (2. * AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  set_dxdX(X, trans);
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  //ucon[1] *= (1. / (r - R0));
  //ucon[1] /= dr_dx(X[1]);
  //ucon[2] *=
  //    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
  //ucon[2] /= dth_dx(X[2]);

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  geom = get_geometry(ii, jj, 0, CENT) ;
  alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  gamma = ucon[0] * alpha;

  beta[1] = alpha * alpha * geom->gcon[0][1];
  beta[2] = alpha * alpha * geom->gcon[0][2];
  beta[3] = alpha * alpha * geom->gcon[0][3];

  Pr[U1] = ucon[1] + beta[1] * gamma / alpha;
  Pr[U2] = ucon[2] + beta[2] * gamma / alpha;
  Pr[U3] = ucon[3] + beta[3] * gamma / alpha;
}

double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}

