/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR EGGBEATER TORUS                                     *
 *   ZERO NET FLUX IN DISK, NET FIELD THREADING INNER REGION                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define CLASSIC_BFIELD (0)
#define NO_BFIELD (0)

// Local functions
void coord_transform(double *Pr, int i, int j);
double lfish_calc(double rmax) ;

static double BHflux;
void set_problem_params()
{
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
  double rho_av, rhomax, umax, beta, bsq_ij, bsq_max, q, rmax, Nloops; 

  // Fishbone-Moncrief parameters
  /*if (N3 == 1) {
    MAD = 1; // SANE initial field not supported in 2D
    rin = 6.;
    rmax = 12.;
    Nloops = 1;
    DTd *= DTf;
    DTf = 1;
  } else {
    rin = 20.;
    rmax = 41;
    Nloops = 4;
  }*/

  rin = 20.;
  rmax = 41.;
  Nloops = 4;

  /*#if CLASSIC_BFIELD
  rin = 6.;
  rmax = 12.;
  #endif*/

  printf("Nloops = %e\n", Nloops);
  l = lfish_calc(rmax);
  kappa = 1.e-3;

  // Plasma beta for initial magnetic field
  beta = 1.e2;

  rhomax = 0.;
  umax = 0.;
  ZSLOOP(-1, N1, -1, N2, -1, N3) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

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
      P[i][j][k][UU] = u * (1. + 4.e-2 * (get_rand() - 0.5));
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

  /*#if NO_BFIELD
  return;
  #endif*/

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


/*
#if CLASSIC_BFIELD
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
#else
    */
  double rstart, rend;
  //if (MAD == 0) {
    int rstart_set = 0, rend_set = 0;
    ZSLOOP(-1, N1, -1, N2, -1, N3) {
      P[i][j][k][RHO] /= rhomax;
      P[i][j][k][UU] /= rhomax;

      // set rstart and rend such that they bound the rho >= 0.2 region at the
      // midplane
      if (j == N2/2+NG && P[i][j][k][RHO] > 0.2 && rstart_set == 0) {
        coord(i, j, k, FACE1, X);
        bl_coord(X, &rstart, &th);
        rstart_set = 1;
        printf("rstart_set!\n");
      }

      coord(i,j,k,CENT,X);
      bl_coord(X, &r, &th);
      if (j == N2/2+NG && r/rstart*P[i][j][k][RHO] < 0.2 && rstart_set == 1 && 
          rend_set == 0) {
        coord(i+1, j, k, FACE1, X);
        bl_coord(X, &rend, &th);
        rend_set = 1;
        printf("rend_set!\n");
      }
    }

    if (mpi_nprocs() > 1) {
      if (mpi_io_proc())
        fprintf(stdout, "MPI nodes > 1: Using hard-coded loop normalization!\n");
      rstart = 2.335894e+01;
      rend = 3.737892e+02;
    } else {
      printf("rstart = %e rend = %e\n", rstart, rend);
    }
    umax /= rhomax;
    rhomax = 1.;
    fixup(P);
    bound_prim(P);

    // first find corner-centered vector potential
    ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
    ZSLOOP(0, N1, 0, N2, 0, 0) {
      rho_av = 0.25*(P[i][j  ][0][RHO] + P[i-1][j  ][0][RHO] +
                     P[i][j-1][0][RHO] + P[i-1][j-1][0][RHO])
        *(1. + 0.0*(get_rand() - 0.5));

      coord(i, j, k, CORN, X);
      bl_coord(X, &r, &th);

      q = r/rstart*rho_av/rhomax - 0.2;
      double lamb = (log(rend) - log(rstart))/(M_PI*Nloops);
      A[i][j] = 0.;
      if (r > rstart && r < rend && q > 0.) {
        if (N3 > 1) {
          A[i][j] = q*sin(M_PI/2. - th)*sin((log(r) - log(rstart))/lamb);
        } else { 
          A[i][j] = q;
        }
      }
    }

    // Differentiate to find cell-centered B, and begin normalization
    bsq_max = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT) ;

      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
          + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1]
               - A[i + 1][j] - A[i + 1][j + 1]) /
               (2. * dx[1] * geom->g);

      P[i][j][k][B3] = 0.;

      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    }
    bsq_max = mpi_max(bsq_max);
 
    // Normalize loops separately
    double ldrloop = (log(rend) - log(rstart))/Nloops;

    // Hardcoded normalization factors (for compatibility with MPI)
    // Does this depend on spin? Rout? Rin? Include some checks...
    double facs[4] = {1.626547e+00,
                   1.902312e+00,
                   4.663646e+00,
                   2.693245e+01};

    for (int n = 0; n < Nloops; n++) {
      double rmin = exp(log(rstart) + n*ldrloop);
      double rmax = exp(log(rstart) + (n+1)*ldrloop);
      double beta_min = 1.e100;
      double bsq_max = 0.;
      double umax = 0.;
      ZLOOP {
        coord(i, j, k, CENT, X);
        bl_coord(X, &r, &th);
        if (r > rmin && r < rmax) {
          // Inside loop
          geom = get_geometry(i, j, k, CENT);
          double bsq = bsq_calc(P[i][j][k], geom);
          if (bsq > bsq_max) bsq_max = bsq;
          if (P[i][j][k][UU] > umax) umax = P[i][j][k][UU];
          double beta_ij = (gam - 1.)*P[i][j][k][UU]/(0.5*bsq);
          if (beta_ij < beta_min) {
            beta_min = beta_ij;
          }
        }
      }

      // Rescale vector potential for this loop
      ZSLOOP(0,N1,0,N2,0,0) {
        coord(i, j, k, CORN, X);
        bl_coord(X, &r, &th);
        if (r > rmin && r < rmax) {
          if (mpi_nprocs() > 1) {
            A[i][j] *= facs[n];
          } else {
            A[i][j] *= sqrt(beta_min/beta);
          }
        }
      }
    }

    // Calculate B again
    ZLOOP {
      geom = get_geometry(i, j, k, CENT) ;

      // Flux-ct
      P[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
          + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1]
               - A[i + 1][j] - A[i + 1][j + 1]) /
               (2. * dx[1] * geom->g);
    }

  // Initialize a net magnetic field inside the initial torus
  ZSLOOP(0, N1, 0, N2, 0, 0) {
    A[i][j] = 0.;
    coord(i,j,k,CORN,X);
    bl_coord(X, &r, &th);
    double x = r*sin(th);
    double z = r*cos(th);
    double a_hyp = 20.;
    double b_hyp = 60.;
    double x_hyp = a_hyp*sqrt(1. + pow(z/b_hyp,2));
    
    q = (pow(x,2) - pow(x_hyp,2))/pow(x_hyp,2);
    if (x < x_hyp) {
      A[i][j] = 10.*q;
    }
  }

  // Evaluate net flux
  double Phi = 0.;
  ISLOOP(5, N1-1) {
    JSLOOP(0, N2-1) {
      int jglobal = j - NG + global_start[2];
      //int j = N2/2+NG;
      int k = NG;
      if (jglobal == N2TOT/2) {
        coord(i, j, k, CENT, X);
        bl_coord(X, &r, &th);
        if (r < rin) {
          double B2net =  (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1])/
                          (2.*dx[1]*ggeom[i][j][CENT].g);
          Phi += ggeom[i][j][CENT].g*2.*M_PI*dx[1]*fabs(B2net)/N3CPU;
        }
      }
    }
  }

  if (global_start[1] == 0) {
    JSLOOP(0, N2/2-1) {
      int i = 5 + NG;
      int k = NG;
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      double B1net = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*ggeom[i][j][CENT].g);
      Phi += ggeom[i][j][CENT].g*dx[2]*2.*M_PI*fabs(B1net)/N3CPU;
    }
  }
  Phi = mpi_reduce(Phi);

  double norm = BHflux/(Phi + SMALL);

  ZLOOP {
    geom = get_geometry(i, j, k, CENT) ;

    // Flux-ct
    P[i][j][k][B1] += -norm*(A[i][j] - A[i][j + 1]
        + A[i + 1][j] - A[i + 1][j + 1]) /
        (2. * dx[2] * geom->g);
    P[i][j][k][B2] += norm*(A[i][j] + A[i][j + 1]
             - A[i + 1][j] - A[i + 1][j + 1]) /
             (2. * dx[1] * geom->g);
  }
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

