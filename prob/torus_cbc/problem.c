/******************************************************************************
 *                                                                            *
 * problem.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Local functions
void   coord_transform(double *Pr, int i, int j);
double lfish_calc(double rmax);

static char   bfield_type[STRLEN];
static int    renormalize_densities;
static double rin;
static double rmax;
static double beta;
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
static double kappa_eos;
#endif
#if EOS == EOS_TYPE_TABLE
static double const_ye; // if > 0, set ye this way
#if !GAMMA_FALLBACK
static double entropy;
static double lrho_guess;
#endif
#endif
#if RADIATION && TRACERS
static int ntracers;
#endif

// Make global so it ends on heap
static double           A[N1 + 2 * NG][N2 + 2 * NG];
static grid_prim_type   PsaveLocal;
static grid_double_type lnh;
static grid_int_type    disk_cell;

void set_problem_params() {
  // supports: none, toroidal, classic
  set_param("bfield", &bfield_type);
  set_param("rin", &rin);
  set_param("rmax", &rmax);
  set_param("beta", &beta);
  set_param("renorm_dens", &renormalize_densities);
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  set_param("kappa_eos", &kappa_eos);
#endif
#if EOS == EOS_TYPE_TABLE
  set_param("const_ye", &const_ye);
#if !GAMMA_FALLBACK
  set_param("entropy", &entropy);
  set_param("lrho_guess", &lrho_guess);
#endif
#endif
#if RADIATION && TRACERS
  set_param("ntracers", &ntracers);
#endif
}

void init_prob() {
  double          r, th, sth, cth, ur, uh, up, u, rho, press, X[NDIM];
  struct of_geom *geom;

  // Disk interior
  double l, expm2chi, up1, DD, AA, SS, thin, sthin, cthin, DDin, AAin;
  double SSin, hm1;

// Diagnostics for entropy
#if EOS == EOS_TYPE_TABLE
  double ent, entmax;
#endif

  // Magnetic field
  double rho_av, rhomax, umax, pressmax, hm1max, bsq_ij, bsq_max, q;

  // total mass
  double mtot = 0.0;

  // scale angle
  double thdsqr = 0.0;

  // min and max
  double rhomin, umin;
  rhomin = RHOMINLIMIT;
  umin   = UUMINLIMIT;
#if EOS == EOS_TYPE_TABLE
  umin = EOS_SC_get_minu(rhomin, const_ye, umin);
#endif

#if NVAR_PASSIVE > 0
  PASSTYPE(PASSIVE_START + 0) = PASSTYPE_NUMBER;
#endif

#if NVAR_PASSIVE > 1
  PASSTYPE(PASSIVE_START + 1) = PASSTYPE_NUMBER;
#endif

#if EOS == EOS_TYPE_TABLE
  double ye, ye_atm;
#if !GAMMA_FALLBACK
  // guesses. Should be function local.
  double lrho0 = lrho_guess;
  // find constant entropy isocontour in table
  struct of_adiabat adiabat;
  int status = EOS_SC_find_adiabat_1d(entropy, const_ye, &adiabat);

  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[Torus]: "
        "Failed to find isocontour for entropy\n"
        "\ts        = %e\n"
        "\tye       = %e\n",
        entropy, const_ye);
    exit(1);
  }

  double hm1_min = EOS_SC_hm1_min_adiabat(&adiabat);
  // DEBUG
  if (mpi_io_proc()) {
    fprintf(stderr,
        "\nADIABAT FOUND!\n"
        "hm1_min     = %e\n"
        "hm1_min_cgs = %e\n"
        "lrho_min    = %e\n"
        "lrho_max    = %e\n"
        "\n",
        hm1_min, hm1_min * CL * CL, adiabat.lrho_min, adiabat.lrho_max);
    sleep(1);
    // EOS_SC_print_adiabat(&adiabat);
    // sleep(1);
    fprintf(stderr, "\n");
  }
// sleep(1);
// exit(1);
#endif
#endif // EOS_TYPE_TABLE

  l = lfish_calc(rmax);

  rhomax   = -INFINITY;
  umax     = -INFINITY;
  pressmax = -INFINITY;
  hm1max   = -INFINITY;
  ZSLOOP(-1, N1, -1, N2, -1, N3) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    sth = sin(th);
    cth = cos(th);

    // Calculate lnh
    DD = r * r - 2. * r + a * a;
    AA = ((r * r + a * a) * (r * r + a * a) - DD * a * a * sth * sth);
    SS = r * r + a * a * cth * cth;

    thin = M_PI / 2.;

    sthin = sin(thin);
    cthin = cos(thin);
    DDin  = rin * rin - 2. * rin + a * a;
    AAin  = ((rin * rin + a * a) * (rin * rin + a * a) -
            DDin * a * a * sthin * sthin);
    SSin  = rin * rin + a * a * cthin * cthin;

    if (r >= rin) {
      lnh[i][j][k] =
          0.5 * log((1. + sqrt(1. + 4. * (l * l * SS * SS) * DD /
                                        (AA * sth * AA * sth))) /
                    (SS * DD / AA)) -
          0.5 * sqrt(1. + 4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth)) -
          2. * a * r * l / AA -
          (0.5 * log((1. + sqrt(1. + 4. * (l * l * SSin * SSin) * DDin /
                                         (AAin * AAin * sthin * sthin))) /
                     (SSin * DDin / AAin)) -
              0.5 * sqrt(1. + 4. * (l * l * SSin * SSin) * DDin /
                                  (AAin * AAin * sthin * sthin)) -
              2. * a * rin * l / AAin);
    } else {
      lnh[i][j][k] = 1.;
    }

#if EOS == EOS_TYPE_TABLE
    ye     = const_ye; // may need to change this eventually
    ye_atm = 0.5;
#endif

    // regions outside torus
    if (lnh[i][j][k] < 0. || r < rin) {
      disk_cell[i][j][k] = 0;

      // Nominal values; real value set by fixup
      rho = rhomin;
      u   = umin;

      ur = 0.;
      uh = 0.;
      up = 0.;

      P[i][j][k][RHO] = rho;
      P[i][j][k][UU]  = u;

      PsaveLocal[i][j][k][RHO] = rho;
      PsaveLocal[i][j][k][UU]  = u;

      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;

#if EOS == EOS_TYPE_TABLE
      PsaveLocal[i][j][k][YE]    = ye_atm;
      PsaveLocal[i][j][k][YE_EM] = ye_atm;
      PsaveLocal[i][j][k][ATM]   = IS_ATM;
      P[i][j][k][YE]             = ye_atm;
      P[i][j][k][YE_EM]          = ye_atm;
      P[i][j][k][ATM]            = IS_ATM;
#endif
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      disk_cell[i][j][k] = 1;

      hm1 = exp(lnh[i][j][k]) - 1.;
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
      rho = pow(hm1 * (gam - 1.) / (kappa_eos * gam), 1. / (gam - 1.));
      u   = kappa_eos * pow(rho, gam) / (gam - 1.);
#elif EOS == EOS_TYPE_TABLE
      hm1 += hm1_min;
      if (lrho0 < EOS_SC_get_min_lrho())
        lrho0 = lrho_guess;
      // if (lrho0 > log10(1e2*RHO_unit)) lrho0 = lrho_guess;
      EOS_SC_isoentropy_hm1(hm1, &adiabat, &lrho0, &rho, &u);
#else
      fprintf(stderr, "[Torus]: Bad EOS chosen.\n");
      exit(1);
#endif

      rho = MY_MAX(rho, rhomin);
      u   = MY_MAX(u, umin);

      ur = 0.;
      uh = 0.;

      // Calculate u^phi
      expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      up1      = sqrt((-1. + sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      up       = (2. * a * r * sqrt(1. + up1 * up1) / sqrt(AA * SS * DD) +
            sqrt(SS / AA) * up1 / sth);

      PsaveLocal[i][j][k][RHO] = rho;
      PsaveLocal[i][j][k][UU]  = u;

      P[i][j][k][RHO] = rho;
      P[i][j][k][UU]  = u;

      P[i][j][k][U1] = ur;
      P[i][j][k][U2] = uh;
      P[i][j][k][U3] = up;

#if EOS == EOS_TYPE_TABLE
      PsaveLocal[i][j][k][YE]    = ye;
      PsaveLocal[i][j][k][YE_EM] = ye;
      PsaveLocal[i][j][k][ATM]   = NOT_ATM;
      P[i][j][k][YE]             = ye;
      P[i][j][k][YE_EM]          = ye;
      P[i][j][k][ATM]            = NOT_ATM;
#endif

      if (rho > rhomax)
        rhomax = rho;
      if (u > umax && r > rin && lnh[i][j][k] >= 0)
        umax = u;
      if (hm1 > hm1max)
        hm1max = hm1;
      P[i][j][k][UU] *= (1. + 4.e-2 * (get_rand() - 0.5));

      // Convert from 4-velocity to 3-velocity
      coord_transform(P[i][j][k], i, j);
    }

    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
  } // ZSLOOP

#if EOS == EOS_TYPE_TABLE && !GAMMA_FALLBACK
  EOS_SC_adiabat_free(&adiabat);
#endif

  // get rhomax, umax globally

  // DEBUG
  /*
  printf("[%d] max: (rho, u, hm1) = (%g, %g, %g)\n",
   mpi_myrank(), rhomax, umax, hm1max);
  */

  umax   = mpi_max(umax);
  rhomax = mpi_max(rhomax);
  hm1max = mpi_max(hm1max);

  if (mpi_io_proc()) {
    fprintf(stdout, "rhomax   = %f\n", rhomax);
    fprintf(stdout, "umax     = %f\n", umax);
    fprintf(stdout, "hm1max   = %f\n", hm1max);
  }
  // Normalize densities
  ZSLOOP(-1, N1, -1, N2, -1, N3) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    if (renormalize_densities) {
      P[i][j][k][RHO] /= rhomax;
      P[i][j][k][UU] /= rhomax;
      PsaveLocal[i][j][k][RHO] /= rhomax;
      PsaveLocal[i][j][k][UU] /= rhomax;
    }

    if (r > rin && lnh[i][j][k] >= 0.) {
#if EOS == EOS_TYPE_TABLE
      EOS_SC_fill(PsaveLocal[i][j][k], extra[i][j][k]);
#endif
      press = EOS_pressure_rho0_u(
          PsaveLocal[i][j][k][RHO], PsaveLocal[i][j][k][UU], extra[i][j][k]);
      if (press > pressmax)
        pressmax = press;
    }
  }

  // print diagnostics
  pressmax = mpi_max(pressmax);
  if (mpi_io_proc()) {
    fprintf(stdout, "After renormalization:\n");
    fprintf(stdout, "\tpressmax = %e\n", pressmax);
  }
  if (renormalize_densities) {
    umax /= rhomax;
    rhomax = 1.;
  }
  if (mpi_io_proc()) {
    fprintf(stdout, "Beginning fixup.\n"); // debug
  }
  fixup(P, extra);

// initialize hot atmosphere
// Broken.
#if EOS == EOS_TYPE_TABLE //&& !COLD_FLOORS
  ZLOOP {
    if (P[i][j][k][ATM] < ATM_THRESH) {
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      P[i][j][k][UU] = P[i][j][k][RHO] / r;
    }
  }
  fixup(P, extra);
#endif

  if (mpi_io_proc()) {
    fprintf(stdout, "Bounding prim.\n"); // debug
  }
  bound_prim(P);
  if (mpi_io_proc()) {
    fprintf(stdout, "Fixup finished.\n"); // debug
  }
  ZLOOP {
    if (disk_cell[i][j][k]) {
      double rho_integrand =
          P[i][j][k][RHO] * ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3];
      mtot += rho_integrand;

      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      thdsqr += rho_integrand * (th - M_PI / 2) * (th - M_PI / 2);
    }
  }
  mtot       = mpi_reduce(mtot);
  thdsqr     = mpi_reduce(thdsqr);
  double thd = sqrtf(thdsqr / mtot);
  if (mpi_io_proc()) {
    printf("TOTAL MASS:\n"
           "\tcode units = %g\n"
           "\tcgs        = %g\n"
           "\tMsun       = %g\n",
        mtot, mtot * M_unit, mtot * M_unit / MSUN);
    printf("Opening angle of initial torus:\n"
           "\tthd        = %g\n",
        thd);
  }
// debug
#if EOS == EOS_TYPE_TABLE
  if (mpi_io_proc()) {
    fprintf(stdout, "Calculating max entropy:\n");
  }
  entmax = -1.0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
    if (r > rin && lnh[i][j][k] >= 0.) {
      EOS_SC_fill(P[i][j][k], extra[i][j][k]);
      ent = EOS_entropy_rho0_u(P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
      if (ent > entmax)
        entmax = ent;
    }
  }
  entmax = mpi_max(entmax);
  if (mpi_io_proc()) {
    fprintf(stdout, "Maximum entropy in disk = %e\n", entmax);
  }
#endif

#if RADIATION && TRACERS
  {
    if (mpi_io_proc()) {
      fprintf(stdout, "Setting up tracers\n");
    }
    sample_all_tracers(ntracers, 0, 0.0, disk_cell, P);
    int    ntracers    = mpi_reduce_int(count_tracers_local());
    double tracer_mass = get_total_tracer_mass();
    if (mpi_io_proc()) {
      fprintf(stdout,
          "TRACERS\n"
          "\tTracer count      = %d\n"
          "\tTracer total mass = %g\n",
          ntracers, tracer_mass);
    }
  }
#endif // TRACERS

  // If we're just testing initial data with no B-fields,
  // we're done.
  if (strcmp(bfield_type, "none") == 0) {
    printf("No magnetic field.\n");
    ZLOOP {
      P[i][j][k][B1] = 0.;
      P[i][j][k][B2] = 0.;
      P[i][j][k][B3] = 0.;
    }
    return;
  }

  // Do magnetic field
  if (strcmp(bfield_type, "classic") == 0) {
    // Classic harm disk.
    // Find vector potential at corners
    ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
    ZSLOOP(0, N1, 0, N2, 0, 0) {
      rho_av = 0.25 * (P[i][j][0][RHO] + P[i - 1][j][0][RHO] +
                          P[i][j - 1][0][RHO] + P[i - 1][j - 1][0][RHO]);

      coord(i, j, k, CORN, X);
      bl_coord(X, &r, &th);

      q = rho_av / rhomax - 0.2;
      if (q > 0.)
        A[i][j] = q;
    }

    // Differentiate to find cell-centered B, and begin normalization
    bsq_max = 0.;
    ZLOOP {
      geom = get_geometry(i, j, k, CENT);

      // Flux-ct
      P[i][j][k][B1] =
          -(A[i][j] - A[i][j + 1] + A[i + 1][j] - A[i + 1][j + 1]) /
          (2. * dx[2] * geom->g);
      P[i][j][k][B2] = (A[i][j] + A[i][j + 1] - A[i + 1][j] - A[i + 1][j + 1]) /
                       (2. * dx[1] * geom->g);
      P[i][j][k][B3] = 0.;

      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max)
        bsq_max = bsq_ij;
    }
    bsq_max = mpi_max(bsq_max);
  } else if (strcmp(bfield_type, "toroidal") == 0) {
    // Magnetic field scales with rho/rhomax.
    // Loop in z is safe b/c initial data is axisymmetric
    bsq_max = 0.;
    ZLOOP {
      P[i][j][k][B1] = P[i][j][k][B2] = 0;

      q = P[i][j][k][RHO] / rhomax - 0.2;

      // P[i][j][k][B3] = (q > 0.) ? q : 0.0;
      P[i][j][k][B3] = (q > 0.) ? 1 : 0.0;

      geom   = get_geometry(i, j, k, CENT);
      bsq_ij = bsq_calc(P[i][j][k], geom);
      if (bsq_ij > bsq_max)
        bsq_max = bsq_ij;
    }
  } else {
    fprintf(stderr, "Unknown bfield configuration: %s\n", bfield_type);
    exit(1);
  }

  // Normalize to set field strength
  double beta_act = pressmax / (0.5 * bsq_max);
  double norm     = sqrt(beta_act / beta);
  ZLOOP {
    P[i][j][k][B1] *= norm;
    P[i][j][k][B2] *= norm;
    P[i][j][k][B3] *= norm;
  }

  // Check normalization
  bsq_max = 0.;
  umax    = 0.;
  ZLOOP {
    geom   = get_geometry(i, j, k, CENT);
    bsq_ij = bsq_calc(P[i][j][k], geom);

    if (bsq_ij > bsq_max)
      bsq_max = bsq_ij;

    if (P[i][j][k][UU] > umax)
      umax = P[i][j][k][UU];
  }
  bsq_max = mpi_max(bsq_max);
  if (mpi_io_proc()) {
    printf("MAX bsq = %e Pmax = %e\n", bsq_max, pressmax);
    printf("FINAL BETA: %e\n", pressmax / (0.5 * bsq_max));
  }
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(double *Pr, int ii, int jj) {
  double          X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double          AA, BB, CC, discr;
  double          alpha, gamma, beta[NDIM];
  struct of_geom *geom, blgeom;

  coord(ii, jj, 0, CENT, X);
  bl_coord(X, &r, &th);
  blgset(ii, jj, &blgeom);

  ucon[1] = Pr[U1];
  ucon[2] = Pr[U2];
  ucon[3] = Pr[U3];

  AA = blgeom.gcov[0][0];
  BB = 2. * (blgeom.gcov[0][1] * ucon[1] + blgeom.gcov[0][2] * ucon[2] +
                blgeom.gcov[0][3] * ucon[3]);
  CC = 1. + blgeom.gcov[1][1] * ucon[1] * ucon[1] +
       blgeom.gcov[2][2] * ucon[2] * ucon[2] +
       blgeom.gcov[3][3] * ucon[3] * ucon[3] +
       2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
                blgeom.gcov[1][3] * ucon[1] * ucon[3] +
                blgeom.gcov[2][3] * ucon[2] * ucon[3]);

  discr   = BB * BB - 4. * AA * CC;
  ucon[0] = (-BB - sqrt(discr)) / (2. * AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16 * sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  for (int mu = 0; mu < NDIM; mu++)
    tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      tmp[mu] += trans[mu][nu] * ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++)
    ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  set_dxdX(X, trans);
  for (int mu = 0; mu < NDIM; mu++)
    tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      tmp[mu] += trans[mu][nu] * ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++)
    ucon[mu] = tmp[mu];

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  geom  = get_geometry(ii, jj, 0, CENT);
  alpha = 1.0 / sqrt(-geom->gcon[0][0]);
  gamma = ucon[0] * alpha;

  beta[1] = alpha * alpha * geom->gcon[0][1];
  beta[2] = alpha * alpha * geom->gcon[0][2];
  beta[3] = alpha * alpha * geom->gcon[0][3];

  Pr[U1] = ucon[1] + beta[1] * gamma / alpha;
  Pr[U2] = ucon[2] + beta[2] * gamma / alpha;
  Pr[U3] = ucon[3] + beta[3] * gamma / alpha;
}

double lfish_calc(double r) {
  return (
      ((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
          ((-2. * a * r * (pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2))) /
                  sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
              ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) * (2. + r))) /
                  sqrt(1 + (2. * a) / pow(r, 1.5) - 3. / r))) /
      (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
          (pow(a, 2) + (-2. + r) * r)));
}
