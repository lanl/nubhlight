/******************************************************************************
 *                                                                            *
 * problem.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NUMERICAL DISK                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Local functions
void   coord_transform(double *Pr, int i, int j);

static char   bfield_type[STRLEN];
static int    renormalize_densities;
static double rin;
static double rmax;
static double beta;
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
static double kappa_eos;
#endif
#if !GAMMA_FALLBACK
static double entropy;
static double lrho_guess;
#endif
#if RADIATION && TRACERS
static int ntracers;
#endif

// Make global so it ends on heap
static double           A[N1 + 2 * NG][N2 + 2 * NG];
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
#if !GAMMA_FALLBACK
  set_param("entropy", &entropy);
  set_param("lrho_guess", &lrho_guess);
#endif
#if RADIATION && TRACERS
  set_param("ntracers", &ntracers);
#endif
}


void init_prob() {
    double r, th, u, rho, press, X[NDIM];
    struct of_geom *geom;
    
    // Disk interior
    double hm1;
    
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
    
    set_prim(P);
    
#if EOS == EOS_TYPE_TABLE
#if !GAMMA_FALLBACK
  // guesses. Should be function local.
  double lrho0 = lrho_guess;
  // find constant entropy isocontour in table
  struct of_adiabat adiabat;
  int status = EOS_SC_find_adiabat_1d(entropy, const_ye, &adiabat);
    // SUDI: how to change const_ye ?

  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[Torus]: "
        "Failed to find isocontour for entropy\n"
        "\ts        = %e\n"
        "\tye       = %e\n",
        entropy, const_ye); // SUDI: how to change const_ye ?
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
    
    rhomax   = -INFINITY;
    umax     = -INFINITY;
    pressmax = -INFINITY;
    hm1max   = -INFINITY;
    ZSLOOP(-1, N1, -1, N2, -1, N3) {
        coord(i, j, k, CENT, X);
        bl_coord(X, &r, &th);
        // regions outside torus
        if (r < rin) {
            disk_cell[i][j][k] = 0;
        }
        // regions inside torus
        else {
            disk_cell[i][j][k] = 1;
            
            hm1 = exp(c[i][j][k]) - 1.;
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
            rho = pow(hm1 * (gam - 1.) / (kappa_eos * gam), 1. / (gam - 1.));
            u   = kappa_eos * pow(rho, gam) / (gam - 1.);
#elif EOS == EOS_TYPE_TABLE
            hm1 += hm1_min;
            if (lrho0 < EOS_SC_get_min_lrho())
                lrho0 = lrho_guess;
            // if (lrho0 > log10(1e2*RHO_unit)) lrho0 = lrho_guess;
            EOS_SC_isoentropy_hm1(hm1, &adiabat, &lrho0, &rho, &u);
            //SUDI: how values of rho and u are feed to this function ?
#else
            fprintf(stderr, "[Torus]: Bad EOS chosen.\n");
            exit(1);
#endif
            rho = P[i][j][k][RHO]
            if (rho > rhomax)
                rhomax = rho;
            u = P[i][j][k][UU]
            if (u > umax)
                umax = u;
            if (hm1 > hm1max)
              hm1max = hm1;
              P[i][j][k][UU] *= (1. + 4.e-2 * (get_rand() - 0.5));
            //SUDI: what is it doing ?
            
            // Convert from 4-velocity to 3-velocity
            coord_transform(P[i][j][k], i, j);
            //SUDI: We need to get 3-vel, right ?
        }
    }//ZSLOOP
        
#if EOS == EOS_TYPE_TABLE && !GAMMA_FALLBACK
    EOS_SC_adiabat_free(&adiabat);
#endif
    

    // get rhomax, umax globally
    // DEBUG
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
      }

      if (r > rin) {
  #if EOS == EOS_TYPE_TABLE
        EOS_SC_fill(P[i][j][k], extra[i][j][k]);
  #endif
        press = EOS_pressure_rho0_u(
            P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
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
    fixup(P, extra); // SUDI: is it required ?
    
//    // initialize hot atmosphere, SUDI: is it required?
//    // Broken.
//    #if EOS == EOS_TYPE_TABLE //&& !COLD_FLOORS
//      ZLOOP {
//        if (P[i][j][k][ATM] < ATM_THRESH) {
//          coord(i, j, k, CENT, X);
//          bl_coord(X, &r, &th);
//          P[i][j][k][UU] = P[i][j][k][RHO] / r;
//        }
//      }
//      fixup(P, extra);
//    #endif
    
// SUDI: is bound_prim required ?
//    if (mpi_io_proc()) {
//      fprintf(stdout, "Bounding prim.\n"); // debug
//    }
//    bound_prim(P);
//    if (mpi_io_proc()) {
//      fprintf(stdout, "Fixup finished.\n"); // debug
//    }
    
    
    ZLOOP {
      if (disk_cell[i][j][k]) {
        double rho_integrand =
            P[i][j][k][RHO] * ggeom[i][j][newk][CENT].g * dx[1] * dx[2] * dx[3];
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
    if (r > rin) {
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

// SUDI: magnetic field config should be tested
    //possibly classic and toroidal won't work in this setup
    //needs to be revised
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
    } else if (strcmp(bfield_type, "uniformz") == 0) {
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

#if METRIC == NUMERICAL
/* Set primitive from Carpet 3d profile */
void set_prim(grid_prim_type P){
    
    double *xp, *yp, *zp;
    double *rho, *press, *ye, *velx, *vely, *velz, *lapse, *betax, *betay, *betaz;
    int np = (N1 + 2. * NG) * (N2 + 2. * NG) * (N3 + 2. * NG);
    double extra[EOS_NUM_EXTRA];
    
    xp = malloc(np * sizeof(*xp));
    yp = malloc(np * sizeof(*yp));
    zp = malloc(np * sizeof(*zp));
    
    rho = malloc(np * sizeof(*rho));
    press = malloc(np * sizeof(*press));
    ye = malloc(np * sizeof(*ye));
    
    velx = malloc(np * sizeof(*velx));
    vely = malloc(np * sizeof(*vely));
    velz = malloc(np * sizeof(*velz));
    
    lapse = malloc(np * sizeof(*lapse));
    betax = malloc(np * sizeof(*betax));
    betay = malloc(np * sizeof(*betay));
    betaz = malloc(np * sizeof(*betaz));
    
    int iflat = 0;
    double X[NDIM];
    ZLOOPALL {
      coord(i, j, k, CENT, X);
        xp[iflat] = X[1];
        yp[iflat] = X[2];
        zp[iflat] = X[3];
        iflat++;
      }

    /* Read metadata */
    cprof3d_file_t * dfile = cprof3d_open_file(carpetprofpath);
    
    /* Open dataset: rho, press, ye, velx, vely, velz,
     lapse, betax,betay, betaz */
    cprof3d_dset_t * dset_rho = cprof3d_read_dset(dfile, "rho");
    cprof3d_dset_t * dset_press = cprof3d_read_dset(dfile, "press");
    cprof3d_dset_t * dset_ye = cprof3d_read_dset(dfile, "Ye");
    cprof3d_dset_t * dset_velx = cprof3d_read_dset(dfile, "velx");
    cprof3d_dset_t * dset_vely = cprof3d_read_dset(dfile, "vely");
    cprof3d_dset_t * dset_velz = cprof3d_read_dset(dfile, "velz");
    cprof3d_dset_t * dset_lapse = cprof3d_read_dset(dfile, "lapse");
    cprof3d_dset_t * dset_betax = cprof3d_read_dset(dfile, "betax");
    cprof3d_dset_t * dset_betay = cprof3d_read_dset(dfile, "betay");
    cprof3d_dset_t * dset_betaz = cprof3d_read_dset(dfile, "betaz");
    

    //Interpolate on the grid
    bool set_all_points_rho = cprof3d_interp(dset_rho,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            rho);
    bool set_all_points_press = cprof3d_interp(dset_press,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            press);
    bool set_all_points_ye = cprof3d_interp(dset_ye,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            ye);
    bool set_all_points_velx = cprof3d_interp(dset_velx,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            velx);
    bool set_all_points_vely = cprof3d_interp(dset_vely,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            vely);
    bool set_all_points_velz = cprof3d_interp(dset_velz,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            velz);
    bool set_all_points_lapse = cprof3d_interp(dset_lapse,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* lapse is a scalar */
            np, xp, yp, zp,
            lapse);
    bool set_all_points_betax = cprof3d_interp(dset_betax,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gzz is a scalar */
            np, xp, yp, zp,
            betax);
    bool set_all_points_betay = cprof3d_interp(dset_betay,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gzz is a scalar */
            np, xp, yp, zp,
            betay);
    bool set_all_points_betaz = cprof3d_interp(dset_betaz,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gzz is a scalar */
            np, xp, yp, zp,
            betaz);

    
    /* Check if interpolation of metric gone wrong*/
    if(!set_all_points_rho) {
        fprintf(stderr, "Something went wrong with the interpolation rho!");
        return;
    }
    if(!set_all_points_press) {
        fprintf(stderr, "Something went wrong with the interpolation press!");
        return;
    }
    if(!set_all_points_ye) {
        fprintf(stderr, "Something went wrong with the interpolation ye!");
        return;
    }
    if(!set_all_points_velx) {
        fprintf(stderr, "Something went wrong with the interpolation velx!");
        return;
    }
    if(!set_all_points_vely) {
        fprintf(stderr, "Something went wrong with the interpolation vely!");
        return;
    }
    if(!set_all_points_velz) {
        fprintf(stderr, "Something went wrong with the interpolation velz!");
        return;
    }
    if(!set_all_points_lapse) {
        fprintf(stderr, "Something went wrong with the interpolation lapse!");
        return;
    }
    if(!set_all_points_betax) {
        fprintf(stderr, "Something went wrong with the interpolation betax!");
        return;
    }
    if(!set_all_points_betay) {
        fprintf(stderr, "Something went wrong with the interpolation betay!");
        return;
    }
    if(!set_all_points_betaz) {
        fprintf(stderr, "Something went wrong with the interpolation betaz!");
        return;
    }
    
    double rho_cactus2cgs = 6.17747e+17; // 1 Cactus unit = 6.17747e+17 gm / cm^3
    double uu_cactus2cgs = 5.30986e+38; // 1 dyn / cm^2 = 1 erg / cm^3 = 5.59378e-55 / (6.6721e-6)^3 Cactus unit
    
    iflat = 0;
    ZLOOPALL{
        double rho_code = rho[iflat] *  rho_cactus2cgs / RHO_unit;
        double press_code = press[iflat] * uu_cactus2cgs / U_unit;
        
        P[i][j][k][RHO] = rho_code; //* (6.17244e+17 / (11357801703.091352)) ; // 1 rest-mass density unit in Cactus= 6.17747e+17 g/cm^3

            #if EOS == EOS_TYPE_TABLE
            {
                P[i][j][k][YE] = ye[iflat];
                PASSTYPE(YE) = PASSTYPE_NUMBER;
                P[i][j][k][YE_EM] = ye[iflat];
                P[i][j][k][ATM] = NOT_ATM;

                extra[EOS_YE] = ye[iflat];

                #if NVAR_PASSIVE >= 4
                P[i][j][k][PASSIVE_START+3] = ye[iflat] * rho_code;
                PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
                #endif
            }
            #endif

        P[i][j][k][UU] = EOS_u_press(press_code, rho_code, extra);
        P[i][j][k][B1] = 0.;
        P[i][j][k][B2] = 0.;
        P[i][j][k][B3] = 0.;
        
        //the velocities in nubhlight are U^i = u^i/u^0 where u^\mu is the contravariant four-velocity
        //v^i = u^i / W + beta^i / alpha; v^i is the Carpet profile velocity
        //W = alpha * u^t,so U^i = u^i / u^0 = alpha * v^i - beta^i.
        
        P[i][j][k][U1] = lapse[iflat] * velx[iflat] - betax[iflat];
        P[i][j][k][U2] = lapse[iflat] * vely[iflat] - betay[iflat];
        P[i][j][k][U3] = lapse[iflat] * velz[iflat] - betaz[iflat];
    
    
        iflat++;
    } // ZLOOPALL end
    
//    printf("%.16lf\n", P[0][0][0][RHO]);
//    printf("%.16lf\n", momendensx[12][14][13]);
//    printf("%lf\n", betax[0]);
//    printf("%lf\n", lapse[0] * velx[0] - betax[0]);
//    printf("%lf\n", P[0][25][25][U1]);
    
    /* Free memory */
    cprof3d_del_dset(dset_rho);
    cprof3d_del_dset(dset_press);
    cprof3d_del_dset(dset_ye);
    cprof3d_del_dset(dset_velx);
    cprof3d_del_dset(dset_vely);
    cprof3d_del_dset(dset_velz);
    
    cprof3d_del_dset(dset_lapse);
    cprof3d_del_dset(dset_betax);
    cprof3d_del_dset(dset_betay);
    cprof3d_del_dset(dset_betaz);
    
    cprof3d_close_file(dfile);
    
    free(rho);
    free(press);
    free(ye);
    free(velx);
    free(vely);
    free(velz);
    
    free(lapse);
    free(betax);
    free(betay);
    free(betaz);
    
    free(zp);
    free(yp);
    free(xp);
} // set_prim end
#endif
