/******************************************************************************
 *                                                                            *
 * problem.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR NUMERICAL DISK                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Local functions
void    set_prim(grid_prim_type P);

static char   bfield_type[STRLEN];
static int    renormalize_densities;
static double beta;
#if RADIATION && TRACERS
static int ntracers;
#endif

// Make global so it ends on heap
static double           A[N1 + 2 * NG][N2 + 2 * NG];

void set_problem_params() {
  // supports: none, toroidal, classic
  set_param("bfield", &bfield_type);
  set_param("beta", &beta);
  set_param("renorm_dens", &renormalize_densities);
#if RADIATION && TRACERS
  set_param("ntracers", &ntracers);
#endif
}


void init_prob() {
    double r, th, u, rho, press, X[NDIM];
    struct of_geom *geom;
    
    // Diagnostics for entropy
    double ent, entmax;
    
    // Magnetic field
    double rho_av, rhomax, umax, pressmax, bsq_ij, bsq_max, q;
    
    //read from 3d profile
    set_prim(P);
    
    rhomax   = -INFINITY;
    umax     = -INFINITY;
    pressmax = -INFINITY;
    ZSLOOP(-1, N1, -1, N2, -1, N3) {
        rho = P[i][j][k][RHO];
        if (rho > rhomax)
            rhomax = rho;
        u = P[i][j][k][UU];
        if (u > umax)
            umax = u;
    }//ZSLOOP

    // get rhomax, umax globally
    // DEBUG
    umax   = mpi_max(umax);
    rhomax = mpi_max(rhomax);

    if (mpi_io_proc()) {
      fprintf(stdout, "rhomax   = %f\n", rhomax);
      fprintf(stdout, "umax     = %f\n", umax);
    }
    
    // Normalize densities
    ZSLOOP(-1, N1, -1, N2, -1, N3) {
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);

      if (renormalize_densities) {
        P[i][j][k][RHO] /= rhomax;
        P[i][j][k][UU] /= rhomax;
      }

        EOS_SC_fill(P[i][j][k], extra[i][j][k]);
        press = EOS_pressure_rho0_u(
            P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
        if (press > pressmax)
          pressmax = press;
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
    if (mpi_io_proc()) {
      fprintf(stdout, "Fixup finished.\n"); // debug
    }
    
    if (mpi_io_proc()) {
      fprintf(stdout, "Bounding prim.\n"); // debug
    }
    bound_prim(P);
    
//// debug
    if (mpi_io_proc()) {
    fprintf(stdout, "Calculating max entropy:\n");
    }
    entmax = -1.0;
    ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    EOS_SC_fill(P[i][j][k], extra[i][j][k]);
    ent = EOS_entropy_rho0_u(P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
      if (ent > entmax)
        entmax = ent;
    }

    entmax = mpi_max(entmax);
    if (mpi_io_proc()) {
    fprintf(stdout, "Maximum entropy in disk = %e\n", entmax);
    }

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


//#if METRIC == NUMERICAL
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

            
        P[i][j][k][YE] = ye[iflat];
        PASSTYPE(YE) = PASSTYPE_NUMBER;
        P[i][j][k][YE_EM] = ye[iflat];
        P[i][j][k][ATM] = NOT_ATM;

        extra[EOS_YE] = ye[iflat];

        #if NVAR_PASSIVE >= 4
        P[i][j][k][PASSIVE_START+3] = ye[iflat] * rho_code;
        PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
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
//#endif
