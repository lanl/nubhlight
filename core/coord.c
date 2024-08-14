// ======================================================================
//  copyright 2020. Triad National Security, LLC. All rights
//  reserved. This program was produced under U.S. Government contract
//  89233218CNA000001 for Los Alamos National Laboratory (LANL), which
//  is operated by Triad National Security, LLC for the U.S. Department
//  of Energy/National Nuclear Security Administration. All rights in
//  the program are reserved by Triad National Security, LLC, and the
//  U.S. Department of Energy/National Nuclear Security
//  Administration. The Government is granted for itself and others
//  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
//  license in this material to reproduce, prepare derivative works,
//  distribute copies to the public, perform publicly and display
//  publicly, and to permit others to do so.
// ======================================================================

/******************************************************************************
 *                                                                            *
 * COORD.C                                                                    *
 *                                                                            *
 * SIMULATION VOLUME COORDINATES                                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
//#include "cprof3d.h"

/*      ASSUMING X3 SYMMETRY IN COORDINATES AND SPACETIME
 *      -- given the indices i,j and location in the cell, return with
 *         the values of X1,X2 there;
 *      -- the locations are defined by :
 *          -----------------------
 *          |                     |
 *          |                     |
 *          |FACE1   CENT         |
 *          |                     |
 *          |CORN    FACE2        |
 *          ----------------------
 *
 */

void coord(int i, int j, int k, int loc, double X[NDIM]) {
  X[0] = t;

  i += global_start[1];
  j += global_start[2];
  k += global_start[3];
  if (loc == FACE1) {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == FACE2) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == FACE3) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k - NG) * dx[3];
  } else if (loc == CENT) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == CORN) {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k - NG) * dx[3];
  } else {
    fprintf(stderr, "Invalid coordinate location!\n");
    exit(1);
  }
}

// Assumes Boyer-Lindquist coordinates
int bl_i_of_r(double r) {
  double X1  = log(r);
  int    ind = (int)((X1 - startx[1]) / dx[1]);
  return ind;
}

double thG_of_X(const double X[NDIM]) {
  return M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
}

#if DEREFINE_POLES
void thJ_of_X(const double X[NDIM], double *y, double *thJ) {
  *y   = 2. * X[2] - 1.;
  *thJ = poly_norm * (*y) *
             (1. + pow((*y) / poly_xt, poly_alpha) / (poly_alpha + 1.)) +
         0.5 * M_PI;
}
#endif // DEREFINE_POLES

double r_of_X(const double X[NDIM]) { return exp(X[1]); }

double th_of_X(const double X[NDIM]) {
  double thG = thG_of_X(X);

#if DEREFINE_POLES
  double y, thJ;
  thJ_of_X(X, &y, &thJ);
  return thG + exp(mks_smooth * (startx[1] - X[1])) * (thJ - thG);
#else
  return thG;
#endif
}

// Boyer-Lindquist coordinate of point X
void bl_coord(const double X[NDIM], double *r, double *th) {
  *r  = r_of_X(X);
  *th = th_of_X(X);

// Avoid singularity at polar axis
#if (COORDSINGFIX)
  if (fabs(*th) < SINGSMALL) {
    if ((*th) >= 0)
      *th = SINGSMALL;
    if ((*th) < 0)
      *th = -SINGSMALL;
  }
  if (fabs(M_PI - (*th)) < SINGSMALL) {
    if ((*th) >= M_PI)
      *th = M_PI + SINGSMALL;
    if ((*th) < M_PI)
      *th = M_PI - SINGSMALL;
  }
#endif
}

// Cartesian coordinate Xcart = {t,x,y,z} of point X
void cart_coord(const double X[NDIM], double Xcart[NDIM]) {
  Xcart[0] = X[0]; // t
#if METRIC == MKS
  {
    double r, th, ph;
    bl_coord(X, &r, &th);
    ph       = X[3];
    Xcart[1] = r * sin(th) * cos(ph); // r
    Xcart[2] = r * sin(th) * sin(ph); // th
    Xcart[3] = r * cos(th);           // ph
  }
#else
  {
    Xcart[1] = X[1];
    Xcart[2] = X[2];
    Xcart[3] = X[3];
  }
#endif // METRIC == MKS
}

// Jacobians.
// Jcon maps contravariant 4-vector from HARM coordinates
// to Boyer-Lindquist coordinates
// Jcov maps covariant 4-vector from HARM coordinates
// to Boyer-Lindquist coordinates
// CONVENTION: rows are top index, columns are bottom index
//             not required for tensor notation.
//             Is required for matrix notation.
void jac_harm_to_bl(
    const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]) {
  double r, th;
  bl_coord(X, &r, &th);

  double drdX1   = exp(X[1]);
  double dthGdX2 = M_PI + M_PI * (1 - hslope) * cos(2 * M_PI * X[2]);
  double dthdX1, dthdX2;
#if DEREFINE_POLES
  {
    double thG = thG_of_X(X);
    double y, thJ;
    thJ_of_X(X, &y, &thJ);
    double dydX2   = 2.;
    double dthJdy  = poly_norm * (1 + pow(y / poly_xt, poly_alpha));
    double dthJdX2 = dthJdy * dydX2;
    dthdX1 = -mks_smooth * (thJ - thG) * exp(mks_smooth * (startx[1] - X[1]));
    dthdX2 =
        dthGdX2 + exp(mks_smooth * (startx[1] - X[1])) * (dthJdX2 - dthGdX2);
  }
#else
  {
    dthdX1 = 0.0;
    dthdX2 = dthGdX2;
  }
#endif

  // First set all values to zero. Then we'll set nonzero values
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Jcon[mu][nu] = 0.;
      Jcov[mu][nu] = 0.;
    }
  }
  // Jcon
  Jcon[0][0] = 1.;     // t
  Jcon[1][1] = drdX1;  // r
  Jcon[2][1] = dthdX1; // th
  Jcon[2][2] = dthdX2;
  Jcon[3][3] = 1.; // phi
  // Jcov
  Jcov[0][0] = 1.;                                 // t
  Jcov[1][1] = 1. / (drdX1 + SMALL);               // r
  Jcov[2][1] = -dthdX1 / (drdX1 * dthdX2 + SMALL); // th
  Jcov[2][2] = 1. / (dthdX2 + SMALL);
  Jcov[3][3] = 1.;
}

// Jcon maps contravariant 4-vector from BL coordinates to Minkowski
// Jcov maps a covariant one from BL to Minkowski
void jac_bl_to_cart(
    const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]) {
  double r, th, ph;
  double Xcart[NDIM];

  bl_coord(X, &r, &th);
  ph = X[3];
  cart_coord(X, Xcart);

  // First set all values to zero. Then we'll set nonzero values
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Jcon[mu][nu] = 0.;
      Jcov[mu][nu] = 0.;
    }
  }
  // Jcon
  Jcon[0][0] = 1.;                // t
  Jcon[1][1] = sin(th) * cos(ph); // x
  Jcon[1][2] = r * cos(th) * cos(ph);
  Jcon[1][3] = -r * sin(th) * sin(ph);
  Jcon[2][1] = sin(th) * sin(ph); // y
  Jcon[2][2] = r * cos(th) * sin(ph);
  Jcon[2][3] = r * sin(th) * cos(ph);
  Jcon[3][1] = cos(th); // z
  Jcon[3][2] = -r * sin(th);
  // Jcov
  Jcov[0][0] = 1.;
  Jcov[1][1] = cos(ph) * sin(th);
  Jcov[1][2] = sin(ph) * sin(th);
  Jcov[1][3] = cos(th);
  Jcov[2][1] = cos(ph) * cos(th) / (r + SMALL);
  Jcov[2][2] = cos(th) * sin(ph) / (r + SMALL);
  Jcov[2][3] = -sin(th) / (r + SMALL);
  Jcov[3][1] = -sin(ph) / (r * sin(th) + SMALL);
  Jcov[3][2] = cos(ph) / (r * sin(th) + SMALL);
}

void jac_harm_to_cart(
    const double X[NDIM], double Jcov[NDIM][NDIM], double Jcon[NDIM][NDIM]) {
#if METRIC == MKS
  {
    double J_h2bl_cov[NDIM][NDIM], J_h2bl_con[NDIM][NDIM];
    double J_bl2c_cov[NDIM][NDIM], J_bl2c_con[NDIM][NDIM];
    jac_harm_to_bl(X, J_h2bl_cov, J_h2bl_con);
    jac_bl_to_cart(X, J_bl2c_cov, J_bl2c_con);

    for (int mupp = 0; mupp < NDIM; mupp++) {
      for (int mu = 0; mu < NDIM; mu++) {
        Jcon[mupp][mu] = 0.0;
        Jcov[mu][mupp] = 0.0;
        for (int mup = 0; mup < NDIM; mup++) {
          Jcon[mupp][mu] += J_bl2c_con[mupp][mup] * J_h2bl_con[mup][mu];
          Jcov[mu][mupp] += J_h2bl_cov[mu][mup] * J_bl2c_cov[mup][mupp];
        }
      }
    }
  }
#else
  { DLOOP2 Jcon[mu][nu] = Jcov[mu][nu] = (mu == nu ? 1 : 0); }
#endif // METRIC
}

// TODO: This is redundant with jac_harm_to_bl.
// We should reduce duplicate code.
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]) {
  memset(dxdX, 0, NDIM * NDIM * sizeof(double));
#if METRIC == MINKOWSKI
  for (int mu = 0; mu < NDIM; mu++) {
    dxdX[mu][mu] = 1.;
  }
#elif METRIC == MKS
  double Jcov[NDIM][NDIM];
  jac_harm_to_bl(X, Jcov, dxdX);
  /*
  // This is the old formula, less readable than the new one
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
#if DEREFINE_POLES
  dxdX[2][1] = -exp(mks_smooth * (startx[1] - X[1])) * mks_smooth *
               (M_PI / 2. - M_PI * X[2] +
                   poly_norm * (2. * X[2] - 1.) *
                       (1 + (pow((-1. + 2 * X[2]) / poly_xt, poly_alpha)) /
                                (1 + poly_alpha)) -
                   1. / 2. * (1. - hslope) * sin(2. * M_PI * X[2]));
  dxdX[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]) +
               exp(mks_smooth * (startx[1] - X[1])) *
                   (-M_PI +
                       2. * poly_norm *
                           (1. + pow((2. * X[2] - 1.) / poly_xt, poly_alpha) /
                                     (poly_alpha + 1.)) +
                       (2. * poly_alpha * poly_norm * (2. * X[2] - 1.) *
                           pow((2. * X[2] - 1.) / poly_xt, poly_alpha - 1.)) /
                           ((1. + poly_alpha) * poly_xt) -
                       (1. - hslope) * M_PI * cos(2. * M_PI * X[2]));
#else
  dxdX[2][2] = M_PI - (hslope - 1.) * M_PI * cos(2. * M_PI * X[2]);
#endif // derefine_poles
  dxdX[3][3] = 1.;
  // debug
  // Use this to check that the old and new formulae agree
  double Jcov[NDIM][NDIM];
  double Jcon[NDIM][NDIM];
  jac_harm_to_bl(X, Jcov, Jcon);
  DLOOP2 {
    if ((fabs(Jcon[mu][nu]) > 1e-10 || fabs(dxdX[mu][nu]) > 1e-10) &&
        (fabs(Jcon[mu][nu] - dxdX[mu][nu]) >=
            1e-8 * 0.5 * (fabs(Jcon[mu][nu]) + fabs(dxdX[mu][nu])))) {
      printf("BAD Jcon! %d %d, %e %e\n", mu, nu, Jcon[mu][nu], dxdX[mu][nu]);
      exit(1);
    }
  }
  */
#endif // metric
}

// Insert metric here
void set_metric(double X[NDIM], struct of_geom *g) {
  set_gcov(X, g->gcov);
  g->g     = gcon_func(g->gcov, g->gcon);
  g->alpha = 1.0 / sqrt(-(g->gcon[0][0]));
   // printf("The set_metric is called \n");
}

#if METRIC == NUMERICAL // we need to add this if statement in order to avoid
// compilation failure for non-numerical metric
/* Set numerical metric */
void num_set_metric(grid_geom_type ggeom){
    
    double *xp, *yp, *zp;
    
    double *gxx, *gxy, *gxz, *gyy, *gyz, *gzz, *lapse, *betax, *betay, *betaz;
    
    double betaxcon,betaycon,betazcon,beta2;
    
    int np = (N1 + 2. * NG) * (N2 + 2. * NG) * (N3 + 2. * NG) * NPG;
    
    //int np = N1 * N2 * N3 * NPG;
    
    xp = malloc(np * sizeof(*xp));
    yp = malloc(np * sizeof(*yp));
    zp = malloc(np * sizeof(*zp));
    
    gxx = malloc(np * sizeof(*gxx));
    gxy = malloc(np * sizeof(*gxy));
    gxz = malloc(np * sizeof(*gxz));

    gyy = malloc(np * sizeof(*gyy));
    gyz = malloc(np * sizeof(*gyz));

    gzz = malloc(np * sizeof(*gzz));
    
    lapse = malloc(np * sizeof(*lapse));
    
    betax = malloc(np * sizeof(*betax));
    betay = malloc(np * sizeof(*betay));
    betaz = malloc(np * sizeof(*betaz));
    
    
    int iflat = 0;
    double X[NDIM];
    ZSLOOP(-NG, N1 - 1 + NG,-NG, N2 - 1 + NG,-NG, N3 - 1 + NG) {
//    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3 - 1) {
      LOCLOOP {
        coord(i, j, k, loc, X);
        //int const iflat = i * (N2 + 2. * NG) * (N3 + 2. * NG) + j * (N3 + 2. * NG) + k;
        xp[iflat] = X[1];
        yp[iflat] = X[2];
        zp[iflat] = X[3];
        //printf("The value of coordinate tuple X[1]: %f, X[2]: %f, X[3] %f", xp[iflat], yp[iflat], zp[iflat]);
        iflat++;
      }
    }
    
    //printf("The value of coordinate tuple X[1]: %f, X[2]: %f, X[3] %f", xp[0], yp[0], zp[0]);

    
    /* Read metadata */
    cprof3d_file_t * dfile = cprof3d_open_file(carpetprofpath); //TODO:pass through runtime arguments
    
    /* Open dataset: gxx,gxy,gxz,gyy,gyz,gzz */
    cprof3d_dset_t * dset_gxx = cprof3d_read_dset(dfile, "gxx");
    cprof3d_dset_t * dset_gxy = cprof3d_read_dset(dfile, "gxy");
    cprof3d_dset_t * dset_gxz = cprof3d_read_dset(dfile, "gxz");
    cprof3d_dset_t * dset_gyy = cprof3d_read_dset(dfile, "gyy");
    cprof3d_dset_t * dset_gyz = cprof3d_read_dset(dfile, "gyz");
    cprof3d_dset_t * dset_gzz = cprof3d_read_dset(dfile, "gzz");
    
    /* Open dataset: lapse */
    cprof3d_dset_t * dset_lapse = cprof3d_read_dset(dfile, "lapse");
    
    /* Open dataset: shift */
    cprof3d_dset_t * dset_betax = cprof3d_read_dset(dfile, "betax");
    cprof3d_dset_t * dset_betay = cprof3d_read_dset(dfile, "betay");
    cprof3d_dset_t * dset_betaz = cprof3d_read_dset(dfile, "betaz");

    /* Interpolate metric on the grid*/
    bool set_all_points_gxx = cprof3d_interp(dset_gxx,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gxx is a scalar */
            np, xp, yp, zp,
            gxx);
    bool set_all_points_gxy = cprof3d_interp(dset_gxy,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gxy is a scalar */
            np, xp, yp, zp,
            gxy);
    bool set_all_points_gxz = cprof3d_interp(dset_gxz,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gxz is a scalar */
            np, xp, yp, zp,
            gxz);
    bool set_all_points_gyy = cprof3d_interp(dset_gyy,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gyy is a scalar */
            np, xp, yp, zp,
            gyy);
    bool set_all_points_gyz = cprof3d_interp(dset_gyz,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gyz is a scalar */
            np, xp, yp, zp,
            gyz);
    bool set_all_points_gzz = cprof3d_interp(dset_gzz,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gzz is a scalar */
            np, xp, yp, zp,
            gzz);
    
    /* Interpolate lapse on the grid*/
    bool set_all_points_lapse = cprof3d_interp(dset_lapse,
            cprof3d_cmap_reflecting_xy, /* default */
            cprof3d_transf_default,     /* gzz is a scalar */
            np, xp, yp, zp,
            lapse);
    
    /* Interpolate shift on the grid*/
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
    if(!set_all_points_gxx) {
        fprintf(stderr, "Something went wrong with the interpolation gxx!");
        return;
    }
    if(!set_all_points_gxy) {
        fprintf(stderr, "Something went wrong with the interpolation gxy!");
        return;
    }
    if(!set_all_points_gxz) {
        fprintf(stderr, "Something went wrong with the interpolation gxz!");
        return;
    }if(!set_all_points_gyy) {
        fprintf(stderr, "Something went wrong with the interpolation gyy!");
        return;
    }
    if(!set_all_points_gyz) {
        fprintf(stderr, "Something went wrong with the interpolation gyz!");
        return;
    }
    if(!set_all_points_gzz) {
        fprintf(stderr, "Something went wrong with the interpolation gzz!");
        return;
    }
    
    /* Check if interpolation of lapse gone wrong*/
    if(!set_all_points_lapse) {
        fprintf(stderr, "Something went wrong with the interpolation lapse!");
        return;
    }
    
    /* Check if interpolation of shift gone wrong*/
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
    
    iflat = 0;
    
    ZSLOOP(-NG, N1 - 1 + NG,-NG, N2 - 1 + NG,-NG, N3 - 1 + NG) {
//    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3 - 1) {
        LOCLOOP {
            struct of_geom *g = &ggeom[i][j][newk][loc];
            memset(g->gcov, 0, NDIM * NDIM * sizeof(double)); // initialization and memory allocation for gcov

            g->gcov[1][1] = gxx[iflat];
            g->gcov[1][2] = gxy[iflat];
            g->gcov[1][3] = gxz[iflat];
              
            g->gcov[2][1] = g->gcov[1][2];
            g->gcov[2][2] = gyy[iflat];
            g->gcov[2][3] = gyz[iflat];

            g->gcov[3][1] = g->gcov[1][3];
            g->gcov[3][2] = g->gcov[2][3];
            g->gcov[3][3] = gzz[iflat];
            
            betaxcon = g->gcov[1][1]*betax[iflat];
            betaycon = g->gcov[2][2]*betay[iflat];
            betazcon = g->gcov[3][3]*betaz[iflat];
            
            g->gcov[0][1] = betaxcon;
            g->gcov[0][2] = betaycon;
            g->gcov[0][3] = betazcon;
            
            g->gcov[1][0] = g->gcov[0][1];
            g->gcov[2][0] = g->gcov[0][2];
            g->gcov[3][0] = g->gcov[0][3];
          
            beta2 =  betaxcon * betax[iflat] + betaycon * betay[iflat] + betazcon * betaz[iflat];

            g->gcov[0][0] = -lapse[iflat]*lapse[iflat] + beta2;
            // setting detg and alpha
            g->g     = gcon_func(g->gcov, g->gcon);
            g->alpha = lapse[iflat]; //1.0 / sqrt(-(g->gcon[0][0]));
            iflat++;
        } // LOCLOOP end
    } // ZSLOOP end
    
    /* Free memory */
    cprof3d_del_dset(dset_gxx);
    cprof3d_del_dset(dset_gxy);
    cprof3d_del_dset(dset_gxz);
    cprof3d_del_dset(dset_gyy);
    cprof3d_del_dset(dset_gyz);
    cprof3d_del_dset(dset_gzz);
    
    cprof3d_del_dset(dset_lapse);
    cprof3d_del_dset(dset_betax);
    cprof3d_del_dset(dset_betay);
    cprof3d_del_dset(dset_betaz);
    
    cprof3d_close_file(dfile);
    
    free(gxx);
    free(gxy);
    free(gxz);
    free(gyy);
    free(gyz);
    free(gzz);
    
    free(lapse);
    free(betax);
    free(betay);
    free(betaz);
    
    free(zp);
    free(yp);
    free(xp);
} // num_set_metric end

#endif

void set_gcov(double X[NDIM], double gcov[NDIM][NDIM]) {
  memset(gcov, 0, NDIM * NDIM * sizeof(double));

#if METRIC == MINKOWSKI
  gcov[0][0] = -1.;
  for (int j = 1; j < NDIM; j++) {
    gcov[j][j] = 1.;
      
  }
#elif METRIC == SPHERICAL
  gcov[0][0] = -1.;
  gcov[1][1] = 1.;
  gcov[2][2] = pow(X[1], 2.);
  gcov[3][3] = pow(X[1] * sin(X[2]), 2.);
#elif METRIC == MKS
  double sth, cth, s2, rho2;
  double r, th;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2   = sth * sth;
  rho2 = r * r + a * a * cth * cth;

  gcov[0][0] = -1. + 2. * r / rho2;
  gcov[0][1] = 2. * r / rho2;
  gcov[0][3] = -2. * a * r * s2 / rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2. * r / rho2;
  gcov[1][3] = -a * s2 * (1. + 2. * r / rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2));

  // Apply coordinate transformation to code coordinates X
  double dxdX[NDIM][NDIM];
  set_dxdX(X, dxdX);

  double gcov_ks[NDIM][NDIM];
  memcpy(gcov_ks, gcov, NDIM * NDIM * sizeof(double));
  memset(gcov, 0, NDIM * NDIM * sizeof(double));

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int lam = 0; lam < NDIM; lam++) {
        for (int kap = 0; kap < NDIM; kap++) {
          gcov[mu][nu] += gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
        }
      }
    }
  }
#endif // METRIC
}

// Establish X coordinates
void set_points() {
#if METRIC == MINKOWSKI || METRIC == SPHERICAL || METRIC == NUMERICAL
  startx[1] = x1Min;
  startx[2] = x2Min;
  startx[3] = x3Min;
  dx[1]     = (x1Max - x1Min) / N1TOT;
  dx[2]     = (x2Max - x2Min) / N2TOT;
  dx[3]     = (x3Max - x3Min) / N3TOT;
#if RADIATION
  startx_rad[1] = startx[1];
  startx_rad[2] = startx[2];
  startx_rad[3] = startx[3];
  stopx_rad[1]  = startx_rad[1] + N1TOT * dx[1];
  stopx_rad[2]  = startx_rad[2] + N2TOT * dx[2];
  stopx_rad[3]  = startx_rad[3] + N3TOT * dx[3];
#endif
#elif METRIC == MKS
  // Calculate some radii determined by the geometry
  Reh = 1. + sqrt(1. - a * a);
  double z1 =
      1 + pow(1 - a * a, 1. / 3.) * (pow(1 + a, 1. / 3.) + pow(1 - a, 1. / 3.));
  double z2 = sqrt(3 * a * a + z1 * z1);
  Risco     = 3 + z2 - sqrt((3 - z1) * (3 + z1 + 2 * z2));

  // Set Rin such that we have 5 zones completely inside the event horizon
  // If xeh = log(Reh), xin = log(Rin), and xout = log(Rout),
  // then we want xeh = xin + 5.5 * (xout - xin) / N1TOT, or solving/replacing:
  Rin = exp((N1TOT * log(Reh) / 5.5 - log(Rout)) / (-1. + N1TOT / 5.5));

  startx[1]     = log(Rin);
  startx[2]     = 0.;
  startx[3]     = 0.;
  dx[1]         = log(Rout / Rin) / N1TOT;
  dx[2]         = 1. / N2TOT;

#if QUADRANT_SYMMETRY
  dx[3]         = (M_PI / 2.) / N3TOT;
#else
  dx[3] = 2. * M_PI / N3TOT;
#endif // QUADRANT_SYMMETRY

#if RADIATION
  startx_rad[1] = log(Reh);
  startx_rad[2] = startx[2];
  startx_rad[3] = startx[3];
  stopx_rad[1]  = log(Rout_rad);
  stopx_rad[2]  = startx[2] + N2TOT * dx[2];
  stopx_rad[3]  = startx[3] + N3TOT * dx[3];
#endif
  poly_norm     = 0.5 * M_PI * 1. /
              (1. + 1. / (poly_alpha + 1.) * 1. / pow(poly_xt, poly_alpha));
// mks_smooth = 0.5;
#endif // MKS
  startx_proc[1] = startx[1] + global_start[1] * dx[1];
  startx_proc[2] = startx[2] + global_start[2] * dx[2];
  startx_proc[3] = startx[3] + global_start[3] * dx[3];
  stopx_proc[1]  = startx_proc[1] + dx[1] * N1;
  stopx_proc[2]  = startx_proc[2] + dx[2] * N2;
  stopx_proc[3]  = startx_proc[3] + dx[3] * N3;

  stopx[1] = startx[1] + dx[1] * N1TOT;
  stopx[2] = startx[2] + dx[2] * N2TOT;
  stopx[3] = startx[3] + dx[3] * N3TOT;
}

void set_grid() {

  // Set up boundaries, steps in coordinate grid
  set_points();
  dV = dx[1] * dx[2] * dx[3];


#if RADIATION
  dt_light_min = 1. / SMALL;
#endif

#if METRIC == NUMERICAL //TODO: make radiation work with numerical
//    ZSLOOP(-NG, N1 - 1 + NG,-NG, N2 - 1 + NG,-NG, N3 - 1 + NG){
//        LOCLOOP{
//            coord(i, j, k, loc, X);
//        }
//    }
    
    num_set_metric(ggeom);
    
    double X[NDIM];
    
    coord(16, 17, 18, CENT, X);
    
//    printf("Test ggeom gcov[3][0] %f %f %f %f\n", ggeom[16][17][18][CENT].gcov[3][0], X[1],X[2],X[3]);
//    printf("Test ggeom gcov[3][1] %f %f %f %f\n", ggeom[16][17][18][CENT].gcov[3][1], X[1],X[2],X[3]);
//    printf("Test ggeom gcov[3][2] %f %f %f %f\n", ggeom[16][17][18][CENT].gcov[3][2], X[1],X[2],X[3]);
//    printf("Test ggeom gcov[3][3] %f %f %f %f\n", ggeom[16][17][18][CENT].gcov[3][3], X[1],X[2],X[3]);
    
    ZSLOOP(-NG, N1 - 1 + NG,-NG, N2 - 1 + NG,-NG, N3 - 1 + NG){ //Sudi: is this correct ??
        
    //ZSLOOP(-NG, N1 - 1,-NG, N2 - 1,-NG, N3 - 1){

        // Only required in zone center /* connection coefficient calculation */
        num_conn_func(ggeom, conn[i][j][newk], i , j , newk);
        
        //printf("Test ggeom %f \n", conn[i][j][newk][1][1][1]);
        
#else
        double X[NDIM];
        for (int mu = 0; mu < NDIM; mu++)
          X[mu] = 0.;
  //ISLOOP(-NG, N1 - 1 + NG) {
    //JSLOOP(-NG, N2 - 1 + NG) {
    ZSLOOP(-NG, N1 - 1 + NG,-NG, N2 - 1 + NG,-NG, N3 - 1 + NG){
        LOCLOOP { // loop variable is loc
            coord(i, j, 0, loc, X);
            set_metric(X, &(ggeom[i][j][newk][loc])); // should I set k = 1 ?
        } // LOCLOOP
        
        // Only required in zone center
        conn_func(X, &ggeom[i][j][newk][CENT], conn[i][j][newk]); // why conn[i][j] ?? should I set k = 1 ?
        
#endif
        
#if RADIATION
        // Set minimum light crossing time for each zone
        dt_light[i][j]           = 1.e30;
        double light_phase_speed = SMALL;
        double dt_light_local    = 0.;
        
        for (int mu = 1; mu < NDIM; mu++) {
            if (pow(ggeom[i][j][newk][CENT].gcon[0][mu], 2.) -
                ggeom[i][j][newk][CENT].gcon[mu][mu] * ggeom[i][j][newk][CENT].gcon[0][0] >=
                0.) {
                double cplus      = fabs((-ggeom[i][j][newk][CENT].gcon[0][mu] +
                                          sqrt(pow(ggeom[i][j][newk][CENT].gcon[0][mu], 2.) -
                                               ggeom[i][j][newk][CENT].gcon[mu][mu] *
                                               ggeom[i][j][newk][CENT].gcon[0][0])) /
                                         (ggeom[i][j][newk][CENT].gcon[0][0]));
                double cminus     = fabs((-ggeom[i][j][newk][CENT].gcon[0][mu] -
                                          sqrt(pow(ggeom[i][j][newk][CENT].gcon[0][mu], 2.) -
                                               ggeom[i][j][newk][CENT].gcon[mu][mu] *
                                               ggeom[i][j][newk][CENT].gcon[0][0])) /
                                         (ggeom[i][j][newk][CENT].gcon[0][0]));
                light_phase_speed = MY_MAX(cplus, cminus);
            } else {
                light_phase_speed = SMALL;
            }
            
            dt_light_local += 1. / (dx[mu] / light_phase_speed);
            
            if (dx[mu] / light_phase_speed < dt_light[i][j]) {
                dt_light[i][j] = dx[mu] / light_phase_speed;
            }
        }
        dt_light_local = 1. / dt_light_local;
        if (dt_light_local < dt_light_min)
            dt_light_min = dt_light_local;
#endif // RADIATION
        //    }  // JSLOOP
        //  }    // ISLOOP
    } // ZSLOOP

#if RADIATION
  /*dt_light_min = 1.e300;
  ISLOOP(0, N1 - 1) {
    JSLOOP(0, N2 - 1) {
      if (dt_light_min > dt_light[i][j]) {
        dt_light_min = dt_light[i][j];
      }
    }
  }*/
  dt_light_min = mpi_min(dt_light_min);
#endif // radiation
}

void cart_to_sph(const double X[NDIM], double *r, double *th, double *ph) {
  double x = X[1];
  double y = X[2];
  double z = X[3];
  *r       = sqrt(x * x + y * y + z * z);
  *th      = acos(z / (*r));
  *ph      = atan2(y, x);
  // atan2 has range [-pi,pi]. I want [0,2*pi].
  if (*ph < 0)
    *ph += 2 * M_PI;
}

void zero_arrays() {
  ZSLOOP(-NG, N1 - 1 + NG, -NG, N2 - 1 + NG, -NG, N3 - 1 + NG) {
    PLOOP {
      P[i][j][k][ip]  = 0.;
      Ph[i][j][k][ip] = 0.;
      F1[i][j][k][ip] = 0.;
      F2[i][j][k][ip] = 0.;
      F3[i][j][k][ip] = 0.;
    }
    pflag[i][j][k]          = 0;
    fail_save[i][j][k]      = 0;
    fixup_required[i][j][k] = 0;
  } // ZSLOOP
}
