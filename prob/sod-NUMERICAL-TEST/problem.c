/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOD SHOCKTUBE                                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "hdf5.h"

static double tscale;
static double pscale;
void set_problem_params()
{
  set_param("tscale", &tscale);
  set_param("pscale", &pscale);
}

void init_prob()
{
//  // Make problem nonrelativistic
//  // double tscale = 1.e-2;
//  tf /= tscale;
//  dt /= tscale;
//  DTd /= tscale;
//  DTl /= tscale;
//
//  #if METRIC != NUMERICAL
//  double extra[EOS_NUM_EXTRA];
//  double X[NDIM];
//  ZLOOP {
//    coord(i, j, k, CENT, X);
//
//    double rhogas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;
//
//    // rescale pressure in relativistic case since
//    // P must produce subluminal speeds
//    double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
//    pgas *= pscale;
//
//    #if NVAR_PASSIVE > 0 || EOS == EOS_TYPE_TABLE
//    double ye = (X[1] < 0.5 || X[1] > 1.5) ? 0.45 : 0.25;
//    #endif
//
//    #if EOS == EOS_TYPE_TABLE
//    {
//      P[i][j][k][YE] = ye;
//      PASSTYPE(YE) = PASSTYPE_NUMBER;
//      P[i][j][k][YE_EM] = ye;
//
//      extra[EOS_YE] = ye;
//
//      #if NVAR_PASSIVE >= 4
//      double yedens = ye*rhogas;
//      P[i][j][k][PASSIVE_START+3] = yedens;
//      PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
//      #endif
//    }
//    #endif
//
//    P[i][j][k][RHO] = rhogas;
//    P[i][j][k][UU] = EOS_u_press(pgas, rhogas, extra);
//    P[i][j][k][U1] = 0.;
//    P[i][j][k][U2] = 0.;
//    P[i][j][k][U3] = 0.;
//    P[i][j][k][B1] = 0.;
//    P[i][j][k][B2] = 0.;
//    P[i][j][k][B3] = 0.;
//
//
//  } // ZLOOP
//#endif

  #if METRIC == NUMERICAL
  set_prim(P); // TODO:: read ye from profile
  #endif

//  // Rescale to make problem nonrelativistic
//  ZLOOP {
//    P[i][j][k][UU] *= tscale*tscale;
//    P[i][j][k][U1] *= tscale;
//    P[i][j][k][U2] *= tscale;
//    P[i][j][k][U3] *= tscale;
//    P[i][j][k][B1] *= tscale;
//    P[i][j][k][B2] *= tscale;
//    P[i][j][k][B3] *= tscale;
//  }
}



/* Set primitive from Carpet 3d profile */
void set_prim(grid_prim_type P){
    
    double *xp, *yp, *zp;
    double *rho, *press, *ye, *velx, *vely, *velz, *lapse, *betax, *betay, *betaz;
    int np = (N1 + 2. * NG) * (N2 + 2. * NG) * (N3 + 2. * NG);
    double extra[EOS_NUM_EXTRA];
    
    double momendensx[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG], momendensy[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG], momendensz[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
    
    
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
    
        momendensx[i][j][k] = P[i][j][k][RHO] * P[i][j][k][U1];
        momendensy[i][j][k] = P[i][j][k][RHO] * P[i][j][k][U2];
        momendensz[i][j][k] = P[i][j][k][RHO] * P[i][j][k][U3];
    
    
        iflat++;
    } // ZLOOPALL end
    
//    printf("%.16lf\n", P[0][0][0][RHO]);
//    printf("%.16lf\n", momendensx[12][14][13]);
//    printf("%lf\n", betax[0]);
//    printf("%lf\n", lapse[0] * velx[0] - betax[0]);
//    printf("%lf\n", P[0][25][25][U1]);
    
    double *p_x = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double)), *p_y = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double)), *p_z = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double));
    double *ux = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double)), *uy = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double)), *uz = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double));
    
    double *rhodens = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double));
    
    int n   = 0;
    ZLOOP {
    
      p_x[n] =  momendensx[i][j][k];
      p_y[n] =  momendensy[i][j][k];
      p_z[n] =  momendensz[i][j][k];
        
      ux[n] = P[i][j][k][U1];
      uy[n] = P[i][j][k][U2];
      uz[n] = P[i][j][k][U3];
        
      rhodens[n] = P[i][j][k][RHO];
    
      n++;
    }

    //printf("%.16lf\n", p_x[12]);
    
    
    // Create an HDF5 file
    hid_t file = H5Fcreate("momendens.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    int np3 = N1 * N2 * N3;

    // Create a dataspace for the dataset
    hsize_t dims[1] = {np3};
    hid_t dataspace = H5Screate_simple(1, dims, NULL);

    // Create dataset for px
    hid_t dataset;
    dataset = H5Dcreate2(file, "px", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_x);
    H5Dclose(dataset);

    // Create dataset for py
    dataset = H5Dcreate2(file, "py", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_y);
    H5Dclose(dataset);
//
    // Create dataset for pz
    dataset = H5Dcreate2(file, "pz", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_z);
    H5Dclose(dataset);
    
    // Create dataset for ux
    dataset = H5Dcreate2(file, "ux", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ux);
    H5Dclose(dataset);

    // Create dataset for uy
    dataset = H5Dcreate2(file, "uy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uy);
    H5Dclose(dataset);
//
    // Create dataset for pz
    dataset = H5Dcreate2(file, "uz", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, uz);
    H5Dclose(dataset);
    
    // Create dataset for rhodens
    dataset = H5Dcreate2(file, "rhodens", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rhodens);
    H5Dclose(dataset);
    
    // Close everything
    H5Sclose(dataspace);
    //H5Gclose(group);  // optional
    H5Fclose(file);
    
    free(p_x);
    free(p_y);
    free(p_z);
    
    free(ux);
    free(uy);
    free(uz);
    
    
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
