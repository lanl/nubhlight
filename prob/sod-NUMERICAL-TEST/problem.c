/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOD SHOCKTUBE                                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
//#include "cprof3d.h"

static double tscale;
static double pscale;
void set_problem_params()
{
  set_param("tscale", &tscale);
  set_param("pscale", &pscale);
}

void init_prob()
{
  // Make problem nonrelativistic
  // double tscale = 1.e-2;
  tf /= tscale;
  dt /= tscale;
  DTd /= tscale;
  DTl /= tscale;

  #if METRIC != NUMERICAL
  double extra[EOS_NUM_EXTRA];
  double X[NDIM];
  ZLOOP {
    coord(i, j, k, CENT, X);
      
    double rhogas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;

    // rescale pressure in relativistic case since
    // P must produce subluminal speeds
    double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
    pgas *= pscale;

    #if NVAR_PASSIVE > 0 || EOS == EOS_TYPE_TABLE
    double ye = (X[1] < 0.5 || X[1] > 1.5) ? 0.45 : 0.25;
    #endif

    #if EOS == EOS_TYPE_TABLE
    {
      P[i][j][k][YE] = ye;
      PASSTYPE(YE) = PASSTYPE_NUMBER;
      P[i][j][k][YE_EM] = ye;

      extra[EOS_YE] = ye;

      #if NVAR_PASSIVE >= 4
      double yedens = ye*rhogas;
      P[i][j][k][PASSIVE_START+3] = yedens;
      PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
      #endif
    }
    #endif

    P[i][j][k][RHO] = rhogas;
    P[i][j][k][UU] = EOS_u_press(pgas, rhogas, extra);
    P[i][j][k][U1] = 0.;
    P[i][j][k][U2] = 0.;
    P[i][j][k][U3] = 0.;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
    

  } // ZLOOP
#endif

  #if METRIC == NUMERICAL
  set_prim(P); // TODO:: read ye from profile
  #endif

  // Rescale to make problem nonrelativistic
  ZLOOP {
    P[i][j][k][UU] *= tscale*tscale;
    P[i][j][k][U1] *= tscale;
    P[i][j][k][U2] *= tscale;
    P[i][j][k][U3] *= tscale;
    P[i][j][k][B1] *= tscale;
    P[i][j][k][B2] *= tscale;
    P[i][j][k][B3] *= tscale;
  }
}



/* Set primitive from Carpet 3d profile */
void set_prim(grid_prim_type P){
    
    double *xp, *yp, *zp;
    double *rho, *eps, *ye, *velx, *vely, *velz;
    int np = (N1 + 2. * NG) * (N2 + 2. * NG) * (N3 + 2. * NG);
    
    xp = malloc(np * sizeof(*xp));
    yp = malloc(np * sizeof(*yp));
    zp = malloc(np * sizeof(*zp));
    
    rho = malloc(np * sizeof(*rho));
    eps = malloc(np * sizeof(*eps));
    ye = malloc(np * sizeof(*ye));
    velx = malloc(np * sizeof(*velx));
    vely = malloc(np * sizeof(*vely));
    velz = malloc(np * sizeof(*velz));
    
    int iflat = 0;
    double X[NDIM];
    ZLOOP {
      coord(i, j, k, CENT, X);
        xp[iflat] = X[1];
        yp[iflat] = X[2];
        zp[iflat] = X[3];
        iflat++;
      }

    /* Read metadata */
    cprof3d_file_t * dfile = cprof3d_open_file("245414.h5"); //TODO:pass through runtime arguments
    
    /* Open dataset: rho, eps, ye, velx, vely, velz */
    cprof3d_dset_t * dset_rho = cprof3d_read_dset(dfile, "rho");
    cprof3d_dset_t * dset_eps = cprof3d_read_dset(dfile, "eps");
    cprof3d_dset_t * dset_ye = cprof3d_read_dset(dfile, "Ye");
    cprof3d_dset_t * dset_velx = cprof3d_read_dset(dfile, "velx");
    cprof3d_dset_t * dset_vely = cprof3d_read_dset(dfile, "vely");
    cprof3d_dset_t * dset_velz = cprof3d_read_dset(dfile, "velz");
    

    /* Interpolate on the grid*/
    bool set_all_points_rho = cprof3d_interp(dset_rho,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            rho);
    bool set_all_points_eps = cprof3d_interp(dset_eps,
            cprof3d_cmap_reflecting_xy,
            cprof3d_transf_default,
            np, xp, yp, zp,
            eps);
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

    
    /* Check if interpolation of metric gone wrong*/
    if(!set_all_points_rho) {
        fprintf(stderr, "Something went wrong with the interpolation rho!");
        return;
    }
    if(!set_all_points_eps) {
        fprintf(stderr, "Something went wrong with the interpolation eps!");
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
    
    iflat = 0;
    ZLOOP{
            P[i][j][k][RHO] = rho[iflat];
            P[i][j][k][UU] = eps[iflat];
        
        
            #if EOS == EOS_TYPE_TABLE
            {
                P[i][j][k][YE] = ye[iflat];
                PASSTYPE(YE) = PASSTYPE_NUMBER;
                P[i][j][k][YE_EM] = ye[iflat];

                //extra[EOS_YE] = ye[iflat];

                #if NVAR_PASSIVE >= 4
                P[i][j][k][PASSIVE_START+3] = ye[iflat] * rho[iflat];
                PASSTYPE(PASSIVE_START+3) = PASSTYPE_INTRINSIC;
                #endif
            }
            #endif
        
            //P[i][j][k][YE] = ye[iflat];
            
            P[i][j][k][U1] = velx[iflat];
            P[i][j][k][U2] = vely[iflat];
            P[i][j][k][U3] = velz[iflat];
            
            P[i][j][k][B1] = 0.;
            P[i][j][k][B2] = 0.;
            P[i][j][k][B3] = 0.;
            iflat++;
    } // ZLOOP end
    
    /* Free memory */
    cprof3d_del_dset(dset_rho);
    cprof3d_del_dset(dset_eps);
    cprof3d_del_dset(dset_ye);
    cprof3d_del_dset(dset_velx);
    cprof3d_del_dset(dset_vely);
    cprof3d_del_dset(dset_velz);
    
    cprof3d_close_file(dfile);
    
    free(rho);
    free(eps);
    free(ye);
    free(velx);
    free(vely);
    free(velz);
    
    free(zp);
    free(yp);
    free(xp);
} // set_prim end

