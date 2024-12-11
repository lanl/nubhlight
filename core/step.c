/******************************************************************************
 *                                                                            *
 * STEP.C                                                                     *
 *                                                                            *
 * ADVANCES FLUID QUANTITIES BY ONE TIMESTEP                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt,
    grid_prim_type Pf, int stage);
double fluxcalc(grid_prim_type Pr);
void   flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
void   lr_to_flux(double p_l[NVAR], double p_r[NVAR], struct of_geom *geom,
      int dir, double Flux[NVAR], double *maxSpeed);
#if RADIATION
void apply_rad_force(grid_prim_type Pr, double Dt);
#endif

// DEBUG
#if RECORD_DT_MIN
double dt_min = INFINITY;
#endif

void step() {
  double ndt;
#if RADIATION
  double dt_cool;
#if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS && \
    RAD_NUM_TYPES >= 4 && NEUTRINO_OSCILLATIONS
  // not used for timestep control. Used to turn oscillations on or
  // off.
  double dt_osc = get_dt_oscillations();
#endif // oscillations
#endif // radiation
  dtsave = dt;

  // Need both P_n and P_n+1 to calculate current
  ZSLOOP(-NG, N1 - 1 + NG, -NG, N2 - 1 + NG, -NG, N3 - 1 + NG) {
    PLOOP { Psave[i][j][k][ip] = P[i][j][k][ip]; }
  }

#if RADIATION && ESTIMATE_THETAE
  estimate_Thetae(P, extra, t, dt);
#endif

  // Predictor step
  ndt = advance(P, P, 0.5 * dt, Ph, 0);
#if ELECTRONS
  heat_electrons(P, P, Ph, 0.5 * dt);
#endif
  fixup(Ph, extra);
  fixup_utoprim(Ph, extra);
#if ELECTRONS
  fixup_electrons(Ph);
#if RADIATION
  coulomb(P, P, Ph, 0.5 * dt);
  fixup_electrons(Ph);
#endif
#endif
  bound_prim(Ph);

// Radiation step
#if RADIATION
  precompute_microphysics();
  make_superphotons(Ph, extra, t, dt);
  // check_nu_type("after make"); // DEBUG
  push_superphotons(P, Ph, dt);
  // check_nu_type("after push"); // DEBUG
  interact(Ph, extra, t, dt);
  // check_nu_type("after interact"); // DEBUG
  bound_superphotons(Ph, t, dt);
  // check_nu_type("after bound"); // DEBUG
#if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS && \
    RAD_NUM_TYPES >= 4 && NEUTRINO_OSCILLATIONS
  int oscillations_active = (dt_osc <= dt);
  if (mpi_io_proc()) {
    printf("\t[Oscillations] Active? %d dt_osc = %.14e\n", oscillations_active,
        dt_osc);
  }
  if (oscillations_active) { // TOOD(JMM): Some safety factor?
    accumulate_local_angles();
    oscillate(local_moments, Gnu);
  }
  // check_nu_type("after oscillate"); // DEBUG
#endif // OSCILLATIONS
#endif

  // Corrector step
  ndt = advance(P, Ph, dt, P, 1);
#if ELECTRONS
  heat_electrons(P, Ph, P, dt);
#endif
  fixup(P, extra);
  fixup_utoprim(P, extra);
#if ELECTRONS
  fixup_electrons(P);
#if RADIATION
  coulomb(P, Ph, P, dt);
  fixup_electrons(P);
#endif
#endif
  bound_prim(P);

#if RADIATION
  // Apply radiation four-force to fluid
  apply_rad_force(P, dt);
  fixup(P, extra);
  fixup_utoprim(P, extra);
#if ELECTRONS
  apply_rad_force_e(Ph, P, radG, dt);
  fixup_electrons(P);
#endif
  bound_prim(P);

  // Reset radG
  memset((void *)&radG[0][0][0][0], 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * (NDIM + NRADCOMP) *
          sizeof(double));

#if TRACERS
  prune_tracers();
#endif // TRACERS
#endif // RADIATION

  // Reset fixup mask
  memset((void *)&fixup_required[0][0][0], 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * sizeof(int));

  // Increment time
  t += dt;

#if RADIATION
  dt_cool = get_min_dt_cool(P, extra);
  ndt     = MY_MIN(cour * dt_light_min, cour_cool * dt_cool);
#if RECORD_DT_MIN
  if (ndt < dt_min)
    dt_min = ndt;
  // DEBUG
  printf("dt_cool       = %g\n"
         "dt_light_min  = %g\n"
         "dt_global_min = %g\n",
      dt_cool, dt_light_min, dt_min);
#endif
#endif

  // Set next timestep
  if (ndt > SAFE * dt) {
    ndt = SAFE * dt;
  }

  // ndt = 0.01; // DEBUG
  dt = ndt;

  // dt = ndt;
  dt = mpi_min(dt);

// if dt is too small, try to advance anyway. Also complain.
#if RADIATION
  if (dt < SMALL) { // || dt_cool/dt_light_min < 1e-2) {
    fprintf(stderr,
        "ERROR: dt too smalL! Trying to advance anyway.\n"
        "\tdt     = %g\n"
        "\ttrying = %g\n",
        dt, SMALL + dtsave / SAFE);
    dt = SMALL + dtsave / SAFE;
  }
#endif

  // Don't step beyond end of run
  if (t + dt > tf) {
    dt = tf - t;
  }
}

double advance(grid_prim_type Pi, grid_prim_type Pb, double Dt,
    grid_prim_type Pf, int stage) {
  double          ndt, U[NVAR], dU[NVAR];
  struct of_state qi;

#pragma omp parallel for collapse(3)
  ZLOOP PLOOP Pf[i][j][k][ip] = Pi[i][j][k][ip];

  timer_start(TIMER_FLUXCALC);
  ndt = fluxcalc(Pb);

#if METRIC == MKS
  fix_flux(F1, F2, F3);
#endif

  flux_ct(F1, F2, F3);
  timer_stop(TIMER_FLUXCALC);

  // Evaluate diagnostics based on fluxes
  timer_start(TIMER_DIAG);
  diag_flux(F1, F2, F3);
  timer_stop(TIMER_DIAG);

#if NO_GRMHD_UPDATE
#pragma omp parallel for collapse(3)
  ZLOOP pflag[i][j][k] = 0;
  return ndt;
#endif

  // Update Pi to Pf
  timer_start(TIMER_UPDATE);
#pragma omp parallel for schedule(guided) private(dU, qi, U) collapse(3)
  ZLOOP {
    source(Pb[i][j][k], &(ggeom[i][j][CENT]), i, j, dU, Dt, extra[i][j][k]);
    get_state(Pi[i][j][k], &(ggeom[i][j][CENT]), &qi);
    primtoflux(Pi[i][j][k], &qi, 0, 0, &(ggeom[i][j][CENT]), U);

    PLOOP {
      U[ip] +=
          Dt * (-(F1[i + 1][j][k][ip] - F1[i][j][k][ip]) / dx[1] -
                   (F2[i][j + 1][k][ip] - F2[i][j][k][ip]) / dx[2] -
                   (F3[i][j][k + 1][ip] - F3[i][j][k][ip]) / dx[3] + dU[ip]);
    }

    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pf[i][j][k]);
    if (pflag[i][j][k])
      fail_save[i][j][k] = 1;
  } // ZLOOP
  timer_stop(TIMER_UPDATE);

  return (ndt);
}

#if RADIATION
void apply_rad_force(grid_prim_type Pr, double Dt) {
  double          U[NVAR];
  struct of_state q;

  timer_start(TIMER_UPDATE);

#pragma omp parallel for schedule(guided) private(q, U) collapse(3)
  ZLOOP {
    // DEBUG
    /*
    double zonevol = dV*L_unit*L_unit*L_unit*ggeom[i][j][CENT].g;
    printf("Dt = %g\n",Dt);
    printf("dU/dt = %g\n",Dt*radG[i][j][k][0]*U_unit);
    printf("dE/dt = %g\n",zonevol*Dt*radG[i][j][k][0]*U_unit);
    */
    // Store primitive variables before cooling for supercooling diagnostics
    PLOOP psupersave[i][j][k][ip] = Pr[i][j][k][ip];

    get_state(Pr[i][j][k], &(ggeom[i][j][CENT]), &q);
    primtoflux(Pr[i][j][k], &q, 0, 0, &(ggeom[i][j][CENT]), U);

    for (int ip = 1; ip < 5; ip++) {
      U[ip] += Dt * radG[i][j][k][ip - 1];
    }

    DLOOP1 { radG_int[i][j][k][mu] += Dt * radG[i][j][k][mu]; }

// TODO generalize this if need be
#if RADIATION == RADTYPE_NEUTRINOS
    {
      U[YE] += Dt * radG[i][j][k][RADG_YE];
      U[YE_EM] += Dt * radG[i][j][k][RADG_YE_EM];
      radG_int[i][j][k][RADG_YE] += Dt * radG[i][j][k][RADG_YE];
      radG_int[i][j][k][RADG_YE_EM] += Dt * radG[i][j][k][RADG_YE_EM];
    }
#endif

    pflag[i][j][k] = Utoprim(U, &(ggeom[i][j][CENT]), Pr[i][j][k]);
    /*if (pflag[i][j][k] == 5) {
      Pr[i][j][k][UU] = 0.;
      pflag[i][j][k] = 0;
      fixup1zone(i, j, k, Pr[i][j][k]);
    }*/

    if (pflag[i][j][k]) {
      fail_save[i][j][k] = 1;
    }
  } // ZLOOP
  timer_stop(TIMER_UPDATE);
}
#endif

double fluxcalc(grid_prim_type Pr) {
  double P_l[NMAX + 2 * NG][NVAR], P_r[NMAX + 2 * NG][NVAR];
  double cij, cmax1, cmax2, cmax3;
  double Ptmp[NMAX + 2 * NG][NVAR];

  cmax1 = cmax2 = cmax3 = 0.;
#pragma omp parallel private(Ptmp, P_l, P_r, cij) reduction(max      \
                                                            : cmax1) \
    reduction(max                                                    \
              : cmax2) reduction(max                                 \
                                 : cmax3)
  {
#pragma omp for collapse(2) nowait
    JSLOOP(-1, N2) {
      KSLOOP(-1, N3) {

        ISLOOP(-NG, N1 - 1 + NG) PLOOP Ptmp[i][ip] = Pr[i][j][k][ip];

        reconstruct(Ptmp, N1, P_l, P_r);

        ISLOOP(0, N1) {
          lr_to_flux(
              P_r[i - 1], P_l[i], &(ggeom[i][j][FACE1]), 1, F1[i][j][k], &cij);
          cmax1 = (cij > cmax1 ? cij : cmax1);
        } // ISLOOP
      }   // KSLOOP
    }     // JSLOOP

#pragma omp for collapse(2) nowait
    ISLOOP(-1, N1) {
      KSLOOP(-1, N3) {

        JSLOOP(-NG, N2 - 1 + NG) PLOOP Ptmp[j][ip] = Pr[i][j][k][ip];

        reconstruct(Ptmp, N2, P_l, P_r);

        JSLOOP(0, N2) {
          lr_to_flux(
              P_r[j - 1], P_l[j], &(ggeom[i][j][FACE2]), 2, F2[i][j][k], &cij);
          cmax2 = (cij > cmax2 ? cij : cmax2);
        } // JSLOOP
      }   // KSLOOP
    }     // ISLOOP

#pragma omp for collapse(2)
    ISLOOP(-1, N1) {
      JSLOOP(-1, N2) {

        KSLOOP(-NG, N3 - 1 + NG) PLOOP Ptmp[k][ip] = Pr[i][j][k][ip];

        reconstruct(Ptmp, N3, P_l, P_r);

        KSLOOP(0, N3) {
          lr_to_flux(
              P_r[k - 1], P_l[k], &(ggeom[i][j][FACE3]), 3, F3[i][j][k], &cij);
          cmax3 = (cij > cmax3 ? cij : cmax3);
        } // KSLOOP
      }   // JSLOOP
    }     // ISLOOP
  }       // omp parallel

  // Otherwise timestep changes with MPI!
  cmax1 = mpi_max(cmax1);
  cmax2 = mpi_max(cmax2);
  cmax3 = mpi_max(cmax3);

  double ndt1 = cour * dx[1] / cmax1;
  double ndt2 = cour * dx[2] / cmax2;
  double ndt3 = cour * dx[3] / cmax3;

  return (1. / (1. / ndt1 + 1. / ndt2 + 1. / ndt3));
}

void flux_ct(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3) {
  static double emf1[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
  static double emf2[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
  static double emf3[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];

// Calculate EMFs via average to corners (Toth approach)
#pragma omp parallel
  {
#pragma omp for collapse(3)
    ZSLOOP(0, N1, 0, N2, 0, N3) {
      emf3[i][j][k] = 0.25 * (F1[i][j][k][B2] + F1[i][j - 1][k][B2] -
                                 F2[i][j][k][B1] - F2[i - 1][j][k][B1]);
      emf2[i][j][k] = -0.25 * (F1[i][j][k][B3] + F1[i][j][k - 1][B3] -
                                  F3[i][j][k][B1] - F3[i - 1][j][k][B1]);
      emf1[i][j][k] = 0.25 * (F2[i][j][k][B3] + F2[i][j][k - 1][B3] -
                                 F3[i][j][k][B2] - F3[i][j - 1][k][B2]);
    }

// Rewrite EMFs as fluxes, after Toth
#pragma omp for collapse(3) nowait
    ZSLOOP(0, N1, 0, N2 - 1, 0, N3 - 1) {
      F1[i][j][k][B1] = 0.;
      F1[i][j][k][B2] = 0.5 * (emf3[i][j][k] + emf3[i][j + 1][k]);
      F1[i][j][k][B3] = -0.5 * (emf2[i][j][k] + emf2[i][j][k + 1]);
    }
#pragma omp for collapse(3) nowait
    ZSLOOP(0, N1 - 1, 0, N2, 0, N3 - 1) {
      F2[i][j][k][B1] = -0.5 * (emf3[i][j][k] + emf3[i + 1][j][k]);
      F2[i][j][k][B2] = 0.;
      F2[i][j][k][B3] = 0.5 * (emf1[i][j][k] + emf1[i][j][k + 1]);
    }
#pragma omp for collapse(3)
    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3) {
      F3[i][j][k][B1] = 0.5 * (emf2[i][j][k] + emf2[i + 1][j][k]);
      F3[i][j][k][B2] = -0.5 * (emf1[i][j][k] + emf1[i][j + 1][k]);
      F3[i][j][k][B3] = 0.;
    }
  } // omp parallel
}

void lr_to_flux(double P_l[NVAR], double P_r[NVAR], struct of_geom *geom,
    int dir, double Flux[NVAR], double *maxSpeed) {
  struct of_state state_l, state_r;
  double          F_l[NVAR], F_r[NVAR], U_l[NVAR], U_r[NVAR];
  double          cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, ctop;

  if (geom->g < SMALL) {
    PLOOP Flux[ip] = 0.;
    *maxSpeed      = 0.;
    return;
  }

  get_state(P_l, geom, &state_l);
  get_state(P_r, geom, &state_r);

  primtoflux(P_l, &state_l, dir, 0, geom, F_l);
  primtoflux(P_r, &state_r, dir, 0, geom, F_r);

  primtoflux(P_l, &state_l, 0, 0, geom, U_l);
  primtoflux(P_r, &state_r, 0, 0, geom, U_r);

  mhd_vchar(P_l, &state_l, geom, dir, &cmax_l, &cmin_l);
  mhd_vchar(P_r, &state_r, geom, dir, &cmax_r, &cmin_r);

  cmax = fabs(MY_MAX(MY_MAX(0., cmax_l), cmax_r));
  cmin = fabs(MY_MAX(MY_MAX(0., -cmin_l), -cmin_r));
  ctop = MY_MAX(cmax, cmin);

  PLOOP Flux[ip] = 0.5 * (F_l[ip] + F_r[ip] - ctop * (U_r[ip] - U_l[ip]));

  *maxSpeed = ctop;
}
