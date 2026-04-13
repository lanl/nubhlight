/******************************************************************************
 *                                                                            *
 * DIAG.C                                                                     *
 *                                                                            *
 * DIAGNOSTIC OUTPUT                                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// static double lnumin, lnumax, dlnu;

void reset_log_variables() {
#if RADIATION
  step_tot = step_lost = step_made = step_abs = step_scatt = step_rec = 0;
  step_sent = step_rcvd = step_fail = 0;
  tracer_tot                        = 0;
#if RADIATION == RADTYPE_NEUTRINOS
  TYPELOOP rad_type_counts[itp] = 0.0;
#endif
#endif
}

void reset_dump_variables() {
#if RADIATION
  memset(nuLnu, 0, (NULNU_IDX0)*NTH * NPHI * NU_BINS_SPEC * sizeof(double));
  memset(Nem, 0, N123G * sizeof(int));
  memset(Nabs, 0, N123G * sizeof(int));
  memset(Nsc, 0, N123G * sizeof(int));
  memset(Nem_phys, 0, RAD_NUM_TYPES * N123G * sizeof(double));
  memset(Nabs_phys, 0, RAD_NUM_TYPES * N123G * sizeof(double));
  memset(Jrad, 0, (MAXNSCATT + 2) * N123G * sizeof(double));
  memset(Esuper, 0, N123G * sizeof(double));
  memset(Nsuper, 0, N123G * sizeof(int));
  memset(radG_int, 0, RAD_NUM_TYPES * N123G * sizeof(double));
  memset(dtau_avg, 0, (RAD_SCATT_TYPES + 1) * N123G * sizeof(double));
  memset(en_int_avg, 0, (RAD_SCATT_TYPES + 1) * N123G * sizeof(double));
#if RZ_HISTOGRAMS
  memset(rz_r_orig_hist, 0, RZ_HISTOGRAMS_N * sizeof(double));
  memset(rz_z_orig_hist, 0, RZ_HISTOGRAMS_N * sizeof(double));
#if NEUTRINO_OSCILLATIONS
  memset(osc_rz_r_orig_hist, 0, RZ_HISTOGRAMS_N * sizeof(double));
  memset(osc_rz_z_orig_hist, 0, RZ_HISTOGRAMS_N * sizeof(double));
#endif // NEUTRINO_OSCILLATIONS
#endif // RZ_HISTOGRAMS
#if NEUTRINO_OSCILLATIONS
  memset(local_osc_count, 0, LOCAL_ANGLES_NX1*LOCAL_ANGLES_NX2*sizeof(double));
#endif // NEUTRINO_OSCILLATIONS
#endif // RADIATION
}

void diag(int call_code) {
  double          U[NVAR], divb;
  double          pp      = 0.;
  double          divbmax = 0.;
  double          rmed    = 0.;
  double          e       = 0.;
  struct of_geom *geom;
  struct of_state q;
  static FILE *   ener_file;

  // Write diagnostics to dump directory
  char log_fnam[STRLEN];
  strcpy(log_fnam, dumpdir);
  strcat(log_fnam, "diag.out");

  if (call_code == DIAG_INIT) {
    // Set things up
    if (mpi_io_proc()) {
      if (is_restart) {
        ener_file = fopen(log_fnam, "a");
      } else {
        ener_file = fopen(log_fnam, "w");
      }
      if (ener_file == NULL) {
        fprintf(stderr, "error opening energy output file\n");
        exit(1);
      }
    }

    // dOmega = sin(theta) dtheta = -d cos(theta)
    /*#if RADIATION
    double dOtot = 0.;
    JSLOOP(0, N2-1) {
      KSLOOP(0, N3-1) {
        int i = N1+NG-1;
        double XL[NDIM], XR[NDIM], thL, thR;
        coord(i, j, k, FACE2, XL);
        coord(i, j+1, k, FACE2, XR);
        thL = th_of_X(XL);
        thR = th_of_X(XR);
        dOmega[j][k] = 2.*M_PI/(N3TOT)*(-cos(thR) + cos(thL));
        dOtot += dOmega[j][k];
      }
    }
    #endif // RADIATION*/
  }

  // Calculate conserved quantities
  if ((call_code == DIAG_INIT || call_code == DIAG_LOG ||
          call_code == DIAG_FINAL) &&
      !failed) {
    pp      = 0.;
    e       = 0.;
    rmed    = 0.;
    divbmax = 0.;
    ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3 - 1) {
      geom = get_geometry(i, j, k, CENT);
      get_state(P[i][j][k], geom, &q);
      primtoflux(P[i][j][k], &q, 0, 0, geom, U);

      rmed += U[RHO] * dV;
      pp += U[U3] * dV;
      e += U[UU] * dV;

      divb = flux_ct_divb(i, j, k);

      if (divb > divbmax) {
        divbmax = divb;
      }
    }
  }

  rmed    = mpi_io_reduce(rmed);
  pp      = mpi_io_reduce(pp);
  e       = mpi_io_reduce(e);
  divbmax = mpi_io_max(divbmax);

#if RADIATION
  set_Rmunu();
#endif

  // Get total mass and energy
  double mass_proc = 0.;
  double egas_proc = 0.;
#if RADIATION
  double erad_proc = 0.;
#endif
  double Phi_proc         = 0.;
  double jet_EM_flux_proc = 0.;
  double lum_eht_proc     = 0.;
  ZLOOP {
#if EOS == EOS_TYPE_TABLE
    EOS_SC_fill(P[i][j][k], extra[i][j][k]);
#endif

    struct of_state q;
    double          U[NVAR];
    get_state(P[i][j][k], &(ggeom[i][j][CENT]), &q);
    primtoflux(P[i][j][k], &q, 0, 0, &(ggeom[i][j][CENT]), U);
    mass_proc += U[0] * dx[1] * dx[2] * dx[3];
    egas_proc += U[1] * dx[1] * dx[2] * dx[3];
    double rho   = P[i][j][k][RHO];
    double ug    = P[i][j][k][UU];
    double Pg    = EOS_pressure_rho0_u(rho, ug, extra[i][j][k]);
    double bsq   = dot(q.bcon, q.bcov);
    double Bmag  = sqrt(bsq);
    double C_eht = 0.2;
    double j_eht = pow(rho, 3.) * pow(Pg, -2.) *
                   exp(-C_eht * pow(rho * rho / (Bmag * Pg * Pg), 1. / 3.));
    lum_eht_proc += j_eht * dx[1] * dx[2] * dx[3] * ggeom[i][j][CENT].g;
    if (global_start[1] == 0 && i == 5 + NG) {
      Phi_proc +=
          0.5 * fabs(P[i][j][k][B1]) * dx[2] * dx[3] * ggeom[i][j][CENT].g;

      double P_EM[NVAR];
      PLOOP  P_EM[ip] = P[i][j][k][ip];
      P_EM[RHO]       = 0.;
      P_EM[UU]        = 0.;
      get_state(P_EM, &(ggeom[i][j][CENT]), &q);
      double sig = dot(q.bcon, q.bcov) / P[i][j][k][RHO];
      if (sig > 1.) {
        primtoflux(P_EM, &q, 1, 1, &(ggeom[i][j][CENT]), U);
        jet_EM_flux_proc += -U[U1] * dx[2] * dx[3];
      }
    }

#if RADIATION
    erad_proc +=
        Rmunu[i][j][k][0][0] * dx[1] * dx[2] * dx[3] * ggeom[i][j][CENT].g;
#endif
  }
  double mdot_all    = mpi_io_reduce(mdot);
  double edot_all    = mpi_io_reduce(edot);
  double ldot_all    = mpi_io_reduce(ldot);
  double mdot_eh_all = mpi_io_reduce(mdot_eh);
  double edot_eh_all = mpi_io_reduce(edot_eh);
  double ldot_eh_all = mpi_io_reduce(ldot_eh);
  double mass        = mpi_reduce(mass_proc);
  double egas        = mpi_reduce(egas_proc);
  double Phi         = mpi_reduce(Phi_proc);
  double phi         = Phi / sqrt(mdot_all + SMALL);
  double jet_EM_flux = mpi_reduce(jet_EM_flux_proc);
  double lum_eht     = mpi_reduce(lum_eht_proc);
#if RADIATION
  double erad     = mpi_reduce(erad_proc);
  double lum_proc = 0.;

  // Get last shell entirely enclosed by stopx_rad[1] (assume r(X^1) independent
  // of X^2, X^3)
  int stopi_rad = -1;
  for (int i = NG + 1; i <= N1 + NG; i++) {
    double X[NDIM], Xprev[NDIM];
    coord(i - 1, NG, NG, FACE1, Xprev);
    coord(i, NG, NG, FACE1, X);
    if (X[1] > stopx_rad[1] && Xprev[1] < stopx_rad[1]) {
      stopi_rad = i - 2;
    }
  }
  if (stopi_rad == -1)
    stopi_rad = N1 + NG - 1;

  if (stopi_rad >= 0) {
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        lum_proc -= Rmunu[stopi_rad][j][k][1][0] * dx[2] * dx[3] *
                    ggeom[stopi_rad][j][CENT].g;
      }
    }
  }
  double lum = mpi_reduce(lum_proc);
  double eff = lum / (mdot + SMALL);

#if ELECTRONS
  int    num_super = 0.;
  double lum_super = 0.;
  ZLOOP {
    num_super += Nsuper[i][j][k];
    lum_super += Esuper[i][j][k];
  }
  num_super = mpi_reduce_int(num_super);
  lum_super = mpi_reduce(lum_super) / DTd;
#endif
#endif // RADIATION

  if (call_code == DIAG_DUMP) {
    dump();
#if RADIATION && TRACERS
    dump_tracers();
#endif
    dump_cnt++;
    reset_dump_variables();
  }

  if (call_code == DIAG_FINAL) {
    dump();
#if RADIATION && TRACERS
    dump_tracers();
#endif
    reset_dump_variables();
  }

  if (call_code == DIAG_INIT || call_code == DIAG_LOG ||
      call_code == DIAG_FINAL) {

    if (mpi_io_proc()) {
#if EOS == EOS_TYPE_TABLE
      EOS_SC_fill(P[N1 / 2][N2 / 2][N3 / 2], extra[N1 / 2][N2 / 2][N3 / 2]);
#endif
      fprintf(stdout, "LOG      t=%g \t divbmax: %g\n", t, divbmax);
#if RADIATION
#if RADIATION == RADTYPE_NEUTRINOS
      fprintf(stdout, "         nleptons = %g \t %%change = %g\n", lepton_tot,
          dlepton_perc);
#endif
#endif
      fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ", t, rmed,
          pp, e,
          EOS_adiabatic_constant(P[N1 / 2][N2 / 2][N3 / 2][RHO],
              P[N1 / 2][N2 / 2][N3 / 2][UU], extra[N1 / 2][N2 / 2][N3 / 2]),
          P[N1 / 2][N2 / 2][N3 / 2][UU]);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot_all, edot_all, ldot_all);
      fprintf(ener_file, "%15.8g %15.8g ", mass, egas);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", Phi, phi, jet_EM_flux);
      fprintf(ener_file, "%15.8g ", divbmax);
#if RADIATION
      printf("step_abs = %i step_abs_tot = %i\n", step_abs, step_abs_all);
      fprintf(ener_file, "%i %i %i %i %i %i %i %i ", step_made, step_abs,
          step_scatt, step_lost, step_rec, step_tot, step_sent, step_rcvd);
      fprintf(ener_file, "%i %i %i %i %i %i %i %i ", step_made_all,
          step_abs_all, step_scatt_all, step_lost_all, step_rec_all,
          step_tot_all, step_sent_all, step_rcvd_all);
      fprintf(ener_file, "%15.8g ", load_imbalance);
      fprintf(ener_file, "%15.8g %15.8g ", tune_emiss, tune_scatt);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", erad, lum, eff);
#if ELECTRONS
      fprintf(ener_file, "%i %15.8g ", num_super, lum_super);
#endif
#if RADIATION == RADTYPE_NEUTRINOS
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", lepton_tot, dlepton_tot,
          dlepton_perc);
#endif
#if TRACERS
      fprintf(ener_file, "%i ", tracer_tot_all);
      fprintf(ener_file, "%i ", step_tot_all - tracer_tot_all);
#endif
#endif
      fprintf(ener_file, "%15.8g ", lum_eht);
      fprintf(ener_file, "%15.8g %15.8g %15.8g ", mdot_eh_all, edot_eh_all,
          ldot_eh_all);
      fprintf(ener_file,
          "%15.8g %15.8g %15.8g %15.8g "
          "%15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g ",
          get_time_per_step(TIMER_UPDATE), get_time_per_step(TIMER_FLUXCALC),
          get_time_per_step(TIMER_FIXUP), get_time_per_step(TIMER_BOUND),
          get_time_per_step(TIMER_DIAG), get_time_per_step(TIMER_OUT),
          get_time_per_step(TIMER_MAKE), get_time_per_step(TIMER_PUSH),
          get_time_per_step(TIMER_INTERACT), get_time_per_step(TIMER_MICRO),
          get_time_per_step(TIMER_ALL));
#if ELECTRONS
      fprintf(ener_file, "%15.8g ", get_time_per_step(TIMER_ELECTRON));
#endif
#if NEUTRINO_OSCILLATIONS || LOCAL_ANGULAR_DISTRIBUTIONS
      fprintf(ener_file, "%15.8g ", get_time_per_step(TIMER_OSCILLATIONS));
#endif
      fprintf(ener_file, "\n");
      fflush(ener_file);
    }

#if RADIATION
    track_ph();
#endif
  }
}

// Diagnostic routines
void fail(int fail_type) {
  failed = 1;

  fprintf(stderr, "\n\nFAIL: [%d %d %d] %d\n", icurr, jcurr, kcurr, fail_type);

  area_map(icurr, jcurr, kcurr, P);

  diag(DIAG_FINAL);

  exit(0);
}

// Map out region around failure point
void area_map(int i, int j, int k, grid_prim_type prim) {
  fprintf(stderr, "*** AREA MAP ***\n");

  PLOOP {
    fprintf(stderr, "variable %d \n", ip);
    fprintf(stderr, "i = \t %12d %12d %12d\n", i - 1, i, i + 1);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j + 1,
        prim[i - 1][j + 1][k][ip], prim[i][j + 1][k][ip],
        prim[i + 1][j + 1][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j,
        prim[i - 1][j][k][ip], prim[i][j][k][ip], prim[i + 1][j][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j - 1,
        prim[i - 1][j - 1][k][ip], prim[i][j - 1][k][ip],
        prim[i + 1][j - 1][k][ip]);
  }

  fprintf(stderr, "****************\n");
}

// Evaluate flux based diagnostics; put results in global variables
void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3) {
  mdot = edot = ldot = 0.;
  mdot_eh = edot_eh = ldot_eh = 0.;
  int iEH                     = NG + 5;
  if (global_start[1] == 0) {
#pragma omp parallel for \
      reduction(+:mdot) reduction(+:edot) reduction(+:ldot) \
      reduction(+:mdot_eh) reduction(+:edot_eh) reduction(+:ldot_eh) \
      collapse(2)
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        mdot += F1[NG][j][k][RHO] * dx[2] * dx[3];
        edot += (F1[NG][j][k][UU] - F1[NG][j][k][RHO]) * dx[2] * dx[3];
        ldot += F1[NG][j][k][U3] * dx[2] * dx[3];
        mdot_eh += F1[iEH][j][k][RHO] * dx[2] * dx[3];
        edot_eh += (F1[iEH][j][k][UU] - F1[iEH][j][k][RHO]) * dx[2] * dx[3];
        ldot_eh += F1[iEH][j][k][U3] * dx[2] * dx[3];
      }
    }
  }
}

double flux_ct_divb(int i, int j, int k) {
  return fabs(
      0.25 *
          (P[i][j][k][B1] * ggeom[i][j][CENT].g +
              P[i][j - 1][k][B1] * ggeom[i][j - 1][CENT].g +
              P[i][j][k - 1][B1] * ggeom[i][j][CENT].g +
              P[i][j - 1][k - 1][B1] * ggeom[i][j - 1][CENT].g -
              P[i - 1][j][k][B1] * ggeom[i - 1][j][CENT].g -
              P[i - 1][j - 1][k][B1] * ggeom[i - 1][j - 1][CENT].g -
              P[i - 1][j][k - 1][B1] * ggeom[i - 1][j][CENT].g -
              P[i - 1][j - 1][k - 1][B1] * ggeom[i - 1][j - 1][CENT].g) /
          dx[1] +
      0.25 *
          (P[i][j][k][B2] * ggeom[i][j][CENT].g +
              P[i - 1][j][k][B2] * ggeom[i - 1][j][CENT].g +
              P[i][j][k - 1][B2] * ggeom[i][j][CENT].g +
              P[i - 1][j][k - 1][B2] * ggeom[i - 1][j][CENT].g -
              P[i][j - 1][k][B2] * ggeom[i][j - 1][CENT].g -
              P[i - 1][j - 1][k][B2] * ggeom[i - 1][j - 1][CENT].g -
              P[i][j - 1][k - 1][B2] * ggeom[i][j - 1][CENT].g -
              P[i - 1][j - 1][k - 1][B2] * ggeom[i - 1][j - 1][CENT].g) /
          dx[2] +
      0.25 *
          (P[i][j][k][B3] * ggeom[i][j][CENT].g +
              P[i][j - 1][k][B3] * ggeom[i][j - 1][CENT].g +
              P[i - 1][j][k][B3] * ggeom[i - 1][j][CENT].g +
              P[i - 1][j - 1][k][B3] * ggeom[i - 1][j - 1][CENT].g -
              P[i][j][k - 1][B3] * ggeom[i][j][CENT].g -
              P[i][j - 1][k - 1][B3] * ggeom[i][j - 1][CENT].g -
              P[i - 1][j][k - 1][B3] * ggeom[i - 1][j][CENT].g -
              P[i - 1][j - 1][k - 1][B3] * ggeom[i - 1][j - 1][CENT].g) /
          dx[3]);
}

#if RADIATION
void record_superphoton(double X[NDIM], struct of_photon *ph) {
  // Do not do this for tracers
  if (ph->type == TYPE_TRACER)
    return;

  double lnumin = log(numin);
  double lnumax = log(numax);
  double dlnu   = (lnumax - lnumin) / NU_BINS_SPEC;
  int    i, j, k;
  Xtoijk(X, &i, &j, &k);

  int nscatt = ph->nscatt;
  nscatt     = MY_MIN(nscatt, MAXNSCATT);

  // Preserve index sanity
  if (j < NG)
    j = N2 + j;
  if (j >= N2 + NG)
    j = j - N2;
  if (k < NG) {
    if (N3CPU == 1)
      k = N3 + k;
    else
      k = NG;
  }
  if (k >= N3 + NG) {
    if (N3CPU == 1)
      k = k - N3;
    else
      k = N3 + NG - 1;
  }

  // Assume X0 symmetry in metric
  double nu = -ph->Kcov[2][0] * ME * CL * CL / HPL;

  int thbin, phibin, nubin = (log(nu) - lnumin) / dlnu;
  get_nuLnu_bin(X, &thbin, &phibin);

  // Store dE / dlognu dOmega dt
  if (nubin >= 0 && nubin < NU_BINS_SPEC) {
#if DIAGNOSTICS_USE_RADTYPES
    {
#pragma omp atomic
      nuLnu[ph->type][thbin][phibin][nubin] -=
          ph->w * ph->Kcov[2][0] * ME * CL * CL / (dlnu * DTd * T_unit);
    }
#else
    {
#pragma omp atomic
      nuLnu[nscatt][thbin][phibin][nubin] -=
          ph->w * ph->Kcov[2][0] * ME * CL * CL / (dlnu * DTd * T_unit);
    }
#endif

#pragma omp atomic
    step_rec++;
  }
}

void report_load_imbalance() {
  if (mpi_io_proc()) {
    fprintf(stdout,
        "\n******** LOAD IMBALANCE *********\n"
        "   good == 0 <= %.2f <= 1 == bad\n"
        "   step_tot_min = %d\n"
        "   step_tot_max = %d\n"
        "*********************************\n\n",
        load_imbalance, step_tot_min, step_tot_max);
  }
}

void        bin_all_superphotons() {
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        bin_superphoton_direction(ph);
      }
      ph = ph->next;
    }
  }
}

#if RADIATION == RADTYPE_NEUTRINOS
void count_leptons(grid_prim_type P, double dt, int nstep) {
  timer_start(TIMER_DIAG);
  lepton_lost_step = mpi_reduce(lepton_lost_local);
  // printf("lepton lost on %d step = %g\n", mpi_myrank(), lepton_lost_local);
  // if (mpi_io_proc()) {
  //   printf("lepton lost this step = %g\n", lepton_lost_step);
  //   printf("lepton lost last = %g\n", lepton_lost);
  //   printf("lepton lost now = %g\n", lepton_lost + lepton_lost_step);
  // }
  lepton_last = lepton_tot;
  lepton_lost += lepton_lost_step;
  lepton_lost_local = 0.0;
  lepton_gas = lepton_rad = 0.0;
  double reference_vol    = dx[1] * dx[2] * dx[3] * pow(L_unit, 3.);
#pragma omp parallel for collapse(3)
  ZLOOP {
    double gamma, alpha, ucon0;
    mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
    alpha = ggeom[i][j][CENT].alpha;
    ucon0 = gamma / alpha;

    double nb_cell     = ucon0 * P[i][j][k][RHO] * RHO_unit / MP;
    double n_cell      = nb_cell * P[i][j][k][YE];
    double cell_vol    = reference_vol * ggeom[i][j][CENT].g;
    double lepton_cell = n_cell * cell_vol;
#pragma omp atomic
    lepton_gas += lepton_cell;
  }
  lepton_gas = mpi_reduce(lepton_gas);
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
#pragma omp atomic
        lepton_rad += (ph->w) * get_lepton_sign(ph);
      }
      ph = ph->next;
    }
  }
  lepton_rad   = mpi_reduce(lepton_rad);
  lepton_tot   = lepton_rad + lepton_gas + lepton_lost;
  dlepton_tot  = (lepton_tot - lepton_last) / dt;
  dlepton_perc = 100. * (lepton_tot - lepton_last) / lepton_tot;
// diagnostic currently only works when mpi is off. I don't know why.
#if COMPLAIN_ON_LEPTON_NONCON
  if (nstep > 1 && fabs(dlepton_perc) > DLEPTON_THRESH && mpi_nprocs() > 1) {
    if (mpi_io_proc()) {
      double dlep_abs = (lepton_tot - lepton_last);
      fprintf(stderr,
          "Lepton number not conserved!\n"
          "\tlepton_gas   = %15.8g\n"
          "\tlepton_rad   = %15.8g\n"
          "\tlepton_lost  = %15.8g\n"
          "\tlepton_tot   = %15.8g\n"
          "\tlepton_last  = %15.8g\n"
          "\tdlep_abs     = %15.8g\n"
          "\tdlepdt       = %15.8g\n"
          "\tdlepton_perc = %15.8g\n",
          lepton_gas, lepton_rad, lepton_lost, lepton_tot, lepton_last,
          dlep_abs, dlepton_tot, dlepton_perc);
    }
    // uncomment this when the code actually works
    exit(1);
  }
#endif
  timer_stop(TIMER_DIAG);
}

#if RZ_HISTOGRAMS
void generate_rz_histograms() {
  timer_start(TIMER_DIAG);

  #pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      double Xorig[NDIM], Xorig_cart[NDIM];
      int i = ph->origin[1];
      int j = ph->origin[2];
      int k = ph->origin[3];
      coord(i, j, k, CENT, Xorig);
      cart_coord(Xorig, Xorig_cart);

      double x = Xorig_cart[1];
      double y = Xorig_cart[2];
      double rcyl = MY_MIN(sqrtf(x*x + y*y), fabs(rz_rmax) - SMALL);
      double z = MY_MIN(fabs(Xorig_cart[3]), fabs(rz_zmax) - SMALL);

      int ir = rcyl / (delta_rcyl + SMALL);
      int iz = z / (delta_z + SMALL);

#pragma omp atomic
      rz_r_orig_hist[ir] += ph->w;
#pragma omp atomic
      rz_z_orig_hist[iz] += ph->w;

#if NEUTRINO_OSCILLATIONS
      if (ph->osc_count) {
#pragma omp atomic
        osc_rz_r_orig_hist[ir] += ph->w;
#pragma omp atomic
        osc_rz_z_orig_hist[iz] += ph->w;
      }
#endif // NEUTRINO_OSCILLATIONS

      ph = ph->next;
    }
  }
  mpi_dbl_allreduce_array((double *)rz_r_orig_hist, RZ_HISTOGRAMS_N);
  mpi_dbl_allreduce_array((double *)rz_z_orig_hist, RZ_HISTOGRAMS_N);
#if NEUTRINO_OSCILLATIONS
  mpi_dbl_allreduce_array((double *)osc_rz_r_orig_hist, RZ_HISTOGRAMS_N);
  mpi_dbl_allreduce_array((double *)osc_rz_z_orig_hist, RZ_HISTOGRAMS_N);
#endif

  timer_stop(TIMER_DIAG);
}
#endif // RZ_HISTOGRAMS

void print_rad_types() {
  timer_start(TIMER_DIAG);
  double   total_count          = 0;
  TYPELOOP rad_type_counts[itp] = 0.0;
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
#pragma omp atomic
        rad_type_counts[ph->type] += ph->w;
      }
      ph = ph->next;
    }
  }
  TYPELOOP rad_type_counts[itp] = mpi_reduce(rad_type_counts[itp]);
  TYPELOOP total_count += rad_type_counts[itp];
  TYPELOOP rad_type_counts[itp] /= total_count;
  TYPELOOP rad_type_counts[itp] *= 100.0;
  if (mpi_io_proc()) {
    fprintf(stdout, "\n********* NEUTRINO TYPES ********\n");
    fprintf(stdout, " TOTAL COUNT = %.4g \n", total_count);
    fprintf(stdout, "*********************************\n");
    fprintf(stdout, " TYPE                  PERCENTAGE\n");
    fprintf(stdout, "*********************************\n");
    fprintf(stdout, " ELECTRON                 %.2f %%\n",
        rad_type_counts[NU_ELECTRON]);
    fprintf(stdout, " ANTI                     %.2f %%\n",
        rad_type_counts[ANTINU_ELECTRON]);
    fprintf(stdout, " X                        %.2f %%\n",
        rad_type_counts[NU_HEAVY]);
#if RAD_NUM_TYPES >= 4
    fprintf(stdout, " ANTIX                    %.2f %%\n",
        rad_type_counts[ANTINU_HEAVY]);
#endif
    fprintf(stdout, "*********************************\n\n");
  }
  timer_stop(TIMER_DIAG);
}
#endif

#endif // RADIATION
