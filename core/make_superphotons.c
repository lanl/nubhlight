/******************************************************************************
 *                                                                            *
 * MAKE_SUPERPHOTONS.C                                                        *
 *                                                                            *
 * EMISSION OF MONTE CARLO SAMPLES                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
static double lnu_min, lnu_max, dlnu, nusamp[NU_BINS + 1], Ns;
static double dnzs[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][RAD_NUM_TYPES];
/*
static double dlepton; // DEBUG
static int navgs = 0;
static double lepfrac= 0;
*/

void sample_photon(int i, int j, int k, double t, double dt, int type,
    double dndlnu[NU_BINS + 1], struct of_photon *tmp, double Econ[NDIM][NDIM],
    double Ecov[NDIM][NDIM], const struct of_microphysics *m,
    double Bcon[NDIM]);
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS + 1],
    int type, const struct of_microphysics *m);

double get_wgt(double nu, double dtau) {
#if EXPTAU_WEIGHTS
  return wgtC / nu * exp(dtau);
#endif
  return wgtC / nu;
}

// Minkowski space (frame of nu) estimate for dtau. Max dtau = 100 to avoid
// numerical errors during exponentiation
double get_dtau(
    double nu, int type, double dt, const struct of_microphysics *m) {
#if EXPTAU_WEIGHTS
  {
    double dtau = alpha_inv_abs(nu, type, m, M_PI / 2.) * L_unit * dt / nu;
    return MY_MIN(dtau, 100.);
  }
#else
  { return 0.; }
#endif
}

void make_superphotons(
    grid_prim_type Prad, grid_eosvar_type extra, double t, double dt) {
#if EMISSION
  timer_start(TIMER_MAKE);
  get_dnz(Prad, extra);

  int step_made_local = 0;
  // dlepton = 0; // DEBUG

#pragma omp parallel reduction(+ : step_made_local)
  {
    struct of_photon *tmp, *head = photon_lists[omp_get_thread_num()];
    // struct of_microphysics m;
    double dndlnu[NU_BINS + 1];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    // double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double X[NDIM];
    int    nz;

    ZLOOP {
      TYPELOOP {
        nz = (int)dnzs[i][j][k][itp];
        if (dnzs[i][j][k][itp] - nz > get_rand())
          nz++;
        if (nz > 0) {
          // Set up zone
          coord(i, j, k, CENT, X);

          make_tetrad(i, j, k, Ucon_grd[i][j][k], Bcon_grd[i][j][k],
              ggeom[i][j][CENT].gcov, Econ, Ecov);
          get_dndlnu(i, j, k, dt, dndlnu, itp, &(m_grd[i][j][k]));

          // Create superphotons in pairs
          for (int n = 0; n < nz; n++) {
            tmp       = safe_malloc(sizeof(struct of_photon));
            tmp->next = safe_malloc(sizeof(struct of_photon));
            sample_photon(i, j, k, t, dt, itp, dndlnu, tmp, Econ, Ecov,
                &(m_grd[i][j][k]), Bcon_grd[i][j][k]);

#if KILL_ALL_PACKETS
            {
              free(tmp->next);
              free(tmp);
            }
#else
            {
              (tmp->next)->next = head;
              head              = tmp;
            }
#endif
          } // n < nz

          //#pragma omp atomic
          step_made_local += 2 * nz;
        } // nz > 0
      }   // TYPELOOP
    }     // ZLOOP

    // Prepend created superphotons to each thread's global list
    photon_lists[omp_get_thread_num()] = head;
  } // omp parallel

  step_made += step_made_local;
#if KILL_ALL_PACKETS
  { step_lost += step_made_local; }
#endif

  // DEBUG
  /*
  double C = 4*M_PI*U_unit*cnu_flat/((numax-numin)*T_unit);
  printf("dlepton = %g\n",dlepton);
  ZLOOP {
    double zoneVol = dV*pow(L_unit,3)*ggeom[i][j][CENT].g;
    double dlep_ana =
  -2*Prad[i][j][k][YE]*log(numax/numin)*C*(1/HPL)*dt*T_unit*zoneVol;
    printf("dlep_ana = %g\n",dlep_ana);
    printf("dlep_ana/dlep = %g\n", dlep_ana/(-dlepton));\
    lepfrac = (lepfrac*navgs + (dlep_ana/(-dlepton)))/(navgs+1);
    navgs++;
    printf("avg_lepfrac = %g\n", lepfrac);
    printf("Ye = %g\n",Prad[i][j][k][YE]);
  }
  */
  timer_stop(TIMER_MAKE);
#endif // EMISSION
}

void sample_photon(int i, int j, int k, double t, double dt, int type,
    double dndlnu[NU_BINS + 1], struct of_photon *ph, double Econ[NDIM][NDIM],
    double Ecov[NDIM][NDIM], const struct of_microphysics *m,
    double Bcon[NDIM]) {
  double            nu, th, cth[2], sth[2], phi, sphi[2], cphi[2];
  double            K_tetrad[NDIM];
  struct of_photon *tmp[2];
  tmp[0]       = ph;
  tmp[1]       = ph->next;
  tmp[1]->next = NULL;

  // Sample emissivity to get frequency
  do {
    nu = exp(get_rand() * (lnu_max - lnu_min) + lnu_min);
  } while (get_rand() > linear_interp_log(nu, dndlnu, lnu_min, dlnu));

  // Get weight from global weight parameter
  double dtau   = get_dtau(nu, type, dt, m);
  double weight = get_wgt(nu, dtau);

  // Sample emissivity in solid angle
  double jmax = jnu(nu, type, m, 0.5 * M_PI);
  do {
    cth[0] = 2. * get_rand() - 1.;
    th     = acos(cth[0]);
  } while (get_rand() > jnu(nu, type, m, th) / jmax);

  sth[0]  = sqrt(1. - cth[0] * cth[0]);
  phi     = 2. * M_PI * get_rand();
  cphi[0] = cos(phi);
  sphi[0] = sin(phi);

  // Second photon antiparallel in fluid frame
  cth[1]  = -cth[0];
  sth[1]  = sth[0];
  cphi[1] = -cphi[0];
  sphi[1] = -sphi[0];

  double E = nu * HPL / (ME * CL * CL);

  /*
  #pragma omp atomic
  dlepton += weight;
  */

  for (int n = 0; n < 2; n++) {
    // Initial zeros
    memset(tmp[n]->X, 0, NSUP * NDIM * sizeof(double));
    memset(tmp[n]->Kcov, 0, NSUP * NDIM * sizeof(double));
    memset(tmp[n]->Kcon, 0, NSUP * NDIM * sizeof(double));

    // Set position
    tmp[n]->X[2][0] = t + dt / 2.;
    coord(i, j, k, CENT, tmp[n]->X[2]);

    // Randomize phi for visualization if in axisymmetry and MKS
    if (N3TOT == 1 && METRIC == MKS)
      tmp[n]->X[2][3] = 2. * M_PI * get_rand();

    // Get coordinate frame wavevector
    K_tetrad[0] = -E;
    K_tetrad[1] = E * cth[n];
    K_tetrad[2] = E * cphi[n] * sth[n];
    K_tetrad[3] = E * sphi[n] * sth[n];

    tetrad_to_coord(Ecov, K_tetrad, tmp[n]->Kcov[2]);

    K_tetrad[0] *= -1.;
    tetrad_to_coord(Econ, K_tetrad, tmp[n]->Kcon[2]);

    // Re-do this to ensure k.k == 0?

    // Set superphoton weight
    tmp[n]->w = 0.5 * weight;

    // Diagnostics
    tmp[n]->nscatt    = 0;
    tmp[n]->origin[0] = nstep;
    tmp[n]->origin[1] = i;
    tmp[n]->origin[2] = j;
    tmp[n]->origin[3] = k;

    // Superphoton type
    tmp[n]->type = type;

    tmp[n]->t0 = t + dt / 2.;

    if (!is_null(tmp[n]->Kcov[2], tmp[n]->Kcon[2], tmp[n]->Kcov[2][0], 0.,
            &(tmp[n]->KdotKprev))) {
      double gamma;
      mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
      fprintf(stderr,
          "Error! K.K != 0 initially!\n"
          "K.K make err [%i %i %i] nu = %e w = %e n = %i K.K = %e\n"
          "K_0 = %e gamma = %e\n",
          i, j, k, nu, weight, n, tmp[n]->KdotKprev, tmp[n]->Kcov[2][0], gamma);
    }

    if (tmp[n]->Kcov[2][0] > 0.) {
      tmp[n]->w = 0.;
    }

    // Record radiation four-force
    for (int mu = 0; mu < NDIM; mu++) {
#pragma omp atomic
      radG[i][j][k][mu] -= 1 / (dt * dx[1] * dx[2] * dx[3]) * kphys_to_num *
                           tmp[n]->w * tmp[n]->Kcov[2][mu];
    }
#if RADIATION == RADTYPE_NEUTRINOS
    {
      double gamma, alpha, ucon0;
      mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
      alpha = ggeom[i][j][CENT].alpha;
      ucon0 = gamma / alpha;

#pragma omp atomic
      radG[i][j][k][RADG_YE] -=
          ((1 / (dt * dx[1] * dx[2] * dx[3])) * ucon0 * tmp[n]->w *
              (MP / M_unit) * get_lepton_sign(tmp[n]));

#pragma omp atomic
      radG[i][j][k][RADG_YE_EM] -=
          ((1 / (dt * dx[1] * dx[2] * dx[3])) * ucon0 * tmp[n]->w *
              (MP / M_unit) * get_lepton_sign(tmp[n]));
    }
#endif // RADTYPE_NEUTRINOS

#pragma omp atomic
    Jrad[0][i][j][k] -= (dt / DTd) * tmp[n]->Kcov[2][0] * kphys_to_num *
                        tmp[n]->w /
                        (ggeom[i][j][CENT].g * dt * dx[1] * dx[2] * dx[3]);

#pragma omp atomic
    Nem[i][j][k] += 1;

#pragma omp atomic
    Nem_phys[i][j][k][tmp[n]->type] += tmp[n]->w;

    if (get_rand() < ((double)Nph_to_track) / (nph_per_proc * mpi_nprocs())) {
      tmp[n]->is_tracked = 1;
    } else {
      tmp[n]->is_tracked = 0;
    }
    tmp[n]->osc_count = 0;
  }
}

#define TINY (1.e-200)
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS + 1],
    int type, const struct of_microphysics *m) {
  for (int n = 0; n < NU_BINS; n++) {
    dndlnu[n] = 0.;
  }

  double dndlnu_max = -1.e100;
  for (int n = 0; n <= NU_BINS; n++) {
    double Jsamp = Jnu(nusamp[n], type, m);
    Jsamp *= dx[1] * dx[2] * dx[3] * pow(L_unit, 3.) * ggeom[i][j][CENT].g;

    double wgt = get_wgt(nusamp[n], get_dtau(nusamp[n], type, dt, m));

    // dndlnu[n] = Jsamp/(wgtC/nusamp[n]*HPL + TINY);
    dndlnu[n] = Jsamp / (wgt * HPL + TINY);

    if (dndlnu[n] > dndlnu_max) {
      dndlnu_max = dndlnu[n];
    }
  }

  for (int n = 0; n <= NU_BINS; n++) {
    dndlnu[n] /= dndlnu_max;
  }
}
#undef TINY

void set_weight(grid_prim_type Prad, grid_eosvar_type extra) {
  double Jtot;
  double zoneVol = dV * L_unit * L_unit * L_unit;

  // Set static variables
  Ns      = tune_emiss / (pow(sim_vol, 1. / 3.) * T_unit / CL);
  lnu_min = log(numin);
  lnu_max = log(numax);
  dlnu    = (lnu_max - lnu_min) / NU_BINS;
  for (int n = 0; n <= NU_BINS; n++) {
    nusamp[n] = exp(n * dlnu + lnu_min);
  }
  Jtot = 0.;

#pragma omp parallel
  {
#pragma omp for collapse(3) reduction(+ : Jtot)
    ZLOOP {
      struct of_microphysics m;
      double                 Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
      get_fluid_zone(i, j, k, Prad, extra, &m, Ucon, Ucov, Bcon, Bcov);
      TYPELOOP {
        for (int n = 0; n <= NU_BINS; n++) {
          Jtot += Jnu(nusamp[n], itp, &m) * zoneVol * ggeom[i][j][CENT].g;
        }
      } // TYPELOOP
    }   // ZLOOP
  }     // omp parallel

  Jtot = mpi_reduce(Jtot);

  wgtC = Jtot / (HPL * Ns) * nusamp[0];
  // printf("wgtC = %g\n",wgtC); // DEBUG
}

// Use gsl to integrate dNs/dlnu in each zone
struct of_params {
  int                     type;
  struct of_microphysics *microphysics;
};
double f(double x, void *params) {
  struct of_params *      p    = (struct of_params *)params;
  struct of_microphysics *m    = p->microphysics;
  int                     type = p->type;

  double nu    = exp(x);
  double Jsamp = Jnu(nu, type, m) * nu;
  double wgt   = get_wgt(nu, get_dtau(nu, type, dt, m));
  if (isinf(Jsamp) || wgt < SMALL) {
    return 0.;
  } else {
    return Jsamp / (nu * wgt);
  }
}

// Calculate number of superphotons to produce per thread in each zone
void        get_dnz(grid_prim_type Prad, grid_eosvar_type extra) {
#pragma omp parallel
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    double                     result, error;
    gsl_function               F;
    F.function = &f;

    double zoneVol = dV * L_unit * L_unit * L_unit;
#pragma omp for collapse(3) schedule(dynamic)
    ZLOOP {
      // struct of_microphysics m;
      // double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];

      // Ignore emission outside region of interest
      double X[NDIM];
      coord(i, j, k, CENT, X);
      if (X[1] < startx_rad[1] || X[1] > stopx_rad[1]) {
        TYPELOOP dnzs[i][j][k][itp] = 0.;
        continue;
      }
// don't emit in the atmosphere
#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK && !GAMMA_FALLBACK
      if (Prad[i][j][k][RHO] < rho_poly_thresh || Prad[i][j][k][UU] < SMALL) {
        TYPELOOP dnzs[i][j][k][itp] = 0.;
        continue;
      }
#if METRIC == MKS
      if (Prad[i][j][k][ATM] < ATM_THRESH) {
        TYPELOOP dnzs[i][j][k][itp] = 0.;
        continue;
      }
#endif // metric == mks
#endif // eos == eos_type_table

      // Get number of superphotons to be emitted
      TYPELOOP {
        struct of_params params;
        params.microphysics = &(m_grd[i][j][k]);
        params.type         = itp;
        F.params            = &params;
        gsl_integration_qags(
            &F, lnu_min, lnu_max, 1.e100, 1.e-4, 1000, w, &result, &error);
        result /= HPL;
        // result /= wgtC;
        result *= zoneVol;
        result *= ggeom[i][j][CENT].g;
        result *= dt * T_unit;

        if (isnan(result / nthreads)) {
          dnzs[i][j][k][itp] = 0.;
        } else {
          dnzs[i][j][k][itp] = result / nthreads;
        }
      } // TYPELOOP
    }   // ZLOOP
    gsl_integration_workspace_free(w);
  } // pragma omp parallel
  // DEBUG
  /*
  ZLOOP {
    struct of_microphysics m;
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double zoneVol = dV*L_unit*L_unit*L_unit*ggeom[i][j][CENT].g;
    get_fluid_zone(i, j, k, Prad, extra, &m, Ucon, Ucov, Bcon, Bcov);
    printf("Ye[%d][%d][%d] = %g\n",i,j,k,Prad[i][j][k][YE]);
    printf("zoneVol = %g\n",zoneVol);
    printf("Jnu_flat/f     =
  %g\n",4*M_PI*U_unit*cnu_flat/((numax-numin)*T_unit)); printf("int*Jnu_flat/f =
  %g\n",Jnu(numax,NU_ELECTRON,&m)*(numax-numin)/(2*m.Ye)); printf("dt =
  %g\n",dt); printf("U_unit = %g\n",U_unit); printf("numax = %g\n",numax);
    printf("numin = %g\n",numin);
    printf("T_unit = %g\n",T_unit);
    double denom = (numax-numin)*T_unit;
    printf("denom = %g\n",denom);
    double num = 4*M_PI*U_unit*cnu_flat;
    printf("num = %g\n",num);
    printf("cnu_flat = %g\n",cnu_flat);
    TYPELOOP {
      double dnz_true = dnzs[i][j][k][itp]*(HPL*wgtC)*nthreads/zoneVol;
      double dnz_expected = Jnu(numax,itp,&m)*(numax-numin)*dt*T_unit;
      printf("\tdnz[%d][%d][%d][%d] =
  %g\n",i,j,k,itp,dnzs[i][j][k][itp]*nthreads); printf("\tdnz*e[%d][%d][%d][%d]
  = %g\n",i,j,k,itp, dnz_true); printf("\t(dnz*e)_expected[%d][%d][%d][%d] =
  %g\n", i,j,k,itp, dnz_expected); if (fabs(dnz_expected) > 0) {
  printf("\t\texpected/true = %g\n", dnz_expected/(dnz_true));
  if (fabs((dnz_expected/dnz_true) - 1.) > 1e-10) exit(1);
      }
    }
  }
  */
}

#endif // RADIATION
