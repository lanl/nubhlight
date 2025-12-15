/******************************************************************************
 *                                                                            *
 * RAD_UTILS.C                                                                *
 *                                                                            *
 * HELPER FUNCTIONS FOR RADIATION INFRASTRUCTURE                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void init_rad(grid_prim_type Prad) {
  set_units();

  ZLOOP {
    sim_vol +=
        ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3] * L_unit * L_unit * L_unit;
  }
  sim_vol = mpi_reduce(sim_vol);

  photon_lists     = safe_malloc(nthreads * sizeof(struct of_photon *));
  photon_mpi_lists = safe_malloc(nthreads * sizeof(struct of_photon *));
#pragma omp parallel
  {
    photon_lists[omp_get_thread_num()]     = NULL;
    photon_mpi_lists[omp_get_thread_num()] = NULL;
  }

  init_emissivity();

  init_superphoton_resolution();

#if SCATTERING
  init_all_hotcross();
#endif

#if RADIATION == RADTYPE_NEUTRINOS
  lepton_tot = lepton_last = lepton_lost = 0.;
#endif
}

void init_superphoton_resolution() {
  made_tune_proc = abs_tune_proc = scatt_tune_proc = 0;

#if METRIC == MKS
  double cross_section;
  dt_tune_emiss = 0.5;
  dt_tune_scatt = Rout_rad;
  if (RADIATION == RADTYPE_NEUTRINOS) {
    cross_section = RAD_NUM_TYPES * RAD_SCATT_TYPES * NUSIGMA0;
  } else {
    cross_section = 16. * pow(10, 2) * THOMSON;
  }
  // goose-tuned to be small compared to actual required value
  // needed so that scattering only turns on after the first adjustment
  tune_scatt = 0.1 / (Rout_rad * L_unit * cross_section * Ne_unit);
  if (mpi_io_proc()) { // DEBUG
    printf("cross_section = %g\ntune_scatt = %g\n", cross_section, tune_scatt);
  }
#else
  dt_tune_emiss = tf;
  dt_tune_scatt = MY_MAX(MY_MAX(N1TOT * dx[1], N2TOT * dx[2]), N3TOT * dx[3]);
  dt_tune_scatt = MY_MAX(dt_tune_scatt, tf);
#endif
  if (t_tune_emiss <= 0.0)
    t_tune_emiss = 2. * dt_tune_emiss;
  if (t_tune_scatt <= 0.0)
    t_tune_scatt = 2. * dt_tune_scatt;
#if ELECTRONS
  if (strcmp(init_from_grmhd, "No") == 0 && fel0 < 0.1) {
    t_tune_emiss = 500.;
    t_tune_scatt = 500.;
  }
#endif
}

void update_superphoton_resolution(
    grid_prim_type Prad, grid_eosvar_type extra) {

#if KILL_ALL_PACKETS
  return;
#else

  double L, real, ideal, correction;
  int    made_tune, abs_tune, scatt_tune;

  made_tune_proc += step_made;
  abs_tune_proc += step_abs;
  scatt_tune_proc += step_scatt;

#if METRIC == MKS
  L = Rout_rad;
#else
  L = MY_MAX(MY_MAX(N1TOT * dx[1], N2TOT * dx[2]), N3TOT * dx[3]);
#endif

  if (t >= t_tune_emiss) {
    made_tune = mpi_reduce_int(made_tune_proc) / mpi_nprocs();
    abs_tune  = mpi_reduce_int(abs_tune_proc) / mpi_nprocs();

    real  = made_tune - abs_tune;
    ideal = MY_MIN(SCATT_BIAS_SAFETY, 1.0) * dt_tune_emiss * nph_per_proc / L;
    correction = ideal / real;

    // Limit strength of correction
    if (correction < 0.) {
      correction = 4. / 3.;
    } else {
      correction = MY_MIN(correction, 4. / 3.);
      correction = MY_MAX(correction, 1. / 2.); // MY_MAX(correction, 3./4.);
    }

    // If no superphotons are being emitted (yet) don't modify emission strength
    // if (real < SMALL) correction = 1.;

    tune_emiss *= correction;

    if (mpi_io_proc()) {
      fprintf(stdout, "Emission correction! tune = %e correction = %e\n",
          tune_emiss, correction);
    }

    t_tune_emiss += dt_tune_emiss;
    set_weight(Prad, extra);
    made_tune_proc = abs_tune_proc = 0;
  }

  scatt_tune = mpi_reduce_int(scatt_tune_proc) / mpi_nprocs();
  real       = scatt_tune;
  ideal      = dt_tune_scatt * nph_per_proc / L;
  correction = ideal / real;
  if (t >= t_tune_scatt || (correction < 0.25 && METRIC == MKS)) {
    // scatt_tune = mpi_reduce_int(scatt_tune_proc)/mpi_nprocs();

    real       = scatt_tune;
    ideal      = dt_tune_scatt * nph_per_proc / L;
    correction = ideal / real;

    // Limit strength of correction
    correction = MY_MIN(correction, 1.25); // MY_MIN(correction, 2.0);
    correction = MY_MAX(correction, 0.25);

    // If no superphotons are being emitted (yet) don't modify emission strength
    if (real < SMALL)
      correction = 1.;

    tune_scatt *= correction;

    if (mpi_io_proc()) {
      fprintf(stdout, "Scattering correction! tune = %e correction = %e\n",
          tune_scatt, correction);
    }

    t_tune_scatt += dt_tune_scatt;
    scatt_tune_proc = 0;
  }
#endif
}

double linear_interp_log(double x, double *table, double lx_min, double dlx) {
  double lx = log(x);
  double dn = (lx - lx_min) / dlx;
  int    n  = (int)dn;
  dn        = dn - n;

  return (1. - dn) * table[n] + dn * table[n + 1];
}

// Remove superphoton from list and release memory
void list_remove(struct of_photon **ph, struct of_photon **ph_head,
    struct of_photon **ph_prev) {
  if (*ph_prev != NULL) {
    (*ph_prev)->next = (*ph)->next;
    free(*ph);
    *ph = (*ph_prev)->next;
  } else {
    *ph_head = (*ph)->next;
    free(*ph);
    *ph = *ph_head;
  }
}

double get_Thetae(double Prad[NVAR]) {
  double Thetae;
#if ELECTRONS
  Thetae = Prad[KEL] * pow(Prad[RHO], game - 1.) * Thetae_unit;
#else
  Thetae = Prad[UU] / Prad[RHO] * Thetae_unit;
#endif

  return MY_MIN(Thetae, thetae_max);
}

// Used for Maxwell distribution. Calculate k_b T/(m*c*c)
// where c is scatterer mass
double scatterer_dimensionless_temp(
    int radtype, int interaction, const struct of_microphysics *m) {
#if RADIATION == RADTYPE_LIGHT
  return m->Thetae;
#elif RADIATION == RADTYPE_NEUTRINOS
  {
    double mass;
#if MULTISCATT_TEST
    {
      // Neutrons are scatterers
      // Temperature is assumed to be in ergs
      // TODO: set temperature to zero for test?
      mass = MP + ME;
      // return m->T/((MP+ME)*CL*CL);
    }
#else // normal neutrino scattering
    {
      if (interaction == RSCATT_TYPE_P)
        mass = MP;
      else if (interaction == RSCATT_TYPE_N)
        mass = MN;
      else if (interaction == RSCATT_TYPE_A)
        mass = (MN + MP) * m->Abar;
      else if (interaction == RSCATT_TYPE_ALPHA)
        mass = 4 * (MN + MP);
      else if (((radtype == NU_ELECTRON) || (radtype == ANTINU_ELECTRON)) &&
               interaction == RSCATT_TYPE_E) {
        mass = ME;
      } else {
        fprintf(stderr, "rad_utils: Unknown interaction type!\n");
        exit(1);
      }
    }
#endif
    return m->T / (mass * CL * CL);
  }
#endif
}

// Used for scattering optical depth. Calculate number
// density of scatter particles.
double scatterer_number_density(
    int radtype, int interaction, const struct of_microphysics *m) {
#if RADIATION == RADTYPE_LIGHT
  { return m->Ne; }
#elif RADIATION == RADTYPE_NEUTRINOS
  {
#if MULTISCATT_TEST
    {
      // use neutron number density
      double rho_n = (m->rho) * (m->Xi[MF_XN]);
      double Nn    = rho_n / (MP + ME);
      return Nn;
    }
#else  // normal neutrino scattering
    {
      int    i;
      double mass;
      // assume there are equal number protons and electrons
      if ((interaction == RSCATT_TYPE_P) ||
          ((interaction == RSCATT_TYPE_E) &&
              ((radtype == NU_ELECTRON) || (radtype == ANTINU_ELECTRON)))) {
        i    = MF_XP;
        mass = MP;
      } else if (interaction == RSCATT_TYPE_N) {
        i    = MF_XN;
        mass = MN;
      } else if (interaction == RSCATT_TYPE_A) {
        i    = MF_XH;
        mass = (MN + MP) * (m->Abar);
      } else if (interaction == RSCATT_TYPE_ALPHA) {
        i    = MF_XA;
        mass = 4 * (MN + MP);
      } else {
        fprintf(stderr, "rad_utils: unknown interaction type!\n");
        exit(1);
      }
      return (m->rho) * (m->Xi[i]) / mass;
    }
#endif // type of neutrino scattering
  }
#endif
}

void precompute_microphysics() {
  timer_start(TIMER_MICRO);
#pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOPALL {
    double X[NDIM];

#if EOS == EOS_TYPE_TABLE && METRIC == MKS
    if (P[i][j][k][ATM] < ATM_THRESH)
      continue;
#endif

    coord(i, j, k, CENT, X);
    if (X[1] <= stopx_rad[1] + (NG + 1) * dx[1]) {
      get_fluid_zone(i, j, k, P, extra, &(m_grd[i][j][k]), Ucon_grd[i][j][k],
          Ucov_grd[i][j][k], Bcon_grd[i][j][k], Bcov_grd[i][j][k]);
#if ESTIMATE_THETAE
      m_grd[i][j][k].Thetae = get_Thetae_est(i, j, k);
#endif
    }
  }
  timer_stop(TIMER_MICRO);
}

void get_fluid_zone(int i, int j, int k, grid_prim_type Prad,
    grid_eosvar_type extra, struct of_microphysics *m, double Ucon[NDIM],
    double Ucov[NDIM], double Bcon[NDIM], double Bcov[NDIM]) {
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;

#if EOS == EOS_TYPE_TABLE
  EOS_SC_fill(Prad[i][j][k], extra[i][j][k]);
#endif

  Bp[1] = Prad[i][j][k][B1] * B_unit;
  Bp[2] = Prad[i][j][k][B2] * B_unit;
  Bp[3] = Prad[i][j][k][B3] * B_unit;

  Vcon[1] = Prad[i][j][k][U1];
  Vcon[2] = Prad[i][j][k][U2];
  Vcon[3] = Prad[i][j][k][U3];

  // Get Ucov
  VdotV = 0.;
  for (int l = 1; l < NDIM; l++) {
    for (int m = 1; m < NDIM; m++) {
      VdotV += ggeom[i][j][CENT].gcov[l][m] * Vcon[l] * Vcon[m];
    }
  }
  Vfac    = sqrt(-1. / ggeom[i][j][CENT].gcon[0][0] * (1. + fabs(VdotV)));
  Ucon[0] = -Vfac * ggeom[i][j][CENT].gcon[0][0];
  for (int l = 1; l < NDIM; l++)
    Ucon[l] = Vcon[l] - Vfac * ggeom[i][j][CENT].gcon[0][l];
  lower(Ucon, ggeom[i][j][CENT].gcov, Ucov);

  // Get Bcon, Bcov, and B
  UdotBp = 0.;
  for (int l = 1; l < NDIM; l++)
    UdotBp += Ucov[l] * Bp[l];
  Bcon[0] = UdotBp;
  for (int l = 1; l < NDIM; l++)
    Bcon[l] = (Bp[l] + Ucon[l] * UdotBp) / Ucon[0];
  lower(Bcon, ggeom[i][j][CENT].gcov, Bcov);

  m->B = sqrt(Bcon[0] * Bcov[0] + Bcon[1] * Bcov[1] + Bcon[2] * Bcov[2] +
              Bcon[3] * Bcov[3]);
#if RADIATION == RADTYPE_NEUTRINOS
  double rho = Prad[i][j][k][RHO];
  double u   = Prad[i][j][k][UU];
  double T   = EOS_temperature(rho, u, extra[i][j][k]);
  m->rho     = rho * RHO_unit;
  m->T       = T * TEMP_unit;
  m->Ye      = Prad[i][j][k][YE];
  EOS_SC_mass_fractions(m->Xi, extra[i][j][k]);
  EOS_SC_avg_ions(&(m->Abar), &(m->Zbar), extra[i][j][k]);
#if BURROWS_OPACITIES
  // Todo: move this to only the spots it matters?
  fill_opac_emis_burrows(m);
#elif HDF5_OPACITIES
  fill_opac_emis_hdf(m);
#endif
#elif RADIATION == RADTYPE_LIGHT
  m->Ne     = Prad[i][j][k][RHO] * Ne_unit;
  m->Thetae = get_Thetae(Prad[i][j][k]);
  // Prevent highly magnetized regions from emitting due to bad internal energy
  double sigma = pow(m->B / B_unit, 2.) / (m->Ne / Ne_unit);
  if (sigma > sigma_max) {
    m->Thetae = SMALL;
  }
#endif
}

int is_null(double Kcov[NDIM], double Kcon[NDIM], double K0, double KdotKprev,
    double *KdotK) {
  *KdotK = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    *KdotK += Kcov[mu] * Kcon[mu];
  }

  double K0sqr = pow(K0, 2.);

  if (fabs(*KdotK - KdotKprev) / K0sqr < kdotk_tol) {
    return 1;
  } else {
    return 0;
  }
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k) {
  *i = (X[1] - startx[1]) / dx[1] + NG - global_start[1];
  *j = (X[2] - startx[2]) / dx[2] + NG - global_start[2];
  *k = (X[3] - startx[3]) / dx[3] + NG - global_start[3];
}

void copy_photon(struct of_photon *ph, struct of_photon *phc) {
  for (int mu = 0; mu < NDIM; mu++) {
    for (int n = 0; n < NSUP; n++) {
      phc->X[n][mu]    = ph->X[n][mu];
      phc->Kcov[n][mu] = ph->Kcov[n][mu];
      phc->Kcon[n][mu] = ph->Kcon[n][mu];
    }
    phc->origin[mu] = ph->origin[mu];
  }
  phc->w          = ph->w;
  phc->KdotKprev  = ph->KdotKprev;
  phc->type       = ph->type;
  phc->nscatt     = ph->nscatt;
  phc->t0         = ph->t0;
  phc->is_tracked = ph->is_tracked;
  phc->osc_count  = ph->osc_count;
}

void print_ph_diag(struct of_photon *ph) {
  printf(" --- PHOTON STATE --- \n");
  for (int n = 0; n < 3; n++) {
    for (int mu = 0; mu < NDIM; mu++) {
      printf("[%i][%i] X = %e Kcov = %e Kcon = %e\n", n, mu, ph->X[n][mu],
          ph->Kcov[n][mu], ph->Kcon[n][mu]);
    }
    double r, th;
    bl_coord(ph->X[n], &r, &th);
    printf("r, th = %e %e\n", r, th);
  }
  printf("origin = %i %i %i %i\n", ph->origin[0], ph->origin[1], ph->origin[2],
      ph->origin[3]);
  printf("w = %e\n", ph->w);
  printf("type = %d\n", ph->type);
  printf("osc_count = %d\n", ph->osc_count);
  printf("K.Kprev = %e\n", ph->KdotKprev);
  printf("nscatt = %i\n", ph->nscatt);
  printf("t0 = %e\n", ph->t0);
}

// Use second-order interpolation to get X^{\mu}, K_{\mu} at time t_interp
int get_X_K_interp(struct of_photon *ph, double t_interp, grid_prim_type P,
    double X[NDIM], double Kcov[NDIM], double Kcon[NDIM]) {
  double *Xa, *Xb, *Kcova, *Kcovb, *Kcona, *Kconb;

  if (t_interp == ph->X[2][0]) {
    for (int mu = 0; mu < NDIM; mu++) {
      X[mu]    = ph->X[2][mu];
      Kcov[mu] = ph->Kcov[2][mu];
      Kcon[mu] = ph->Kcon[2][mu];
    }
    return SPH_INTERP_SUCCESS;
  }

  double KdotKprev;
  if (t_interp < (1. - 1.e-50) * ph->X[1][0]) {
    Xa        = ph->X[0];
    Xb        = ph->X[1];
    Kcova     = ph->Kcov[0];
    Kcovb     = ph->Kcov[1];
    Kcona     = ph->Kcon[0];
    Kconb     = ph->Kcon[1];
    KdotKprev = dot(ph->Kcov[0], ph->Kcon[0]);
  } else if (t_interp <= (1. + 1.e-50) * ph->X[2][0]) {
    Xa        = ph->X[1];
    Xb        = ph->X[2];
    Kcova     = ph->Kcov[1];
    Kcovb     = ph->Kcov[2];
    Kcona     = ph->Kcon[1];
    Kconb     = ph->Kcon[2];
    KdotKprev = dot(ph->Kcov[1], ph->Kcon[1]);
  } else {
    printf("BAD INTERPOLATION\n");
    for (int mu = 0; mu < NDIM; mu++) {
      printf(
          "X[][%i] = %e %e %e\n", mu, ph->X[0][mu], ph->X[1][mu], ph->X[2][mu]);
    }
    printf("startx[] = %e %e %e\n", startx[1], startx[2], startx[3]);
    for (int nint = 0; nint < 3; nint++) {
      double r, th;
      bl_coord(ph->X[nint], &r, &th);
      printf("r, th = %e %e\n", r, th);
    }
    printf("nscatt = %i\n", ph->nscatt);
    printf("t_interp = %e t = %e\n", t_interp, t);
    printf("origin: %i %i %i %i\n", ph->origin[0], ph->origin[1], ph->origin[2],
        ph->origin[3]);
    exit(-1);
  }

  double fac = (t_interp - Xa[0]) / (Xb[0] - Xa[0]);

  for (int mu = 0; mu < NDIM; mu++) {
    X[mu]    = fac * Xb[mu] + (1. - fac) * Xa[mu];
    Kcov[mu] = fac * Kcovb[mu] + (1. - fac) * Kcova[mu];
    Kcon[mu] = fac * Kconb[mu] + (1. - fac) * Kcona[mu];
  }

  /*if (Kcov[0] > 0. && ph->w > 0.) {
    printf("Kcov[0] > 0 in interp!\n");
    printf("t_interp = %e\n", t_interp);
    for (int mu = 0; mu < 3; mu++) {
      printf("ph->X[%i][] = %e %e %e %e\n", mu, ph->X[mu][0], ph->X[mu][1],
  ph->X[mu][2], ph->X[mu][3]); printf("ph->Kcov[%i][] = %e %e %e %e\n", mu,
  ph->Kcov[mu][0], ph->Kcov[mu][1], ph->Kcov[mu][2], ph->Kcov[mu][3]);
    }
    printf("nscatt = %i, w = %e\n", ph->nscatt, ph->w);
  }*/

  double kdotk = dot(Kcon, Kcov);
  int    status;
  // if (1) {
  if (fabs((kdotk - KdotKprev) / (Kcov[0] * Kcov[0])) > 100. * kdotk_tol) {
    if (t_interp < (1. - 1.e-50) * ph->X[1][0]) {
      for (int mu = 0; mu < NDIM; mu++) {
        X[mu]    = ph->X[0][mu];
        Kcov[mu] = ph->Kcov[0][mu];
        Kcon[mu] = ph->Kcon[0][mu];
      }
      double KdotKprev = dot(Kcov, Kcon);
      // Should not use Phalf for this push
      status = push_X_K(
          X, Kcov, Kcon, P, P, KdotKprev, ph->type, t_interp - ph->X[0][0]);
    } else {
      for (int mu = 0; mu < NDIM; mu++) {
        X[mu]    = ph->X[1][mu];
        Kcov[mu] = ph->Kcov[1][mu];
        Kcon[mu] = ph->Kcon[1][mu];
      }
      double KdotKprev = dot(Kcov, Kcon);
      // Should not use Phalf for this push
      status = push_X_K(
          X, Kcov, Kcon, P, P, KdotKprev, ph->type, t_interp - ph->X[1][0]);
    }
    if (status == PUSH_FAIL) {
      fprintf(stderr, "get_X_K_interp failed!\n");
      fprintf(stderr, "X[] = %e %e %e %e\n", X[0], X[1], X[2], X[3]);
      return SPH_INTERP_FAIL;
    }
  }

  return SPH_INTERP_SUCCESS;
}

// Does superphoton need to be pushed from step t to t + dt?
int to_be_pushed(double t, double dt, struct of_photon *ph) {
  if (ph->type == TYPE_TRACER)
    return 1;

  int i, j, k;
  Xtoijk(ph->X[2], &i, &j, &k);

  if (ph->X[2][0] < t + dt) {
    return 1;
  } else {
    return 0;
  }
}

// How big of a superphoton push?
double get_dtpush(struct of_photon *ph, double dt) {
  if (ph->type == TYPE_TRACER)
    return dt;

  int i, j, k;
  Xtoijk(ph->X[2], &i, &j, &k);
  return MY_MAX(cour * dt_light[i][j], dt);
}

// Move photon from donor list to head of recipient list; advance donor list
void swap_ph(struct of_photon **donor, struct of_photon **recipient) {
  struct of_photon *tmp;
  if (*recipient == NULL) {
    *recipient         = *donor;
    *donor             = (*donor)->next;
    (*recipient)->next = NULL;
  } else {
    tmp        = *donor;
    *donor     = (*donor)->next;
    tmp->next  = *recipient;
    *recipient = tmp;
  }
}

void set_Rmunu() {
  memset((void *)Rmunu, 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * NDIM * NDIM *
          sizeof(double));
  memset((void *)Nsph, 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * sizeof(int));
  memset((void *)nph, 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * sizeof(double));
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    double            X[NDIM], Kcov[NDIM], Kcon[NDIM];
    while (ph != NULL) {

      if (ph->type == TYPE_TRACER) {
        ph = ph->next;
        continue;
      }

      int i, j, k;
      get_X_K_interp(ph, t, P, X, Kcov, Kcon);
      Xtoijk(X, &i, &j, &k);

      double volume = ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3];

#pragma omp atomic
      Nsph[i][j][k] += 1;

#pragma omp atomic
      nph[i][j][k] += ph->w / (volume * pow(L_unit, 3));

      DLOOP2 {
#pragma omp atomic
        Rmunu[i][j][k][mu][nu] +=
            kphys_to_num * Kcon[mu] * Kcov[nu] * ph->w / (Kcon[0] * volume);
      }
      ph = ph->next;
    } // photon loop
  }   // omp parallel
}

void get_nuLnu_bin(double X[NDIM], int *thbin, int *phibin) {
  double r, th, phi;
  bl_coord(X, &r, &th);
  phi = fmod(X[3], 2. * M_PI);
  // phi = X[3] % (2.*M_PI);

  double dth  = M_PI / NTH;
  double dphi = 2. * M_PI / NPHI;

  *thbin  = (int)(th / dth);
  *phibin = (int)(phi / dphi);
}

void bin_superphoton_direction(const struct of_photon *ph) {
  // DO NOT DO THIS FOR TRACERS
  if (ph->type == TYPE_TRACER)
    return;

  const int klevel = 2;
  double    lnumin = log(numin);
  double    lnumax = log(numax);
  double    dlnu   = (lnumax - lnumin) / NU_BINS_SPEC;
  double    dth    = M_PI / NTH;
  double    dphi   = 2. * M_PI / NPHI;
  double    Kcov[NDIM];

  // Assume X0 symmetry in metric
  double E   = -ph->Kcov[klevel][0] * ME * CL * CL;
  double nu  = E / HPL;
  double lnu = log(nu);
  DLOOP1 {
    // put wavevector on unit sphere. Not necessary,
    // but convenient for debugging
    Kcov[mu] = -ph->Kcov[klevel][mu] / ph->Kcov[klevel][0];
  }

  // Uses HARM coordinates for angles. Not physical angles.
  // use with caution!
  double r, th, phi;
  cart_to_sph(Kcov, &r, &th, &phi);

  int phibin = (int)(phi / dphi);
  int thbin  = (int)(th / dth);
  int nubin  = (int)((lnu - lnumin) / dlnu);

  if (nubin >= 0 && nubin < NU_BINS_SPEC) {
// Unlike for camera BC, we don't want a light curve. We just want
// a probability distribution.
#if DIAGNOSTICS_USE_RADTYPES
    {
#pragma omp atomic
      nuLnu[ph->type][thbin][phibin][nubin] += ph->w;
    }
#else
    {
      int nscatt = MY_MIN(ph->nscatt, MAXNSCATT);

#pragma omp atomic
      nuLnu[nscatt][thbin][phibin][nubin] += ph->w;
    }
#endif
  }
}

double get_dt_cool_rad_zone(double u, struct of_microphysics *m) {
  double u_cgs       = fabs(u * U_unit);
  double J           = fabs(get_J(m));
  double dt_cool_cgs = (u_cgs + SMALL) / (J + SMALL);
  double dt_cool     = dt_cool_cgs / T_unit;
  return dt_cool;
}

void set_cooling_time(
    grid_double_type tau_cool, grid_prim_type P, grid_eosvar_type extra) {
#pragma omp parallel for collapse(3)
  ZLOOP {
    double                 X[NDIM];
    struct of_microphysics m;
    double                 Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];

    coord(i, j, k, CENT, X);
    if (X[1] < startx_rad[1] || X[1] > stopx_rad[1]) {
      tau_cool[i][j][k] = INFINITY;
    }
#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK && !GAMMA_FALLBACK
    else if (P[i][j][k][RHO] < rho_poly_thresh || P[i][j][k][UU] < SMALL) {
      tau_cool[i][j][k] = INFINITY;
    }
#endif
    else {
      get_fluid_zone(i, j, k, P, extra, &m, Ucon, Ucov, Bcon, Bcov);
#if ESTIMATE_THETAE
      m.Thetae = get_Thetae_est(i, j, k);
#endif

      tau_cool[i][j][k] = get_dt_cool_rad_zone(P[i][j][k][UU], &m);
    }
  }
}

double get_min_dt_cool(grid_prim_type P, grid_eosvar_type extra) {
  double dt_cool = INFINITY;
#pragma omp parallel for reduction(min : dt_cool) collapse(3)
  ZLOOP {
    double                 X[NDIM];
    struct of_microphysics m;
    double                 Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];

    // ignore outside of emissivity region
    coord(i, j, k, CENT, X);
    if (X[1] < startx_rad[1] || X[1] > stopx_rad[1])
      continue;
// ignore atmosphere
#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK && !GAMMA_FALLBACK
    if (P[i][j][k][RHO] < rho_poly_thresh || P[i][j][k][UU] < SMALL)
      continue;
#if METRIC == MKS
    if (P[i][j][k][ATM] < ATM_THRESH)
      continue;
#endif
#endif

    get_fluid_zone(i, j, k, P, extra, &m, Ucon, Ucov, Bcon, Bcov);
#if ESTIMATE_THETAE
    m.Thetae = get_Thetae_est(i, j, k);
#endif

    double dtc_zone = get_dt_cool_rad_zone(P[i][j][k][UU], &m);
    if (dtc_zone < dt_cool)
      dt_cool = dtc_zone;
  }
  dt_cool = mpi_min(dt_cool);
  return dt_cool;
}

#if RADIATION == RADTYPE_NEUTRINOS
void        record_lepton_flux(const struct of_photon *ph) {
#pragma omp atomic
  lepton_lost_local += (ph->w) * get_lepton_sign(ph);
}

int get_lepton_sign(const struct of_photon *ph) {
  if (ph->type == NU_ELECTRON)
    return 1;
  if (ph->type == ANTINU_ELECTRON)
    return -1;
  return 0;
}
int nu_is_heavy(const int radtype) {
  return ((radtype == NU_HEAVY) || (radtype == ANTINU_HEAVY));
}

// for debugging
void        check_nu_type(const char *location) {
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type < TYPE_TRACER || ph->type > RAD_NUM_TYPES) {
        fprintf(stderr,
            "[%s] Photon has bad type!\n"
            "\tw = %g\n"
            "\tKdotKprev = %g\n"
            "\ttype = %d\n"
            "\tnscatt = %d\n"
            "\t[X][0]    = [%g, %g, %g %g]\n"
            "\t[X][1]    = [%g, %g, %g %g]\n"
            "\t[X][2]    = [%g, %g, %g %g]\n"
            "\t[Kcon][0] = [%g, %g, %g %g]\n"
            "\t[Kcon][1] = [%g, %g, %g %g]\n"
            "\t[Kcon][2] = [%g, %g, %g %g]\n"
            "\t[Kcov][0] = [%g, %g, %g %g]\n"
            "\t[Kcov][1] = [%g, %g, %g %g]\n"
            "\t[Kcov][2] = [%g, %g, %g %g]\n"
            "\t[origin] = [%d, %d, %d %d]\n"
            "\tt0 = %g\n"
            "\tis_tracked = %d\n",
            location, ph->w, ph->KdotKprev, ph->type, ph->nscatt, ph->X[0][0],
            ph->X[0][1], ph->X[0][2], ph->X[0][3], ph->X[1][0], ph->X[1][1],
            ph->X[1][2], ph->X[1][3], ph->X[2][0], ph->X[2][1], ph->X[2][2],
            ph->X[2][3], ph->Kcon[0][0], ph->Kcon[0][1], ph->Kcon[0][2],
            ph->Kcon[0][3], ph->Kcon[1][0], ph->Kcon[1][1], ph->Kcon[1][2],
            ph->Kcon[1][3], ph->Kcon[2][0], ph->Kcon[2][1], ph->Kcon[2][2],
            ph->Kcon[2][3], ph->Kcov[0][0], ph->Kcov[0][1], ph->Kcov[0][2],
            ph->Kcov[0][3], ph->Kcov[1][0], ph->Kcov[1][1], ph->Kcov[1][2],
            ph->Kcov[1][3], ph->Kcov[2][0], ph->Kcov[2][1], ph->Kcov[2][2],
            ph->Kcov[2][3], ph->origin[0], ph->origin[1], ph->origin[2],
            ph->origin[3], ph->t0, ph->is_tracked);
        exit(1);
      }
      ph = ph->next;
    }
  }
}

#endif // NEUTRINOS

unsigned long int count_particles_local() {
  unsigned long int count = 0;

#pragma omp parallel reduction(+ : count)
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->w > SMALL)
        count++;
      ph = ph->next;
    }
  }
  return count;
}

#endif // RADIATION
