/******************************************************************************
 *                                                                            *
 * ESTIMATE_THETAE.C                                                          *
 *                                                                            *
 * ESTIMATE THETAE AT END OF TIMESTEP DUE TO RADIATION SOURCES                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
#if ESTIMATE_THETAE

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define NTHREADS (24)

// ESTIMATE THETAE BEFORE PUSH TO AVOID ISSUES WITH MPI COMMUNICATION

double Thetae_est[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
double Thetae_old[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];

double Ucon[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NDIM];
double Ucov[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NDIM];
double Bcov[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NDIM];
double Ne[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
double Bmag[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];

int        Nsph_zone[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
int        Nsph_cntd[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
double *   w[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
double *   nu[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
double *   dlam[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
double *   theta[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG][NTHREADS];
static int type        = 0; // TODO: remove me when type loops are in place.
static int interaction = 0; // TODO: remove me when type loops are in place.

double get_Thetae_est(int i, int j, int k) { return Thetae_est[i][j][k]; }

// Rootfinding parameters and function
struct of_params {
  int    i, j, k;
  double rho, Ne, Bmag, Thetaei, Ucon0, dt;
  double extra[EOS_NUM_EXTRA];
};

double dEdt(double Thetae, void *params) {
  struct of_params *p     = (struct of_params *)params;
  int               i     = p->i;
  int               j     = p->j;
  int               k     = p->k;
  double            rho   = p->rho;
  double            dt    = p->dt;
  double *          extra = p->extra;

  // TODO:  encapsulate electrons in EOS framework.
  double uei, uef;
#if ELECTRONS && EOS == EOS_TYPE_GAMMA
  uei = Ne[i][j][k] * Thetae_old[i][j][k] * ME * CL * CL / (game - 1.);
  uef = Ne[i][j][k] * Thetae * ME * CL * CL / (game - 1.);
#else
  uei = EOS_u_N_Theta(rho, Ne[i][j][k], Thetae_old[i][j][k], extra);
  uef = EOS_u_N_Theta(rho, Ne[i][j][k], Thetae, extra);
#endif

  struct of_microphysics micro;
  micro.Thetae = Thetae;
  micro.Ne     = Ne[i][j][k];
  micro.B      = Bmag[i][j][k];

  double J   = get_J(&m);
  double vol = ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3] * dt;

  // Loop over superphotons
  double udotG_abs   = 0.;
  double udotG_scatt = 0.;
  for (int n = 0; n < nthreads; n++) {
    for (int m = 0; m < Nsph_zone[i][j][k][n]; m++) {
      double udotk = HPL * nu[i][j][k][n][m] / (ME * CL * CL);

      // Absorption
      double alpha_inv_a =
          alpha_inv_abs(nu[i][j][k][n][m], type, &micro, theta[i][j][k][n][m]);
      double dtau_a =
          alpha_inv_a * L_unit * HPL / (ME * CL * CL) * dlam[i][j][k][n][m];
      double dw_a = w[i][j][k][n][m] * (1. - exp(-dtau_a));
      udotG_abs += kphys_to_num * dw_a * udotk / vol;

      // dEtau_abs += HPL*nu[i][j][k][n][m]*w[i][j][k][n][m]*(1. -
      // exp(-dtau_a));

      // Scattering (assuming h \nu << k_B T_e)
      double alpha_inv_s =
          alpha_inv_scatt(nu[i][j][k][n][m], type, interaction, &micro);
      double dtau_s =
          alpha_inv_s * L_unit * HPL / (ME * CL * CL) * dlam[i][j][k][n][m];
      double dw_s = w[i][j][k][n][m] * (1. - exp(-dtau_s));
      double amp =
          1. + 4. * Thetae - 2. * pow(Thetae, 3. / 2.) + 16. * pow(Thetae, 2.);
      udotG_scatt += kphys_to_num * (1. - amp) * dw_s * udotk / vol;
      // dEdtau_scatt -= HPL*nu[i][j][k][n][m]*;
    }
  }
  // dEdtau_abs *= Ucon[i][j][k][0]/(/* d3xi! */dt*T_unit);
  // dEdtau_scatt *= Ucon[i][j][k][0]/(dt*T_unit);
  udotG_abs *= U_unit / T_unit;
  udotG_scatt *= U_unit / T_unit;

  // if (udotG_abs != 0. || udotG_scatt != 0.) {
  /*if (i == 80 && j == 64) {
    printf("Thetae Ne Bmag = %e %e %e\n", Thetae, Ne[i][j][k], Bmag[i][j][k]);
    printf("%e %e %e %e\n", Ucon[i][j][k][0]*(uef - uei)/dt,
      J, udotG_abs, udotG_scatt);
  }*/

  // Solve entropy equation in cgs:
  //   d u_e / d \tau = -emission + absorption - upscattering
  double resid = Ucon[i][j][k][0] * (uef - uei) / (dt * T_unit);
  resid += J;
  resid += udotG_abs;
  resid += udotG_scatt;
  return resid;
}

void estimate_Thetae(
    grid_prim_type P, grid_eosvar_type extra, double t, double dt) {
  if (NTHREADS != nthreads) {
    fprintf(stderr, "NTHREADS = %i nthreads = %i! Exiting...\n", NTHREADS,
        nthreads);
    exit(-1);
  }
  struct of_microphysics micro;

// Count superphotons in each zone
#pragma omp parallel
  {
    int               n  = omp_get_thread_num();
    struct of_photon *ph = photon_lists[n];
    while (ph != NULL) {
      int    i, j, k;
      double X[NDIM], Kcov[NDIM], Kcon[NDIM];
      get_X_K_interp(ph, t, P, X, Kcov, Kcon);
      Xtoijk(X, &i, &j, &k);
      Nsph_zone[i][j][k][n]++;
      ph = ph->next;
    }
  } // omp parallel

// malloc required memory and store ucon for convenience
#pragma omp parallel for collapse(3)
  ZLOOP {
    double Bcon[NDIM];
    get_fluid_zone(i, j, k, P, extra, &micro, Ucon[i][j][k], Ucov[i][j][k],
        Bcon, Bcov[i][j][k]);
    Ne[i][j][k]         = micro.Ne;
    Thetae_old[i][j][k] = micro.Thetae;
    Bmag[i][j][k]       = micro.B;
    for (int n = 0; n < nthreads; n++) {
      w[i][j][k][n]     = safe_malloc(Nsph_zone[i][j][k][n] * sizeof(double));
      nu[i][j][k][n]    = safe_malloc(Nsph_zone[i][j][k][n] * sizeof(double));
      dlam[i][j][k][n]  = safe_malloc(Nsph_zone[i][j][k][n] * sizeof(double));
      theta[i][j][k][n] = safe_malloc(Nsph_zone[i][j][k][n] * sizeof(double));
      Nsph_cntd[i][j][k][n] = 0;
    }
  } // omp parallel

// Create per-zone lists of w, nu, dlam
#pragma omp parallel
  {
    int               n  = omp_get_thread_num();
    struct of_photon *ph = photon_lists[n];
    while (ph != NULL) {
      int    i, j, k;
      double X[NDIM], Kcov[NDIM], Kcon[NDIM];
      get_X_K_interp(ph, t, P, X, Kcov, Kcon);
      Xtoijk(X, &i, &j, &k);

      if (i < NG || i > NG + N1 - 1 || j < NG || j > NG + N2 - 1 || k < NG ||
          k > NG + N3 - 1) {
        printf("BAD PH????\n");
        printf("[%i] %i %i %i X[] = %e %e %e %e\n", n, i, j, k, X[0], X[1],
            X[2], X[3]);
        for (int mu = 0; mu < 3; mu++) {
          printf("X[%i][] = %e %e %e %e\n", mu, ph->X[mu][0], ph->X[mu][1],
              ph->X[mu][2], ph->X[mu][3]);
          printf("Kcov[%i][] = %e %e %e %e\n", mu, ph->Kcov[mu][0],
              ph->Kcov[mu][1], ph->Kcov[mu][2], ph->Kcov[mu][3]);
          printf("Kcon[%i][] = %e %e %e %e\n", mu, ph->Kcon[mu][0],
              ph->Kcon[mu][1], ph->Kcon[mu][2], ph->Kcon[mu][3]);
        }
        printf("origin %i %i %i %i\n", ph->origin[0], ph->origin[1],
            ph->origin[2], ph->origin[3]);
        printf("nscatt = %i\n", ph->nscatt);
        exit(-1);
      }

      double freq = 0.;
      for (int mu = 0; mu < NDIM; mu++) {
        freq -= Ucon[i][j][k][mu] * Kcov[mu];
      }
      freq *= ME * CL * CL / HPL;

      w[i][j][k][n][Nsph_cntd[i][j][k][n]]    = ph->w;
      nu[i][j][k][n][Nsph_cntd[i][j][k][n]]   = freq;
      dlam[i][j][k][n][Nsph_cntd[i][j][k][n]] = dt / Kcon[0];
      theta[i][j][k][n][Nsph_cntd[i][j][k][n]] =
          get_bk_angle(X, Kcon, Ucov[i][j][k], Bcov[i][j][k], Bmag[i][j][k]);

      Nsph_cntd[i][j][k][n]++;
      ph = ph->next;
    }
  } // omp parallel

// In each zone, rootfind to Thetae that matches sources. If failure, simply
// return current Thetae
#pragma omp parallel
  {
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *           s;
    gsl_function                 F;
    F.function = &dEdt;
    T          = gsl_root_fsolver_brent;
    s          = gsl_root_fsolver_alloc(T);

    static double extra[EOS_NUM_EXTRA];
#pragma omp threadprivate(extra)
#if EOS == EOS_TYPE_GAMMA
    EOS_ELOOP { extra[e] = 0.0; }
#endif

#pragma omp for collapse(3)
    ZLOOP {
#if EOS == EOS_TYPE_TABLE
      extra[EOS_YE] = P[i][j][k][YE];
#endif

      // fill parameters
      struct of_params params;
      params.i       = i;
      params.j       = j;
      params.k       = k;
      params.rho     = P[i][j][k][RHO];
      params.Ne      = Ne[i][j][k];
      params.Bmag    = Bmag[i][j][k];
      params.Thetaei = Thetae_old[i][j][k];
      params.Ucon0   = Ucon[i][j][k][0];
      params.dt      = dt;
      // hopefully this doesn't break openmp
      EOS_ELOOP { params.extra[e] = extra[e]; }

      F.params         = &params;
      double r         = 0.;
      double Thetae_lo = 0.5 * Thetae_old[i][j][k];
      double Thetae_hi = 2. * Thetae_old[i][j][k];

      // Test interval for sanity
      double rmin = dEdt(Thetae_lo, &params);
      double rmax = dEdt(Thetae_hi, &params);
      if (rmin * rmax > 0.) {
        fprintf(stderr, "[%i %i %i] Root not bracketed!\n", i, j, k);
        Thetae_est[i][j][k] = Thetae_old[i][j][k];
        continue;
      }

      gsl_root_fsolver_set(s, &F, Thetae_lo, Thetae_hi);

      int iter = 0, max_iter = 100;
      int status;
      do {
        status    = gsl_root_fsolver_iterate(s);
        r         = gsl_root_fsolver_root(s);
        Thetae_lo = gsl_root_fsolver_x_lower(s);
        Thetae_hi = gsl_root_fsolver_x_upper(s);
        status    = gsl_root_test_interval(Thetae_lo, Thetae_hi, 0, 0.001);
      } while (status == GSL_CONTINUE && iter < max_iter);

      if (status != GSL_SUCCESS) {
        Thetae_est[i][j][k] = Thetae_old[i][j][k];
      } else {
        Thetae_est[i][j][k] = r;
      }
    }

    gsl_root_fsolver_free(s);
  } // omp parallel

// Clean up mallocs and reset counters
#pragma omp parallel for collapse(3)
  ZLOOP {
    for (int n = 0; n < nthreads; n++) {
      free(w[i][j][k][n]);
      free(nu[i][j][k][n]);
      free(dlam[i][j][k][n]);
      free(theta[i][j][k][n]);
      Nsph_zone[i][j][k][n] = 0.;
    }
  } // omp parallel
}

#endif // ESTIMATE_THETAE
#endif // RADIATION
