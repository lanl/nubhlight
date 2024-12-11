/******************************************************************************
 *                                                                            *
 * SCATTERING.C                                                               *
 *                                                                            *
 * RELATIVISTIC SCATTERING KERNEL                                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION

// RADIATION HOT CROSS SECTION
#define NW 220
#define NT 90
#define MINW (1.e-16) // w = h nu/m_e c^2
// #define MAXW (1.e10)
#define MINT (0.001) // Theta = k_b T / m c^2, m = mass of scatterer
#define MAXT (1.e11)
#define MAXGAMMA (12.)
#define DMUE (0.05)
#define DGAMMAE (0.05)
#define N_HC_TABLES (1)
#define TRUNCATE_HOTCROSS (1)
#if RADIATION == RADTYPE_NEUTRINOS
#if MULTISCATT_TEST
const char *HOTCROSS[N_HC_TABLES] = {"hotcross_multiscatt.dat"};
#else
const char *HOTCROSS[N_HC_TABLES] = {"hotcross_quad.dat"};
#endif
#else
const char *HOTCROSS[N_HC_TABLES] = {"hotcross.dat"};
#endif
double table[N_HC_TABLES][NW + 1][NT + 1];
double MAXW, dlw, dlT, lminw, lmint;

void   sample_gas_particle(double Ktetrad[NDIM], double Pelectron[NDIM],
      const struct of_microphysics *m, int type, int interaction);
void   sample_beta(double Thetae, double *gamma_e, double *beta_e);
double sample_y(double Thetae);
double sample_mu(double beta_e);
double get_total_cross_section(double k, const struct of_microphysics *m,
    int type, int interaction, int normalized);
void   sample_scattered_rad(double k[NDIM], double p[NDIM], double kp[NDIM],
      const struct of_microphysics *m, int type, int interaction);
void   boost(double v[NDIM], double u[NDIM], double vp[NDIM]);
void   sample_cross_section(double k, double *k0p, double *cth,
      const struct of_microphysics *m, int type, int interaction);
double total_cross_num(hc_ftype f, double w, double Thetae);
double dNdgammae(double thetae, double gammae);
double boostcross(hc_ftype f, double w, double mue, double gammae);
void   init_hc_table(
      hc_ftype f, double table[NW + 1][NT + 1], const char *hc_name);
double interpolate_hc_table(double w, double theta, hc_ftype f,
    double table[NW + 1][NT + 1], const char *hc_name);

// rejection sampling
void rejection_sample(dsdom_ftype f, double xmax, double k, double *k0p,
    double *cth, const struct of_microphysics *m, int type, int interaction);
// compton scattering
double sample_klein_nishina(double k0);
double klein_nishina(double a, double ap);
double hc_klein_nishina(double we, double mue);
// Neutrino physics
double hc_quad(double we, double mue);
double hc_quad_max();
double nu_cross_factor(
    double sigma, int type, int interaction, const struct of_microphysics *m);
double total_cross_ions(double sigma_hc, double A, double Z);
double nu_cross_delta(int type, int interaction);
double nu_cross(double w, double mu, const struct of_microphysics *m, int type,
    int interaction);
double nu_cross_max(const struct of_microphysics *m, int type, int interaction);
// multiscatt test
#if MULTISCATT_TEST
double hc_flat(double we, double mue);
double ms_flat(double w, double mu, const struct of_microphysics *m, int type,
    int interaction);
double ms_flat_max(const struct of_microphysics *m, int type, int interaction);
#endif

// Scattering temperature too small
int scatt_temp_too_small(const struct of_microphysics *m) {
#if MULTISCATT_TEST
  return 0;
#else // normal scattering
#if RADIATION == RADTYPE_LIGHT
  { return m->Thetae < 10. * SMALL; }
#elif RADIATION == RADTYPE_NEUTRINOS
  { return (KBOL * m->T) / (ME * CL * CL) < 10. * SMALL; }
#endif // neutrinos
#endif // Not multiscatt test
}

int scatter_superphoton(grid_prim_type P, grid_eosvar_type extra,
    struct of_photon *ph, double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
    int interaction) {
  const int fail    = 0;
  const int success = 1;

  double Pelectron[NDIM], gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];

  double Ktetrad[NDIM], Ktetrad_scatt[NDIM];

  int i, j, k;
  int bad_scatter = 0;
  Xtoijk(X, &i, &j, &k);

  set_gcov(X, gcov);
  gcon_func(gcov, gcon);

  normalize_null(gcov, Kcon);
  normalize_null_cov(gcon, Kcov);

  // Quality control
  if (ph->type == TYPE_TRACER)
    return fail;
  if (Kcon[0] < 0. || Kcov[0] > 0.)
    return fail;
  DLOOP1 {
    if (is_practically_nan(Kcon[mu]) || is_practically_nan(Kcov[mu])) {
      return fail;
    }
  }

  double                 Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
  struct of_microphysics m;
  get_fluid_zone(i, j, k, P, extra, &m, Ucon, Ucov, Bcon, Bcov);

  if (scatt_temp_too_small(&m))
    return fail;

  make_tetrad(i, j, k, Ucon, Bcon, gcov, Econ, Ecov);

  coord_to_tetrad(Ecov, Kcon, Ktetrad);

  DLOOP1 {
    if (is_practically_nan(Ktetrad[mu])) {
      bad_scatter = 1;
    }
  }
  if (bad_scatter) {
    fprintf(stderr,
        "Bad tetrad!\n"
        "\tUcon    = [%e %e %e %e]\n"
        "\tUcov    = [%e %e %e %e]\n"
        "\tBcon    = [%e %e %e %e]\n"
        "\tBcov    = [%e %e %e %e]\n"
        "\tKtetrad = [%e %e %e %e]\n"
        "\tKcoord  = [%e %e %e %e]\n"
        "\n"
        "\tEcon    = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n"
        "\n"
        "\tEcov    = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n"
        "\t        = [%e %e %e %e]\n",
        Ucon[0], Ucon[1], Ucon[2], Ucon[3], Ucov[0], Ucov[1], Ucov[2], Ucov[3],
        Bcon[0], Bcon[1], Bcon[2], Bcon[3], Bcov[0], Bcov[1], Bcov[2], Bcov[3],
        Ktetrad[0], Ktetrad[1], Ktetrad[2], Ktetrad[3], Kcon[0], Kcon[1],
        Kcon[2], Kcon[3], Econ[0][0], Econ[0][1], Econ[0][2], Econ[0][3],
        Econ[1][0], Econ[1][1], Econ[1][2], Econ[1][3], Econ[2][0], Econ[2][1],
        Econ[2][2], Econ[2][3], Econ[3][0], Econ[3][1], Econ[3][2], Econ[3][3],
        Ecov[0][0], Ecov[0][1], Ecov[0][2], Ecov[0][3], Ecov[1][0], Ecov[1][1],
        Ecov[1][2], Ecov[1][3], Ecov[2][0], Ecov[2][1], Ecov[2][2], Ecov[2][3],
        Ecov[3][0], Ecov[3][1], Ecov[3][2], Ecov[3][3]);
    return fail;
  }

  sample_gas_particle(Ktetrad, Pelectron, &m, ph->type, interaction);

  DLOOP1 {
    if (is_practically_nan(Pelectron[mu])) {
#if RADIATION == RADTYPE_LIGHT
      printf("m.Thetae = %e m.Ne = %e m.B = %e\n", m.Thetae, m.Ne, m.B);
#endif
      printf("Pelectron[%i] = %e!\n", mu, Pelectron[mu]);
      printf(
          "K = %e %e %e %e\n", Ktetrad[0], Ktetrad[1], Ktetrad[2], Ktetrad[3]);
      printf("k.k = %e\n", -Ktetrad[0] * Ktetrad[0] + Ktetrad[1] * Ktetrad[1] +
                               Ktetrad[2] * Ktetrad[2] +
                               Ktetrad[3] * Ktetrad[3]);
      return fail;
    }
  }

  sample_scattered_rad(
      Ktetrad, Pelectron, Ktetrad_scatt, &m, ph->type, interaction);

  DLOOP1 {
    if (is_practically_nan(Ktetrad_scatt[mu])) {
      // printf("Ktetrad_scatt[%i] = %e!\n", mu, Ktetrad_scatt[mu]);
      bad_scatter = 1;
    }
  }

  // If scatter fails, return failure and complain
  if (bad_scatter) {
    fprintf(stderr,
        "Bad scatter after sample_scattered_rad:\n"
        "\tph->Kcon[2] is bad\n"
        "\tUcon    = [%e %e %e %e]\n"
        "\tUcov    = [%e %e %e %e]\n"
        "\tBcon    = [%e %e %e %e]\n"
        "\tBcov    = [%e %e %e %e]\n"
        "\tKtetrad = [%e %e %e %e]\n"
        "\tKcon    = [%e %e %e %e]\n",
        Ucon[0], Ucon[1], Ucon[2], Ucon[3], Ucov[0], Ucov[1], Ucov[2], Ucov[3],
        Bcon[0], Bcon[1], Bcon[2], Bcon[3], Bcov[0], Bcov[1], Bcov[2], Bcov[3],
        Ktetrad[0], Ktetrad[1], Ktetrad[2], Ktetrad[3], ph->Kcon[2][0],
        ph->Kcon[2][1], ph->Kcon[2][2], ph->Kcon[2][3]);
    return fail;
  }

  // Check NAN after each of these

  tetrad_to_coord(Econ, Ktetrad_scatt, ph->Kcon[2]);

  DLOOP1 {
    if (is_practically_nan(ph->Kcon[2][mu])) {
      // fprintf(stderr,"ph->Kcon[2][%i] = %e!\n", mu, ph->Kcon[2][mu]);
      bad_scatter = 1;
    }
  }

  if (bad_scatter) {
    fprintf(stderr,
        "Bad scatter after tetrad_to_coord:\n"
        "\tph->Kcon[2] is bad\n"
        "\ttype        = %d\n"
        "\tinteraction = %d\n"
        "\tUcon    = [%e %e %e %e]\n"
        "\tUcov    = [%e %e %e %e]\n"
        "\tBcon    = [%e %e %e %e]\n"
        "\tBcov    = [%e %e %e %e]\n"
        "\tKtetrad = [%e %e %e %e]\n"
        "\tKcon    = [%e %e %e %e]\n",
        ph->type, interaction, Ucon[0], Ucon[1], Ucon[2], Ucon[3], Ucov[0],
        Ucov[1], Ucov[2], Ucov[3], Bcon[0], Bcon[1], Bcon[2], Bcon[3], Bcov[0],
        Bcov[1], Bcov[2], Bcov[3], Ktetrad[0], Ktetrad[1], Ktetrad[2],
        Ktetrad[3], ph->Kcon[2][0], ph->Kcon[2][1], ph->Kcon[2][2],
        ph->Kcon[2][3]);
    return fail;
  }

  // Ensure scattered superphoton is sane
  normalize_null(gcov, ph->Kcon[2]);

  DLOOP1 {
    if (is_practically_nan(ph->Kcon[2][mu])) {
      // fprintf(stderr,"after norm ph->Kcon[2][%i] = %e!\n", mu,
      // ph->Kcon[2][mu]);
      bad_scatter = 1;
    }
  }

  if (bad_scatter) {
    fprintf(stderr,
        "Bad scatter after normalize_Null:\n"
        "\tph-Kcon[2] is bad\n"
        "\ttype        = %d\n"
        "\tinteraction = %d\n"
        "\tUcon    = [%e %e %e %e]\n"
        "\tUcov    = [%e %e %e %e]\n"
        "\tBcon    = [%e %e %e %e]\n"
        "\tBcov    = [%e %e %e %e]\n"
        "\tKtetrad = [%e %e %e %e]\n"
        "\tKcon    = [%e %e %e %e]\n",
        ph->type, interaction, Ucon[0], Ucon[1], Ucon[2], Ucon[3], Ucov[0],
        Ucov[1], Ucov[2], Ucov[3], Bcon[0], Bcon[1], Bcon[2], Bcon[3], Bcov[0],
        Bcov[1], Bcov[2], Bcov[3], Ktetrad[0], Ktetrad[1], Ktetrad[2],
        Ktetrad[3], ph->Kcon[2][0], ph->Kcon[2][1], ph->Kcon[2][2],
        ph->Kcon[2][3]);
    return fail;
  }

  lower(ph->Kcon[2], gcov, ph->Kcov[2]);

  DLOOP1 {
    if (is_practically_nan(ph->Kcov[2][mu])) {
      // fprintf(stderr,"after lower ph->Kcov[2][%i] = %e!\n", mu,
      // ph->Kcov[2][mu]);
      bad_scatter = 1;
    }
  }

  if (bad_scatter) {
    fprintf(stderr,
        "Bad scatter after lower Kcon:\n"
        "\tph-Kcov[2] is bad\n"
        "\ttype        = %d\n"
        "\tinteraction = %d\n"
        "\tUcon    = [%e %e %e %e]\n"
        "\tUcov    = [%e %e %e %e]\n"
        "\tBcon    = [%e %e %e %e]\n"
        "\tBcov    = [%e %e %e %e]\n"
        "\tKcon    = [%e %e %e %e]\n"
        "\tKtetrad = [%e %e %e %e]\n"
        "\tKcov    = [%e %e %e %e]\n",
        ph->type, interaction, Ucon[0], Ucon[1], Ucon[2], Ucon[3], Ucov[0],
        Ucov[1], Ucov[2], Ucov[3], Bcon[0], Bcon[1], Bcon[2], Bcon[3], Bcov[0],
        Bcov[1], Bcov[2], Bcov[3], Ktetrad[0], Ktetrad[1], Ktetrad[2],
        Ktetrad[3], ph->Kcon[2][0], ph->Kcon[2][1], ph->Kcon[2][2],
        ph->Kcon[2][3], ph->Kcov[2][0], ph->Kcov[2][1], ph->Kcov[2][2],
        ph->Kcov[2][3]);
    return fail;
  }

  // Ensure scattered superphoton is sane
  // normalize_null(gcov, ph->Kcon[2]);
  // normalize_null_cov(gcon, ph->Kcov[2]);

  /*if (ph->Kcov[2][0] > 0.) {
    printf("Kcov[0] > 0 after scattering!\n");
    printf("ph->X[] = %e %e %e %e\n", ph->X[2][0], ph->X[2][1], ph->X[2][2],
      ph->X[2][3]);
    printf("ph->Kcon[] = %e %e %e %e\n", ph->Kcon[2][0], ph->Kcon[2][1],
      ph->Kcon[2][2], ph->Kcon[2][3]);
    printf("ph->Kcov[] = %e %e %e %e\n", ph->Kcov[2][0], ph->Kcov[2][1],
      ph->Kcov[2][2], ph->Kcov[2][3]);
  }*/

  return success;
}

// Procedure from Canfield et al. 1987
void sample_gas_particle(double k[NDIM], double p[NDIM],
    const struct of_microphysics *m, int type, int interaction) {
  double beta_e, mu, phi, cphi, sphi, gamma_e, sigma;
  double K, sth, cth, x1, n0dotv0, v0, v1;
  double n0x, n0y, n0z;
  double v0x, v0y, v0z;
  double v1x, v1y, v1z;
  double v2x, v2y, v2z;
  int    sample_cnt = 0;
  double Thetae     = scatterer_dimensionless_temp(type, interaction, m);
  double factor     = 1.0;

  do {
    sample_beta(Thetae, &gamma_e, &beta_e);
    mu = sample_mu(beta_e);

    // Sometimes |mu| > 1 from roundoff error. Fix it
    if (mu > 1.)
      mu = 1.;
    else if (mu < -1.)
      mu = -1;

    // Frequency in electron rest frame
    K = gamma_e * (1. - beta_e * mu) * k[0];

    sigma = get_total_cross_section(K, m, type, interaction, 1);

    x1 = factor * get_rand();

    sample_cnt++;

    if (sample_cnt > 1000000) {
      fprintf(stderr,
          "in sample_gas_particle:\n"
          "\t type, int,  mu, gamma_e, K, sigma, x1, factor:\n"
          "\t%d %d %g %g %g %g %g %g %g\n",
          type, interaction, Thetae, mu, gamma_e, K, sigma, x1, factor);

      // Kluge to prevent stalling for large values of \Theta_e
      // TODO: does this work?
      // Thetae *= 0.5 ;
      factor *= 0.5;
      sample_cnt = 0;
    }
  } while (x1 >= sigma);

  // First unit vector for coordinate system
  v0x = k[1];
  v0y = k[2];
  v0z = k[3];
  v0  = sqrt(v0x * v0x + v0y * v0y + v0z * v0z);
  v0x /= v0;
  v0y /= v0;
  v0z /= v0;

  // Pick zero-angle for coordinate system
  get_ran_dir_3d(&n0x, &n0y, &n0z);
  n0dotv0 = v0x * n0x + v0y * n0y + v0z * n0z;

  // Second unit vector
  v1x = n0x - (n0dotv0)*v0x;
  v1y = n0y - (n0dotv0)*v0y;
  v1z = n0z - (n0dotv0)*v0z;

  // Normalize
  v1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
  v1x /= v1;
  v1y /= v1;
  v1z /= v1;

  // Find one more unit vector using cross product; automatically normalized
  v2x = v0y * v1z - v0z * v1y;
  v2y = v0z * v1x - v0x * v1z;
  v2z = v0x * v1y - v0y * v1x;

  // Resolve new momentum vector along unit vectors and create a four-vector p
  phi  = get_rand() * 2. * M_PI; // uniform orientation
  sphi = sin(phi);
  cphi = cos(phi);
  cth  = mu;
  sth  = sqrt(1. - mu * mu);

  p[0] = gamma_e;
  p[1] = gamma_e * beta_e * (cth * v0x + sth * (cphi * v1x + sphi * v2x));
  p[2] = gamma_e * beta_e * (cth * v0y + sth * (cphi * v1y + sphi * v2y));
  p[3] = gamma_e * beta_e * (cth * v0z + sth * (cphi * v1z + sphi * v2z));

  if (beta_e < 0) {
    fprintf(stderr, "betae error: %g %g %g %g\n", p[0], p[1], p[2], p[3]);
  }
}

void sample_beta(double Thetae, double *gamma_e, double *beta_e) {
  double y = sample_y(Thetae);

  *gamma_e = y * y * Thetae + 1.;
  *beta_e  = sqrt(1. - 1. / ((*gamma_e) * (*gamma_e)));
  *beta_e += SMALL; // to prevent numerical problems for zero temperature
}

double sample_y(double Thetae) {
  double S_3, pi_3, pi_4, pi_5, pi_6, y, x1, x2, x, prob, num, den;

  pi_3 = sqrt(M_PI) / 4.;
  pi_4 = sqrt(0.5 * Thetae) / 2.;
  pi_5 = 3. * sqrt(M_PI) * Thetae / 8.;
  pi_6 = Thetae * sqrt(0.5 * Thetae);

  S_3 = pi_3 + pi_4 + pi_5 + pi_6;

  pi_3 /= S_3;
  pi_4 /= S_3;
  pi_5 /= S_3;
  pi_6 /= S_3;

  int max_samp = 100000;
  int n        = 0;
  do {
    n++;
    x1 = get_rand();

    if (x1 < pi_3) {
      x = get_chisq(3);
    } else if (x1 < pi_3 + pi_4) {
      x = get_chisq(4);
    } else if (x1 < pi_3 + pi_4 + pi_5) {
      x = get_chisq(5);
    } else {
      x = get_chisq(6);
    }

    // Translate between Canfield et al. and standard chisq distribution
    y = sqrt(x / 2);

    x2  = get_rand();
    num = sqrt(1. + 0.5 * Thetae * y * y);
    den = 1. + y * sqrt(0.5 * Thetae);

    prob = num / den;

  } while (x2 >= prob && n < max_samp);

  if (n >= max_samp) {
    fprintf(stderr, "FAILED TO SAMPLE Y! Thetae = %e\n", Thetae);
    exit(-1);
  }

  return y;
}

double sample_mu(double beta_e) {
  double mu, x1;

  x1 = get_rand();
  mu = (1. - sqrt(1. + 2. * beta_e + beta_e * beta_e - 4. * beta_e * x1)) /
       beta_e;
  return mu;
}

/*
 * True total cross section
 * in rest frame of the scattering particle
 * Distinct from hot cross section,
 * which is looked up in total_cross_lkup
 */
double get_total_cross_section(double k, const struct of_microphysics *m,
    int type, int interaction, int normalized) {
#if RADIATION == RADTYPE_LIGHT
  {
    double sigma = hc_klein_nishina(k, 0);
    if (!normalized)
      sigma *= THOMSON;
    return sigma;
  }
#else // RADIATION == RADTYPE_NEUTRINOS
  {
#if MULTISCATT_TEST
    {
      double sigma = hc_flat(k, 0);
      if (!normalized) {
        int i = interaction;
        sigma *= (2 * 2 * i + 1) * NUSIGMA0 * pow(ms_theta_nu0, 2.0) / (4);
      }
      return sigma;
    }
#else  // Normal neutrino scattering
    { // TODO: doesn't work with electrons
      double sigma = hc_quad(k, 0);
      if (normalized) {
        return sigma / hc_quad_max();
      }
      return nu_cross_factor(sigma, type, interaction, m);
    }
#endif // multiscatt test?
  }
#endif // RADTYPE_NEUTRINOS
}

void sample_scattered_rad(double k[NDIM], double p[NDIM], double kp[NDIM],
    const struct of_microphysics *m, int type, int interaction) {
  double ke[4], kpe[4];
  double k0p;
  double n0x, n0y, n0z, n0dotv0, v0x, v0y, v0z, v1x, v1y, v1z, v2x, v2y, v2z,
      v1, dir1, dir2, dir3;
  double cth, sth, phi, cphi, sphi;

  boost(k, p, ke);

  sample_cross_section(ke[0], &k0p, &cth, m, type, interaction);
  sth = sqrt(fabs(1. - cth * cth));

  // Unit vector 1 for scattering coordinate system is oriented along initial
  // photon wavevector

  // Explicitly compute kemag instead of using ke[0] to ensure that photon is
  // created normalized and doesn't inherit the light cone errors from the
  // original photon
  double kemag = sqrt(ke[1] * ke[1] + ke[2] * ke[2] + ke[3] * ke[3]);
  v0x          = ke[1] / kemag;
  v0y          = ke[2] / kemag;
  v0z          = ke[3] / kemag;

  // Randomly pick zero-angle for scattering coordinate system.
  get_ran_dir_3d(&n0x, &n0y, &n0z);
  n0dotv0 = v0x * n0x + v0y * n0y + v0z * n0z;

  // Unit vector 2
  v1x = n0x - (n0dotv0)*v0x;
  v1y = n0y - (n0dotv0)*v0y;
  v1z = n0z - (n0dotv0)*v0z;
  v1  = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
  v1x /= v1;
  v1y /= v1;
  v1z /= v1;

  // Find one more unit vector using cross product; automatically normalized
  v2x = v0y * v1z - v0z * v1y;
  v2y = v0z * v1x - v0x * v1z;
  v2z = v0x * v1y - v0y * v1x;

  // Resolve new momentum vector along unit vectors
  // Create a four-vector p
  // Solve for orientation of scattered photon

  // Find phi for new photon
  phi  = 2. * M_PI * get_rand();
  sphi = sin(phi);
  cphi = cos(phi);

  p[1] *= -1.;
  p[2] *= -1.;
  p[3] *= -1.;

  dir1 = cth * v0x + sth * (cphi * v1x + sphi * v2x);
  dir2 = cth * v0y + sth * (cphi * v1y + sphi * v2y);
  dir3 = cth * v0z + sth * (cphi * v1z + sphi * v2z);

  kpe[0] = k0p;
  kpe[1] = k0p * dir1;
  kpe[2] = k0p * dir2;
  kpe[3] = k0p * dir3;

  // Transform k back to lab frame
  boost(kpe, p, kp);

  if (kp[0] < 0 || isnan(kp[0])) {
    fprintf(stderr, "in sample_scattered_photon:\n");
    fprintf(stderr, "kp[0], kpe[0]: %g %g\n", kp[0], kpe[0]);
    fprintf(stderr, "kpe: %g %g %g %g\n", kpe[0], kpe[1], kpe[2], kpe[3]);
    fprintf(stderr, "k:  %g %g %g %g\n", k[0], k[1], k[2], k[3]);
    fprintf(stderr, "p:   %g %g %g %g\n", p[0], p[1], p[2], p[3]);
    fprintf(stderr, "kp:  %g %g %g %g\n", kp[0], kp[1], kp[2], kp[3]);
  }
}

void boost(double v[NDIM], double u[NDIM], double vp[NDIM]) {
  double g, V, n1, n2, n3, gm1;

  g   = u[0];
  V   = sqrt(fabs(1. - 1. / (g * g)));
  n1  = u[1] / (g * V + SMALL);
  n2  = u[2] / (g * V + SMALL);
  n3  = u[3] / (g * V + SMALL);
  gm1 = g - 1.;

  // Lorentz boost into frame u from lab frame
  vp[0] = u[0] * v[0] - (u[1]) * v[1] - (u[2]) * v[2] - (u[3]) * v[3];
  vp[1] = -u[1] * v[0] + (1. + n1 * n1 * gm1) * v[1] + (n1 * n2 * gm1) * v[2] +
          (n1 * n3 * gm1) * v[3];
  vp[2] = -u[2] * v[0] + (n2 * n1 * gm1) * v[1] + (1. + n2 * n2 * gm1) * v[2] +
          (n2 * n3 * gm1) * v[3];
  vp[3] = -u[3] * v[0] + (n3 * n1 * gm1) * v[1] + (n3 * n2 * gm1) * v[2] +
          (1. + n3 * n3 * gm1) * v[3];
}

void sample_cross_section(double k, double *k0p, double *cth,
    const struct of_microphysics *m, int type, int interaction) {
#if RADIATION == RADTYPE_LIGHT
  {
    *k0p = sample_klein_nishina(k);
    *cth = 1. - 1. / (*k0p) + 1. / k;
  }
#elif RADIATION == RADTYPE_NEUTRINOS
  {
#if MULTISCATT_TEST
    {
      double xmax = ms_flat_max(m, type, interaction);
      rejection_sample(ms_flat, xmax, k, k0p, cth, m, type, interaction);
    }
#else  // Normal neutrino scattering
    {
      // TODO: doesn't work for electrons
      double xmax = nu_cross_max(m, type, interaction);
      rejection_sample(nu_cross, xmax, k, k0p, cth, m, type, interaction);
    }
#endif // neutrino scattering
  }
#endif // radiation type
}

double sample_klein_nishina(double k0) {
  double k0pmin, k0pmax, k0p_tent, x1;
  int    n = 0;

  // A low efficiency sampling algorithm, particularly for large k0. Limiting
  // efficiency is log(2 k0)/(2 k0)
  k0pmin = k0 / (1. + 2. * k0); // at theta = Pi
  k0pmax = k0;                  // at theta = 0

  do {
    // Tentative value
    k0p_tent = k0pmin + (k0pmax - k0pmin) * get_rand();

    // Rejection sample in box of height = kn(kmin)
    x1 = 2. * (1. + 2. * k0 + 2. * k0 * k0) / (k0 * k0 * (1. + 2. * k0));
    x1 *= get_rand();

    n++;
  } while (x1 >= klein_nishina(k0, k0p_tent));

  return k0p_tent;
}

double klein_nishina(double a, double ap) {
  double ch = 1. + 1. / a - 1. / ap;
  double kn = (a / ap + ap / a - 1. + ch * ch) / (a * a);

  return kn;
}

void init_all_hotcross() {
#if RADIATION == RADTYPE_LIGHT
  { init_hc_table(hc_klein_nishina, table[0], HOTCROSS[0]); }
#elif RADIATION == RADTYPE_NEUTRINOS
  {
#if MULTISCATT_TEST
    { init_hc_table(hc_flat, table[0], HOTCROSS[0]); }
#else  // Normal neutrino scattering
    { init_hc_table(hc_quad, table[0], HOTCROSS[0]); }
#endif // neutrino scattering type
  }
#endif // RADIATION TYPE
}

void init_hc_table(
    hc_ftype f, double table[NW + 1][NT + 1], const char *hc_name) {
  int    nread;
  double lw, lT;
  FILE * fp;

  MAXW  = 100. * (HPL * numax) / (ME * CL * CL);
  dlw   = log10(MAXW / MINW) / NW;
  dlT   = log10(MAXT / MINT) / NT;
  lminw = log10(MINW);
  lmint = log10(MINT);

  // Create file if needed using IO proc
  if (mpi_io_proc()) {
    fp = fopen(hc_name, "r");
    if (fp == NULL) {
      fprintf(stdout, "Making lookup table for %s cross section...\n", hc_name);
#pragma omp parallel for collapse(2)
      for (int i = 0; i <= NW; i++) {
        for (int j = 0; j <= NT; j++) {
          lw          = lminw + i * dlw;
          lT          = lmint + j * dlT;
          table[i][j] = log10(total_cross_num(f, pow(10., lw), pow(10., lT)));
          if (isnan(table[i][j])) {
            fprintf(stderr, "NAN for %s cross section: %d %d %g %g\n", hc_name,
                i, j, lw, lT);
            exit(0);
          }
        }
      }
      fprintf(stdout, "Lookup table created.\n\n");

      fprintf(stdout, "Writing lookup table to file...\n");
      fp = fopen(hc_name, "w");
      if (fp == NULL) {
        fprintf(stderr, "Couldn't write to file %s\n", hc_name);
        exit(0);
      }
      for (int i = 0; i <= NW; i++) {
        for (int j = 0; j <= NT; j++) {
          lw = lminw + i * dlw;
          lT = lmint + j * dlT;
          fprintf(fp, "%d %d %g %g %15.10g\n", i, j, lw, lT, table[i][j]);
        }
      }
      fprintf(stderr, "Lookup table written.\n\n");
    } // fp == NULL
    fclose(fp);
  } // mpi_io_proc()

  mpi_barrier();

  // Read lookup table with every MPI processor
  fp = fopen(hc_name, "r");
  if (fp == NULL) {
    fprintf(stderr, "rank %i: file %s not found.\n", mpi_myrank(), hc_name);
    exit(-1);
  }
  for (int i = 0; i <= NW; i++) {
    for (int j = 0; j <= NT; j++) {
      nread = fscanf(fp, "%*d %*d %*f %*f %lf\n", &table[i][j]);
      if (isnan(table[i][j]) || nread != 1) {
        fprintf(stderr, "Error on table %s read: %d %d\n", hc_name, i, j);
        exit(0);
      }
    }
  }
  fclose(fp);
}

double interpolate_hc_table(double w, double thetae, hc_ftype f,
    double table[NW + 1][NT + 1], const char *hc_name) {
  int    i, j;
  double lw, lT, di, dj, lcross;

// DEBUG
// double cross_section = RAD_SCATT_TYPES*NUSIGMA0;
// return 0.1/(Rout_rad*L_unit*cross_section*Ne_unit);

// if hotcross takes too long!
#if TRUNCATE_HOTCROSS
  if (w <= MINW)
    w = MINW + fabs(0.01 * MINW);
  if (w >= MAXW)
    w = MAXW - fabs(0.01 * MAXW);
  if (thetae <= MINT)
    thetae = MINT + fabs(0.01 * MINT);
  if (thetae >= MAXT)
    thetae = MAXT - fabs(0.01 * MAXT);
#endif

  // In-bounds for table
  if ((w > MINW && w < MAXW) && (thetae > MINT && thetae < MAXT)) {
    lw = log10(w);
    lT = log10(thetae);
    i  = (int)((lw - lminw) / dlw);
    j  = (int)((lT - lmint) / dlT);
    di = (lw - lminw) / dlw - i;
    dj = (lT - lmint) / dlT - j;

    lcross = (1. - di) * (1. - dj) * table[i][j] +
             di * (1. - dj) * table[i + 1][j] +
             (1. - di) * dj * table[i][j + 1] + di * dj * table[i + 1][j + 1];

    if (isnan(lcross)) {
      fprintf(stderr, "NAN cross section %s: %g %g %d %d %g %g\n", hc_name, lw,
          lT, i, j, di, dj);
    }

    return pow(10., lcross);
  }

  fprintf(stderr, "Cross section %s out of bounds: %g [%g,%g] %g [%g,%g]\n",
      hc_name, w, MINW, MAXW, thetae, MINT, MAXT);

  return total_cross_num(f, w, thetae);
}

double total_cross_lkup(
    double w, int type, int interaction, const struct of_microphysics *m) {
  double thetae = scatterer_dimensionless_temp(type, interaction, m);
#if RADIATION == RADTYPE_LIGHT
  {
    // Cold/low-energy: Use Thomson cross section
    if (w * thetae < 1.e-6)
      return (THOMSON);

    // Cold, but possible high-energy photon: use Klein-Nishina
    if (thetae < MINT)
      return (hc_klein_nishina(w, 0) * THOMSON);

    double sigma = interpolate_hc_table(
        w, thetae, hc_klein_nishina, table[0], HOTCROSS[0]);
    return THOMSON * sigma;
  }
#elif RADIATION == RADTYPE_NEUTRINOS
  {
    if (thetae < MINT) { // easy zero-temperature limit
      return get_total_cross_section(w, m, type, interaction, 0);
    }
#if MULTISCATT_TEST
    {
      double sigma =
          interpolate_hc_table(w, thetae, hc_flat, table[0], HOTCROSS[0]);
      return 4 * M_PI * sigma * ms_flat(w, 0, m, type, interaction) /
             hc_flat(w, 0);
    }
#else  // Normal neutrino scattering
    {
      if (nu_is_heavy(type) && interaction == RSCATT_TYPE_E) {
        return 0.0; // heavy cannot scatter off of electrons
      }
      double sigma =
          interpolate_hc_table(w, thetae, hc_quad, table[0], HOTCROSS[0]);
      return nu_cross_factor(sigma, type, interaction, m);
    }
#endif // neutrino scattering type
  }
#endif // radiation type
}

double total_cross_num(hc_ftype f, double w, double thetae) {
  double dmue, dgammae, mue, gammae, maxwell, cross;

  if (isnan(w)) {
    fprintf(stderr, "NAN cross section: %g %g\n", w, thetae);
    return 0.;
  }

// Check for easy limits
#if RADIATION == RADTYPE_LIGHT
  {
    if (thetae < MINT && w < MINW)
      return 1.;
    if (thetae < MINT)
      return hc_klein_nishina(w, 0);
  }
#endif

  dmue    = DMUE;
  dgammae = thetae * DGAMMAE;

  // Integrate over mu_e and gamma_e, where mu_e is the cosine of the angle
  // between K and U_e, and the angle k is assumed to lie, wlog, along the z
  // z axis
  cross = 0.;
  for (mue = -1. + 0.5 * dmue; mue < 1.; mue += dmue)
    for (gammae = 1. + 0.5 * dgammae; gammae < 1. + MAXGAMMA * thetae;
         gammae += dgammae) {
      maxwell = 0.5 * dNdgammae(thetae, gammae);

      cross += dmue * dgammae * boostcross(f, w, mue, gammae) * maxwell;

      if (isnan(cross)) {
        fprintf(stderr, "NAN cross section: %g %g %g %g %g %g\n", w, thetae,
            mue, gammae, dNdgammae(thetae, gammae),
            boostcross(f, w, mue, gammae));
      }
    }

  return cross;
}

// Normalized (per unit proper electron number density) electron distribution
double dNdgammae(double thetae, double gammae) {
  double K2f;

  if (thetae > 1.e-2) {
    K2f = gsl_sf_bessel_Kn(2, 1. / thetae) * exp(1. / thetae);
  } else {
    K2f = sqrt(M_PI * thetae / 2.) +
          15. / 8. * sqrt(M_PI / 2.) * pow(thetae, 1.5) +
          105. / 128. * sqrt(M_PI / 2.) * pow(thetae, 2.5) -
          315. / 1024. * sqrt(M_PI / 2.) * pow(thetae, 3.5);
  }

  return (gammae * sqrt(gammae * gammae - 1.) / (thetae * K2f)) *
         exp(-(gammae - 1.) / thetae);
}

double boostcross(hc_ftype f, double w, double mue, double gammae) {
  double we, boostcross, v;

  // Energy in electron rest frame
  v  = sqrt(gammae * gammae - 1.) / gammae;
  we = w * gammae * (1. - mue * v);

  boostcross = f(we, mue) * (1. - mue * v);

#if RADIATION == RADTYPE_LIGHT
  if (boostcross > 2) {
    fprintf(stderr,
        "w, mue, gammae: %g %g %g\n"
        "v, we, boostcross: %g %g %g\n"
        "kn: %g %g %g\n",
        w, mue, gammae, v, we, boostcross, v, we, boostcross);
    exit(1);
  }
#endif

  if (isnan(boostcross)) {
    fprintf(stderr, "isnan: %g %g %g\n", w, mue, gammae);
    return 0.;
  }

  return boostcross;
}

/*
 * TODO: for elastic scattering, we could
 *       inverse-transform sample for significant speedup
 */
void rejection_sample(dsdom_ftype f, double xmax, double k, double *k0p,
    double *cth, const struct of_microphysics *m, int type, int interaction) {
  double mu, x, sigma;
  do {
    mu    = 2 * get_rand() - 1.;
    x     = get_rand() * xmax;
    sigma = f(k, mu, m, type, interaction);
  } while (x >= sigma);

  *k0p = k; // because scattering is elastic
  *cth = mu;
}

double hc_klein_nishina(double we, double mue) {
  double sigma;

  if (we < 1.e-3) {
    sigma = 1. - 2. * we + 5.2 * we * we - 13.3 * we * we * we +
            1144 * we * we * we * we / 35.;
  } else {
    sigma = (3. / 4.) * (2. / (we * we) +
                            (1. / (2. * we) - (1. + we) / (we * we * we)) *
                                log(1. + 2. * we) +
                            (1. + we) / ((1. + 2. * we) * (1. + 2. * we)));
  }

  return sigma;
}

#if RADIATION == RADTYPE_NEUTRINOS
#if MULTISCATT_TEST
double hc_flat(double we, double mue) { return 1.0; }

double ms_flat(double w, double mu, const struct of_microphysics *m, int type,
    int interaction) {
  int    i     = interaction;
  double sigma = hc_flat(w, 0);
  sigma *= (2 * 2 * i + 1) * NUSIGMA0 * pow(ms_theta_nu0, 2.0) / (16 * M_PI);
  sigma *= (1 + ms_delta0 * pow(mu, 2 * i + 1));
  return sigma;
}

double ms_flat_max(const struct of_microphysics *m, int type, int interaction) {
  if (ms_delta0 > 0)
    return ms_flat(1, 1, m, type, interaction);
  if (ms_delta0 < 0)
    return ms_flat(1, -1, m, type, interaction);
  return ms_flat(1, 0, m, type, interaction);
}

#else  // not multiscatt test

double nu_cross(double w, double mu, const struct of_microphysics *m, int type,
    int interaction) {
  // TODO: doesn't work for electrons
  double sigma = hc_quad(w, mu);
  double delta = nu_cross_delta(type, interaction);
  sigma *= nu_cross_factor(sigma, type, interaction, m);
  sigma /= 4. * M_PI;
  sigma *= 1 + delta * mu;
  return sigma;
}

double nu_cross_max(
    const struct of_microphysics *m, int type, int interaction) {
  // TODO: doesn't work for electrons
  double sigma = hc_quad_max();
  double delta = nu_cross_delta(type, interaction);
  sigma        = nu_cross_factor(sigma, type, interaction, m);
  sigma /= 4. * M_PI;
  sigma *= 1 + fabs(delta);
  return sigma;
}

// Burrows, Reddy, Thomson, arXiv:astro-ph/0404432
double total_cross_ions(double sigma_hc, double A, double Z) {
  const double CFF  = 1.; // form factor
  const double CLOS = 1.; // electron polarization
  const double Sion = 1.; // ion screening
  double       W    = 1. - 2. * (Z / A) * (1. - 2. * S2THW);
  double out = ((1. / 16.) * NUSIGMA0 * sigma_hc * A * A * (W * CFF + CLOS) *
                (W * CFF + CLOS) * Sion);
  return out;
}

// Burrows, Reddy, Thomson, arXiv:astro-ph/0404432
double nu_cross_delta(int type, int interaction) {
  // heavy cannot scatter off of electrons
  if (nu_is_heavy(type) && interaction == RSCATT_TYPE_E)
    return 0.0;
  if (interaction == RSCATT_TYPE_P) {
    double Cpv = 0.5 + 2.0 * S2THW;
    double Cpa = 0.5;
    double num = (Cpv - 1) * (Cpv - 1) - GA2 * (Cpa - 1) * (Cpa - 1);
    double den = (Cpv - 1) * (Cpv - 1) + 3 * GA2 * (Cpa - 1) * (Cpa - 1);
    return num / den;
  } else if (interaction == RSCATT_TYPE_N) {
    return (1 - GA2) / (1 + 3 * GA2);
  } else if (interaction == RSCATT_TYPE_A || interaction == RSCATT_TYPE_ALPHA) {
    return 1.0;
  } else {
    fprintf(stderr, "total_cross_lkup: cross not implemented\n");
    exit(1);
  }
}

// Burrows, Reddy, Thomson, arXiv:astro-ph/0404432
double nu_cross_factor(
    double sigma, int type, int interaction, const struct of_microphysics *m) {
  // heavy cannot scatter off of electrons
  if (nu_is_heavy(type) && interaction == RSCATT_TYPE_E)
    return 0.0;
  if (interaction == RSCATT_TYPE_P) {
    return 0.25 * NUSIGMA0 * sigma *
           (4. * S4THW - 2. * S2THW + 0.25 * (1 + 3. * GA2));
  } else if (interaction == RSCATT_TYPE_N) {
    return (NUSIGMA0 / 4.) * sigma * (1 + 3 * GA2) / 4.;
  } else if (interaction == RSCATT_TYPE_A) {
    return total_cross_ions(sigma, m->Abar, m->Zbar);
  } else if (interaction == RSCATT_TYPE_ALPHA) {
    return total_cross_ions(sigma, 4.0, 2.0);
  } else { // TODO: add neutrino-electron scattering
    fprintf(stderr, "total_cross_lkup: cross not implemented\n");
    exit(1);
  }
}

double hc_quad(double we, double mue) { // we in units of hnu/Mec^2
  double nu_maxw = (HPL * numax) / (ME * CL * CL);
  if (we >= nu_maxw)
    we = nu_maxw;
  return we * we;
}

double hc_quad_max() {
  double nu_maxw = (HPL * numax) / (ME * CL * CL);
  return hc_quad(nu_maxw, 0);
}
#endif // MULTISCATT_TEST
#endif // Neutrinos

#endif // RADIATION
