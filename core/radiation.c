/******************************************************************************
 *                                                                            *
 * RADIATION.C                                                                *
 *                                                                            *
 * MODEL-INDEPENDENT RADIATION QUANTITIES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
double Bnu_inv(double nu,
    const struct of_microphysics
        *m) { // unused for neutrinos. Use tables by Burrows et al.
#if RADIATION == RADTYPE_LIGHT
  double x;
  x = HPL * nu / (ME * CL * CL * m->Thetae);

  if (x < 1.e-3) // Taylor-expand small arguments for numerical accuracy
    return (
        (2. * HPL / (CL * CL)) / (x / 24. * (24. + x * (12. + x * (4. + x)))));
  else
    return ((2. * HPL / (CL * CL)) / (exp(x) - 1.));
#endif
#if RADIATION == RADTYPE_NEUTRINOS
  return 0.; // STUB
#endif
}

double jnu_inv(
    double nu, int type, const struct of_microphysics *m, double theta) {
  double j;

  j = jnu(nu, type, m, theta);

  return (j / (nu * nu));
}

// Invariant scattering opacity
double alpha_inv_scatt(
    double nu, int type, int interaction, const struct of_microphysics *m) {
  double Eg    = HPL * nu / (ME * CL * CL);
  double ne    = scatterer_number_density(type, interaction, m);
  double sigma = total_cross_lkup(Eg, type, interaction, m);

  return nu * sigma * ne;
}

// Invariant absorption opacity
double alpha_inv_abs(
    double nu, int type, const struct of_microphysics *m, double theta) {
#if RADIATION == RADTYPE_NEUTRINOS
  {
#if BURROWS_OPACITIES
    {
      double alpha = alpha_nu_burrows(nu, type, m);
      return alpha * nu; // units of cm^-1 s^-1
    }
#elif HDF5_OPACITIES
    {
      double alpha = alpha_nu_hdf(nu, type, m);
      return alpha * nu; // units of cm^-1 s^-1
    }
#else
    {
      return 0.; // STUB
    }
#endif
  }
#else
  {
    double jnu_invariant = jnu_inv(nu, type, m, theta);

    // Otherwise synchrotron emissivity will fail in this region
    if (isnan(jnu_invariant)) {
      return 0.;
    } else {
      return (jnu_invariant / (Bnu_inv(nu, m) + 1.e-100));
    }
  }
#endif
}

// Get frequency in fluid frame, in Hz
double get_fluid_nu(double X[4], double Kcov[4], double Ucon[NDIM]) {
  double ener, nu;

  // Energy in electron rest-mass units
  ener = -(Kcov[0] * Ucon[0] + Kcov[1] * Ucon[1] + Kcov[2] * Ucon[2] +
           Kcov[3] * Ucon[3]);

  nu = ener * ME * CL * CL / HPL;

  return nu;
}

// Return angle between magnetic field and wavevector
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
    double Bcov[NDIM], double B) {
  double k, mu;

  if (B == 0.)
    return (M_PI / 2.);

  k = fabs(K[0] * Ucov[0] + K[1] * Ucov[1] + K[2] * Ucov[2] + K[3] * Ucov[3]);

  mu = (K[0] * Bcov[0] + K[1] * Bcov[1] + K[2] * Bcov[2] + K[3] * Bcov[3]) /
       fabs(k * B);

#define SMALLEPS (1.e-10)

  if (fabs(mu) > 1.) {
    mu /= fabs(mu);

    // We need |mu| < 1 or else acos(mu) will return nan.
    if (mu > 0.)
      mu -= SMALLEPS;
    else
      mu += SMALLEPS;
  }

  return (acos(mu));
}
#endif // RADIATION
