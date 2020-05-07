/******************************************************************************
 *                                                                            *
 * OPAC_EMIS_NEUTRINO.C                                                       *
 *                                                                            *
 * Reader for Tabulated Neutrino Opacities                                    *
 * First published in Burrows, Reddy, Thompson, arXiv:astro-ph/0404432, 2004  *
 *                                                                            *
 ******************************************************************************/

#include "constants.h"
#include "decs.h"
#include <hdf5.h>

#if RADIATION == RADTYPE_NEUTRINOS && BURROWS_OPACITIES
#include <hdf5_hl.h>
#ifdef __INTEL_COMPILER
#define FORT_OPAC_CALL(name) opacity_table_module_mp_##name##_
#else
#ifdef __GNUC__
#define FORT_OPAC_CALL(name) __opacity_table_module_MOD_##name
#endif
#endif
#define FORT_CALL(name) name##_

// needs following global variables defined
// opac_param_file
// opac_file

// These are set once and then left unchanged
static double lnu_min, lnu_max, dlnu;
static double mev_min, mev_max;
static double lnugroup[RAD_NUM_TYPES]
                      [NU_BINS + 1]; // frequency of groups in table
static double nugroup[RAD_NUM_TYPES][NU_BINS + 1];
static double egroup[RAD_NUM_TYPES][NU_BINS + 1]; // energy of nugroup in MeV
static int    ngrp[RAD_NUM_TYPES] = {NU_BINS + 1, NU_BINS + 1, NU_BINS + 1};
static int    ngroups             = 3 * (NU_BINS + 1);

void opac_emis(double rho, double T, double Ye, double *kappa, double *jg,
    double *sc, double *delta);
void test_opac_emis();

void FORT_OPAC_CALL(get_opacity_emissivity)(double *ab, double *sc,
    double *delta, double *eta, double *dab, double *dsc, double *ddelta,
    double *deta, double *rho, double *ye, double *temp, double *freq, int *ngr,
    int *ngroups, int *get_k, int *get_s, int *get_d, int *get_j, int *get_dk,
    int *get_ds, int *get_dd, int *get_de);
void FORT_CALL(init_opacity_table)(double *energy, int *ngrp, int *ntot,
    char *param, char *table, size_t plen, size_t tlen);

void init_opac_emis_burrows() {
  // initialize static variables and arrays
  lnu_min = log(numin);
  lnu_max = log(numax);
  mev_min = HPL * numin / MEV;
  mev_max = HPL * numax / MEV;
  dlnu    = (lnu_max - lnu_min) / NU_BINS;
  TYPELOOP {
    for (int n = 0; n < NU_BINS + 1; n++) {
      lnugroup[itp][n] = n * dlnu + lnu_min;
      nugroup[itp][n]  = exp(lnugroup[itp][n]);
      egroup[itp][n]   = HPL * nugroup[itp][n] / MEV;
    }
  }

  FORT_CALL(init_opacity_table)
  (&(egroup[0][0]), ngrp, &ngroups, opac_param_file, opac_file,
      strlen(opac_param_file), strlen(opac_file));

  return;
}

/*
 * given microphysics fills it with opacity and emisivity data
 */
void fill_opac_emis_burrows(struct of_microphysics *m) {
  const double rho = m->rho;       // should be in g/cm^3
  const double T   = (m->T) / MEV; // Should be in MeV
  const double Ye  = m->Ye;

  opac_emis(rho, T, Ye, &(m->alphagrp[0][0]), &(m->jgrp[0][0]), NULL, NULL);
}

/*
 * Given microphysics, a frequency (Hz) and a neutrino species,
 * returns emissivity in units of
 * (1/4pi) ergs s^-1 cm^-3 Hz^-1 = (1/4pi) erg / cm^3
 * corresponding to j_nu = dE/(dx^3 dt dOmega dnu)
 */
double jnu_burrows(double nu, int type, const struct of_microphysics *m) {
  if (type < 0) {
    fprintf(stderr, "[jnu_burrows]: Bad photon type!\n");
    exit(1);
  }

  const double lnu = log(fabs(nu));

  double j = interp_1d(lnu, lnu_min, lnu_max, 0, NU_BINS + 1,
      &(lnugroup[type][0]), &(m->jgrp[type][0]));

  return j;
}

double Jnu_burrows(double nu, int type, const struct of_microphysics *m) {
  return 4 * M_PI * jnu_burrows(nu, type, m);
}

/* Integrates emissivity over solid angle and frequency to return
 * total emissivity per cell
 */
double int_jnudnudOmega_burrows(const struct of_microphysics *m) {
  double Jtot = 0.0;
  TYPELOOP {
    for (int n = 0; n < NU_BINS + 1; n++) {
      Jtot += (m->jgrp[itp][n]) * nugroup[itp][n];
    }
  }
  Jtot *= 4 * M_PI * dlnu;
  return Jtot;
}

/*
 * Given microphysics, a frequency (Hz) and a neutrino species,
 * returns absorption coefficient in units of cm^-1
 */
double alpha_nu_burrows(double nu, int type, const struct of_microphysics *m) {
  if (type < 0) {
    fprintf(stderr, "[alpha_nu_burrows]: bad photon type!\n");
    exit(1);
  }
  const double lnu = log(fabs(nu));

  double alpha = interp_1d(lnu, lnu_min, lnu_max, 0, NU_BINS + 1,
      &(lnugroup[type][0]), &(m->alphagrp[type][0]));

  return alpha;
}

/*
 * Given rho, T, Ye,
 * fills opacities kappa, emissivities, jg,
 * scattering opacities sc, and scattering anisotropy delta
 * in arrays that are of shape
 * RAD_NUM_TYPES x NU_BINS+1
 * Where the electron neutrinos are first, then antis, then heavies
 */
void opac_emis(double rho, double T, double Ye,
    double *kappa, // opacities    (1/cm)
    double *jg,    // emissivities (erg s^{-1} cm^{-3} Hz^{-1} / 4pi)
    double *sc,    // scattering opacities
    double *delta  // scattering anisotropy
) {

  int get_kappa = (kappa == NULL ? 0 : 1);
  int get_emis  = (jg == NULL ? 0 : 1);
  int get_scat  = (sc == NULL ? 0 : 1);
  int get_delta = (delta == NULL ? 0 : 1);
  int get_dk    = 0;
  int get_de    = 0;
  int get_ds    = 0;
  int get_dd    = 0;

  double *dlkdlnu = NULL;
  double *dljdlnu = NULL;
  double *dlsdlnu = NULL;
  double *dlddlnu = NULL;

  FORT_OPAC_CALL(get_opacity_emissivity)
  (kappa, sc, delta, jg, dlkdlnu, dlsdlnu, dlddlnu, dljdlnu, &rho, &Ye, &T,
      &(egroup[0][0]), ngrp, &ngroups, &get_kappa, &get_scat, &get_delta,
      &get_emis, &get_dk, &get_ds, &get_dd, &get_de);
  if (get_emis) {
    // convert from
    // (erg s^{-1} cm^{-3} MeV^{-1} / 4pi)
    // to
    // (erg s^{-1} cm^{-3} Hz^{-1} / 4pi)
    for (int n = 0; n < ngroups; n++) {
      jg[n] *= HPL / MEV;
    }
  }

  return;
}

#if EOS == EOS_TYPE_TABLE
void opac_emis_to_hdf(const char *name, const struct of_tablebounds *b) {

  // Derivatives
  const double dlrho = (b->lrho_max - b->lrho_min) / (b->Nrho - 1);
  const double dlT   = (b->lT_max - b->lT_min) / (b->NT - 1);
  const double dYe   = (b->Ye_max - b->Ye_min) / (b->NYe - 1);

  int     size = (b->Nrho) * (b->NT) * (b->NYe) * RAD_NUM_TYPES * (NU_BINS + 1);
  double *tab_lrho = safe_malloc((b->Nrho) * sizeof(double));
  double *tab_lT   = safe_malloc((b->NT) * sizeof(double));
  double *tab_Ye   = safe_malloc((b->NYe) * sizeof(double));
  if (mpi_io_proc()) {
    printf("lrho, lt, ye malloced\n");
  }
  double *jg = safe_malloc(size * sizeof(double));
  double *kg = safe_malloc(size * sizeof(double));
  if (mpi_io_proc()) {
    printf("jg and kg malloced\n");
  }
  for (int irho = 0; irho < b->Nrho; irho++) {
    double lrho    = b->lrho_min + irho * (dlrho);
    double rho     = pow(10., lrho);
    tab_lrho[irho] = lrho;
    for (int iT = 0; iT < b->NT; iT++) {
      double lT  = b->lT_min + iT * (dlT);
      double T   = pow(10., lT);
      tab_lT[iT] = lT;
      for (int iY = 0; iY < b->NYe; iY++) {
        double Ye  = b->Ye_min + iY * (dYe);
        tab_Ye[iY] = Ye;
        int e      = (b->Nrho * ((iY)*b->NT + (iT)) + (irho)) * RAD_NUM_TYPES *
                (NU_BINS + 1);
        opac_emis(rho, T, Ye, &(kg[e]), &(jg[e]), NULL, NULL);
      }
    }
  }
  if (mpi_io_proc()) {
    hid_t   file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[]  = {b->NYe, b->NT, b->Nrho, RAD_NUM_TYPES, NU_BINS + 1};

    {
      hid_t group =
          H5Gcreate(file_id, "units", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5LTset_attribute_string(file_id, "units", "density", "log_{10}(g/cm^3)");
      H5LTset_attribute_string(
          file_id, "units", "temperature", "log_{10}(MeV)");
      H5LTset_attribute_string(file_id, "units", "Ye", "N/A");
      H5LTset_attribute_string(file_id, "units", "frequency", "log_e(Hz)");
      H5LTset_attribute_string(file_id, "units", "emissivity", "cgs");
      H5LTset_attribute_string(file_id, "units", "opacity", "cgs");
      H5Gclose(group);
    }

    {
      hid_t group = H5Gcreate(
          file_id, "dimensions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      int rad_num_types = RAD_NUM_TYPES;
      int num_nu        = NU_BINS + 1;
      H5LTset_attribute_int(file_id, "dimensions", "numRho", &b->Nrho, 1);
      H5LTset_attribute_int(file_id, "dimensions", "numT", &b->NT, 1);
      H5LTset_attribute_int(file_id, "dimensions", "numYe", &b->NYe, 1);
      H5LTset_attribute_int(
          file_id, "dimensions", "numRadTypes", &rad_num_types, 1);
      H5LTset_attribute_int(file_id, "dimensions", "numNu", &num_nu, 1);
      H5LTset_attribute_string(
          file_id, "dimensions", "index order", "rho,T,Ye,type,nu");
      H5Gclose(group);
    }

    H5LTmake_dataset_double(file_id, "lrho", 1, &(dims[2]), tab_lrho);
    H5LTmake_dataset_double(file_id, "lT", 1, &(dims[1]), tab_lT);
    H5LTmake_dataset_double(file_id, "Ye", 1, &(dims[0]), tab_Ye);
    H5LTmake_dataset_double(file_id, "lnu", 1, &(dims[4]), &(lnugroup[0][0]));
    H5LTmake_dataset_double(file_id, "emis", 5, dims, jg);
    H5LTmake_dataset_double(file_id, "opac", 5, dims, kg);
    H5Fclose(file_id);
  }
  free(jg);
  free(kg);
  free(tab_lrho);
  free(tab_lT);
  free(tab_Ye);
}
#endif

void test_opac_emis() {

  int     g, one, zero;
  int     nr1 = NU_BINS + 1;
  double *jg, *kappa;

  one  = 1;
  zero = 0;

  double rho0 = 8.e13;
  double T0   = 8.;
  double ymin = 0.035;
  double ymax = 0.2;
  int    ny   = (0.2 - 0.035) * 29. / (0.56 - 0.035);

  jg    = safe_malloc(ngroups * sizeof(double));
  kappa = safe_malloc(ngroups * sizeof(double));

  fprintf(stdout, "\n\n");
  for (int g = 0; g < nr1; g++) {
    fprintf(stdout, "%g ", egroup[0][g]);
  }
  fprintf(stdout, "\n\n");

  for (int iy = 0; iy <= ny; iy++) {
    double ye = ymin + iy * (ymax - ymin) / ny;
    FORT_OPAC_CALL(get_opacity_emissivity)
    (kappa, NULL, NULL, jg, NULL, NULL, NULL, NULL, &rho0, &ye, &T0,
        &(egroup[0][0]), ngrp, &ngroups, &one, &zero, &zero, &one, &zero, &zero,
        &zero, &zero);
    fprintf(stdout, "%g\n", ye);
    for (g = 0; g < nr1; g++) {
      fprintf(stdout, "%g ", jg[g] / (kappa[g] * pow(egroup[0][g], 3)));
    }
    fprintf(stdout, "\n\n");
    sleep(1);
    for (g = 1; g < nr1; g++) {
      fprintf(stdout, "%g ", jg[g] - jg[g - 1]);
    }
    fprintf(stdout, "\n\n");
    sleep(1);
    for (g = 1; g < nr1; g++) {
      fprintf(stdout, "%g ", kappa[g] - kappa[g - 1]);
    }
    fprintf(stdout, "\n\n");
    sleep(1);
  }

  free(jg);
  free(kappa);
  return;
}

#endif // BURROWS_OPACITIES
