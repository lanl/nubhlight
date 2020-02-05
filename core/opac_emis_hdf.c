/******************************************************************************
 *                                                                            *
 * OPAC_EMIS_HDF.C                                                            *
 *                                                                            *
 * Reader for Tabulated Neutrino Opacities                                    *
 * Tabulated format designed by Jonah Miller                                  *
 *                                                                            *
 ******************************************************************************/

// TODO
// Should this be re-designed to match nulib?

#include "decs.h"

#if (RADIATION == RADTYPE_NEUTRINOS) && HDF5_OPACITIES

#include "constants.h"
#include <hdf5.h>
#include <hdf5_hl.h>

// energy/group info
static double lnu_min, lnu_max, dlnu;
static double mev_min, mev_max;
static double lnugroup[RAD_NUM_TYPES]
                      [NU_BINS + 1]; // frequency of groups in table
static double nugroup[RAD_NUM_TYPES]
                     [NU_BINS + 1]; // frequency of groups in table
static double egroup[RAD_NUM_TYPES][NU_BINS + 1]; // energy of nugroup in MeV

// element access macros
#define ELEM3D(irho, iT, iY) (Nrho * ((iY)*NT + (iT)) + (irho))
#define ELEM4D(irho, iT, iY, it) \
  (RAD_NUM_TYPES * (Nrho * ((iY)*NT + (iT)) + (irho)) + (it))
#define ELEM5D(irho, iT, iY, it, iNu)                                   \
  ((NU_BINS + 1) *                                                      \
          (RAD_NUM_TYPES * (Nrho * ((iY)*NT + (iT)) + (irho)) + (it)) + \
      iNu)
#define ELEM5DGRID(irho, iT, iY, it, iNu) \
  (NNu * (RAD_NUM_TYPES * (Nrho * ((iY)*NT + (iT)) + (irho)) + (it)) + iNu)
#define TABLE_TOL 1e-10

// HDF table info
static int     Nrho, NT, NYe, NNu;
static int     tab_size_hdf, tab_size;
static double  tab_dlrho, tab_dlT, tab_dYe;
static double  tab_lrho_min, tab_lT_min, tab_Ye_min;
static double  tab_lrho_max, tab_lT_max, tab_Ye_max;
static double  tab_rho_min, tab_T_min;
static double  tab_rho_max, tab_T_max;
static double *tab_lrho;
static double *tab_lT;
static double *tab_Ye;
static double *tab_lNu;
static double *tab_J_hdfgrid, *tab_J;
static double *tab_K_hdfgrid, *tab_K;

// private functions
static void   interp_3d_to_1d(const double rho, const double T, double Ye,
      double *restrict tab_in, double tab_out[RAD_NUM_TYPES][NU_BINS + 1]);
static double rho2lrho(const double rho);
static double T2lT(const double temp);
static double catch_ye(const double ye);

void init_opac_emis_hdf(char *name) {
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

  // load HDF
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(name, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // size parameters
  H5LTget_attribute_int(file_id, "dimensions", "numRho", &Nrho);
  H5LTget_attribute_int(file_id, "dimensions", "numT", &NT);
  H5LTget_attribute_int(file_id, "dimensions", "numYe", &NYe);
  H5LTget_attribute_int(file_id, "dimensions", "numNu", &NNu);
  tab_size_hdf = Nrho * NT * NYe * RAD_NUM_TYPES * NNu;
  tab_size     = Nrho * NT * NYe * RAD_NUM_TYPES * (NU_BINS + 1);

  // allocate buffers
  tab_lrho      = safe_malloc(Nrho * sizeof(double));
  tab_lT        = safe_malloc(NT * sizeof(double));
  tab_Ye        = safe_malloc(NYe * sizeof(double));
  tab_lNu       = safe_malloc(NNu * sizeof(double));
  tab_J_hdfgrid = safe_malloc(tab_size_hdf * sizeof(double));
  tab_K_hdfgrid = safe_malloc(tab_size_hdf * sizeof(double));
  tab_J         = safe_malloc(tab_size * sizeof(double));
  tab_K         = safe_malloc(tab_size * sizeof(double));

  // load into buffers
  H5LTread_dataset_double(file_id, "lrho", tab_lrho);
  H5LTread_dataset_double(file_id, "lT", tab_lT);
  H5LTread_dataset_double(file_id, "Ye", tab_Ye);
  H5LTread_dataset_double(file_id, "lnu", tab_lNu);
  H5LTread_dataset_double(file_id, "emis", tab_J_hdfgrid);
  H5LTread_dataset_double(file_id, "opac", tab_K_hdfgrid);

  // done with HDF5
  H5Fclose(file_id);

  // metadata
  tab_lrho_min = tab_lrho[0];
  tab_lT_min   = tab_lT[0];
  tab_Ye_min   = tab_Ye[0];
  tab_lrho_max = tab_lrho[Nrho - 1];
  tab_lT_max   = tab_lT[NT - 1];
  tab_Ye_max   = tab_Ye[NYe - 1];
  tab_rho_min  = pow(10., tab_lrho_min);
  tab_rho_max  = pow(10., tab_lrho_max);
  tab_T_min    = pow(10., tab_lT_min);
  tab_T_max    = pow(10., tab_lT_max);
  tab_dlrho    = tab_lrho[1] - tab_lrho[0];
  tab_dlT      = tab_lT[1] - tab_lT[0];
  tab_dYe      = tab_Ye[1] - tab_Ye[0];

  // interpolate to internal NUBINS+1 grid
  for (int irho = 0; irho < Nrho; irho++) {
    for (int iT = 0; iT < NT; iT++) {
      for (int iY = 0; iY < NYe; iY++) {
        TYPELOOP {
          int elem4d = ELEM5DGRID(irho, iT, iY, itp, 0);
          for (int inu = 0; inu < NU_BINS + 1; inu++) {
            int elem5d    = ELEM5D(irho, iT, iY, itp, inu);
            tab_J[elem5d] = interp_1d(lnugroup[itp][inu], tab_lNu[0],
                tab_lNu[NNu - 1], 0, NNu, tab_lNu, &(tab_J_hdfgrid[elem4d]));
            tab_K[elem5d] = interp_1d(lnugroup[itp][inu], tab_lNu[0],
                tab_lNu[NNu - 1], 0, NNu, tab_lNu, &(tab_K_hdfgrid[elem4d]));
          }
        }
      }
    }
  }

  return;
}

void fill_opac_emis_hdf(struct of_microphysics *m) {
  const double rho = m->rho;       // should be in g/cm^3
  const double T   = (m->T) / MEV; // Should be in MeV
  const double Ye  = m->Ye;

  interp_3d_to_1d(rho, T, Ye, tab_J, m->jgrp);
  interp_3d_to_1d(rho, T, Ye, tab_K, m->alphagrp);
}

/*
 * Given microphysics, a frequency (Hz) and a neutrino species,
 * returns emissivity in units of
 * (1/4pi) ergs s^-1 cm^-3 Hz^-1 = (1/4pi) erg / cm^3
 * corresponding to j_nu = dE/(dx^3 dt dOmega dnu)
 */
double jnu_hdf(double nu, int type, const struct of_microphysics *m) {
  if (type < 0) {
    fprintf(stderr, "[jnu_hdf]: Bad photon type!\n");
    exit(1);
  }

  const double lnu = log(fabs(nu));

  double j = interp_1d(lnu, lnu_min, lnu_max, 0, NU_BINS + 1,
      &(lnugroup[type][0]), &(m->jgrp[type][0]));

  return j;
}

double Jnu_hdf(double nu, int type, const struct of_microphysics *m) {
  return 4 * M_PI * jnu_hdf(nu, type, m);
}

/* Integrates emissivity over solid angle and frequency to return
 * total emissivity per cell
 */
double int_jnudnudOmega_hdf(const struct of_microphysics *m) {
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
double alpha_nu_hdf(double nu, int type, const struct of_microphysics *m) {
  if (type < 0) {
    fprintf(stderr, "[alpha_nu_hdf]: bad photon type!\n");
    exit(1);
  }
  const double lnu = log(fabs(nu));

  double alpha = interp_1d(lnu, lnu_min, lnu_max, 0, NU_BINS + 1,
      &(lnugroup[type][0]), &(m->alphagrp[type][0]));

  return alpha;
}

static void interp_3d_to_1d(const double rho, const double T, double Ye,
    double *restrict tab_in, double tab_out[RAD_NUM_TYPES][NU_BINS + 1]) {

  // convert to log space and put things on the table
  const double lrho = rho2lrho(rho); // log (g/cm^3)
  const double lT   = T2lT(T);       // log (MeV)
  Ye                = catch_ye(Ye);  // unitless
  // Indices
  const int irho = (lrho - tab_lrho_min) / tab_dlrho;
  const int iT   = (lT - tab_lT_min) / tab_dlT;
  const int iY   = (Ye - tab_Ye_min) / tab_dYe;

  // Offsets
  const double delrho = (lrho - tab_lrho[irho]) / tab_dlrho;
  const double delT   = (lT - tab_lT[iT]) / tab_dlT;
  const double delY   = (Ye - tab_Ye[iY]) / tab_dYe;

  TYPELOOP {
    NULOOP {
      tab_out[itp][inu] =
          ((1 - delrho) * (1 - delT) * (1 - delY) *
                  tab_in[ELEM5D(irho, iT, iY, itp, inu)] +
              delrho * (1 - delT) * (1 - delY) *
                  tab_in[ELEM5D(irho + 1, iT, iY, itp, inu)] +
              (1 - delrho) * delT * (1 - delY) *
                  tab_in[ELEM5D(irho, iT + 1, iY, itp, inu)] +
              delrho * delT * (1 - delY) *
                  tab_in[ELEM5D(irho + 1, iT + 1, iY, itp, inu)] +
              (1 - delrho) * (1 - delT) * delY *
                  tab_in[ELEM5D(irho, iT, iY + 1, itp, inu)] +
              delrho * (1 - delT) * delY *
                  tab_in[ELEM5D(irho + 1, iT, iY + 1, itp, inu)] +
              (1 - delrho) * delT * delY *
                  tab_in[ELEM5D(irho, iT + 1, iY + 1, itp, inu)] +
              delrho * delT * delY *
                  tab_in[ELEM5D(irho + 1, iT + 1, iY + 1, itp, inu)]);
      // printf("\t%d, %d: %e\n",itp,inu,tab_out[itp][inu]); // DEBUG
    }
  }
}

static double rho2lrho(const double rho) {
  // dont' fall off the table
  if (rho < tab_rho_min) {
    return tab_lrho_min;
  } else if (rho >= tab_rho_max) {
    return tab_lrho_max - TABLE_TOL / 100.;
  } else {
    return log10(rho);
  }
}

static double catch_ye(const double ye) {
  if (ye < tab_Ye_min) {
    return tab_Ye_min;
  } else if (ye >= tab_Ye_max) {
    return tab_Ye_max - TABLE_TOL / 100.;
  } else {
    return ye;
  }
}

static double T2lT(const double temp) {
  if (temp < tab_T_min) {
    return tab_lT_min;
  } else if (temp >= tab_T_max) {
    return tab_lT_max - TABLE_TOL / 100.;
  } else {
    return log10(temp);
  }
}

#undef TABLE_TOL
#undef ELEM3D
#undef ELEM4D
#undef ELEM5D
#undef ELEM5DGRID

#endif // RADIATION == RADTYPE_NEUTRINOS && HDF_OPACITIES
