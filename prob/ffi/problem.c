/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR UNIT TEST FOR MULTISCATT                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double Nsph_tot;
static double Nph_tot;

static double fracs[RAD_NUM_TYPES];
static double tot_frac;

static double E_MEV;
static double rho0;
static double uu0;
static int    direction;

// Problem-specific variables to set at runtime
void set_problem_params() {
  set_param("Nsph_tot", &Nsph_tot);
  set_param("Nph_tot", &Nph_tot);

  // fraction of each species
  set_param("frac_e", &fracs[NU_ELECTRON]);
  set_param("frac_ebar", &fracs[ANTINU_ELECTRON]);
  set_param("frac_x", &fracs[NU_HEAVY]);
  set_param("frac_xbar", &fracs[ANTINU_HEAVY]);


  // neutrino energy
  set_param("E_MEV", &E_MEV);

  // background gas
  set_param("rho0", &rho0);
  set_param("uu0", &uu0);

  set_param("direction", &direction);
}

double distmu(int type, double x) {
  if (type == NU_ELECTRON) {
    return 1.0 - (3. / 4.) * x * x;
  } else if (type == ANTINU_ELECTRON) {
    return (1. / 4.) + (3. / 4.) * x * x;
  } else {
    return 1.0;
  }
}

// Initialize dynamical variables
void init_prob() {
  if ((RADIATION != RADTYPE_NEUTRINOS) || !LOCAL_ANGULAR_DISTRIBUTIONS ||
      !NEUTRINO_OSCILLATIONS || (RAD_NUM_TYPES != 4)) {
    if (mpi_io_proc()) {
      printf("Test needs the following physics active!\n"
             "\tRADIATION == RADTYPE_NEUTRINOS: %d\n"
             "\tLOCAL_ANGULAR_DISTRIBUTIONS: %d\n"
             "\tNEUTRINO_OSCILLATIONS: %d\n"
             "\tRAD_NUM_TYPES == 4: %d\n",
          RADIATION, LOCAL_ANGULAR_DISTRIBUTIONS, NEUTRINO_OSCILLATIONS,
          RAD_NUM_TYPES);
    }
    exit(1);
  }

  tot_frac = 0;
  TYPELOOP tot_frac += fracs[itp];
  if (fabs(tot_frac - 1) > 1e-5) {
    if (mpi_io_proc()) {
      printf("Fractions don't add up to 1!\n"
             "%.14e %.14e %.14e %.14e\n",
             fracs[0], fracs[1], fracs[2], fracs[3]);
    }
    exit(1);
  }

  if ((direction < 1) || (direction > 2)) {
    if (mpi_io_proc()) {
      printf("direction may only be 1 or 2!\n");
    }
    exit(1);
  }

  double Nsph_per_cell = Nsph_tot / (N1TOT * N2TOT * N3TOT);
  int Nsph_per_thread = (int)ceil(Nsph_per_cell/nthreads);
  double nu_cgs        = E_MEV * MEV / HPL;
  double lnu           = log10(nu_cgs);
  double nu            = pow(10., lnu);

  if (mpi_io_proc()) {
    printf("E_MEV                    = %g\n"
           "lnu                      = %g\n"
           "nu                       = %g\n"
           "rho0                     = %g\n"
           "uu0                      = %g\n"
           "Nphys                    = %g\n"
           "Nnum                     = %g\n"
           "Nsph_per_cell            = %g\n"
           "Nsph_per_cell_per_thread = %d\n"
           "Fracs                    = [%g, %g, %g, %g]\n"
           "direction                = %d\n",
           E_MEV, lnu, nu_cgs, rho0, uu0,
           Nph_tot, Nsph_tot, Nsph_per_cell, Nsph_per_thread,
           fracs[0], fracs[1], fracs[2], fracs[3],
           direction);
  }

#pragma omp parallel for collapse(3)
  ZLOOP {
    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU]  = uu0;
    P[i][j][k][U1]  = 0.0;
    P[i][j][k][U2]  = 0.0;
    P[i][j][k][U3]  = 0.0;
    P[i][j][k][B1]  = 0.0;
    P[i][j][k][B2]  = 0.0;
    P[i][j][k][B3]  = 0.0;
#if EOS == EOS_TYPE_TABLE
    P[i][j][k][YE]    = 0.5;
    P[i][j][k][YE_EM] = 0.5;
#endif
  } // ZLOOP

#pragma omp parallel
  {
    double            X[NDIM], K_tetrad[NDIM];
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    ZLOOP {
      coord(i, j, k, CENT, X);

      for (int n = 0; n < Nsph_per_thread; n++) {
        struct of_photon *phadd = safe_malloc(sizeof(struct of_photon));

        // photon position centered on cell
        phadd->X[2][0] = 0.;
        for (int mu = 1; mu < NDIM; mu++) {
          phadd->X[2][mu] = X[mu];
        }

        // sample particle type
        double x       = get_rand();
        double compare = 0;
        TYPELOOP {
          compare += fracs[itp];
          if (x <= compare) {
            phadd->type = itp;
            break;
          }
        }

        // compute mu, phi
        double mu, f, fmax = 1.;
        do { // rejection sample for mu
          mu = 2 * get_rand() - 1;
          x  = get_rand() * fmax;
          f  = distmu(phadd->type, mu);
        } while (x >= f);
        // phi
        double phi = 2 * M_PI * get_rand();
        // angles
        double costh = mu;
        double sinth = sqrt(1 - mu*mu);

        double E    = HPL * nu / (ME * CL * CL);
        K_tetrad[0] = -E;
        if (direction == 1) {
          K_tetrad[1] = E * costh;
          K_tetrad[2] = E * sinth;
        } else { // direction == 2
          K_tetrad[1] = E * sinth;
          K_tetrad[2] = E * costh;
        }
        K_tetrad[3] = 0;

        // Can just use non-moving fluid frame
        DLOOP1 {
          phadd->Kcov[2][mu] = K_tetrad[mu];
          phadd->Kcon[2][mu] = K_tetrad[mu];
        }
        // Minkowski background so Kcov/Kcon only differ by this
        phadd->Kcon[2][0] *= -1.;

        phadd->t0        = 0.;
        phadd->origin[0] = nstep;
        phadd->origin[1] = i;
        phadd->origin[2] = j;
        phadd->origin[3] = k;
        phadd->w         = Nph_tot / Nsph_tot;
        phadd->nscatt    = 0;

        phadd->next = ph;
        ph          = phadd;
      }
    }
    photon_lists[omp_get_thread_num()] = ph;
  } // omp parallel
}
