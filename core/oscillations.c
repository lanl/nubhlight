/******************************************************************************
 *                                                                            *
 * OSCILLATIONS.C                                                             *
 *                                                                            *
 * Routines for neutrino oscillations                                         *
 *                                                                            *
 ******************************************************************************/
#include "decs.h"

#if RADIATION == RADTYPE_NEUTRINOS
#if LOCAL_ANGULAR_DISTRIBUTIONS

double get_dt_ffi() {
  timer_start(TIMER_OSCILLATIONS);
  set_Rmunu(); // So we have Nsph and nph
  double nph_max = 0;
#pragma omp parallel for reduction(max : nph_max) collapse(3)
  ZLOOP {
    nph_max = MY_MAX(nph_max, nph[i][j][k]); // 1/cm^3
  }
  // seconds
  double dt_osc = 1. / (NUFERM * nph_max + SMALL);
  dt_osc /= T_unit; // code units
  dt_osc = mpi_min(dt_osc);
  timer_stop(TIMER_OSCILLATIONS);
  return dt_osc;
}

void accumulate_local_angles() {
  timer_start(TIMER_OSCILLATIONS);
  static const int LOCAL_ANGLES_SIZE = LOCAL_NUM_BASES * LOCAL_ANGLES_NX1 *
                                       LOCAL_ANGLES_NX2 * LOCAL_ANGLES_NMU *
                                       RAD_NUM_TYPES;
  static const int LOCAL_STDDEV_SIZE =
      LOCAL_NUM_BASES * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 * LOCAL_ANGLES_NMU;

  memset(local_angles, 0, LOCAL_ANGLES_SIZE * sizeof(double));
  memset(local_Ns, 0, LOCAL_STDDEV_SIZE * sizeof(double));
  memset(local_wsqr, 0, LOCAL_STDDEV_SIZE * sizeof(double));

#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        int ix1, ix2, icosth[LOCAL_NUM_BASES];
        get_local_angle_bins(ph, &ix1, &ix2, &icosth[0], &icosth[1]);
        for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
#pragma omp atomic
          local_angles[b][ix1][ix2][ph->type][icosth[b]] += ph->w;
#pragma omp atomic
          local_Ns[b][ix1][ix2][icosth[b]] += 1.;
#pragma omp atomic
          local_wsqr[b][ix1][ix2][icosth[b]] += (ph->w) * (ph->w);
        }
      }
      ph = ph->next;
    }
  } // omp parallel

  mpi_dbl_allreduce_array((double *)local_angles, LOCAL_ANGLES_SIZE);
  mpi_dbl_allreduce_array((double *)local_Ns, LOCAL_STDDEV_SIZE);
  mpi_dbl_allreduce_array((double *)local_wsqr, LOCAL_STDDEV_SIZE);

  // Gnu, local_moments are global
#if RAD_NUM_TYPES >= 4
  compute_local_gnu(local_angles, local_Ns, local_wsqr, Gnu);
  compute_local_moments(Gnu, local_moments);
#endif // RAD_NUM_TYPES >= 4
  timer_stop(TIMER_OSCILLATIONS);
}

void get_local_angle_bins(
    struct of_photon *ph, int *pi, int *pj, int *pmu1, int *pmu2) {
  double X[NDIM];
  double Kcov[NDIM];
  double Kcon[NDIM];
  int    k;
  get_X_K_interp(ph, t, P, X, Kcov, Kcon);
  Xtoijk(X, pi, pj, &k);

  // JMM: If we use more complicated bases this is more complicated
  // change this for other basis vectors
  double X1norm      = sqrt(fabs(ggeom[*pi][*pj][CENT].gcov[1][1]));
  double X2norm      = sqrt(fabs(ggeom[*pi][*pj][CENT].gcov[2][2]));
  double X1vec[NDIM] = {0, 1. / (X1norm + SMALL), 0, 0};
  double X2vec[NDIM] = {0, 0, 1. / (X2norm + SMALL), 0};

  double costh1 = 0;
  double costh2 = 0;
  SDLOOP {
    costh1 += X1vec[mu] * Kcov[mu]; // X1^a K_a
    costh2 += X2vec[mu] * Kcov[mu]; // X2^a K_a
  }
  // cos(th) = e^a K_a / |e| |K|
  // |e| = 1 by construction, but K must be normalized
  // note we want to normalize the SPATIAL part of K
  double knorm = 0;
  for (int mu = 1; mu < NDIM; ++mu) {
    for (int nu = 1; nu < NDIM; ++nu) {
      knorm += fabs((ggeom[*pi][*pj][CENT].gcov[mu][nu]) * Kcon[mu] * Kcon[mu]);
    }
  }
  // sqrt the inner product and lets go
  knorm = sqrt(knorm);
  knorm = 1. / (fabs(knorm) + SMALL);
  costh1 *= knorm;
  costh2 *= knorm;

  *pi = MY_MAX(
      0, MY_MIN(LOCAL_ANGLES_NX1 - 1, (X[1] - startx_rad[1]) / local_dx1_rad));
  *pj = MY_MAX(
      0, MY_MIN(LOCAL_ANGLES_NX2 - 1, (X[2] - startx_rad[2]) / local_dx2_rad));
  *pmu1 =
      MY_MAX(0, MY_MIN(LOCAL_ANGLES_NMU - 1, (costh1 + 1) / local_dx_costh));
  *pmu2 =
      MY_MAX(0, MY_MIN(LOCAL_ANGLES_NMU - 1, (costh2 + 1) / local_dx_costh));
}

#if RAD_NUM_TYPES >= 4
void compute_local_gnu(grid_local_angles_type f, grid_Gnu_type local_Ns,
                       grid_Gnu_type local_wsqr, grid_Gnu_type gnu) {
#pragma omp parallel for collapse(4)
  for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
    LOCALXMULOOP {
      const double Ns = local_Ns[b][i][j][imu];
      const double w2 = local_wsqr[b][i][j][imu];

      double   wmean = 0;
      TYPELOOP wmean += fabs(f[b][i][j][itp][imu]);
      wmean /= (Ns + SMALL);

      const double wb2N = Ns * wmean * wmean;
      // should have units of sum(w)/sqrt(N) ~ sqrt(N) wmean
      // middle term scales b/c
      // sqrt(sum w^2) ~ sqrt(N wmean^2) ~ sqrt(N) wmean
      const double stddev =
          sqrt(wb2N + Ns * (w2 - wb2N) / (fabs(Ns - 1) + SMALL));

      // TODO(JMM): Generalize this for six species?
      const double ELN =
          (f[b][i][j][NU_ELECTRON][imu] - f[b][i][j][ANTINU_ELECTRON][imu]);
      const double XLN =
          (f[b][i][j][NU_HEAVY][imu] - f[b][i][j][ANTINU_HEAVY][imu]);

      double g_temp     = ELN - 0.5 * XLN;
      gnu[b][i][j][imu] = (fabs(g_temp) > stddev) * g_temp;
    }
  }
}

// JMM: We can also compute, e.g., the average bin momentum if we need
// to, e.g., compute higher moment integrands
void compute_local_moments(grid_Gnu_type gnu,
                           grid_local_moment_type moments) {
  // We are reducing over mu, but if we just parallel loop over b,i,j,
  // there is no danger of index collisions.
  LOCALMULOOP {
#pragma omp parallel for collapse(3)
    for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
      LOCALXLOOP {
        if (gnu[b][i][j][imu] < 0) {
          // TODO(JMM): Pretty sure this atomic isn't needed
#pragma omp atomic
          moments[b][MOMENTS_A][i][j] += fabs(gnu[b][i][j][imu]);
        } else if (gnu[b][i][j][imu] > 0) {
          // TODO(JMM): Pretty sure this atomic isn't needed
#pragma omp atomic
          moments[b][MOMENTS_B][i][j] += fabs(gnu[b][i][j][imu]);
        } // else nada. We don't care about == 0.
      }
    }
  } // MU LOOP

#pragma omp parallel for collapse(2)
  LOCALXLOOP {
    for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
      moments[b][MOMENTS_DIFF][i][j] =
          fabs(moments[b][MOMENTS_B][i][j] - moments[b][MOMENTS_A][i][j]);
    }
    local_b_osc[i][j] = (local_moments[1][MOMENTS_DIFF][i][j] >
                         local_moments[0][MOMENTS_DIFF][i][j]);
  }
}

void oscillate_ffi(grid_local_angles_type f,
                   grid_local_moment_type local_moments, grid_Gnu_type gnu) {
  timer_start(TIMER_OSCILLATIONS);
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        int ix1, ix2, icosth[LOCAL_NUM_BASES];
        double X[NDIM], Kcov[NDIM], Kcon[NDIM];
        get_local_angle_bins(ph, &ix1, &ix2, &icosth[0], &icosth[1]);

        int    b_osc = local_b_osc[ix1][ix2];
        int    imu   = icosth[b_osc];
        double A     = local_moments[b_osc][MOMENTS_A][ix1][ix2];
        double B     = local_moments[b_osc][MOMENTS_B][ix1][ix2];
        double G     = gnu[b_osc][ix1][ix2][imu];

        // gnu == 0 when we activated stddev trigger. Don't oscillate.
        if (((G != 0) || FORCE_EQUIPARTITION) && (A != 0) && (B != 0)) {
          // finite volumes cell
          int i, j, k;
          get_X_K_interp(ph, t, P, X, Kcov, Kcon);
          Xtoijk(X, &i, &j, &k);

          // compute crossing depth
          double n_tot_supercell = 0;
          TYPELOOP {
            LOCALMULOOP {
              // total neutrinos of all flavors
              n_tot_supercell += fabs(f[0][ix1][ix2][itp][imu]);
            }
          }
          // crossing depth in units of number/cm^3
          // sqrt(AB)/n_tot_supercell is normalized unitless crossing
          // depth but units should be in number/cm^3.
          double crossing_depth = sqrt(fabs(A * B))*nph[i][j][k] / (n_tot_supercell + SMALL);

          // find local FFI time scale nph already computed because of
          // Rmunu This crossing depth formula is slightly more
          // sophisticated than just 1 / G_F n
          double tau = 1. / (T_unit * NUFERM * crossing_depth  + SMALL);

          // If A == B then which region we treat as shallow is
          // unimportant. Psurvive = 1/3 for both regions.
          int    A_is_shallow = A < B;
          int    B_is_shallow = !(A_is_shallow);
          double shallow      = A_is_shallow ? A : B;
          double deep         = A_is_shallow ? B : A;

          int g_in_A = G < 0;
          int g_in_B = G > 0;

          int in_shallow = (A_is_shallow && g_in_A) || (B_is_shallow && g_in_B);

          double peq = nu_is_heavy(ph->type) ? (2. / 3.) : (1. / 3.);
#if FORCE_EQUIPARTITION
          double p_survival = peq;
#else
          double p_survival =
              in_shallow ? peq : (1 - (1 - peq) * shallow / (deep + SMALL));
#endif // FORCE_EQUIPARTITION
          double p_osc = MY_MIN(1.0, dt / (tau + SMALL)) * (1. - p_survival);
          if (get_rand() < p_osc) {
            // JMM:
            // Type order is NUE, NUEBAR, NUX, NUXBAR
            // adding 2 on ring 0, 1, 2, 3
            // moves through without changing to antiparticle.
            ph->type = (ph->type + (RAD_NUM_TYPES / 2)) % RAD_NUM_TYPES;
            ph->osc_count += 1;

            #pragma omp atomic
            local_osc_count[ix1][ix2] += ph->w;
          }
        }
      }
      ph = ph->next;
    }
  } // omp parallel
  timer_stop(TIMER_OSCILLATIONS);
}

#endif // RAD_NUM_TYPES >= 4
#endif // LOCAL_ANGULAR_DISTRIBUTIONS

#if NEUTRINO_OSCILLATIONS_CFI
void compute_cfi_symmetrized_avgs(grid_symm_radtype_type nph_flavor,
                                  grid_symm_radtype_type kappa_avg,
                                  grid_CFI_Gamma_type cfi_gamma) {
  timer_start(TIMER_OSCILLATIONS);

  static const int SYMM_SIZE = (N1 + 2*NG)*(N2 + 2*NG);

  memset(nph_flavor, 0, SYMM_SIZE * RAD_NUM_TYPES * sizeof(double));
  memset(kappa_avg, 0, SYMM_SIZE * RAD_NUM_TYPES * sizeof(double));
  memset(cfi_Gamma, 0, 2 * SYMM_SIZE * sizeof(double));

#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        double X[NDIM], Kcov[NDIM], Kcon[NDIM];
        double kappa[RAD_NUM_TYPES];
        int i, j, k;
        get_X_K_interp(ph, t, P, X, Kcov, Kcon);
        Xtoijk(X, &i, &j, &k);
        
        #pragma omp atomic
        nph_flavor[i][j][ph->type] += ph->w;

        double theta = get_bk_angle(
            X, Kcon, Ucov_grd[i][j][k],
            Bcov_grd[i][j][k], m_grd[i][j][k].B);
        double nu = get_fluid_nu(X, Kcov, Ucon_grd[i][j][k]);
        TYPELOOP {
          kappa[itp] = alpha_inv_abs(nu, itp,
                                     &(m_grd[i][j][k]), theta);
          kappa[itp] /= (nu + SMALL);
        }

        #pragma omp atomic
        kappa_avg[i][j][ph->type] += kappa[ph->type] * ph->w;

        // electrons and antis positive, heavies and antis negative
        int gamma_sign = (ph->type < NU_HEAVY) ? 1 : -1;

        // Check for the need to divide neutrino count for heavies by 2
        double heavy_dupfac = ((RAD_NUM_TYPES < 4) && (ph->type >= NU_HEAVY))? 0.5 : 1.0;

        // check for antiparticles
        int isanti = is_antiparticle(ph);

        int kappa_e_idx = NU_ELECTRON + isanti; // kappa_e for es, kappa_ebar for ebars
        // kappa_x for xs, kappa_xbar for xbars. But account for x = xbar
        int kappa_x_idx = NU_HEAVY + (RAD_NUM_TYPES == 4)*isanti; 

        // es contribute + to this sum, xs minus. Need to divide by 2
        // if heavies and their antis are bundled together.
        double weight = gamma_sign * heavy_dupfac * ph->w;

        // gamma_0 = (1/2) (kappa_e + kappa_x)
        // gamma_1 = (1/2) (kappa_ebar + kappa_xbar)
        //
        // Method A
        // avg gamma_0 = \int E^2 dE (kappa_e + kappa_x) (f_e - f_x) / \int E^2 / \int E^2 dE (f_e - f_x)
        // avg gamma_1 = \int E^2 dE (kappa_ebar + kappa_xbar) (f_ebar - f_xbar) / \int E^2 / \int E^2 dE (f_ebar - f_xbar)
        #pragma omp atomic
        cfi_Gamma[i][j][isanti] += 0.5*(kappa[kappa_e_idx] + kappa[kappa_x_idx])*weight;
      }
      ph = ph->next;
    }
  }

  mpi_dbl_allreduce_array((double *)nph_flavor, SYMM_SIZE * RAD_NUM_TYPES);
  mpi_dbl_allreduce_array((double *)kappa_avg, SYMM_SIZE * RAD_NUM_TYPES);
  mpi_dbl_allreduce_array((double *)cfi_Gamma, 2 * SYMM_SIZE);

  #pragma omp parallel for collapse(2)
  ILOOP {
    JLOOP {
      cfi_Gamma[i][j][0] /= (nph_flavor[i][j][NU_ELECTRON] - nph_flavor[i][j][NU_HEAVY]);
      cfi_Gamma[i][j][1] /= (nph_flavor[i][j][ANTINU_ELECTRON] - nph_flavor[i][j][ANTINU_HEAVY]);
    }
  }

  timer_stop(TIMER_OSCILLATIONS);
}
#endif // NEUTRINO_OSCILLATIONS_CFI

#endif // RADIATION
