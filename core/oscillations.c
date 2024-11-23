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

double get_dt_oscillations() {
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
  return dt_osc;
}

void accumulate_local_angles() {
  static const int LOCAL_ANGLES_SIZE = LOCAL_NUM_BASES * LOCAL_ANGLES_NX1 *
                                       LOCAL_ANGLES_NX2 * LOCAL_ANGLES_NMU *
                                       RAD_NUM_TYPES;
  static const int LOCAL_STDDEV_SIZE =
      LOCAL_NUM_BASES * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 * LOCAL_ANGLES_NMU;

  memset(local_angles, 0, LOCAL_ANGLES_SIZE * sizeof(double));
  memset(local_stddev, 0, LOCAL_STDDEV_SIZE * sizeof(double));

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
          local_stddev[b][ix1][ix2][icosth[b]] += 1.;
        }
      }
      ph = ph->next;
    }
  } // omp parallel

  mpi_dbl_allreduce_array((double *)local_angles, LOCAL_ANGLES_SIZE);
  mpi_dbl_allreduce_array((double *)local_stddev, LOCAL_STDDEV_SIZE);

#pragma omp parallel for collapse(4)
  for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
    LOCALXMULOOP {
      local_stddev[b][i][j][imu] = sqrt(fabs(local_stddev[b][i][j][imu]));
    }
  }

  // Gnu, local_moments are global
#if RAD_NUM_TYPES >= 4
  compute_local_gnu(local_angles, local_stddev, Gnu);
  compute_local_moments(Gnu, local_moments);
#endif // RAD_NUM_TYPES >= 4
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
void compute_local_gnu(
    grid_local_angles_type f, grid_Gnu_type stddev, grid_Gnu_type gnu) {
#pragma omp parallel for collapse(4)
  for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
    LOCALXMULOOP {
      // TODO(JMM): Generalize this for six species?
      double ELN =
          (f[b][i][j][NU_ELECTRON][imu] - f[b][i][j][ANTINU_ELECTRON][imu]);
      double XLN = (f[b][i][j][NU_HEAVY][imu] - f[b][i][j][ANTINU_HEAVY][imu]);

      double   tot = 0;
      TYPELOOP tot += f[b][i][j][itp][imu];
      double   ebar = tot / (stddev[b][i][j][imu] + SMALL);

      double g_temp     = ELN - XLN;
      gnu[b][i][j][imu] = (fabs(g_temp) > ebar) * g_temp;
    }
  }
}

// JMM: We can also compute, e.g., the average bin momentum if we need
// to, e.g., compute higher moment integrands
void compute_local_moments(grid_Gnu_type gnu, grid_local_moment_type moments) {
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

void        oscillate(grid_local_moment_type local_moments, grid_Gnu_type gnu) {
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        int ix1, ix2, icosth1, icosth2;
        get_local_angle_bins(ph, &ix1, &ix2, &icosth1, &icosth2);

        int    b_osc = local_b_osc[ix1][ix2];
        double A     = local_moments[b_osc][MOMENTS_A][ix1][ix2];
        double B     = local_moments[b_osc][MOMENTS_B][ix1][ix2];
        int    imu   = b_osc ? icosth2 : icosth1;

        // gnu == 0 when we activated stddev trigger. Don't oscillate.
        if (!((gnu[b_osc][ix1][ix2][imu] == 0) || ((A == 0) && (B == 0)))) {
          // If A == B then which region we treat as shallow is
          // unimportant. Psurviv = 1/3 for both regions.
          int    A_is_shallow = A < B;
          int    B_is_shallow = !(A_is_shallow);
          double shallow      = A_is_shallow ? A : B;
          double deep         = A_is_shallow ? B : A;

          int g_in_A = gnu[b_osc][ix1][ix2][imu] < 0;
          int g_in_B = gnu[b_osc][ix1][ix2][imu] > 0;

          int in_shallow = (A_is_shallow && g_in_A) || (B_is_shallow && g_in_B);
          int in_deep    = !(in_shallow);

          double p_survival =
              in_shallow ? PEQ : (1 - (1 - PEQ) * B / (A + SMALL));
          double p_osc = 1. - p_survival;
          if (get_rand() < p_osc) {
            // JMM:
            // Type order is NUE, NUEBAR, NUX, NUXBAR
            // adding 2 on ring 0, 1, 2, 3
            // moves through without changing to antiparticle.
            ph->type = (ph->type + (RAD_NUM_TYPES / 2)) % RAD_NUM_TYPES;
          }
        }
      }
      ph = ph->next;
    }
  } // omp parallel
}

#endif // RAD_NUM_TYPES >= 4

#endif // LOCAL_ANGULAR_DISTRIBUTIONS
#endif // RADIATION
