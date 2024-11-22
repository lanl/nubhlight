/******************************************************************************
 *                                                                            *
 * OSCILLATIONS.C                                                             *
 *                                                                            *
 * Routines for neutrino oscillations                                         *
 *                                                                            *
 ******************************************************************************/
#include "decs.h"

#if RADIATION
#if LOCAL_ANGULAR_DISTRIBUTIONS

void accumulate_local_angles() {
  static const int LOCAL_ANGLES_SIZE = 2 * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 *
                                       LOCAL_ANGLES_NMU * RAD_NUM_TYPES;
  memset(local_angles, 0, LOCAL_ANGLES_SIZE * sizeof(double));

#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        int ix1, ix2, icosth1, icosth2;
        get_local_angle_bins(ph, &ix1, &ix2, &icosth1, &icosth2);
#pragma omp atomic
        local_angles[0][ix1][ix2][ph->type][icosth1] += ph->w;
#pragma omp atomic
        // local_angles global
        local_angles[1][ix1][ix2][ph->type][icosth2] += ph->w;
      }
      ph = ph->next;
    }
  } // omp parallel

  mpi_dbl_allreduce_array((double *)local_angles, LOCAL_ANGLES_SIZE);

  // Gnu, local_moments are global
#if RAD_NUM_TYPES >= 4
  compute_local_gnu(local_angles, Gnu);
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
void compute_local_gnu(grid_local_angles_type f, grid_Gnu_type gnu) {
#pragma omp parallel for collapse(4)
  for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
    for (int i = 0; i < LOCAL_ANGLES_NX1; ++i) {
      for (int j = 0; j < LOCAL_ANGLES_NX2; ++j) {
        for (int imu = 0; imu < LOCAL_ANGLES_NMU; ++imu) {
          // TODO(JMM): Generalize this for six species?
          double ELN =
              (f[b][i][j][NU_ELECTRON][imu] - f[b][i][j][ANTINU_ELECTRON][imu]);
          double XLN =
              (f[b][i][j][NU_HEAVY][imu] - f[b][i][j][ANTINU_HEAVY][imu]);
          gnu[b][i][j][imu] = ELN - XLN;
        }
      }
    }
  }
}

// JMM: We can also compute, e.g., the average bin momentum if we need
// to, e.g., compute higher moment integrands
void compute_local_moments(grid_Gnu_type gnu, grid_local_moment_type moments) {
  // We are reducing over mu, but if we just parallel loop over b,i,j,
  // there is no danger of index collisions.
  for (int imu = 0; imu < LOCAL_ANGLES_NMU; ++imu) {
#pragma omp parallel for collapse (3)
    for (int b = 0; b < LOCAL_NUM_BASES; ++b) {
      for (int i = 0; i < LOCAL_ANGLES_NX1; ++i) {
        for (int j = 0; j < LOCAL_ANGLES_NX2; ++j) {
          if (gnu[b][i][j][imu] < 0) {
            // TODO(JMM): Pretty sure this atomic isn't needed
#pragma omp atomic
            moments[b][MOMENTS_A][i][j] += gnu[b][i][j][imu];
          } else if (gnu[b][i][j][imu] > 0) {
            // TODO(JMM): Pretty sure this atomic isn't needed
#pragma omp atomic
            moments[b][MOMENTS_B][i][j] += gnu[b][i][j][imu];
          } // else nada. We don't care about == 0.
        }
      }
    }
  }
}
#endif // RAD_NUM_TYPES >= 4

#endif // LOCAL_ANGULAR_DISTRIBUTIONS
#endif // RADIATION
