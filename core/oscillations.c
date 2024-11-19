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

void local_accum_superph(double X[NDIM],
                         double Kcov[NDIM], double Kcon[NDIM],
                         double w, int type,
                         struct of_geom *pgeom,
                         grid_local_angles_type local_angles);

void accumulate_local_angles() {
  static const int LOCAL_ANGLES_SIZE = 2 * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 *
                                       LOCAL_ANGLES_NMU * RAD_NUM_TYPES;
  memset(local_angles, 0, LOCAL_ANGLES_SIZE * sizeof(double));

#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type != TYPE_TRACER) {
        double X[NDIM];
        double Kcov[NDIM];
        double Kcon[NDIM];
        int i, j, k;
        get_X_K_interp(ph, t, P, X, Kcov, Kcon);
        Xtoijk(X, &i, &j, &k);
        local_accum_superph(X, Kcov, Kcon, ph->w, ph->type,
                            &ggeom[i][j][CENT], local_angles);
      }
      ph = ph->next;
    }
  }

  mpi_dbl_allreduce_array((double*)local_angles, LOCAL_ANGLES_SIZE);
}

void local_accum_superph(double X[NDIM],
                         double Kcov[NDIM], double Kcon[NDIM],
                         double w, int type,
                         struct of_geom *pgeom,
                         grid_local_angles_type local_angles) {
  if (type == TYPE_TRACER) return;
  // JMM: If we use more complicated bases this is more complicated
  // change this for other basis vectors
  double X1norm = sqrt(fabs(pgeom->gcov[1][1]));
  double X2norm = sqrt(fabs(pgeom->gcov[2][2]));
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
      knorm += fabs((pgeom->gcov[mu][nu])*Kcon[mu]*Kcon[mu]);
    }
  }
  // sqrt the inner product and lets go
  knorm = sqrt(knorm);
  knorm = 1./(fabs(knorm) + SMALL);
  costh1 *= knorm;
  costh2 *= knorm;

  int ix1     = MY_MAX(0, MY_MIN(LOCAL_ANGLES_NX1 - 1, (X[1] - startx_rad[1]) / local_dx1_rad));
  int ix2     = MY_MAX(0, MY_MIN(LOCAL_ANGLES_NX2 - 1, (X[2] - startx_rad[2]) / local_dx2_rad));
  int icosth1 = MY_MAX(0, MY_MIN(LOCAL_ANGLES_NMU - 1, (costh1 + 1) / local_dx_costh));
  int icosth2 = MY_MAX(0, MY_MIN(LOCAL_ANGLES_NMU - 1, (costh2 + 1) / local_dx_costh));

  #pragma omp atomic
  local_angles[0][ix1][ix2][type][icosth1] += w;
  #pragma omp atomic
  local_angles[1][ix1][ix2][type][icosth2] += w;
}

#endif // LOCAL_ANGULAR_DISTRIBUTIONS
#endif // RADIATION
