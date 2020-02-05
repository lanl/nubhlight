/******************************************************************************
 *                                                                            *
 * TETRADS.C                                                                  *
 *                                                                            *
 * CONSTRUCT AND APPLY TRANSFORMATION MATRICES BETWEEN FLUID AND COORDINATE   *
 * FRAMES                                                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void normalize(double *Vcon, double Gcov[4][4]);
void project_out(double *Vcona, double *Vconb, double Gcov[4][4]);

void coord_to_tetrad(
    double Ecov[NDIM][NDIM], double Kcoord[NDIM], double Ktetrad[NDIM]) {
  for (int mu = 0; mu < NDIM; mu++) {
    Ktetrad[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      Ktetrad[mu] += Ecov[mu][nu] * Kcoord[nu];
    }
  }
}

void tetrad_to_coord(
    double Econ[NDIM][NDIM], double Ktetrad[NDIM], double Kcoord[NDIM]) {
  for (int mu = 0; mu < NDIM; mu++) {
    Kcoord[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      Kcoord[mu] += Econ[nu][mu] * Ktetrad[nu];
    }
  }
}

// Remove this i, j, k!
#define SMALL_VECTOR (1.e-30) // TODO: change this size?
void make_tetrad(int i, int j, int k, double Ucon[NDIM], double trial[NDIM],
    double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]) {
  // closeness of trial to other three vectors
  double X1ness = 0., X2ness = 0., X3ness = 0.;
  DLOOP1 X1ness += Gcov[1][mu] * trial[mu] / sqrt(fabs(Gcov[1][1]));
  DLOOP1 X2ness += Gcov[2][mu] * trial[mu] / sqrt(fabs(Gcov[2][2]));
  DLOOP1 X3ness += Gcov[3][mu] * trial[mu] / sqrt(fabs(Gcov[3][3]));
  X1ness = fabs(X1ness);
  X2ness = fabs(X2ness);
  X3ness = fabs(X3ness);

  // Normalize trial vector
  double norm = 0.;
  DLOOP2 { norm += trial[mu] * trial[nu] * Gcov[mu][nu]; }

  // Time component parallel to u^{\mu}
  DLOOP1 Econ[0][mu] = Ucon[mu];

  // Use X1, X2, X3. Then, whichever is closest to trial, overwrite.
  DLOOP1 Econ[1][mu] = delta(mu, 1);
  DLOOP1 Econ[2][mu] = delta(mu, 2);
  DLOOP1 Econ[3][mu] = delta(mu, 3);

  if (norm > SMALL_VECTOR) {
    // We can use the trial vector
    if (X1ness > X2ness && X1ness > X3ness) {
      // Trial vector is closest to X1. Overwrite
      DLOOP1 Econ[1][mu] = trial[mu];
    } else if (X2ness >= X1ness && X2ness > X3ness) {
      // Trial vector is closest to X2. Overwrite
      DLOOP1 Econ[2][mu] = trial[mu];
    } else { // Trial vector is closest X3. Overwrite
      DLOOP1 Econ[3][mu] = trial[mu];
    }
  }

  // Gram-schmidt and normalization
  normalize(Econ[0], Gcov);
  project_out(Econ[1], Econ[0], Gcov);
  normalize(Econ[1], Gcov);
  project_out(Econ[2], Econ[0], Gcov);
  project_out(Econ[2], Econ[1], Gcov);
  normalize(Econ[2], Gcov);
  project_out(Econ[3], Econ[0], Gcov);
  project_out(Econ[3], Econ[1], Gcov);
  project_out(Econ[3], Econ[2], Gcov);
  normalize(Econ[3], Gcov);

  // Make covariant version
  for (int mu = 0; mu < NDIM; mu++) {
    // lower(Econ[mu], ggeom[i][j][CENT].gcov, Ecov[mu]);
    lower(Econ[mu], Gcov, Ecov[mu]);
  }
  for (int mu = 0; mu < NDIM; mu++) {
    Ecov[0][mu] *= -1.;
  }
}

void normalize(double *Vcon, double Gcov[NDIM][NDIM]) {
  double norm = 0.;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      norm += Vcon[mu] * Vcon[nu] * Gcov[mu][nu];
    }
  }

  norm = sqrt(fabs(norm));
  for (int mu = 0; mu < NDIM; mu++) {
    Vcon[mu] /= norm;
  }
}

void project_out(double *Vcona, double *Vconb, double Gcov[NDIM][NDIM]) {
  double Vconb_sq = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Vconb_sq += Vconb[mu] * Vconb[nu] * Gcov[mu][nu];
    }
  }

  double adotb = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      adotb += Vcona[mu] * Vconb[nu] * Gcov[mu][nu];
    }
  }

  for (int mu = 0; mu < NDIM; mu++) {
    Vcona[mu] -= Vconb[mu] * adotb / Vconb_sq;
  }
}

// Enforce K.K = 0
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]) {
  double A = Gcov[0][0];
  double B = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    B += 2. * Gcov[mu][0] * K[mu];
  }
  double C = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    for (int nu = 1; nu < NDIM; nu++) {
      C += Gcov[mu][nu] * K[mu] * K[nu];
    }
  }

  K[0] = (-B - sqrt(fabs(B * B - 4. * A * C))) / (2. * A);
}

void normalize_null_cov(double Gcon[NDIM][NDIM], double K[NDIM]) {
  double A = Gcon[0][0];
  double B = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    B += 2. * Gcon[mu][0] * K[mu];
  }
  double C = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    for (int nu = 1; nu < NDIM; nu++) {
      C += Gcon[mu][nu] * K[mu] * K[nu];
    }
  }

  K[0] = (-B + sqrt(fabs(B * B - 4. * A * C))) / (2. * A);
}
#undef SMALL_VECTOR
#endif // RADIATION
