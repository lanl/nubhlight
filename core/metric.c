/******************************************************************************
 *                                                                            *
 * METRIC.C                                                                   *
 *                                                                            *
 * HELPER FUNCTIONS FOR METRIC TENSORS                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double gcon_func(double gcov[][NDIM], double gcon[][NDIM]) {
  double gdet = invert(&gcov[0][0], &gcon[0][0]);
  return sqrt(fabs(gdet));
}

// Set the spatial discretization in numerical derivatives
#define DELTA 1.e-5

// Calculate connection coefficient \Gamma^{i}_{j,k} = conn[..][i][j][k]
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM]) {
  double tmp[NDIM][NDIM][NDIM];
  double Xh[NDIM], Xl[NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];

  for (int k = 0; k < NDIM; k++) {
    for (int l = 0; l < NDIM; l++) {
      Xh[l] = X[l];
    }
    for (int l = 0; l < NDIM; l++) {
      Xl[l] = X[l];
    }
    Xh[k] += DELTA;
    Xl[k] -= DELTA;
    set_gcov(Xh, gh);
    set_gcov(Xl, gl);

    for (int i = 0; i < NDIM; i++) {
      for (int j = 0; j < NDIM; j++) {
        conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  } // for k

  // Rearrange to find \Gamma_{ijk}
  for (int i = 0; i < NDIM; i++) {
    for (int j = 0; j < NDIM; j++) {
      for (int k = 0; k < NDIM; k++) {
        tmp[i][j][k] = 0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);
      }
    }
  }

  // Raise index to get \Gamma^i_{jk}
  for (int i = 0; i < NDIM; i++) {
    for (int j = 0; j < NDIM; j++) {
      for (int k = 0; k < NDIM; k++) {
        conn[i][j][k] = 0.;
        for (int l = 0; l < NDIM; l++)
          conn[i][j][k] += geom->gcon[i][l] * tmp[l][j][k];
      }
    }
  }
}

// Lower a contravariant rank-1 tensor to a covariant one
void lower(double ucon[NDIM], double gcov[NDIM][NDIM], double ucov[NDIM]) {
  ucov[0] = gcov[0][0] * ucon[0] + gcov[0][1] * ucon[1] + gcov[0][2] * ucon[2] +
            gcov[0][3] * ucon[3];
  ucov[1] = gcov[1][0] * ucon[0] + gcov[1][1] * ucon[1] + gcov[1][2] * ucon[2] +
            gcov[1][3] * ucon[3];
  ucov[2] = gcov[2][0] * ucon[0] + gcov[2][1] * ucon[1] + gcov[2][2] * ucon[2] +
            gcov[2][3] * ucon[3];
  ucov[3] = gcov[3][0] * ucon[0] + gcov[3][1] * ucon[1] + gcov[3][2] * ucon[2] +
            gcov[3][3] * ucon[3];
}

// Raise a covariant rank-1 tensor to a contravariant one
void raise(double ucov[NDIM], double gcon[NDIM][NDIM], double ucon[NDIM]) {
  ucon[0] = gcon[0][0] * ucov[0] + gcon[0][1] * ucov[1] + gcon[0][2] * ucov[2] +
            gcon[0][3] * ucov[3];
  ucon[1] = gcon[1][0] * ucov[0] + gcon[1][1] * ucov[1] + gcon[1][2] * ucov[2] +
            gcon[1][3] * ucov[3];
  ucon[2] = gcon[2][0] * ucov[0] + gcon[2][1] * ucov[1] + gcon[2][2] * ucov[2] +
            gcon[2][3] * ucov[3];
  ucon[3] = gcon[3][0] * ucov[0] + gcon[3][1] * ucov[1] + gcon[3][2] * ucov[2] +
            gcon[3][3] * ucov[3];
}

// Load local geometry into structure geom
struct of_geom *get_geometry(int ii, int jj, int kk, int loc) {
  // Store current spatial indices in case of errors
  icurr = ii;
  jcurr = jj;
  kcurr = kk;

  return (&(ggeom[ii][jj][loc]));
}
#undef DELTA

// Boyer-Lindquist metric functions
void blgset(int i, int j, struct of_geom *geom) {
  double r, th, X[NDIM];

  coord(i, j, 0, CENT, X);
  bl_coord(X, &r, &th);

  if (th < 0)
    th *= -1.;
  if (th > M_PI)
    th = 2. * M_PI - th;

  geom->g = bl_gdet_func(r, th);
  bl_gcov_func(r, th, geom->gcov);
  bl_gcon_func(r, th, geom->gcon);
}

double bl_gdet_func(double r, double th) {
  double a2, r2;

  a2 = a * a;
  r2 = r * r;
  return (r * r * fabs(sin(th)) * (1. + 0.5 * (a2 / r2) * (1. + cos(2. * th))));
}

void bl_gcov_func(double r, double th, double gcov[][NDIM]) {
  double sth, cth, s2, a2, r2, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcov[mu][nu] = 0.;
    }
  }

  sth = fabs(sin(th));
  s2  = sth * sth;
  cth = cos(th);
  a2  = a * a;
  r2  = r * r;
  DD  = 1. - 2. / r + a2 / r2;
  mu  = 1. + a2 * cth * cth / r2;

  gcov[0][0] = -(1. - 2. / (r * mu));
  gcov[0][3] = -2. * a * s2 / (r * mu);
  gcov[3][0] = gcov[0][3];
  gcov[1][1] = mu / DD;
  gcov[2][2] = r2 * mu;
  gcov[3][3] = r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));
}

void bl_gcon_func(double r, double th, double gcon[][NDIM]) {
  double sth, cth, a2, r2, r3, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcon[mu][nu] = 0.;
    }
  }

  sth = sin(th);
  cth = cos(th);

#if (COORDSINGFIX)
  if (fabs(sth) < SINGSMALL) {
    if (sth >= 0)
      sth = SINGSMALL;
    if (sth < 0)
      sth = -SINGSMALL;
  }
#endif

  a2 = a * a;
  r2 = r * r;
  r3 = r2 * r;
  DD = 1. - 2. / r + a2 / r2;
  mu = 1. + a2 * cth * cth / r2;

  gcon[0][0] = -1. - 2. * (1. + a2 / r2) / (r * DD * mu);
  gcon[0][3] = -2. * a / (r3 * DD * mu);
  gcon[3][0] = gcon[0][3];
  gcon[1][1] = DD / mu;
  gcon[2][2] = 1. / (r2 * mu);
  gcon[3][3] = (1. - 2. / (r * mu)) / (r2 * sth * sth * DD);
}

double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2) {
  return m[4 * r0 + c0] * (m[4 * r1 + c1] * m[4 * r2 + c2] -
                              m[4 * r2 + c1] * m[4 * r1 + c2]) -
         m[4 * r0 + c1] * (m[4 * r1 + c0] * m[4 * r2 + c2] -
                              m[4 * r2 + c0] * m[4 * r1 + c2]) +
         m[4 * r0 + c2] * (m[4 * r1 + c0] * m[4 * r2 + c1] -
                              m[4 * r2 + c0] * m[4 * r1 + c1]);
}

void adjoint(double m[16], double adjOut[16]) {
  adjOut[0] = MINOR(m, 1, 2, 3, 1, 2, 3);
  adjOut[1] = -MINOR(m, 0, 2, 3, 1, 2, 3);
  adjOut[2] = MINOR(m, 0, 1, 3, 1, 2, 3);
  adjOut[3] = -MINOR(m, 0, 1, 2, 1, 2, 3);

  adjOut[4] = -MINOR(m, 1, 2, 3, 0, 2, 3);
  adjOut[5] = MINOR(m, 0, 2, 3, 0, 2, 3);
  adjOut[6] = -MINOR(m, 0, 1, 3, 0, 2, 3);
  adjOut[7] = MINOR(m, 0, 1, 2, 0, 2, 3);

  adjOut[8]  = MINOR(m, 1, 2, 3, 0, 1, 3);
  adjOut[9]  = -MINOR(m, 0, 2, 3, 0, 1, 3);
  adjOut[10] = MINOR(m, 0, 1, 3, 0, 1, 3);
  adjOut[11] = -MINOR(m, 0, 1, 2, 0, 1, 3);

  adjOut[12] = -MINOR(m, 1, 2, 3, 0, 1, 2);
  adjOut[13] = MINOR(m, 0, 2, 3, 0, 1, 2);
  adjOut[14] = -MINOR(m, 0, 1, 3, 0, 1, 2);
  adjOut[15] = MINOR(m, 0, 1, 2, 0, 1, 2);
}

double determinant(double m[16]) {
  return m[0] * MINOR(m, 1, 2, 3, 1, 2, 3) - m[1] * MINOR(m, 1, 2, 3, 0, 2, 3) +
         m[2] * MINOR(m, 1, 2, 3, 0, 1, 3) - m[3] * MINOR(m, 1, 2, 3, 0, 1, 2);
}

double invert(double *m, double *invOut) {
  adjoint(m, invOut);

  double det     = determinant(m);
  double inv_det = 1. / det;
  for (int i = 0; i < 16; ++i) {
    invOut[i] = invOut[i] * inv_det;
  }

  return det;
}
