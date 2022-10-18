/******************************************************************************
 *                                                                            *
 * CURRENT.C                                                                  *
 *                                                                            *
 * CALCULATE CURRENT FROM FLUID VARIABLES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double Fcon_calc(double *prim, int mu, int nu, int i, int j, int k);
int    antisym(int a, int b, int c, int d);
int    pp(int n, int *P);

static grid_prim_type Pa;

void current_calc() {
  double gF0p[NDIM], gF0m[NDIM], gF1p[NDIM], gF1m[NDIM], gF2p[NDIM], gF2m[NDIM];
  double gF3p[NDIM], gF3m[NDIM];

  if (nstep == 0) {
    ZLOOP {
      for (int mu = 0; mu < NDIM; mu++)
        jcon[i][j][k][mu] = 0.;
    }
    return;
  }

  // Calculate time-centered P
  ZSLOOP(-1, N1 - 1 + NG, -1, N2 - 1 + NG, -1, N3 - 1 + NG) {
    PLOOP { Pa[i][j][k][ip] = 0.5 * (P[i][j][k][ip] + Psave[i][j][k][ip]); }
  }

  // Calculate j^{\mu} using centered differences for active zones

  ZLOOP {
    for (int mu = 0; mu < NDIM; mu++)
      jcon[i][j][k][mu] = 0.;

    // Get sqrt{-g}*F^{mu nu} at neighboring points

    // X0
    for (int mu = 0; mu < NDIM; mu++) {
      gF0p[mu] = Fcon_calc(P[i][j][k], 0, mu, i, j, k);
      gF0m[mu] = Fcon_calc(Psave[i][j][k], 0, mu, i, j, k);
    }

    // X1
    for (int mu = 0; mu < NDIM; mu++) {
      gF1p[mu] = Fcon_calc(Pa[i + 1][j][k], 1, mu, i + 1, j, k);
      gF1m[mu] = Fcon_calc(Pa[i - 1][j][k], 1, mu, i - 1, j, k);
    }

    // X2
    for (int mu = 0; mu < NDIM; mu++) {
      gF2p[mu] = Fcon_calc(Pa[i][j + 1][k], 2, mu, i, j + 1, k);
      gF2m[mu] = Fcon_calc(Pa[i][j - 1][k], 2, mu, i, j - 1, k);
    }

    // X3
    for (int mu = 0; mu < NDIM; mu++) {
      gF3p[mu] = Fcon_calc(Pa[i][j][k + 1], 3, mu, i, j, k + 1);
      gF3m[mu] = Fcon_calc(Pa[i][j][k - 1], 3, mu, i, j, k - 1);
    }

    // Difference: D_mu F^{mu nu} = 4 \pi j^nu
    for (int mu = 0; mu < NDIM; mu++) {
      jcon[i][j][k][mu] = (1. / (4. * M_PI * ggeom[i][j][CENT].g)) *
                          ((gF0p[mu] - gF0m[mu]) / dtsave +
                              (gF1p[mu] - gF1m[mu]) / (2. * dx[1]) +
                              (gF2p[mu] - gF2m[mu]) / (2. * dx[2]) +
                              (gF3p[mu] - gF3m[mu]) / (2. * dx[3]));
    }
  }
}

// Return mu, nu component of contravarient Maxwell tensor at grid zone i, j, k
double Fcon_calc(double *prim, int mu, int nu, int i, int j, int k) {
  double ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
  double Fcon, gFcon, dFcon;

  if (mu == nu)
    return 0.;

  ucon_calc(prim, &(ggeom[i][j][CENT]), ucon);
  lower(ucon, ggeom[i][j][CENT].gcov, ucov);
  bcon_calc(prim, ucon, ucov, bcon);
  lower(bcon, ggeom[i][j][CENT].gcov, bcov);

  Fcon = 0.;
  for (int kap = 0; kap < NDIM; kap++) {
    for (int lam = 0; lam < NDIM; lam++) {
      dFcon = (-1. / ggeom[i][j][CENT].g) * antisym(mu, nu, kap, lam) *
              ucov[kap] * bcov[lam];
      Fcon += dFcon;
    }
  }

  gFcon = Fcon * ggeom[i][j][CENT].g;

  return gFcon;
}

// Completely antisymmetric 4D symbol
int antisym(int a, int b, int c, int d) {
  // Check for valid permutation
  if (a < 0 || a > 3)
    return 100;
  if (b < 0 || b > 3)
    return 100;
  if (c < 0 || c > 3)
    return 100;
  if (d < 0 || d > 3)
    return 100;

  // Entries different?
  if (a == b)
    return 0;
  if (a == c)
    return 0;
  if (a == d)
    return 0;
  if (b == c)
    return 0;
  if (b == d)
    return 0;
  if (c == d)
    return 0;

  // Determine parity of permutation
  int p[4] = {a, b, c, d};

  return pp(4, p);
}

// Due to Norm Hardy; good for general n
int pp(int n, int *P) {
  int x;
  int p = 0;
  int v[n];

  for (int j = 0; j < n; j++)
    v[j] = 0;

  for (int j = 0; j < n; j++) {
    if (v[j]) {
      p++;
    } else {
      x = j;
      do {
        x    = P[x];
        v[x] = 1;
      } while (x != j);
    }
  }

  if (p % 2 == 0) {
    return 1;
  } else {
    return -1;
  }
}
