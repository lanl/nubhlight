/******************************************************************************
 *                                                                            *
 * METRIC.C                                                                   *
 *                                                                            *
 * HELPER FUNCTIONS FOR METRIC TENSORS                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Calculate fluxes in direction dir
void primtoflux(double *Pr, struct of_state *q, int dir, int magnetic,
    struct of_geom *geom, double *flux) {
  double mhd[NDIM];

  // Locally static variable because
  // this is used at cell faces, not just cell centers
  // and the global cache is only for cell centers
  static double extra[EOS_NUM_EXTRA];
#pragma omp threadprivate(extra)

  // Particle number flux
  flux[RHO] = Pr[RHO] * q->ucon[dir];

  mhd_calc(Pr, dir, magnetic, q, extra, mhd);

  // MHD stress-energy tensor w/ first index up, second index down
  flux[UU] = mhd[0] + flux[RHO];
  flux[U1] = mhd[1];
  flux[U2] = mhd[2];
  flux[U3] = mhd[3];

  // Dual of Maxwell tensor
  flux[B1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  flux[B2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  flux[B3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];

#if ELECTRONS
  flux[KEL]  = flux[RHO] * Pr[KEL];
  flux[KTOT] = flux[RHO] * Pr[KTOT];
#endif // ELECTRONS

// Passive scalars
#if NVAR_PASSIVE > 0
  PASSLOOP {
    if (PASSTYPE(ipass) == PASSTYPE_NUMBER) {
      flux[ipass] = flux[RHO] * Pr[ipass];
    } else { // PASSTYPE_INTRINSIC
      flux[ipass] = Pr[ipass] * q->ucon[dir];
    }
  }
#endif // passive scalars

  PLOOP flux[ip] *= geom->g;
}

// calculate magnetic field four-vector
void bcon_calc(double *Pr, double *ucon, double *ucov, double *bcon) {
  bcon[0] = Pr[B1] * ucov[1] + Pr[B2] * ucov[2] + Pr[B3] * ucov[3];
  for (int j = 1; j < 4; j++)
    bcon[j] = (Pr[B1 - 1 + j] + bcon[0] * ucon[j]) / ucon[0];
}

// MHD stress-energy tensor with first index up, second index down. A factor of
// sqrt(4 pi) is absorbed into the definition of b.
void mhd_calc(double *Pr, int dir, int magnetic, struct of_state *q,
    double *extra, double *mhd) {
  double r, u, pres, w, bsq, eta, ptot;
  r = fabs(Pr[RHO]) + SMALL;
  u = Pr[UU];

#if EOS == EOS_TYPE_TABLE
  if (!(magnetic)) {
    EOS_SC_fill(Pr, extra);
  }
#endif

  if (magnetic) {
    pres = 0.0;
  } else {
    pres = EOS_pressure_rho0_u(r, u, extra);
  }
  w    = pres + r + u;
  bsq  = dot(q->bcon, q->bcov);
  eta  = w + bsq;
  ptot = pres + 0.5 * bsq;

  for (int mu = 0; mu < NDIM; mu++) {
    mhd[mu] = eta * q->ucon[dir] * q->ucov[mu] + ptot * delta(dir, mu) -
              q->bcon[dir] * q->bcov[mu];
  }
}

// Source terms for equations of motion
void source(double *Ph, struct of_geom *geom, int ii, int jj, double *dU,
    double Dt, double *extra) {
  double          mhd[NDIM][NDIM];
  struct of_state q;

  get_state(Ph, geom, &q);
  mhd_calc(Ph, 0, 0, &q, extra, mhd[0]);
  mhd_calc(Ph, 1, 0, &q, extra, mhd[1]);
  mhd_calc(Ph, 2, 0, &q, extra, mhd[2]);
  mhd_calc(Ph, 3, 0, &q, extra, mhd[3]);

  // Contract mhd stress tensor with connection
  PLOOP dU[ip] = 0.;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      dU[UU] += mhd[mu][nu] * conn[ii][jj][nu][0][mu];
      dU[U1] += mhd[mu][nu] * conn[ii][jj][nu][1][mu];
      dU[U2] += mhd[mu][nu] * conn[ii][jj][nu][2][mu];
      dU[U3] += mhd[mu][nu] * conn[ii][jj][nu][3][mu];
    }
  }

  PLOOP dU[ip] *= geom->g;
}

// Returns b.b (twice magnetic pressure)
double bsq_calc(double *Pr, struct of_geom *geom) {
  struct of_state q;

  get_state(Pr, geom, &q);
  return (dot(q.bcon, q.bcov));
}

// Calculate ucon, ucov, bcon, bcov from primitive variables
void get_state(double *Pr, struct of_geom *geom, struct of_state *q) {
  ucon_calc(Pr, geom, q->ucon);
  lower(q->ucon, geom->gcov, q->ucov);
  bcon_calc(Pr, q->ucon, q->ucov, q->bcon);
  lower(q->bcon, geom->gcov, q->bcov);
}

// Find contravariant four-velocity
void ucon_calc(double *Pr, struct of_geom *geom, double *ucon) {
  double alpha, gamma;

  if (mhd_gamma_calc(Pr, geom, &gamma)) {
    fprintf(stderr, "\nucon_calc(): gamma failure\n");
    fail(FAIL_GAMMA);
  }

  alpha   = geom->alpha;
  ucon[0] = gamma / alpha;
  for (int j = 1; j < NDIM; j++) {
    ucon[j] = Pr[U1 + j - 1] - gamma * alpha * geom->gcon[0][j];
  }
}

// Find gamma-factor wrt normal observer
int mhd_gamma_calc(double *Pr, struct of_geom *geom, double *gamma) {
  double qsq;

  qsq = geom->gcov[1][1] * Pr[U1] * Pr[U1] +
        geom->gcov[2][2] * Pr[U2] * Pr[U2] +
        geom->gcov[3][3] * Pr[U3] * Pr[U3] +
        2. * (geom->gcov[1][2] * Pr[U1] * Pr[U2] +
                 geom->gcov[1][3] * Pr[U1] * Pr[U3] +
                 geom->gcov[2][3] * Pr[U2] * Pr[U3]);

  if (qsq < 0.) {
    if (fabs(qsq) > 1.E-10) { // Then assume not just machine precision
      fprintf(stderr, "gamma_calc():  failed: i,j,qsq = %d %d %28.18e \n",
          icurr, jcurr, qsq);
      fprintf(stderr, "v[1-3] = %28.18e %28.18e %28.18e  \n", Pr[U1], Pr[U2],
          Pr[U3]);
      *gamma = 1.;
      return (1);
    } else {
      qsq = 1.E-10; // Set floor
    }
  }

  *gamma = sqrt(1. + qsq);

  return (0);
}

// Calculate components of magnetosonic velocity from primitive variables
void mhd_vchar(double *Pr, struct of_state *q, struct of_geom *geom, int js,
    double *vmax, double *vmin) {
  double discr, vp, vm, bsq, ee, ef, va2, cs, cs2, cms2, rho, u;
  double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;

  // Locally static variable because
  // this is used at cell faces, not just cell centers
  // and the global cache is only for cell centers
  static double extra[EOS_NUM_EXTRA];
#pragma omp threadprivate(extra)
#if EOS == EOS_TYPE_TABLE
  EOS_SC_fill(Pr, extra);
#endif

  for (int mu = 0; mu < NDIM; mu++) {
    Acov[mu] = 0.;
  }
  Acov[js] = 1.;
  raise(Acov, geom->gcon, Acon);

  for (int mu = 0; mu < NDIM; mu++) {
    Bcov[mu] = 0.;
  }
  Bcov[0] = 1.;
  raise(Bcov, geom->gcon, Bcon);

  // Find fast magnetosonic speed
  bsq = dot(q->bcon, q->bcov);
  rho = fabs(Pr[RHO] + SMALL);
  u   = Pr[UU];
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  u = fabs(u + SMALL);
#endif
  ef  = EOS_enthalpy_rho0_u(rho, u, extra);
  cs  = EOS_sound_speed_rho0_u(rho, u, extra);
  cs2 = cs * cs;
  // alfven velocity must be positive and < 1
  ee  = bsq + fabs(ef) + SMALL;
  va2 = MY_MAX(0.0, bsq / ee);

  cms2 = cs2 + va2 - cs2 * va2;
  if (fabs(cms2) < SMALL)
    cms2 = SMALL;

  // Sanity checks
  if (cs2 < 0. || isnan(cs2) || cms2 < 0. || isnan(cms2)) {
    double temp = EOS_temperature(Pr[RHO], Pr[UU], extra);
    fprintf(stderr,
        "\n\ncms2 too small:\n"
        "\tcs2  = %g\n"
        "\tcms2 = %g\n"
        "\tva2  = %g\n"
        "\trho  = %g\n"
        "\tu    = %g\n"
        "\te    = %g\n"
        "\tef   = %g\n"
        "\tTemp = %g\n"
        "\n",
        cs2, cms2, va2, rho, u, u / rho, ef, temp);
#if POLYTROPE_FALLBACK
    double K, Gam;
    double lrho  = extra[EOS_LRHO];
    double lT    = extra[EOS_LT];
    double ye    = extra[EOS_YE];
    double press = EOS_pressure_rho0_u(rho, u, extra);
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    fprintf(stderr,
        "Polytrope fallback:\n"
        "\trho_poly_thresh = %e\n"
        "\tK               = %e\n"
        "\tGam             = %e\n"
        "\tPres            = %e\n",
        rho_poly_thresh, K, Gam, press);
#endif
#if GAMMA_FALLBACK
    fprintf(stderr, "cs2 w/ gam = %e\n", EOS_Gamma_sound_speed_rho0_u(rho, u));
#endif
    fail(FAIL_COEFF_NEG);
    cms2 = SMALL;
  }
  if (cms2 > 1. || cs2 > 1.) {
    double temp = EOS_temperature(Pr[RHO], Pr[UU], extra);
    fprintf(stderr,
        "\n\ncms2 too large:\n"
        "\tcs2  = %g\n"
        "\tcms2 = %g\n"
        "\tva2  = %g\n"
        "\trho  = %g\n"
        "\tu    = %g\n"
        "\te    = %g\n"
        "\tef   = %g\n"
        "\tTemp = %g\n"
        "\n",
        cs2, cms2, va2, rho, u, u / rho, ef, temp);
#if POLYTROPE_FALLBACK
    double K, Gam;
    double lrho  = extra[EOS_LRHO];
    double lT    = extra[EOS_LT];
    double ye    = extra[EOS_YE];
    double press = EOS_pressure_rho0_u(rho, u, extra);
    EOS_SC_get_polytrope(lrho, lT, ye, &K, &Gam);
    fprintf(stderr,
        "Polytrope fallback:\n"
        "\trho_poly_thresh = %e\n"
        "\tK               = %e\n"
        "\tGam             = %e\n"
        "\tPres            = %e\n",
        rho_poly_thresh, K, Gam, press);
#endif
#if GAMMA_FALLBACK
    fprintf(stderr, "cs2 w/ gam = %e\n", EOS_Gamma_sound_speed_rho0_u(rho, u));
#endif
    fail(FAIL_COEFF_SUP);
    cms2 = 1.;
  }

  // Require that speed of wave measured by observer q->ucon is cms2
  Asq  = dot(Acon, Acov);
  Bsq  = dot(Bcon, Bcov);
  Au   = dot(Acov, q->ucon);
  Bu   = dot(Bcov, q->ucon);
  AB   = dot(Acon, Bcov);
  Au2  = Au * Au;
  Bu2  = Bu * Bu;
  AuBu = Au * Bu;

  A = Bu2 - (Bsq + Bu2) * cms2;
  B = 2. * (AuBu - (AB + AuBu) * cms2);
  C = Au2 - (Asq + Au2) * cms2;

  discr = B * B - 4. * A * C;
  if ((discr < 0.0) && (discr > -1.e-10)) {
    discr = 0.0;
  } else if (discr < -1.e-10) {
    fprintf(stderr,
        "\nBad discriminant!\n"
        "\t [A, B, C, discr, cms2]:\n"
        "\t\t %g %g %g %g %g\n"
        "\t q->ucon: %g %g %g %g\n"
        "\t q->bcon: %g %g %g %g\n"
        "\t Acon: %g %g %g %g\n"
        "\t Bcon: %g %g %g %g\n",
        A, B, C, discr, cms2, q->ucon[0], q->ucon[1], q->ucon[2], q->ucon[3],
        q->bcon[0], q->bcon[1], q->bcon[2], q->bcon[3], Acon[0], Acon[1],
        Acon[2], Acon[3], Bcon[0], Bcon[1], Bcon[2], Bcon[3]);
    fail(FAIL_VCHAR_DISCR);
    discr = 0.;
  }

  discr = sqrt(discr);
  vp    = -(-B + discr) / (2. * A);
  vm    = -(-B - discr) / (2. * A);

  if (vp > vm) {
    *vmax = vp;
    *vmin = vm;
  } else {
    *vmax = vm;
    *vmin = vp;
  }
}
