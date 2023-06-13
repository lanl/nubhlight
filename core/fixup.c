/******************************************************************************
 *                                                                            *
 * FIXUP.C                                                                    *
 *                                                                            *
 * REPAIR INTEGRATION FAILURES                                                *
 * DRIFT FRAME DENSITY AND INTERNAL ENERGY FLOORS FOLLOWING A. TCHEKHOVSKOY
 *   SEE RESSLER ET AL. 2016
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#pragma GCC diagnostic ignored "-Wstringop-overflow"

// Apply floors to density, internal energy
void fixup(grid_prim_type Pv, grid_eosvar_type extra) {
  timer_start(TIMER_FIXUP);
#pragma omp parallel for collapse(3) // schedule(dynamic)
  ZLOOP fixup1zone(i, j, k, Pv[i][j][k], extra[i][j][k]);
  timer_stop(TIMER_FIXUP);
}

void ucon_to_utcon(
    double ucon[NDIM], struct of_geom *geom, double utcon[NDIM]) {
  double alpha, beta[NDIM], gamma;

  alpha = 1. / sqrt(-geom->gcon[0][0]);
  for (int i = 1; i < NDIM; i++) {
    beta[i] = geom->gcon[0][i] * alpha * alpha;
  }
  gamma = alpha * ucon[0];

  utcon[0] = 0.;
  for (int i = 1; i < NDIM; i++) {
    utcon[i] = ucon[i] + gamma * beta[i] / alpha;
  }
}

void ut_calc_3vel(double vcon[NDIM], struct of_geom *geom, double *ut) {
  double AA, BB, CC, DD, one_over_alpha_sq;

  AA = geom->gcov[0][0];
  BB = 2. * (geom->gcov[0][1] * vcon[1] + geom->gcov[0][2] * vcon[2] +
                geom->gcov[0][3] * vcon[3]);
  CC = geom->gcov[1][1] * vcon[1] * vcon[1] +
       geom->gcov[2][2] * vcon[2] * vcon[2] +
       geom->gcov[3][3] * vcon[3] * vcon[3] +
       2. * (geom->gcov[1][2] * vcon[1] * vcon[2] +
                geom->gcov[1][3] * vcon[1] * vcon[3] +
                geom->gcov[2][3] * vcon[2] * vcon[3]);
  DD = 1. / (AA + BB + CC);

  one_over_alpha_sq = -geom->gcon[0][0];

  if (DD < one_over_alpha_sq) {
    DD = one_over_alpha_sq;
  }

  *ut = sqrt(DD);
}

// ORIGINAL FLOORS
/*
double rhoscal, uscal;
rhoscal = pow(r,-1.5);
uscal = rhoscal/r;
rhoflr = RHOMIN*rhoscal;
uflr = UUMIN*uscal;
if (rhoflr < RHOMINLIMIT) rhoflr = RHOMINLIMIT;
if (uflr < UUMINLIMIT) uflr = UUMINLIMIT;
pv[RHO] = MY_MAX(rhoflr, pv[RHO]);
pv[UU] = MY_MAX(uflr, pv[UU]);
if (mhd_gamma_calc(pv, geom, &gamma)) {
  pflag[i][j][k] = -333;
} else {
  if (gamma > GAMMAMAX) {
    f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
    pv[U1] *= f;
    pv[U2] *= f;
    pv[U3] *= f;
  }
}
return;*/

double get_scale(int i, int j, int k) {
#if METRIC == MINKOWSKI
  { return 1.e-2; }
#elif METRIC == MKS
  {
    double scale, r, th, X[NDIM];
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    if (r <= FLR_R0) {
      scale = pow(r, -FLR_POWER1);
    } else {
      scale = pow(FLR_R0, FLR_POWER2 - FLR_POWER1) * pow(r, -FLR_POWER2);
    }
    return scale;
  }
#else
  {
    fprintf(stderr, "[fixup1zone]: Unknown metric!\n");
    exit(1);
  }
#endif // METRIC
}

void fixup1zone(
    int i, int j, int k, double pv[NVAR], double extra[EOS_NUM_EXTRA]) {
  double          rhoflr, uflr, f, gamma;
  struct of_geom *geom;
  struct of_state q;
  double          bsq;
  double          pv_prefloor[NVAR];
  PLOOP           pv_prefloor[ip] = pv[ip];

  double scale = get_scale(i, j, k);

#if EOS == EOS_TYPE_TABLE
  EOS_SC_fill(pv, extra);
#endif

  // Enhance floors in case of large magnetic energy density
  geom = get_geometry(i, j, k, CENT);
  get_state(pv, geom, &q);
  bsq = dot(q.bcon, q.bcov);
  EOS_set_floors(scale, pv[RHO], pv[UU], bsq, &rhoflr, &uflr, extra);

#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK && COLD_FLOORS
  double rhosave = pv[RHO];
#endif // POLYTROPE_FALLBACK

  if (rhoflr > pv[RHO] || uflr > pv[UU]) { // Apply floors
    fixup_required[i][j][k] = 1;
    double trans =
        10. * bsq / (MY_MIN(fabs(pv[RHO]), fabs(pv[UU])) + SMALL) - 1.;
    if (trans > 0.) { // Strongly magnetized region; use drift frame floors
      pv[RHO] = MY_MAX(pv[RHO], rhoflr);
      pv[UU]  = MY_MAX(pv[UU], uflr);

      double betapar, betasqmax, betasq, gamma, ucondr[NDIM], Bcon[NDIM];
      double Bcov[NDIM], udotB, Bsq, B, wold, QdotB, wnew, x, vpar;
      double one_over_ucondr_t, vcon[NDIM], ucon[NDIM], ut, utcon[NDIM];
      trans = MY_MIN(trans, 1.);

      // Set velocity to drift velocity
      betapar   = -q.bcon[0] / ((bsq + SMALL) * q.ucon[0]);
      betasq    = betapar * betapar * bsq;
      betasqmax = 1. - 1. / (GAMMAMAX * GAMMAMAX);
      betasq    = MY_MIN(betasq, betasqmax);

      gamma = 1. / sqrt(1. - betasq);

      DLOOP1 ucondr[mu] = gamma * (q.ucon[mu] + betapar * q.bcon[mu]);

      Bcon[0] = 0.;
      for (int i = 1; i < NDIM; i++) {
        Bcon[i] = pv[B1 - 1 + i];
      }
      lower(Bcon, geom->gcov, Bcov);
      udotB = dot(q.ucon, Bcov);
      Bsq   = dot(Bcon, Bcov);
      B     = sqrt(Bsq);

// Enthalpy before floors are applied
#if EOS == EOS_TYPE_TABLE
      EOS_SC_fill(pv_prefloor, extra);
#endif
      wold  = EOS_enthalpy_rho0_u(pv_prefloor[RHO], pv_prefloor[UU], extra);
      QdotB = udotB * wold * q.ucon[0];

// Apply floors to enthalpy and recompute parallel velocity
#if EOS == EOS_TYPE_TABLE
      EOS_SC_fill(pv, extra);
#endif
      wnew = EOS_enthalpy_rho0_u(pv[RHO], pv[UU], extra);
      x    = 2. * QdotB / (B * wnew * ucondr[0] + SMALL);
      vpar = x / (ucondr[0] * (1. + sqrt(1. + x * x)));

      one_over_ucondr_t = 1. / ucondr[0];

      vcon[0] = 1.;
      for (int i = 1; i < NDIM; i++) {
        vcon[i] = vpar * Bcon[i] / (B + SMALL) + ucondr[i] * one_over_ucondr_t;
      }

      ut_calc_3vel(vcon, geom, &ut);

      DLOOP1 ucon[mu] = ut * vcon[mu];

      ucon_to_utcon(ucon, geom, utcon);

      // Convert 3-velocity to relative 4-velocity and store in primitives
      for (int i = 1; i < NDIM; i++) {
        pv[i + UU] = utcon[i] * trans + pv_prefloor[i + UU] * (1. - trans);
      }

    } else { // Weakly magnetized region; use normal observer frame floors
      double Padd[NVAR], Uadd[NVAR];
      PLOOP  Padd[ip] = 0.;
      PLOOP  Uadd[ip] = 0.;
      Padd[RHO]       = MY_MAX(0.0, rhoflr - pv[RHO]);
      Padd[UU]        = MY_MAX(0.0, uflr - pv[UU]);

      get_state(Padd, &ggeom[i][j][CENT], &q);

      primtoflux(Padd, &q, 0, 0, &ggeom[i][j][CENT], Uadd);

      double Utot[NVAR];
      get_state(pv, &ggeom[i][j][CENT], &q);
      primtoflux(pv, &q, 0, 0, &ggeom[i][j][CENT], Utot);

      PLOOP Utot[ip] += Uadd[ip];

      PLOOP pv[ip] += Padd[ip];

      // Record fails here?
      Utoprim(Utot, &ggeom[i][j][CENT], pv);
    }
  }

#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK && COLD_FLOORS
  // set to zero temperature anywhere in the floor region
  if (rhosave < RHOEPS * rho_poly_thresh) {
    pv[UU] = EOS_SC_get_minu(pv[RHO], pv[YE], scale);
  }
#endif

#if NVAR_PASSIVE > 0
  fixup_passive(i, j, k, pv, pv_prefloor);
#endif

#if ELECTRONS && EOS == EOS_TYPE_GAMMA
  // Reset entropy after floors
  pv[KTOT] = EOS_Gamma_entropy_rho0_u(pv[RHO], pv[UU]);

  // Set KTOTMAX to 3 by controlling u, to avoid anomalous cooling from funnel
  // wall
  double KTOTMAX = 3.;
  if (pv[KTOT] > KTOTMAX) {
    pv[UU]   = KTOTMAX * pow(pv[RHO], gam) / (gam - 1.);
    pv[KTOT] = KTOTMAX;
  }
#endif // ELECTRONS

  // Limit gamma with respect to normal observer
  if (mhd_gamma_calc(pv, geom, &gamma)) {
    pflag[i][j][k] = -333;
  } else {
    if (gamma > GAMMAMAX) {
      f = sqrt((GAMMAMAX * GAMMAMAX - 1.) / (gamma * gamma - 1.));
      pv[U1] *= f;
      pv[U2] *= f;
      pv[U3] *= f;
    }
  }
}

static grid_prim_type Pv_tmp;
static grid_int_type  pflag_tmp, pflag_save;

// Replace bad points with values interpolated from neighbors
#define FLOOP for (int ip = 0; ip < B1; ip++)
void fixup_utoprim(grid_prim_type Pv, grid_eosvar_type extra) {
  timer_start(TIMER_FIXUP);
  int    bad;
  double sum[B1], wsum;

// Flip the logic of the pflag[] so that it now indicates which cells are good
#pragma omp parallel for collapse(3)
  ZSLOOP(-NG, (N1 - 1 + NG), -NG, (N2 - 1 + NG), -NG, (N3 - 1 + NG)) {
    pflag_save[i][j][k] = pflag[i][j][k];
    pflag[i][j][k]      = !pflag[i][j][k];
  }

  // Make sure we are not using ill defined corner regions
  for (int i = 0; i < NG; i++) {
    for (int j = 0; j < NG; j++) {
      for (int k = 0; k < NG; k++) {
        pflag[i][j][k]                     = 0;
        pflag[i + N1 + NG][j][k]           = 0;
        pflag[i][j + N2 + NG][k]           = 0;
        pflag[i][j][k + N3 + NG]           = 0;
        pflag[i + N1 + NG][j + N2 + NG][k] = 0;
        pflag[i + N1 + NG][j][k + N3 + NG] = 0;
        pflag[i][j + N2 + NG][k + N3 + NG] = 0;
        // pflag[i+N1+NG][j+N2+NG][k+N3-1+NG] = 0;
        pflag[i + N1 + NG][j + N2 + NG][k + N3 + NG] = 0;
      }
    }
  }

  // Fix the interior points first
  int fail_consec = 0;
  do {
    bad = 0;

#pragma omp parallel for collapse(3)
    ZSLOOP(-NG, N1 + NG - 1, -NG, N2 + NG - 1, -NG, N3 + NG - 1)
    FLOOP Pv_tmp[i][j][k][ip] = Pv[i][j][k][ip];
#pragma omp parallel for collapse(3)
    ZSLOOP(-NG, N1 + NG - 1, -NG, N2 + NG - 1, -NG, N3 + NG - 1)
    pflag_tmp[i][j][k] = pflag[i][j][k];

    //#pragma omp parallel for collapse(3) reduction(+:bad)
    ZSLOOP(0, (N1 - 1), 0, (N2 - 1), 0, (N3 - 1)) {
      if (pflag_tmp[i][j][k] == 0) {
        wsum          = 0.;
        FLOOP sum[ip] = 0.;
        for (int l = -1; l < 2; l++) {
          for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
              double w = 1. / (abs(l) + abs(m) + abs(n) + 1) *
                         pflag_tmp[i + l][j + m][k + n];
              wsum += w;
              FLOOP sum[ip] += w * Pv_tmp[i + l][j + m][k + n][ip];
            }
          }
        }
        if (wsum < 1.e-10) {
          // No usable neighbors. Average over all neighbors.
          fail_consec++;
          fprintf(stderr,
              "[%i][istart=%i] fixup_utoprim problem: No usable neighbors!\n",
              mpi_myrank(), global_start[1]);
          fprintf(stderr, "i j k = %i %i %i pflag = %d wsum = %e\n", i, j, k,
              pflag_save[i][j][k], wsum);
          // exit(-1); // DEBUG

          for (int l = -1; l < 2; l++) {
            for (int m = -1; m < 2; m++) {
              for (int n = -1; n < 2; n++) {
                double w = 1. / (abs(l) + abs(m) + abs(n) + 1);
                FLOOP  sum[ip] += w * Pv_tmp[i + l][j + m][k + n][ip];
              }
            }
          }
          Pv_tmp[i][j][k][U1] = 0.;
          Pv_tmp[i][j][k][U2] = 0.;
          Pv_tmp[i][j][k][U3] = 0.;
          bad++;
          continue;
        }
        fail_consec           = 0;
        FLOOP Pv[i][j][k][ip] = sum[ip] / wsum;

        // Cell is fixed, can now use for other interpolations
        pflag[i][j][k] = 1;

        fixup1zone(i, j, k, Pv[i][j][k], extra[i][j][k]);
      }
    }
  } while (bad > 0 && fail_consec < N1 * N2 * N3);

  if (fail_consec == N1 * N2 * N3)
    fixup(Pv, extra);

  timer_stop(TIMER_FIXUP);
}
#undef FLOOP
