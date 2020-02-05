/******************************************************************************
 *                                                                            *
 * PUSH_SUPERPHOTONS.C                                                        *
 *                                                                            *
 * INTEGRATE SUPERPHOTON GEODESICS                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
// void push(double dt, struct of_photon *ph);
void push(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM], double dt);
void get_X_source(double Kcon[NDIM], double src[NDIM]);
void get_K_source(
    double X[NDIM], double Kcov[NDIM], double Kcon[NDIM], double src[NDIM]);

#define MAX_SUBDIV (10)

int push_X_K(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
    grid_prim_type P, grid_prim_type Prad, double KdotKprev, int type,
    double dtpush) {
  int    nsubdiv = 0;
  double Xprev[NDIM], Kcovprev[NDIM], Kconprev[NDIM], KdotK = 0.;
  int    good_push;

  // Store initial X and K
  for (int mu = 0; mu < NDIM; mu++) {
    Xprev[mu]    = X[mu];
    Kcovprev[mu] = Kcov[mu];
    Kconprev[mu] = Kcon[mu];
  }

  do {
    for (int mu = 0; mu < NDIM; mu++) {
      X[mu]    = Xprev[mu];
      Kcov[mu] = Kcovprev[mu];
      Kcon[mu] = Kconprev[mu];
    }

    for (int n = 0; n < pow(2, nsubdiv); n++) {
      if (type == TYPE_TRACER) {
#if TRACERS
        push_tracers(X, Kcov, Kcon, P, Prad, dtpush / pow(2, nsubdiv));
#else
        fprintf(stderr,
            "[push_X_K]: Bad photon type!\n"
            "\ttype = %d\n"
            "\tX    = [%g %g %g %g]\n"
            "\tKcov = [%g %g %g %g]\n"
            "\tKcon = [%g %g %g %g]\n",
            type, X[0], X[1], X[2], X[3], Kcov[0], Kcov[1], Kcov[2], Kcov[3],
            Kcon[0], Kcon[1], Kcon[2], Kcon[3]);
        exit(-1);
#endif
      } else {
        push(X, Kcov, Kcon, dtpush / pow(2, nsubdiv));
      }
    }
    if (type == TYPE_TRACER) {
      good_push = 1; // TODO: is this right?
    } else {
      good_push = is_null(Kcov, Kcon, Kcovprev[0], KdotKprev, &KdotK);
    }
    nsubdiv++;
  } while (!good_push && nsubdiv < MAX_SUBDIV);

  if (nsubdiv == MAX_SUBDIV) {
    // fprintf(stderr, "Geodesic failure!\n");
    // fprintf(stderr, "X[] = %e %e %e %e\n", X[0], X[1], X[2], X[3]);
    return PUSH_FAIL;
  }

  return PUSH_SUCCESS;
}

int push_superphoton(struct of_photon *ph, grid_prim_type P,
    grid_prim_type Prad, double dtpush) {
  // Store X and K at step n
  for (int n = 0; n < NSUP - 1; n++) {
    for (int mu = 0; mu < NDIM; mu++) {
      ph->X[n][mu]    = ph->X[n + 1][mu];
      ph->Kcov[n][mu] = ph->Kcov[n + 1][mu];
      ph->Kcon[n][mu] = ph->Kcon[n + 1][mu];
    }
  }

  int status    = push_X_K(ph->X[2], ph->Kcov[2], ph->Kcon[2], P, Prad,
      ph->KdotKprev, ph->type, dtpush);
  ph->KdotKprev = dot(ph->Kcov[2], ph->Kcon[2]);
  return status;
}

void push_superphotons(grid_prim_type P, grid_prim_type Prad, double dt) {
  timer_start(TIMER_PUSH);

  int step_lost_local  = 0;
  int step_tot_local   = 0;
  int tracer_tot_local = 0;

#pragma omp parallel reduction(+:step_lost_local) reduction(+:step_tot_local) \
    reduction(+:tracer_tot_local)
  {
    int               push_status;
    struct of_photon *ph   = photon_lists[omp_get_thread_num()];
    struct of_photon *prev = NULL;
    struct of_photon *head = ph;

    // Push each photon from n to n+1
    while (ph != NULL) {
      if (to_be_pushed(t, dt, ph)) {
        double dtpush = get_dtpush(ph, dt);
        push_status   = push_superphoton(ph, P, Prad, dtpush);

        if (push_status == PUSH_FAIL) {
// fprintf(stderr, "Geodesic push failed!\n");
// fprintf(stderr, "X[] = %e %e %e %e\n", ph->X[2][0], ph->X[2][1],
//  ph->X[2][2], ph->X[2][3]);
#if RADIATION == RADTYPE_NEUTRINOS
          record_lepton_flux(ph);
#endif
          list_remove(&ph, &head, &prev);
          step_lost_local++;
          continue;
        }
      }

      if (ph->type == TYPE_TRACER)
        tracer_tot_local++;
      step_tot_local++;

      prev = ph;
      ph   = ph->next;

    } // ph != NULL

    photon_lists[omp_get_thread_num()] = head;
  } // omp parallel

  step_lost += step_lost_local;
  step_tot += step_tot_local;
  tracer_tot += tracer_tot_local;

  timer_stop(TIMER_PUSH);
}
#undef MAX_SUBDIV

// Second order update to X^{\mu}, K_{\mu} from t to t + dt
// void push(double dt, struct of_photon *ph)
void push(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM], double dt) {
  // Heun's method:
  //   x_n+1 = x_n + dt*(0.5*c1 + 0.5*c2)
  //     c1  = dxdt(t_n, x_n)
  //     c2  = dxdt(t_n + dt, x_n + dt*c1)
  //   y_n+1 = y_n + dt*(0.5*d1 + 0.5*d2)
  //     d1  = dydt(t_n, y_n)
  //     d2  = dydt(t_n + dt, y_n + dt*d1)
  //
  //   where
  //     x = X^{\mu}, \mu = [1, 2, 3] (X[0] known)
  //     y = K_{\mu}, \mu = [1, 2]    (K_0, K_3 conserved)
  //     dydt = -1/(2*k^{\nu})*k_b*k_c*(d g^{bc} / dx^{\mu})
  //     dxdt = k^{\mu} / k^0

  double c1[NDIM], c2[NDIM], d1[NDIM], d2[NDIM];
  double Xtmp[NDIM], Kcontmp[NDIM], Kcovtmp[NDIM];
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
  memset(c1, 0, NDIM * sizeof(double));
  memset(c2, 0, NDIM * sizeof(double));
  memset(d1, 0, NDIM * sizeof(double));
  memset(d2, 0, NDIM * sizeof(double));

  // First stage
  set_gcov(X, gcov);
  gcon_func(gcov, gcon);
  raise(Kcov, gcon, Kcontmp);

  get_X_source(Kcontmp, c1);
  get_K_source(X, Kcov, Kcontmp, d1);

  for (int mu = 0; mu < NDIM; mu++) {
    Xtmp[mu]    = X[mu] + dt * c1[mu];
    Kcovtmp[mu] = Kcov[mu] + dt * d1[mu];
  }

  // Second stage
  set_gcov(Xtmp, gcov);
  gcon_func(gcov, gcon);
  raise(Kcovtmp, gcon, Kcontmp);

  get_X_source(Kcontmp, c2);
  get_K_source(Xtmp, Kcovtmp, Kcontmp, d2);

  X[0] += dt;
  for (int mu = 1; mu < NDIM; mu++) {
    X[mu] += 0.5 * dt * (c1[mu] + c2[mu]);
    Kcov[mu] += 0.5 * dt * (d1[mu] + d2[mu]);
  }

  // Also provide Kcon
  set_gcov(X, gcov);
  gcon_func(gcov, gcon);
  raise(Kcov, gcon, Kcon);
}

void get_X_source(double Kcon[NDIM], double src[NDIM]) {
  for (int mu = 0; mu < NDIM; mu++) {
    src[mu] = Kcon[mu] / Kcon[0];
  }
}

#define DELTA (1.e-6)
void get_K_source(
    double X[NDIM], double Kcov[NDIM], double Kcon[NDIM], double src[NDIM]) {
// Don't do work when sources are trivial
#if METRIC == MINKOWSKI
  static int killing[] = {1, 1, 1, 1};
#elif METRIC == MKS
  static int killing[] = {1, 0, 0, 1};
#endif

  for (int mu = 0; mu < NDIM; mu++) {
    if (killing[mu] == 1) {
      src[mu] = 0.;
    } else {
      src[mu] = 0.;

      // Numerically evaluate d g^{\nu \kap} / dx^{\mu}
      double Xm[NDIM], Xp[NDIM];
      double gcovm[NDIM][NDIM], gcovp[NDIM][NDIM];
      double gconm[NDIM][NDIM], gconp[NDIM][NDIM], dG;
      for (int nu = 0; nu < NDIM; nu++) {
        Xm[nu] = X[nu];
        Xp[nu] = X[nu];
      }
      Xm[mu] -= DELTA;
      Xp[mu] += DELTA;
      set_gcov(Xm, gcovm);
      set_gcov(Xp, gcovp);
      gcon_func(gcovm, gconm);
      gcon_func(gcovp, gconp);

      for (int nu = 0; nu < NDIM; nu++) {
        for (int kap = 0; kap < NDIM; kap++) {
          dG = (gconp[nu][kap] - gconm[nu][kap]) / (Xp[mu] - Xm[mu]);
          src[mu] += Kcov[nu] * Kcov[kap] * dG;
        }
      }

      src[mu] *= -1. / (2. * Kcon[0]);
    }
  }
}
#undef DELTA
#endif // RADIATION
