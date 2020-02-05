/******************************************************************************
 *                                                                            *
 * ELECTRONS.C                                                                *
 *                                                                            *
 * ELECTRON THERMODYNAMICS                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <gsl/gsl_sf_bessel.h>

/* TODO: Encapsulate electron equation of state in EOS framework.
 *       For now, I'v ebroken encapsulation but made the code
 *	 complain/break if EOS_TYPE is not EOS_TYPE_GAMMA
 *	 when ELECTRONS are active.
 *
 *       The best strategy is probably to make a two-temperature
 *       electron EOS one of the available EOS's.
 *
 *       Until then, I'm worried that we might hit a maintainability
 *       problem becasue if the implementation in ELECTRONS changes
 *       the implementation in EOS_TYPE_GAMMA needs to change too.
 *       ~JMM
 */

#if ELECTRONS
#if EOS == EOS_TYPE_GAMMA
void heat_electrons_zone(int i, int j, int k, double Pi[NVAR], double Ps[NVAR],
    double Pf[NVAR], double Dt);
void fixup_electrons_1zone(double P[NVAR]);

void init_electrons() {
  double uel;

  ZSLOOP(-NG, N1 + NG - 1, -NG, NG + N2 - 1, -NG, NG + N3 - 1) {
    // Set electron internal energy to constant fraction of internal energy
    uel = fel0 * P[i][j][k][UU];

    // Initialize entropies
    P[i][j][k][KTOT] = (gam - 1.) * P[i][j][k][UU] * pow(P[i][j][k][RHO], -gam);
    P[i][j][k][KEL]  = (game - 1.) * uel * pow(P[i][j][k][RHO], -game);
  }

  bound_prim(P);
}

void heat_electrons(
    grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt) {
  timer_start(TIMER_ELECTRON);
#pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    heat_electrons_zone(i, j, k, Pi[i][j][k], Ps[i][j][k], Pf[i][j][k], Dt);
  }
  timer_stop(TIMER_ELECTRON);
}

void heat_electrons_zone(int i, int j, int k, double Pi[NVAR], double Ps[NVAR],
    double Pf[NVAR], double Dt) {
  double ktotharm, ktotadv, fel;

  // Calculated and advected entropy at final time
  ktotharm = (gam - 1.) * Pf[UU] / pow(Pf[RHO], gam);
  ktotadv  = Pf[KTOT];

  // Electron heating fraction
  fel = get_fel(i, j, k, Ps);

  // Update final electron entropy according to Ressler+ 2015 Eqn. 27:
  Pf[KEL] += (game - 1.) / (gam - 1.) * pow(Ps[RHO], gam - game) * fel *
             (ktotharm - ktotadv);

  // Diagnostics
  struct of_geom *geom = &ggeom[i][j][CENT];
  struct of_state q;
  get_state(Ps, geom, &q);
  double uadv = ktotadv / (gam - 1.) * pow(Pf[RHO], gam);
  double Qud  = q.ucon[0] * q.ucov[0] * (Pf[UU] - uadv) *
               pow(Ps[RHO] / Pf[RHO], gam) / Dt;
  Qvisc[i][j][k] = fel * Qud;

  // Reset total entropy
  Pf[KTOT] = ktotharm;
}

double get_fel(int i, int j, int k, double P[NVAR]) {
#if BETA_HEAT == 0
  return fel0;
#endif

  struct of_geom geom;
  double         beta, fel, c1, Trat, Tpr, mbeta, qrat;
  double         pres, bsq;
  double         c2, c3, c22, c32;

  c1 = 0.92;

  Tpr        = (gam - 1.) * P[UU] / P[RHO];
  double uel = 1. / (game - 1.) * P[KEL] * pow(P[RHO], game);
  double Tel = (game - 1.) * uel / P[RHO];
  if (Tel <= 0.)
    Tel = SMALL;
  if (Tpr <= 0.)
    Tpr = SMALL;
  Trat = fabs(Tpr / Tel);

  pres = P[RHO] * Tpr; // Proton pressure

  geom = *get_geometry(i, j, k, CENT);
  bsq  = bsq_calc(P, &geom);

  beta = pres / bsq * 2.;
  if (beta > 1.e20)
    beta = 1.e20;
  mbeta = 2. - 0.2 * log10(Trat);

  if (Trat <= 1.) {
    c2 = 1.6 / Trat;
    c3 = 18. + 5. * log10(Trat);
  } else {
    c2 = 1.2 / Trat;
    c3 = 18.;
  }
  c22 = pow(c2, 2.);
  c32 = pow(c3, 2.);

  qrat = c1 * (c22 + pow(beta, mbeta)) / (c32 + pow(beta, mbeta)) *
         exp(-1. / beta) * pow(MP / ME * Trat, .5);
  fel = 1. / (1. + qrat);

  return fel;
}

// Modified Bessel function of second kind with safe inputs
double safe_Kn(int n, double x) {
  if (x > 100.) {
    return exp(-x) * sqrt(M_PI / (2. * x));
  } else {
    return gsl_sf_bessel_Kn(n, x);
  }
}

#if RADIATION
void coulomb(
    grid_prim_type Pi, grid_prim_type Ps, grid_prim_type Pf, double Dt) {
  timer_start(TIMER_ELECTRON);

#pragma omp parallel for collapse(3)
  ZLOOP {
    double          rho     = Ps[i][j][k][RHO];
    double          thetae  = MP / ME * Ps[i][j][k][KEL] * pow(rho, game - 1.);
    double          ue      = Ps[i][j][k][KEL] * pow(rho, game) / (game - 1.);
    double          up      = Ps[i][j][k][UU] - ue;
    double          n       = rho * Ne_unit;
    double          Ti      = up * U_unit * (gamp - 1.) / (n * KBOL);
    double          thetai  = KBOL * Ti / (MP * CL * CL);
    double          thetam  = 1. / (1. / thetae + 1. / thetai);
    double          logCoul = COULOMB_LOG;
    double          Te      = thetae * ME * CL * CL / KBOL;
    struct of_geom *geom;
    struct of_state q;

    // Sanity checks, although electron fixup routine should catch these
    if (!isnan(Te) && !isnan(Ti) && Te > 0. && Ti > 0.) {
      double Qc, term1, term2;

      // Get Coulomb heating rate.
      // Need to handle cases where Thetai < 1e-2, Thetae < 1e-2, and both
      // Thetae and Thetai < 1e-2 separately due to Bessel functions exploding
      double prefac =
          3. / 2. * ME / MP * n * n * logCoul * CL * KBOL * THOMSON * (Ti - Te);
      double thetaCrit = 1.e-2;
      if (thetae < thetaCrit && thetai < thetaCrit) {
        term1 = sqrt(thetam / (M_PI * thetae * thetai / 2.));
        term2 = sqrt(thetam / (M_PI * thetae * thetai / 2.));
      } else if (thetae < thetaCrit) {
        term1 =
            exp(-1. / thetai) / safe_Kn(2, 1. / thetai) * sqrt(thetam / thetae);
        term2 =
            exp(-1. / thetai) / safe_Kn(2, 1. / thetai) * sqrt(thetam / thetae);
      } else if (thetai < thetaCrit) {
        term1 =
            exp(-1. / thetae) / safe_Kn(2, 1. / thetae) * sqrt(thetam / thetai);
        term2 =
            exp(-1. / thetae) / safe_Kn(2, 1. / thetae) * sqrt(thetam / thetai);
      } else {
        term1 = safe_Kn(1, 1. / thetam) /
                (safe_Kn(2, 1. / thetae) * safe_Kn(2, 1. / thetai));
        term2 = safe_Kn(0, 1. / thetam) /
                (safe_Kn(2, 1. / thetae) * safe_Kn(2, 1. / thetai));
      }
      term1 *= (2. * pow(thetae + thetai, 2) + 1.) / (thetae + thetai);
      term2 *= 2.;
      Qc = prefac * (term1 + term2);

      // Convert to code units
      Qc *= T_unit / U_unit;

      // Update electron internal energy
      geom = &ggeom[i][j][CENT];
      get_state(Ps[i][j][k], geom, &q);

      double ue_f =
          Pf[i][j][k][KEL] * pow(Pf[i][j][k][RHO], game) / (game - 1.);
      ue_f += Qc * Dt / q.ucon[0];

      // Record diagnostic
      Qcoul[i][j][k] = q.ucov[0] * Qc;

      // Update electron entropy
      Pf[i][j][k][KEL] = (game - 1.) * ue_f * pow(Pf[i][j][k][RHO], -game);
    }
  } // ZLOOP

  timer_stop(TIMER_ELECTRON);
}

void apply_rad_force_e(grid_prim_type Prh, grid_prim_type Pr,
    grid_fourvector_type radG, double Dt) {
// Apply only to active zones for this proc -- ghost zone four-force
// depositions already communicated over MPI
#pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP {
    struct of_geom *geom = &ggeom[i][j][CENT];
    struct of_state q;
    double          Uel, Urho;
    double          C = 0.;

    // Get fluid state at n + 1/2 where radiation four-force is centered
    get_state(Prh[i][j][k], geom, &q);

    for (int mu = 0; mu < NDIM; mu++) {
      C += -q.ucon[mu] * radG[i][j][k][mu];
    }

    // Get fluid state at n+1 for full update
    get_state(Pr[i][j][k], geom, &q);

    // Remove \sqrt{-g} from radG
    C = C / geom->g;

    Urho = Pr[i][j][k][RHO] * q.ucon[0];
    Uel  = Pr[i][j][k][KEL] * Urho;

    Uel += Dt * (C * (game - 1.) * pow(Prh[i][j][k][RHO], 1. - game));

    // Supercooling diagnostics
    if (Uel < 0.) {
      double          U_1[NVAR], prim_2[NVAR], U_2[NVAR];
      struct of_state q_1, q_2;
      Nsuper[i][j][k]++;

      // (1) Record total energy density after cooling
      get_state(Pr[i][j][k], &ggeom[i][j][CENT], &q_1);
      primtoflux(Pr[i][j][k], &q_1, 0, 0, &ggeom[i][j][CENT], U_1);

      // (2) Calculate total energy density with zero electron energy density
      PLOOP  prim_2[ip] = psupersave[i][j][k][ip];
      double ue         = prim_2[KEL] * pow(prim_2[RHO], game) / (game - 1.);
      prim_2[UU] -= ue;
      get_state(prim_2, &ggeom[i][j][CENT], &q_2);
      primtoflux(prim_2, &q_2, 0, 0, &ggeom[i][j][CENT], U_2);

      // Subtract (2) from (1); integrated over volume, this is fake energy
      Esuper[i][j][k] += fabs((U_1[UU] - U_2[UU]) * dx[1] * dx[2] * dx[3]);
    } // Uel < 0

    Pr[i][j][k][KEL] = Uel / Urho;

    // Reset total entropy
    Pr[i][j][k][KTOT] =
        (gam - 1.) * Pr[i][j][k][UU] * pow(Pr[i][j][k][RHO], -gam);
  } // ZSLOOP
}
#endif // RADIATION

void fixup_electrons(grid_prim_type P) {
  timer_start(TIMER_ELECTRON);
#pragma omp parallel for collapse(3) schedule(dynamic)
  ZLOOP { fixup_electrons_1zone(P[i][j][k]); }
  timer_stop(TIMER_ELECTRON);
}

void fixup_electrons_1zone(double P[NVAR]) {
  double kelmax =
      P[KTOT] * pow(P[RHO], gam - game) / (tptemin + (gam - 1.) / (game - 1.));
  double kelmin =
      P[KTOT] * pow(P[RHO], gam - game) / (tptemax + (gam - 1.) / (game - 1.));

  // Replace NANs with cold electrons
  if (isnan(P[KEL]))
    P[KEL] = kelmin;

  // Enforce maximum Tp/Te
  P[KEL] = MY_MAX(P[KEL], kelmin);

  // Enforce minimum Tp/Te
  P[KEL] = MY_MIN(P[KEL], kelmax);
}
#endif // EOS == EOS_TYPE_GAMMA
#endif // ELECTRONS
