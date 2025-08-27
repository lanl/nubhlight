/******************************************************************************
 *                                                                            *
 * INTERACT.C                                                                 *
 *                                                                            *
 * PROCESS ABSORPTION AND SCATTERING                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#if RADIATION

#if SCATTERING
#ifndef REBALANCE
#define REBALANCE (1)
#endif // rebalance
#endif // scattering

#if RADIATION == RAD_TYPE_LIGHT
#define CONSERVATIVE_BOUND (1)
#else
#define CONSERVATIVE_BOUND (0)
#endif

#if SCATTERING
double conservative_dtau_est(
    double nu, int type, int interaction, const struct of_microphysics *m) {

  double dl        = Rout_rad * L_unit;
  double alpha_inv = alpha_inv_scatt(nu, type, interaction, m);
  double dtau      = (alpha_inv / nu) * dl;

  if (fabs(alpha_inv) <= SMALL * SMALL || alpha_inv < 0)
    dtau = 1.0;

  return dtau;
}

double bound_bias(double bias, double nu, int type, int interaction,
    const struct of_microphysics *m, double uph) {
  // Another BOUND_BIAS method. Assumes large hotspots, may work fine instead
  // assuming ~GM/c^2 length scale for hot spots.
  double dtau = conservative_dtau_est(nu, type, interaction, m);
  bias        = MY_MIN(bias, 1. / (dtau * RAD_NUM_TYPES * RAD_SCATT_TYPES));
  bias        = MY_MAX(bias, 1.);

  return bias;
}

double get_scatt_bias(double nu, int type, int interaction,
    const struct of_microphysics *m, double uph) {
  double Thetae = scatterer_dimensionless_temp(type, interaction, m);
  double amp =
      1. + 4. * Thetae - 2. * pow(Thetae, 3. / 2.) + 16. * pow(Thetae, 2.);
  double bias = tune_scatt * amp;
#if CONSERVATIVE_BOUND
  bias = bound_bias(bias, nu, type, interaction, m, uph);
#endif
  return bias;
}

void rebalance_biases(
    double dtau_scatt[RAD_SCATT_TYPES], double bias_scatt[RAD_SCATT_TYPES]) {
  if (RAD_SCATT_TYPES == 1)
    return;
  double max_prod = -INFINITY;
  double prod;
  SCATTLOOP {
    prod = bias_scatt[iscatt] * dtau_scatt[iscatt];
    if (prod > max_prod)
      max_prod = prod;
  }
  SCATTLOOP {
    // dtau can underflow. If it does, set the bias to 1.
    // This interaction is very subdominant and should be ignored.
    if (dtau_scatt[iscatt] < SMALL) {
      bias_scatt[iscatt] = 1;
    } else {
      bias_scatt[iscatt] = max_prod / (dtau_scatt[iscatt] + SMALL);
    }
  }
}

/*
 * Ensure biases obey relation
 *
 * 0 < sum_{interactions i} ( b_i dtau_i ) = SCATT_BIAS_SAFETY < 1
 *
 * and that each bias >= 1.
 */
void bound_all_biases(double dtau_scatt[RAD_SCATT_TYPES],
    double bias_scatt[RAD_SCATT_TYPES], double nu, double type,
    const struct of_microphysics *m) {

  double sum = 0;
  SCATTLOOP {
    // ensure dtau and bias are physical
    if (bias_scatt[iscatt] < 1.)
      bias_scatt[iscatt] = 1.;
    double dtau =
        dtau_scatt[iscatt]; // conservative_dtau_est(nu,type,iscatt,m);
    double prod = dtau * bias_scatt[iscatt];
    sum += prod;
  }
  double ratio = SCATT_BIAS_SAFETY / sum;
  // printf("Sum = %g\nRatio = %g\n",sum,ratio);
  if (ratio < 1) { // Shrink all biases so sum is safe
    SCATTLOOP bias_scatt[iscatt] *= ratio;
  }
  // Make sure all biases >= 1.
  SCATTLOOP {
    if (isnan(bias_scatt[iscatt]) || bias_scatt[iscatt] < 1 ||
        isinf(bias_scatt[iscatt])) {
      bias_scatt[iscatt] = 1.0;
    }
  }
}
#endif // SCATTERING

#define MAX_INTERACTIONS 100
void interact(grid_prim_type P, grid_eosvar_type extra, double t, double dt) {
  timer_start(TIMER_INTERACT);
#if KILL_ALL_PACKETS
  {
#pragma omp parallel
    {
      struct of_photon *ph = photon_lists[omp_get_thread_num()];
      while (ph != NULL) {
        if (ph->type != TYPE_TRACER)
          ph->w = 0;
        ph = ph->next;
      }
    }
    return;
  }
#endif

#if ABSORPTION || SCATTERING
  // printf("Entering interact.\n"); // DEBUG
  const double d3x = dx[1] * dx[2] * dx[3];

#pragma omp parallel
  {
    int               i, j, k;
    struct of_photon *ph   = photon_lists[omp_get_thread_num()];
    struct of_photon *prev = NULL;
    struct of_photon *head = ph;
    // struct of_microphysics m;
    double X[NDIM], Kcon[NDIM], Kcov[NDIM];
    // double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double dlam, nu;
    double dtau_abs = 0.;
    double dtau_scatt[RAD_SCATT_TYPES];
    double bias_scatt[RAD_SCATT_TYPES];
    double xabs;
    double xscatt[RAD_SCATT_TYPES];
    double dlfrac;
    int    abs_possible, scatt_possible, scatt_to_do;
    int    nint;

    SCATTLOOP {
      dtau_scatt[iscatt] = 0.;
      bias_scatt[iscatt] = 0.;
    }

    while (ph != NULL) {
      if (ph->w < SMALL || ph->type == TYPE_TRACER) {
        prev = ph;
        ph   = ph->next;
        continue;
      }
      // double tmin = t;
      double tmin = MY_MAX(t, ph->t0);
      // Loop with decreasing d \lambda to account for multiple interactions
      // per superphoton per timestep (scattering only)
      nint = 0;
      while (nint < MAX_INTERACTIONS) {
        nint++;

        double dtint  = t + dt - tmin;
        int    status = get_X_K_interp(ph, t + dt / 2., P, X, Kcov, Kcon);
        if (status == SPH_INTERP_FAIL) {
          prev = ph;
          ph   = ph->next;
          break;
        }

        dlam = dtint / Kcon[0];

        Xtoijk(X, &i, &j, &k);

        if (i < 0 || i > N1 + 2 * NG - 1 || j < 0 || j > N2 + 2 * NG - 1 ||
            k < 0 || k > N3 + 2 * NG - 1) {
          fprintf(stderr, "PHOTON IN BAD PLACE! %i %i %i\n", i, j, k);
          print_ph_diag(ph);
#if RADIATION == RADTYPE_NEUTRINOS
          record_lepton_flux(ph);
#endif
          list_remove(&ph, &head, &prev);
#pragma omp atomic
          step_tot--;
          break;
        }

        // Get quantities relating radiation and fluid
// TODO: this copy is unnecessary
// get_fluid_zone(i, j, k, P, extra, &m, Ucon, Ucov, Bcon, Bcov);
/*
DLOOP1 Ucon[mu] = Ucon_grd[i][j][k][mu];
DLOOP1 Ucov[mu] = Ucon_grd[i][j][k][mu];
// DLOOP1 Bcon[mu] = Bcon_grd[i][j][k][mu];
DLOOP1 Bcov[mu] = Bcov_grd[i][j][k][mu];
memcpy((void*)(&m), (void*)(&(m_grd[i][j][k])),
       sizeof(struct of_microphysics));
*/
/*
      #if ESTIMATE_THETAE
      m.Thetae = get_Thetae_est(i, j, k);
      #endif
*/

// skip if in atmosphere
#if EOS == EOS_TYPE_TABLE && METRIC == MKS
        if (P[i][j][k][ATM] < ATM_THRESH) {
          prev = ph;
          ph   = ph->next;
          break;
        }
#endif

        if (scatt_temp_too_small(&(m_grd[i][j][k]))) {
          prev = ph;
          ph   = ph->next;
          break;
        }
        /*double sigma = pow(m.B/B_unit,2.)/(m.Ne/Ne_unit);
        if (sigma > sigma_max) {
          prev = ph;
          ph = ph->next;
          break;
        }*/

        nu = get_fluid_nu(X, Kcov, Ucon_grd[i][j][k]);

        if (nu <= 0 || is_practically_nan(nu)) {
          double gamma;
          mhd_gamma_calc(P[i][j][k], &(ggeom[i][j][CENT]), &gamma);
          fprintf(stderr,
              "Bad NU in interact [%i %i %i]\n"
              "\tNU    = %e\n"
              "\tgamma = %e\n"
              "\tX[]     = [%e %e %e %e]\n"
              "\tKcov[]  = [%e %e %e %e]\n"
              "\tUcon[]  = [%e %e %e %e]\n",
              i, j, k, nu, gamma, X[0], X[1], X[2], X[3], Kcov[0], Kcov[1],
              Kcov[2], Kcov[3], Ucon_grd[i][j][k][0], Ucon_grd[i][j][k][1],
              Ucon_grd[i][j][k][2], Ucon_grd[i][j][k][3]);
#if RADIATION == RADTYPE_NEUTRINOS
          record_lepton_flux(ph);
#endif
          list_remove(&ph, &head, &prev);
#pragma omp atomic
          step_tot--;
          break;
        }

        // Superphoton type
        // DEBUGGING
        /*
        #if RADIATION == RADTYPE_NEUTRINOS
        if (ph->type < 0 || ph->type > RAD_NUM_TYPES) {
          fprintf(stderr,
            "[interact] PHOTON HAS BAD TYPE!\n"
            "\tw = %g\n"
            "\tKdotKprev = %g\n"
            "\ttype = %d\n"
            "\tnscatt = %d\n"
            "\t[X][0]    = [%g, %g, %g %g]\n"
            "\t[X][1]    = [%g, %g, %g %g]\n"
            "\t[X][2]    = [%g, %g, %g %g]\n"
            "\t[Kcon][0] = [%g, %g, %g %g]\n"
            "\t[Kcon][1] = [%g, %g, %g %g]\n"
            "\t[Kcon][2] = [%g, %g, %g %g]\n"
            "\t[Kcov][0] = [%g, %g, %g %g]\n"
            "\t[Kcov][1] = [%g, %g, %g %g]\n"
            "\t[Kcov][2] = [%g, %g, %g %g]\n"
            "\t[origin] = [%d, %d, %d %d]\n"
            "\tt0 = %g\n"
            "\tis_tracked = %d\n",
            location,
            ph->w,
            ph->KdotKprev,
            ph->type,
            ph->nscatt,
            ph->X[0][0],ph->X[0][1],ph->X[0][2],ph->X[0][3],
            ph->X[1][0],ph->X[1][1],ph->X[1][2],ph->X[1][3],
            ph->X[2][0],ph->X[2][1],ph->X[2][2],ph->X[2][3],
            ph->Kcon[0][0],ph->Kcon[0][1],ph->Kcon[0][2],ph->Kcon[0][3],
            ph->Kcon[1][0],ph->Kcon[1][1],ph->Kcon[1][2],ph->Kcon[1][3],
            ph->Kcon[2][0],ph->Kcon[2][1],ph->Kcon[2][2],ph->Kcon[2][3],
            ph->Kcov[0][0],ph->Kcov[0][1],ph->Kcov[0][2],ph->Kcov[0][3],
            ph->Kcov[1][0],ph->Kcov[1][1],ph->Kcov[1][2],ph->Kcov[1][3],
            ph->Kcov[2][0],ph->Kcov[2][1],ph->Kcov[2][2],ph->Kcov[2][3],
            ph->origin[0],ph->origin[1],ph->origin[2],ph->origin[3],
            ph->t0,
            ph->is_tracked
            );
          record_lepton_flux(ph);
                list_remove(&ph, &head, &prev);
                #pragma omp atomic
                step_tot--;
                // break;
          exit(1);
        }
        #endif
        */

// Calculate and sample optical depths along step
#if ABSORPTION
        double theta = get_bk_angle(
            X, Kcon, Ucov_grd[i][j][k], Bcov_grd[i][j][k], m_grd[i][j][k].B);
        dtau_abs = ((HPL * L_unit / (ME * CL * CL)) * dlam *
                    alpha_inv_abs(nu, ph->type, &(m_grd[i][j][k]), theta));
        // DEBUG
        /*
        if (dtau_abs > 1e2) {
          printf("dtau_abs = %g\n"
           "nu       = %g\n"
           "type     = %d\n"
           "alpha    = %g\n"
           "rho      = %g\n"
           "T        = %g\n"
           "Ye       = %g\n"
           "[i,j,k]  = [%d, %d, %d]\n",
           dtau_abs,nu,ph->type,
           alpha_inv_abs(nu, ph->type,
                   &(m_grd[i][j][k]),
                   theta)/nu,
           m_grd[i][j][k].rho,
           m_grd[i][j][k].T/MEV,
           m_grd[i][j][k].Ye,
           i,j,k);
        }
        */
#endif

        /* Strategy:
         * 1. Calculate all biases based on global bias and bound them
         * 2. calculate dtau*bias for all interactions
         * 3. Set every bias so that dtau*bias is equal to largest dtau*bias
         *    for all interactions
         * 4. The total probability of interaction cannot be greater than 1.
         *    Therefore we demand
         *             sum_{interactions i} (b_i dtau_i) < 1
         *    We check this after rebalancing the biases and enforce it
         *    conservatively.
         *
         * NOTE: The idea is that weakly/undersampled interactions will be
         *       enhanced while maintaining stability.
         *       There's a finite "bias" budget we can spend and we want to
         *       spend it enhancing the interactions that need it.
         * NOTE: there are probably wasted FLOPS in this implementation
         * TODO: check to see if these inefficiencies matter
         * ~JMM
         */
#if SCATTERING
        {
          double uph =
              HPL * nu * ph->w / (ggeom[i][j][CENT].g * d3x * pow(L_unit, 3));
          SCATTLOOP {
            dtau_scatt[iscatt] =
                ((HPL * L_unit / (ME * CL * CL)) * dlam *
                    alpha_inv_scatt(nu, ph->type, iscatt, &(m_grd[i][j][k])));
            bias_scatt[iscatt] =
                get_scatt_bias(nu, ph->type, iscatt, &(m_grd[i][j][k]), uph);
          }
#if REBALANCE
          {
            /*
            printf("Raw:\n" // DEBUG
             "\tprod = [ %g %g %g %g ]\n"
             "\tbias = [ %g %g %g %g ]\n"
             "\tdtau = [ %g %g %g %g ]\n",
             bias_scatt[0]*dtau_scatt[0],
             bias_scatt[1]*dtau_scatt[1],
             bias_scatt[2]*dtau_scatt[2],
             bias_scatt[3]*dtau_scatt[3],
             bias_scatt[0],bias_scatt[1],
             bias_scatt[2],bias_scatt[3],
             dtau_scatt[0],dtau_scatt[1],
             dtau_scatt[2],dtau_scatt[3]);
            */
            rebalance_biases(dtau_scatt, bias_scatt);
            /*
            printf("Rebalanced:\n" // DEBUG
             "\tprod = [ %g %g %g %g ]\n"
             "\tbias = [ %g %g %g %g ]\n"
             "\tdtau = [ %g %g %g %g ]\n",
             bias_scatt[0]*dtau_scatt[0],
             bias_scatt[1]*dtau_scatt[1],
             bias_scatt[2]*dtau_scatt[2],
             bias_scatt[3]*dtau_scatt[3],
             bias_scatt[0],bias_scatt[1],
             bias_scatt[2],bias_scatt[3],
             dtau_scatt[0],dtau_scatt[1],
             dtau_scatt[2],dtau_scatt[3]);
            */
            bound_all_biases(
                dtau_scatt, bias_scatt, nu, ph->type, &(m_grd[i][j][k]));
#if CONSERVATIVE_BOUND
            SCATTLOOP {
              bias_scatt[iscatt] = bound_bias(bias_scatt[iscatt], nu, ph->type,
                  iscatt, &(m_grd[i][j][k]), uph);
            }
#endif
            /*
            printf("Bounded:\n" // DEBUG
             "\tprod = [ %g %g %g %g ]\n"
             "\tbias = [ %g %g %g %g ]\n"
             "\tdtau = [ %g %g %g %g ]\n",
             bias_scatt[0]*dtau_scatt[0],
             bias_scatt[1]*dtau_scatt[1],
             bias_scatt[2]*dtau_scatt[2],
             bias_scatt[3]*dtau_scatt[3],
             bias_scatt[0],bias_scatt[1],
             bias_scatt[2],bias_scatt[3],
             dtau_scatt[0],dtau_scatt[1],
             dtau_scatt[2],dtau_scatt[3]);
            */
          }
#else  // NO REBALANCE
          {
            SCATTLOOP {
              bias_scatt[iscatt] = bound_bias(bias_scatt[iscatt], nu, ph->type,
                  iscatt, &(m_grd[i][j][k]), uph);
            }
          }
#endif // REBALANCE
        }
#endif // SCATTERING

        // random variables
        xabs = -log(get_rand());
        SCATTLOOP { xscatt[iscatt] = -log(get_rand()) / bias_scatt[iscatt]; }

        // are we absorbing and/or scattering?
        abs_possible   = (xabs <= dtau_abs) && ABSORPTION;
        scatt_possible = SCATTERING;
        if (scatt_possible) {
          SCATTLOOP {
            scatt_possible =
                (scatt_possible && (xscatt[iscatt] <= dtau_scatt[iscatt]));
          }
        }

        // No interaction
        if (!(abs_possible || scatt_possible)) {
          prev = ph;
          ph   = ph->next;
          break;
        }

// Absorption
#if ABSORPTION
        int do_abs = abs_possible;
        if (abs_possible) {
          double absfrac = xabs / (dtau_abs + SMALL);
          SCATTLOOP {
            double scattfrac = xscatt[iscatt] / (dtau_scatt[iscatt] + SMALL);
            do_abs           = do_abs && (absfrac < scattfrac);
          }
        }
        if (do_abs) {
          dlfrac        = xabs / dtau_abs;
          double tabs   = tmin + dlfrac * dtint;
          int    status = get_X_K_interp(ph, tabs, P, X, Kcov, Kcon);
          if (status == SPH_INTERP_FAIL) {
            prev = ph;
            ph   = ph->next;
            break;
          }

          Xtoijk(X, &i, &j, &k);

          // Boundary transport cannot use MPI with one zone
          if (N1 == 1)
            i = NG;
          if (N2 == 1)
            j = NG;
          if (N3 == 1)
            k = NG;

          for (int mu = 0; mu < NDIM; mu++) {
#pragma omp atomic
            radG[i][j][k][mu] +=
                1. / (dt * d3x) * ph->w * kphys_to_num * Kcov[mu];
          }
#if RADIATION == RADTYPE_NEUTRINOS
          {
#pragma omp atomic
            radG[i][j][k][RADG_YE] +=
                ((1. / (dt * d3x)) * Ucon_grd[i][j][k][0] * ph->w *
                    (MP / M_unit) * get_lepton_sign(ph));
          }
#endif

          Jrad[1][i][j][k] += (dt / DTd) * ph->Kcov[2][0] * kphys_to_num *
                              ph->w / (ggeom[i][j][CENT].g * dt * d3x);

#pragma omp atomic
          step_abs++;

#pragma omp atomic
          Nabs[i][j][k]++;

#pragma omp atomic
          Nabs_phys[i][j][k][ph->type] += ph->w;

          if (dtau_abs < 100) {
#pragma omp atomic
            dtau_avg[0][i][j][k] +=
                dtau_abs * (-ph->Kcov[2][0] * ME * CL * CL) * ph->w;

#pragma omp atomic
            en_int_avg[0][i][j][k] += (-ph->Kcov[2][0] * ME * CL * CL) * ph->w;
          }

          ph->w = 0.;

          tmin = tabs;

          prev = ph;
          ph   = ph->next;
          break;
        }
#endif

        else { // scattering
          // DEBUG
          /*
          printf("\tScattering superphoton\n"
           "\tBounded:\n"
           "\t\tprod = [ %g %g %g %g ]\n"
           "\t\tbias = [ %g %g %g %g ]\n"
           "\t\tdtau = [ %g %g %g %g ]\n",
           bias_scatt[0]*dtau_scatt[0],
           bias_scatt[1]*dtau_scatt[1],
           bias_scatt[2]*dtau_scatt[2],
           bias_scatt[3]*dtau_scatt[3],
           bias_scatt[0],bias_scatt[1],
           bias_scatt[2],bias_scatt[3],
           dtau_scatt[0],dtau_scatt[1],
           dtau_scatt[2],dtau_scatt[3]);
          */

          // Of all possible scattering interactions,
          // we perform the one with the smallest fraction x/dtau
          scatt_to_do = 0;
          dlfrac      = xscatt[0] / (dtau_scatt[0] + SMALL);
          SCATTLOOP {
            double scattfrac = xscatt[iscatt] / (dtau_scatt[iscatt] + SMALL);
            if (scattfrac < dlfrac) {
              scatt_to_do = iscatt;
              dlfrac      = scattfrac;
            }
          }
          // printf("scatt_to_do = %d\n",scatt_to_do); // DEBUG

          double tscatt = tmin + dlfrac * dtint;
          int    status;
          status = get_X_K_interp(ph, tscatt, P, X, Kcov, Kcon);
          if (status == SPH_INTERP_FAIL) {
            prev = ph;
            ph   = ph->next;
            break;
          }

          struct of_photon *phscatt = safe_malloc(sizeof(struct of_photon));

          if (get_rand() < Nph_to_track / (nph_per_proc * mpi_nprocs())) {
            phscatt->is_tracked = 1;
          } else {
            phscatt->is_tracked = 0;
          }

          Xtoijk(X, &i, &j, &k);

          // Initialize scattered photon at position of scattering event
          for (int mu = 0; mu < NDIM; mu++) {
            phscatt->X[2][mu]    = X[mu];
            phscatt->Kcov[2][mu] = Kcov[mu];
            phscatt->Kcon[2][mu] = Kcon[mu];
            for (int n = 0; n < 2; n++) {
              phscatt->X[n][mu]    = 0.;
              phscatt->Kcov[n][mu] = 0.;
              phscatt->Kcon[n][mu] = 0.;
            }
          }
          phscatt->t0        = tscatt;
          phscatt->w         = ph->w / bias_scatt[scatt_to_do];
          phscatt->type      = ph->type;
          phscatt->nscatt    = ph->nscatt + 1;
          phscatt->origin[0] = nstep;
          phscatt->origin[1] = i;
          phscatt->origin[2] = j;
          phscatt->origin[3] = k;
          phscatt->osc_count = ph->osc_count;

          double wsave = ph->w;
          ph->w        = (1. - 1. / bias_scatt[scatt_to_do]) * ph->w;

          int success = scatter_superphoton(
              P, extra, phscatt, X, Kcov, Kcon, scatt_to_do);

          if (!success) {
#pragma omp atomic
            step_fail++;
            free(phscatt);
            prev = ph;
            ph   = ph->next;
            break;
          }

          // Need to reset K.K
          phscatt->KdotKprev = 0.;
          for (int mu = 0; mu < NDIM; mu++) {
            phscatt->KdotKprev += phscatt->Kcov[2][mu] * phscatt->Kcon[2][mu];
          }

          // Boundary transport cannot use MPI with one zone
          if (N1 == 1)
            i = NG;
          if (N2 == 1)
            j = NG;
          if (N3 == 1)
            k = NG;

          // Apply four-force at interaction site
          for (int mu = 0; mu < NDIM; mu++) {
#pragma omp atomic
            radG[i][j][k][mu] += 1. / (dt * d3x) * phscatt->w * kphys_to_num *
                                 (Kcov[mu] - phscatt->Kcov[2][mu]);
          } // TODO: lepton number unchanged for elastic scattering?

          int nscatt = MY_MIN(ph->nscatt, MAXNSCATT - 1);

#pragma omp atomic
          Jrad[nscatt + 2][i][j][k] -=
              (dt / DTd) *
              ((phscatt->Kcov[2][0] - Kcov[0]) * kphys_to_num * phscatt->w) /
              (ggeom[i][j][CENT].g * dt * d3x);

#pragma omp atomic
          dtau_avg[scatt_to_do + 1][i][j][k] +=
              dtau_scatt[scatt_to_do] * (-ph->Kcov[2][0] * ME * CL * CL) *
              wsave;

#pragma omp atomic
          en_int_avg[scatt_to_do + 1][i][j][k] +=
              (-ph->Kcov[2][0] * ME * CL * CL) * wsave;

          // push phscatt as far as possible
          // should not use separate step and half-step for this push
          status = push_superphoton(phscatt, P, P, cour * dt_light[i][j]);
          if (status == PUSH_FAIL) {
            free(phscatt);
            prev = ph;
            ph   = ph->next;
#pragma omp atomic
            step_fail++;
            break;
          }

#pragma omp atomic
          step_scatt++;

#pragma omp atomic
          step_tot++;

#pragma omp atomic
          Nsc[i][j][k]++;

          // Add scattered superphoton to list
          phscatt->next = ph->next;
          ph->next      = phscatt;

          // Allow for additional scatterings
          tmin = tscatt;
          continue;
        }
      }
    } // ph != NULL
    photon_lists[omp_get_thread_num()] = head;

  }    // omp parallel
#endif // ABSORPTION || SCATTERING

  timer_stop(TIMER_INTERACT);
}
#endif // RADIATION
