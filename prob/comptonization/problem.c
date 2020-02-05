/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR COMPTONIZATION                                      *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Problem-specific variables to set at runtime
void set_problem_params()
{
}

// Initialize dynamical variables
void init_prob()
{
  double Te0 = 5.e7; // K
  double ne0 = 2.5e17; // cm^-3
  double nr0 = 2.536851e18; // cm^-3
  double ur0 = 4.731013e8; // erg cm^-3
  double nu0 = ur0/(HPL*nr0);
  double Nr0 = 30000;
  double u0 = ne0*KBOL*Te0/(gam-1.);
  double Tf = (HPL*nr0*nu0 + 2.*ne0*KBOL*Te0/(gam-1.))/(3.*nr0*KBOL + 2.*ne0*KBOL/(gam-1.));
  printf("nu0 = %e cm^-3\n", nu0);
  printf("Tf = %e K\n", Tf);

  ZLOOP {
    P[i][j][k][RHO] = 1.e0;
    P[i][j][k][UU]  = 2.*u0/U_unit;
    P[i][j][k][U1]  = 0.;
    P[i][j][k][U2]  = 0.;
    P[i][j][k][U3]  = 0.;
    P[i][j][k][B1]  = 0.;
    P[i][j][k][B2]  = 0.;
    P[i][j][k][B3]  = 0.;
  }

  #pragma omp parallel
  {
    double X[NDIM], K_tetrad[NDIM];
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    struct of_microphysics micro;
    ZLOOP {
      coord(i, j, k, CENT, X);
      get_fluid_zone(i, j, k, P, extra, &micro, Ucon, Ucov, Bcon, 
        Bcov);
      make_tetrad(i, j, k, Ucon, Bcon, ggeom[i][j][CENT].gcov, Econ, Ecov);

      for (int n = 0; n < Nr0/(N1*N2*N3*nthreads); n++) {
        struct of_photon *phadd = safe_malloc(sizeof(struct of_photon));

        phadd->X[2][0] = 0.;
        for (int mu = 1; mu < NDIM; mu++) {
          phadd->X[2][mu] = X[mu];
        }

        // Monoenergetic, uniform in solid angle
        double cth = 2.*get_rand() - 1.;
        double th = acos(cth);
        double sth = sin(th);
        double phi = 2.*M_PI*get_rand();
        double cphi = cos(phi);
        double sphi = sin(phi);

        double E = HPL*nu0/(ME*CL*CL);

        K_tetrad[0] = -E;
        K_tetrad[1] = E*cth;
        K_tetrad[2] = E*cphi*sth;
        K_tetrad[3] = E*sphi*sth;

        tetrad_to_coord(Ecov, K_tetrad, phadd->Kcov[2]);
        K_tetrad[0] *= -1.;
        tetrad_to_coord(Econ, K_tetrad, phadd->Kcon[2]);

        phadd->t0 = 0.;
        phadd->origin[0] = nstep;
        phadd->origin[1] = i;
        phadd->origin[2] = j;
        phadd->origin[3] = k;
        phadd->w = nr0*dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit/(Nr0/(N1*N2*N3));
	phadd->nscatt = 0;
	phadd->type = 0;

        phadd->next = ph;
        ph = phadd;
      }
    }
    photon_lists[omp_get_thread_num()] = ph;
  } // omp parallel
}

