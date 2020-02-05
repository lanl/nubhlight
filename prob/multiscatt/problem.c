/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR UNIT TEST FOR MULTISCATT                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static int Nr0;
static double E_MEV;
static double rho0;
static double uu0;
static int direction;

// Problem-specific variables to set at runtime
void set_problem_params()
{
  set_param("Nr0", &Nr0);
  set_param("E_MEV", &E_MEV);
  set_param("rho0", &rho0);
  set_param("uu0", &uu0);
  set_param("direction", &direction);
}

// Initialize dynamical variables
void init_prob()
{
  int nr_per_cell = Nr0/(N1TOT*N2TOT*N3TOT);
  double nu_cgs = E_MEV*MEV/HPL;
  double lnu = log10(nu_cgs);
  double nr0 = 1e10; // cm^-3

  if (mpi_io_proc()) {
    printf("E_MEV    = %g\n"
	   "lnu      = %g\n"
	   "nu       = %g\n"
	   "Nr0      = %d\n"
	   "rho0     = %g\n"
	   "uu0      = %g\n",
	   E_MEV, lnu, nu_cgs, Nr0, rho0, uu0);
    printf("SCATTERING      = %d\n"
	   "MULTISCATT_TEST = %d\n"
	   "RADIATION       = %d\n"
	   "RAD_SCATT_TYPES = %d\n",
	   SCATTERING,MULTISCATT_TEST,
	   RADIATION,RAD_SCATT_TYPES);
  }

  ZLOOP {
    P[i][j][k][RHO] = rho0;
    P[i][j][k][UU]  = uu0;
    P[i][j][k][U1]  = 0.0; // 1.e-6;
    P[i][j][k][U2]  = 0.0; // 1.e-6;
    P[i][j][k][U3]  = 0.0; // -1.e-3;
    P[i][j][k][B1]  = 0.0; // 1e-2;
    P[i][j][k][B2]  = 0.0; // 1e-2;
    P[i][j][k][B3]  = 0.0; // 1e-2;
    #if EOS == EOS_TYPE_TABLE
    P[i][j][k][YE]    = 0.5;
    P[i][j][k][YE_EM] = 0.5;
    #endif
  }

  #pragma omp parallel
  {
    double X[NDIM], K_tetrad[NDIM];
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    ZLOOP {
      coord(i, j, k, CENT, X);
      
      for (int n = 0; n < nr_per_cell/nthreads; n++) {
        struct of_photon *phadd = safe_malloc(sizeof(struct of_photon));
	
	// photon position centered on cell
        phadd->X[2][0] = 0.;
        for (int mu = 1; mu < NDIM; mu++) {
          phadd->X[2][mu] = X[mu];
        }
	
	double nu = pow(10.,lnu);
	double E = HPL*nu/(ME*CL*CL);
	K_tetrad[0] = -E;
	K_tetrad[1] = K_tetrad[2] = K_tetrad[3] = 0.;
	K_tetrad[direction] = E;

	DLOOP1 { // Can just use non-moving fluid frame
	  phadd->Kcov[2][mu] = K_tetrad[mu];
	  phadd->Kcon[2][mu] = K_tetrad[mu];
	}
	phadd->Kcon[2][0] *= -1.;
	
        phadd->t0 = 0.;
        phadd->origin[0] = nstep;
        phadd->origin[1] = i;
        phadd->origin[2] = j;
        phadd->origin[3] = k;
        phadd->w = nr0*dx[1]*dx[2]*dx[3]*L_unit*L_unit*L_unit/nr_per_cell;
	phadd->nscatt = 0;
	phadd->type = 0;
	
	phadd->next = ph;
        ph = phadd;
      }
    }
    photon_lists[omp_get_thread_num()] = ph;
  } // omp parallel
}

