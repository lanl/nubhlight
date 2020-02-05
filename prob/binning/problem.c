/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR UNIT TEST FOR SUPERPHOTON BINNING                   *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define ISOTROPIC (0)
#define X1DIR (1)
#define X2DIR (2)
#define X3DIR (3)
#define LINEAR_ANISOTROPY (4)
static int direction;
static int Nr0;

// Problem-specific variables to set at runtime
void set_problem_params()
{
  set_param("direction", &direction);
  set_param("Nr0", &Nr0);
}

// Initialize dynamical variables
void init_prob()
{
  int nr_per_cell = Nr0/(N1TOT*N2TOT*N3TOT);
  double lnumax = log(numax);
  double lnumin = log(numin);
  double lnumu = 0.5*(lnumax + lnumin);
  double lnusigma = (lnumax - lnumin)/(2.*4);
  double nr0 = 1e10; // cm^-3

  if (mpi_io_proc()) {
    printf("lnumin   = %g\n"
	   "lnumax   = %g\n"
	   "lnumu    = %g\n"
	   "lnusigma = %g\n",
	   lnumin,lnumax,
	   lnumu,lnusigma);
  }

  ZLOOP {
    P[i][j][k][RHO] = 1.e0;
    P[i][j][k][UU]  = 1.e0;
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
	
	// Energetics: Gaussian distributed in lognu
	// but throw out photons outside range.
	double lnu = lnumin - 1.0;
	while (lnu < lnumin || lnu > lnumax) {
	  lnu = get_gaussian(lnumu, lnusigma);
	}
	double nu = exp(lnu);
	double E = HPL*nu/(ME*CL*CL);
	K_tetrad[0] = -E;
	if (direction == LINEAR_ANISOTROPY) {
	  double y,cth,pdf;
	  do {
	    y = 2.*get_rand();
	    cth = 2.*get_rand() - 1;
	    pdf = 1 + cth;
	  } while (y >= pdf);
	  double sth = sqrt(fabs(1 - cth*cth));
	  double phi = 2*M_PI*get_rand();
	  double cphi = cos(phi);
	  double sphi = sin(phi);
	  K_tetrad[1] = E*cphi*sth;
	  K_tetrad[2] = E*sphi*sth;
	  K_tetrad[3] = E*cth;
	}
	else if (direction == ISOTROPIC) {
	  // uniform in solid angle
	  double cth = 2.*get_rand() - 1.;
	  double th = acos(cth);
	  double sth = sin(th);
	  double phi = 2.*M_PI*get_rand();
	  double cphi = cos(phi);
	  double sphi = sin(phi);
	  K_tetrad[1] = E*cphi*sth;
	  K_tetrad[2] = E*sphi*sth;
	  K_tetrad[3] = E*cth;
	} else if (X1DIR <= direction && direction <= X3DIR) {
	  K_tetrad[1] = K_tetrad[2] = K_tetrad[3] = 0;
	  K_tetrad[direction] = E;
	} else { // SOMETHING WRONG
	  fprintf(stderr, "Unsupported photon direction.\n");
	  exit(-1);
	}

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
	
	phadd->next = ph;
        ph = phadd;
      }
    }
    photon_lists[omp_get_thread_num()] = ph;
  } // omp parallel
}

