/******************************************************************************
 *                                                                            *
 * TRACERS.C                                                                  *
 *                                                                            *
 * TRACER PARTICLES                                                           *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION && TRACERS

// Lagrange interpolation
#define LA_ORDER (2)
#define LA_STRIDE (LA_ORDER / 2)
#define LA_NP (2 * LA_STRIDE + 1)

void get_nearby_coords(
    int i0, int j0, int k0, int stride, double Xlocal[NDIM][LA_NP]);
double lagrange_basis(int j, double x, const double x_values[], int num_x);
void   tracer_rhs(double X[NDIM], double rhs[NDIM], grid_prim_type P,
      struct of_geom *g, struct of_state *q);

/*
 * Gets X1, X2, X3 values
 * for stride grid points in every direction
 * fill X1, X2, X3 in 1D arrays
 */
void get_nearby_coords(
    int i0, int j0, int k0, int stride, double Xlocal[NDIM][LA_NP]) {
  double np = 2 * stride + 1;
  double X[NDIM];

  for (int i = 0; i < np; i++) {
    int ir = i0 - stride + i;
    coord(ir, j0, k0, CENT, X);
    Xlocal[1][i] = X[1];
  }
  for (int j = 0; j < np; j++) {
    int jr = j0 - stride + j;
    coord(i0, jr, k0, CENT, X);
    Xlocal[2][j] = X[2];
  }
  for (int k = 0; k < np; k++) {
    int kr = k0 - stride + k;
    coord(i0, j0, kr, CENT, X);
    Xlocal[3][k] = X[3];
  }
}

/*
 * Evaluates the j'th Lagrange basis function l_j at point x. The
 * basis function is constructed from the num_x interpolation points
 * and is thus of order k = num_x-1.
 */
double lagrange_basis(int j, double x, const double x_values[], int num_x) {
  double l_j = 1;
  double r;
  if (num_x == 0 || num_x == 1)
    return 1.;
  for (int m = 0; m < num_x; m++) {
    if (m != j) {
      r = (x - x_values[m]) / (x_values[j] - x_values[m]);
      l_j *= r;
    }
  }
  return l_j;
}

/*
 * Interpolates the primitive vector to a point X in a cell using
 * 2nd-order Lagrange interpolation.
 * TODO: Use conservative or monotone interpolation like ENO/WENO
 *       OR reconstruct cell faces and then use trilinear interp
 *
 * NOTE: Assumes v is a contiguous array of shape
 *       (N1+2*NG)x(N2+2*NG)x(N3+2*NG)x(size)
 */
void lagrange_interp_grid_3d(
    double v_interp[], double X[NDIM], double *v, int size) {
  // int si = N1+2*NG;
  int sj = N2 + 2 * NG;
  int sk = N3 + 2 * NG;
  int i0, j0, k0;
  Xtoijk(X, &i0, &j0, &k0);

  double Xlocal[NDIM][LA_NP];
  get_nearby_coords(i0, j0, k0, LA_STRIDE, Xlocal);

  for (int ip = 0; ip < size; ip++)
    v_interp[ip] = 0;
  for (int i = 0; i < LA_NP; i++) {
    int    ir = i0 - LA_STRIDE + i;
    double lx = lagrange_basis(i, X[1], Xlocal[1], LA_NP);
    for (int j = 0; j < LA_NP; j++) {
      int    jr = j0 - LA_STRIDE + j;
      double ly = lagrange_basis(j, X[2], Xlocal[2], LA_NP);
      for (int k = 0; k < LA_NP; k++) {
        int    kr = k0 - LA_STRIDE + k;
        double lz = lagrange_basis(k, X[3], Xlocal[3], LA_NP);
        for (int ip = 0; ip < size; ip++) {
          // Relies on contiguous, row-major ordering!
          int    idx = ip + size * (kr + sk * (jr + sj * ir));
          double f   = v[idx];
          v_interp[ip] += f * lx * ly * lz;
        }
      }
    }
  }
}

void lagrange_interp_prim_3d(
    double P_interp[NVAR], double X[NDIM], grid_prim_type P) {
  lagrange_interp_grid_3d(P_interp, X, &(P[0][0][0][0]), NVAR);
}
#undef LA_ORDER
#undef LA_STRIDE
#undef LA_NP

/*
 * Second-order update to X^mu and K_mu from t to t + dt
 * for tracer particles.
 */
void push_tracers(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM],
    grid_prim_type P, grid_prim_type Ph, double dt) {
  struct of_geom  g;
  struct of_state q;
  double          k[NDIM], rhs[NDIM];

  // RK2
  // x_n+1 = x_n + dt*rhs(k)
  // k     = x_n + 0.5*dt rhs(x_n)

  // predictor
  tracer_rhs(X, rhs, P, &g, &q);
  SDLOOP k[mu] = X[mu] + 0.5 * dt * rhs[mu];

  // corrector
  tracer_rhs(k, rhs, Ph, &g, &q);
  SDLOOP X[mu] = X[mu] + dt * rhs[mu];

  // fill Kcov, Kcon with ucov, ucon at new position
  // used only for diagnostics
  tracer_rhs(X, rhs, Ph, &g, &q);
  DLOOP1 Kcov[mu] = q.ucov[mu];
  DLOOP1 Kcon[mu] = q.ucon[mu];

  // update time
  X[0] += dt;
}

/*
 * Tracer Right-hand-side
 * Tracers obey the equation:
 *
 * (d x^i/dt) = alpha v^i - beta^i = vtilde^i = u^i/u^0
 * where:
 * v^i = u^i/(alpha u^0) + beta^i/alpha
 * alpha = lapse
 * beta = shift
 * u^mu = contravariant 4-velocity
 *
 * NOTE: We choose to calculate ucon and ucov FIRST and use these
 *       to calculate tracer velocity. We do this to make tracers
 *       undependent of changes in primitive vs. conserved variables.
 *
 * NOTE: We calculate and save g and q.
 */
void tracer_rhs(double X[NDIM], double rhs[NDIM], grid_prim_type P,
    struct of_geom *g, struct of_state *q) {
  // Interpolate fluid to X
  double P_interp[NVAR];
  lagrange_interp_prim_3d(P_interp, X, P);

  // Extract metric/covariant quantities
  set_metric(X, g);
  get_state(P_interp, g, q);

  // 3+1 quantities
  SDLOOP rhs[mu] = q->ucon[mu] / q->ucon[0];
}

int count_tracers_local() {
  int count = 0;

#pragma omp parallel reduction(+ : count)
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->type == TYPE_TRACER)
        count++;
      ph = ph->next;
    }
  }
  return count;
}

int tracer_get_id(struct of_photon *ph) {
  if (ph->type == TYPE_TRACER) {
    return ph->nscatt;
  } else {
    fprintf(stderr, "[tracer_get_id]: Not a tracer! type = %d\n", ph->type);
    exit(1);
  }
}

double tracer_get_mass(struct of_photon *ph) {
  if (ph->type == TYPE_TRACER) {
    return ph->w;
  } else {
    fprintf(stderr, "[tracer_get_mass]: Not a tracer!\n");
    exit(1);
  }
}

void tracer_get_X_u(struct of_photon *ph, double t_interp, double X[NDIM],
    double ucov[NDIM], double ucon[NDIM]) {
  if (ph->type != TYPE_TRACER) {
    fprintf(stderr, "[tracer_get_X_u]: Not a tracer!\n");
    exit(1);
  }
  SDLOOP X[mu]    = ph->X[2][mu];
  SDLOOP ucov[mu] = ph->Kcov[2][mu];
  SDLOOP ucon[mu] = ph->Kcon[2][mu];
}

/*
 * Adds a tracer particle to photon list head
 * A tracer particle is positioned (initially) at X
 * at time t and nsteps nsteps.
 * It symbolizes a Lagrangian fluid packet with mass.
 * It has a unique id used to track it through timesteps.
 * Tracers are assumed to be created once and never destroyed.
 * Prims P are passed in to fill tracer velocity data.
 *
 * In designing tracers, we re-use several quantities in the photon
 * struct:
 * - nscatt is used to save the id
 * - w is used to save the mass
 * - Kcov and Kcon are used to save fluid 4-velocities
 */
void set_tracer(struct of_photon *ph, int id, int nstep, double t,
    double X[NDIM], double mass, grid_prim_type P) {
  double P_interp[NVAR];
  lagrange_interp_prim_3d(P_interp, X, P);

  struct of_geom g;
  set_metric(X, &g);

  struct of_state q;
  get_state(P_interp, &g, &q);

  int i, j, k;
  Xtoijk(X, &i, &j, &k);

  DLOOP1 {
    ph->X[2][mu]    = X[mu];
    ph->Kcov[2][mu] = q.ucov[mu];
    ph->Kcon[2][mu] = q.ucon[mu];
  }
  ph->X[2][0] = t;

  ph->t0        = t;
  ph->origin[0] = nstep;
  ph->origin[1] = i;
  ph->origin[2] = j;
  ph->origin[3] = k;

  ph->w      = mass;
  ph->nscatt = id;
  ph->type   = TYPE_TRACER;
}

void make_tracer(struct of_photon **head, int id, int nstep, double t,
    double X[NDIM], double mass, grid_prim_type P) {
  struct of_photon *phadd = safe_malloc(sizeof(struct of_photon));
  set_tracer(phadd, id, nstep, t, X, mass, P);
  phadd->next = *head;
  *head       = phadd;
}

/*
 * Sampling scheme for tracers in a cell
 * returns current last_id
 *
 * Current sampling scheme: randomly uniformly
 * distribute within the logical cell
 *
 * Saves sampled tracers to photon_list head
 *
 * updates photon list
 * updates last_id
 *
 * TODO: evenly distribute within the logical cell
 */
void sample_tracers_in_cell(struct of_photon **head, int nstep, int i, int j,
    int k, int *last_id, int ntracers_per_cell, double t, grid_prim_type P) {
  int    id;
  double r[NDIM];
  double Xcell[NDIM], Xpart[NDIM];
  double mass_cell, mass_part;
  mass_cell = P[i][j][k][RHO] * ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3];
  mass_part = mass_cell / ntracers_per_cell;

  for (int tcr = 0; tcr < ntracers_per_cell; tcr++) {
    id = *last_id + tcr;
    coord(i, j, k, CORN, Xcell);
    SDLOOP r[mu]     = get_rand();
    SDLOOP Xpart[mu] = Xcell[mu] + r[mu] * dx[mu];
    make_tracer(head, id, nstep, t, Xpart, mass_part, P);
  }
  *last_id += ntracers_per_cell;
  // return head;
}

/*
 * Serialized tracer function
 * produces approximately ntracers_tot
 */
grid_int_type    ntracers_p_cell;
grid_double_type tracer_mass;
void sample_all_tracers(long unsigned int ntracers_tot, int nstep, double t,
    grid_int_type tcrs_in_cell, grid_prim_type P) {
  int tcr_total = tracer_max_id(ntracers_tot);
  if (tcr_total < N1TOT * N2TOT * N3TOT) {
    fprintf(stderr,
        "Less than one tracer per cell! Please add more tracers.\n"
        "\tntracers = %d\n"
        "\tNTOT     = %d\n",
        tcr_total, N1TOT * N2TOT * N3TOT);
    exit(1);
  }

  int ntracers_p_rank = tcr_total / mpi_nprocs();
  // buffer id space in case of id collisions
  int id_global = (int)(1.5 * ntracers_p_rank * mpi_myrank());
  int last_id   = id_global;

  // Number of tracers per cell
  int ncells = 0;
  ZLOOP if (tcrs_in_cell[i][j][k]) ncells++;
  ncells          = mpi_reduce_int(ncells);
  int tcrs_p_cell = MY_MAX(tcr_total / ncells, 1);
  while (tcrs_p_cell * ncells < tcr_total)
    tcrs_p_cell++;

  // initalize tracers
  struct of_photon *ph = NULL;
  ZLOOP {
    if (tcrs_in_cell[i][j][k]) {
      sample_tracers_in_cell(&ph, nstep, i, j, k, &last_id, tcrs_p_cell, t, P);
    }
  }

  // distribute generated tracers equitably between threads
  int thread_start = (int)(get_rand() * nthreads);
  while (ph != NULL) {
    for (int n = thread_start; n < thread_start + nthreads; n++) {
      if (ph == NULL)
        break; // in case nphotons % nthreads != 0
      int index = n % nthreads;
      swap_ph(&ph, &(photon_lists[index]));
    }
  }
}

double get_total_tracer_mass() {
  double mass = 0.0;
  for (int thread = 0; thread < nthreads; thread++) {
    struct of_photon *ph = photon_lists[thread];
    while (ph != NULL) {
      if (ph->type == TYPE_TRACER) {
        mass += ph->w;
      }
      ph = ph->next;
    }
  }
  mass = mpi_reduce(mass);
  return mass;
}

int tracer_max_id(long unsigned int ntracers) {
  if (ntracers > INT_MAX) {
    if (mpi_io_proc()) {
      fprintf(stderr,
          "[tracer_max_id]: last_id > INT_MAX!\n"
          "\tntracers = %lu\n"
          "\tINT_MAX  = %d\n"
          "\tPlease reduce number of tracers.\n",
          ntracers, INT_MAX);
      exit(1);
    }
  }
  return (int)ntracers;
}

/*
  Deletes tracers in cells that have been fixed up
 */
void        prune_tracers() {
#pragma omp parallel
  {
    double            X[NDIM], Kcov[NDIM], Kcon[NDIM];
    int               i, j, k;
    struct of_photon *ph   = photon_lists[omp_get_thread_num()];
    struct of_photon *prev = NULL;
    struct of_photon *head = ph;
    while (ph != NULL) {
      if (ph->type == TYPE_TRACER) {
        tracer_get_X_u(ph, t, X, Kcov, Kcon);
        Xtoijk(X, &i, &j, &k);
        if (fixup_required[i][j][k]) {
          list_remove(&ph, &head, &prev);
          continue;
        }
      }
      prev = ph;
      ph   = ph->next;
    } // ph != NULL
    photon_lists[omp_get_thread_num()] = head;
  }
}

#endif // RADIATION && TRACERS
