/******************************************************************************
 *                                                                            *
 * MPI.C                                                                      *
 *                                                                            *
 * HANDLES COMMUNICATION ACROSS MPI NODES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static MPI_Comm     comm, comm_phi;
static MPI_Datatype prim_face_subtype[4], pflag_face_subtype[4];
static int          rank, numprocs, periodic[3], neighbors[3][3][3];
#if RADIATION
static MPI_Datatype mpi_photon_type, mpi_pointer_type;
static MPI_Datatype radG_face_subtype[4];
static MPI_Datatype Jrad_face_subtype[4];
static MPI_Datatype grid_radtype_face_subtype[4];
void mpi_photon_sendrecv(int is, int js, int ks, int ir, int jr, int kr);
#define NGR (1)
#endif

void init_mpi() {
  int coord[3], cpudims[3] = {N1CPU, N2CPU, N3CPU};
  numprocs = N1CPU * N2CPU * N3CPU;

#if RADIATION
#if ((X1L_GAS_BOUND == BC_PERIODIC && X1R_GAS_BOUND == BC_PERIODIC) !=     \
     (X1L_RAD_BOUND == BC_PERIODIC && X1R_RAD_BOUND == BC_PERIODIC)) ||    \
    ((X2L_GAS_BOUND == BC_PERIODIC && X2R_GAS_BOUND == BC_PERIODIC) !=     \
        (X2L_RAD_BOUND == BC_PERIODIC && X2R_RAD_BOUND == BC_PERIODIC)) || \
    ((X3L_GAS_BOUND == BC_PERIODIC && X3R_GAS_BOUND == BC_PERIODIC) !=     \
        (X3L_RAD_BOUND == BC_PERIODIC && X3R_RAD_BOUND == BC_PERIODIC))
  fprintf(stderr, "GAS and RAD must have same topology!\n");
  exit(-1);
#endif
#endif

  // Make MPI communication periodic if required
  periodic[0] = X1L_GAS_BOUND == BC_PERIODIC && X1R_GAS_BOUND == BC_PERIODIC;
  periodic[1] = X2L_GAS_BOUND == BC_PERIODIC && X2R_GAS_BOUND == BC_PERIODIC;
  periodic[2] = X3L_GAS_BOUND == BC_PERIODIC && X3R_GAS_BOUND == BC_PERIODIC;

  // Set up communicator for Cartesian processor topology
  MPI_Cart_create(MPI_COMM_WORLD, 3, cpudims, periodic, 1, &comm);
  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 3, coord);

  // Set up communicator for reductions over X^3
  MPI_Comm_split(comm, coord[2], rank, &comm_phi);

  // Find the ranks of neighbors, including edge/corner neighbors
  int n[3];
  for (int i = -1; i < 2; i++) {
    n[0] = coord[0] + i;
    for (int j = -1; j < 2; j++) {
      n[1] = coord[1] + j;
      for (int k = -1; k < 2; k++) {
        n[2] = coord[2] + k;

        if ((n[0] < 0 || n[0] >= N1CPU) && periodic[0] == 0) {
          neighbors[i + 1][j + 1][k + 1] = MPI_PROC_NULL;
        } else if ((n[1] < 0 || n[1] >= N2CPU) && periodic[1] == 0) {
          neighbors[i + 1][j + 1][k + 1] = MPI_PROC_NULL;
        } else if ((n[2] < 0 || n[2] >= N3CPU) && periodic[2] == 0) {
          neighbors[i + 1][j + 1][k + 1] = MPI_PROC_NULL;
        } else {
          MPI_Cart_rank(comm, n, &neighbors[i + 1][j + 1][k + 1]);
        }
      }
    }
  }

  // Start and stop in global index space
  int sdims[3] = {N1TOT, N2TOT, N3TOT};
  for (int d = 0; d < 3; d++) {
    global_start[d + 1] = coord[d] * sdims[d] / cpudims[d];
    global_stop[d + 1]  = (coord[d] + 1) * sdims[d] / cpudims[d];
  }

  MPI_Datatype prim_type      = MPI_DOUBLE;
  int          prim_sizes[4]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NVAR};
  int          prim_starts[4] = {0, 0, 0, 0};

  int prim_X1_subsizes[4] = {NG, N2, N3, NVAR};
  MPI_Type_create_subarray(4, prim_sizes, prim_X1_subsizes, prim_starts,
      MPI_ORDER_C, prim_type, &prim_face_subtype[1]);
  MPI_Type_commit(&prim_face_subtype[1]);

  int prim_X2_subsizes[4] = {N1 + 2 * NG, NG, N3, NVAR};
  MPI_Type_create_subarray(4, prim_sizes, prim_X2_subsizes, prim_starts,
      MPI_ORDER_C, prim_type, &prim_face_subtype[2]);
  MPI_Type_commit(&prim_face_subtype[2]);

  int prim_X3_subsizes[4] = {N1 + 2 * NG, N2 + 2 * NG, NG, NVAR};
  MPI_Type_create_subarray(4, prim_sizes, prim_X3_subsizes, prim_starts,
      MPI_ORDER_C, prim_type, &prim_face_subtype[3]);
  MPI_Type_commit(&prim_face_subtype[3]);

  MPI_Datatype pflag_type      = MPI_INT;
  int          pflag_sizes[3]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
  int          pflag_starts[3] = {0, 0, 0};

  int pflag_X1_subsizes[3] = {NG, N2, N3};
  MPI_Type_create_subarray(3, pflag_sizes, pflag_X1_subsizes, pflag_starts,
      MPI_ORDER_C, pflag_type, &pflag_face_subtype[1]);
  MPI_Type_commit(&pflag_face_subtype[1]);

  int pflag_X2_subsizes[3] = {N1 + 2 * NG, NG, N3};
  MPI_Type_create_subarray(3, pflag_sizes, pflag_X2_subsizes, pflag_starts,
      MPI_ORDER_C, pflag_type, &pflag_face_subtype[2]);
  MPI_Type_commit(&pflag_face_subtype[2]);

  int pflag_X3_subsizes[3] = {N1 + 2 * NG, N2 + 2 * NG, NG};
  MPI_Type_create_subarray(3, pflag_sizes, pflag_X3_subsizes, pflag_starts,
      MPI_ORDER_C, pflag_type, &pflag_face_subtype[3]);
  MPI_Type_commit(&pflag_face_subtype[3]);

#if RADIATION
  MPI_Datatype radG_type = MPI_DOUBLE;
  int radG_sizes[4]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NDIM + NRADCOMP};
  int radG_starts[4] = {0, 0, 0, 0};

  int radG_X1_subsizes[4] = {NGR, N2 + 2 * NGR, N3 + 2 * NGR, NDIM + NRADCOMP};
  MPI_Type_create_subarray(4, radG_sizes, radG_X1_subsizes, radG_starts,
      MPI_ORDER_C, radG_type, &radG_face_subtype[1]);
  MPI_Type_commit(&radG_face_subtype[1]);

  int radG_X2_subsizes[4] = {N1 + 2 * NGR, NGR, N3 + 2 * NGR, NDIM + NRADCOMP};
  MPI_Type_create_subarray(4, radG_sizes, radG_X2_subsizes, radG_starts,
      MPI_ORDER_C, radG_type, &radG_face_subtype[2]);
  MPI_Type_commit(&radG_face_subtype[2]);

  int radG_X3_subsizes[4] = {N1 + 2 * NGR, N2 + 2 * NGR, NGR, NDIM + NRADCOMP};
  MPI_Type_create_subarray(4, radG_sizes, radG_X3_subsizes, radG_starts,
      MPI_ORDER_C, radG_type, &radG_face_subtype[3]);
  MPI_Type_commit(&radG_face_subtype[3]);

  MPI_Datatype Jrad_type = MPI_DOUBLE;
  int Jrad_sizes[4]  = {MAXNSCATT + 2, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
  int Jrad_starts[4] = {0, 0, 0, 0};

  int Jrad_X1_subsizes[4] = {MAXNSCATT + 2, NGR, N2 + 2 * NGR, N3 + 2 * NGR};
  MPI_Type_create_subarray(4, Jrad_sizes, Jrad_X1_subsizes, Jrad_starts,
      MPI_ORDER_C, Jrad_type, &Jrad_face_subtype[1]);
  MPI_Type_commit(&Jrad_face_subtype[1]);

  int Jrad_X2_subsizes[4] = {MAXNSCATT + 2, N1 + 2 * NGR, NGR, N3 + 2 * NGR};
  MPI_Type_create_subarray(4, Jrad_sizes, Jrad_X2_subsizes, Jrad_starts,
      MPI_ORDER_C, Jrad_type, &Jrad_face_subtype[2]);
  MPI_Type_commit(&Jrad_face_subtype[2]);

  int Jrad_X3_subsizes[4] = {MAXNSCATT + 2, N1 + 2 * NGR, N2 + 2 * NGR, NGR};
  MPI_Type_create_subarray(4, Jrad_sizes, Jrad_X3_subsizes, Jrad_starts,
      MPI_ORDER_C, Jrad_type, &Jrad_face_subtype[3]);
  MPI_Type_commit(&Jrad_face_subtype[3]);

  MPI_Datatype grid_radtype_type = MPI_DOUBLE;
  int radtype_sizes[4] = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, RAD_NUM_TYPES};
  int radtype_starts[4] = {0, 0, 0, 0};

  int radtype_X1_subsizes[4] = {NGR, N2 + 2 * NGR, N3 + 2 * NGR, RAD_NUM_TYPES};
  MPI_Type_create_subarray(4, radtype_sizes, radtype_X1_subsizes,
      radtype_starts, MPI_ORDER_C, grid_radtype_type,
      &grid_radtype_face_subtype[1]);
  MPI_Type_commit(&grid_radtype_face_subtype[1]);

  int radtype_X2_subsizes[4] = {N1 + 2 * NGR, NGR, N3 + 2 * NGR, RAD_NUM_TYPES};
  MPI_Type_create_subarray(4, radtype_sizes, radtype_X2_subsizes,
      radtype_starts, MPI_ORDER_C, grid_radtype_type,
      &grid_radtype_face_subtype[2]);
  MPI_Type_commit(&grid_radtype_face_subtype[2]);

  int radtype_X3_subsizes[4] = {N1 + 2 * NGR, N2 + 2 * NGR, NGR, RAD_NUM_TYPES};
  MPI_Type_create_subarray(4, radtype_sizes, radtype_X3_subsizes,
      radtype_starts, MPI_ORDER_C, grid_radtype_type,
      &grid_radtype_face_subtype[3]);
  MPI_Type_commit(&grid_radtype_face_subtype[3]);

  // Create MPI photon type
  MPI_Type_contiguous(sizeof(struct of_photon *), MPI_BYTE, &mpi_pointer_type);
  MPI_Type_commit(&mpi_pointer_type);
  MPI_Datatype type[PH_ELEM] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
      MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT,
      mpi_pointer_type};
  int          blocklen[PH_ELEM] = {
      NDIM * NSUP, NDIM * NSUP, NDIM * NSUP, 1, 1, 1, 1, NDIM, 1, 1, 1, 1};
  MPI_Aint         disp[PH_ELEM];
  struct of_photon tmp;
  MPI_Get_address(&(tmp.X[0][0]), &(disp[0]));
  MPI_Get_address(&(tmp.Kcov[0][0]), &(disp[1]));
  MPI_Get_address(&(tmp.Kcon[0][0]), &(disp[2]));
  MPI_Get_address(&(tmp.w), &(disp[3]));
  MPI_Get_address(&(tmp.KdotKprev), &(disp[4]));
  MPI_Get_address(&(tmp.type), &(disp[5]));
  MPI_Get_address(&(tmp.nscatt), &(disp[6]));
  MPI_Get_address(tmp.origin, &(disp[7]));
  MPI_Get_address(&(tmp.t0), &(disp[8]));
  MPI_Get_address(&(tmp.is_tracked), &(disp[9]));
  MPI_Get_address(&(tmp.osc_count), &(disp[10]));
  MPI_Get_address(&(tmp.next), &(disp[11]));
  MPI_Aint base;
  MPI_Get_address(&tmp, &base);
  for (int n = 0; n < PH_ELEM; n++) {
    disp[n] -= base;
  }
  MPI_Type_create_struct(PH_ELEM, blocklen, disp, type, &mpi_photon_type);
  MPI_Type_commit(&mpi_photon_type);
#endif // RADIATION

  if (mpi_io_proc()) {
    fprintf(stdout, "MPI initialized on %i processors\n", numprocs);
  }
}

void sync_mpi_boundaries_X1L(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[N1][NG][NG][0], 1, prim_face_subtype[1], neighbors[2][1][1],
      0, &Pr[0][NG][NG][0], 1, prim_face_subtype[1], neighbors[0][1][1], 0,
      comm, &status);
  MPI_Sendrecv(&pflag[N1][NG][NG], 1, pflag_face_subtype[1], neighbors[2][1][1],
      0, &pflag[0][NG][NG], 1, pflag_face_subtype[1], neighbors[0][1][1], 0,
      comm, &status);
}

void sync_mpi_boundaries_X1R(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[NG][NG][NG][0], 1, prim_face_subtype[1], neighbors[0][1][1],
      1, &Pr[N1 + NG][NG][NG][0], 1, prim_face_subtype[1], neighbors[2][1][1],
      1, comm, &status);
  MPI_Sendrecv(&pflag[NG][NG][NG], 1, pflag_face_subtype[1], neighbors[0][1][1],
      1, &pflag[N1 + NG][NG][NG], 1, pflag_face_subtype[1], neighbors[2][1][1],
      1, comm, &status);
}

void sync_mpi_boundaries_X2L(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[0][N2][NG][0], 1, prim_face_subtype[2], neighbors[1][2][1],
      2, &Pr[0][0][NG][0], 1, prim_face_subtype[2], neighbors[1][0][1], 2, comm,
      &status);
  MPI_Sendrecv(&pflag[0][N2][NG], 1, pflag_face_subtype[2], neighbors[1][2][1],
      2, &pflag[0][0][NG], 1, pflag_face_subtype[2], neighbors[1][0][1], 2,
      comm, &status);
}

void sync_mpi_boundaries_X2R(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[0][NG][NG][0], 1, prim_face_subtype[2], neighbors[1][0][1],
      3, &Pr[0][N2 + NG][NG][0], 1, prim_face_subtype[2], neighbors[1][2][1], 3,
      comm, &status);
  MPI_Sendrecv(&pflag[0][NG][NG], 1, pflag_face_subtype[2], neighbors[1][0][1],
      3, &pflag[0][N2 + NG][NG], 1, pflag_face_subtype[2], neighbors[1][2][1],
      3, comm, &status);
}

void sync_mpi_boundaries_X3L(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[0][0][N3][0], 1, prim_face_subtype[3], neighbors[1][1][2], 4,
      &Pr[0][0][0][0], 1, prim_face_subtype[3], neighbors[1][1][0], 4, comm,
      &status);
  MPI_Sendrecv(&pflag[0][0][N3], 1, pflag_face_subtype[3], neighbors[1][1][2],
      4, &pflag[0][0][0], 1, pflag_face_subtype[3], neighbors[1][1][0], 4, comm,
      &status);
}

void sync_mpi_boundaries_X3R(grid_prim_type Pr) {
  MPI_Status status;
  MPI_Sendrecv(&Pr[0][0][NG][0], 1, prim_face_subtype[3], neighbors[1][1][0], 5,
      &Pr[0][0][N3 + NG][0], 1, prim_face_subtype[3], neighbors[1][1][2], 5,
      comm, &status);
  MPI_Sendrecv(&pflag[0][0][NG], 1, pflag_face_subtype[3], neighbors[1][1][0],
      5, &pflag[0][0][N3 + NG], 1, pflag_face_subtype[3], neighbors[1][1][2], 5,
      comm, &status);
}

#if RADIATION
void reduce_radG_buf() {
  ZSLOOP(-NGR, N1 + NGR - 1, -NGR, N2 + NGR - 1, -NGR, N3 + NGR - 1) {
    for (int mu = 0; mu < NDIM + NRADCOMP; mu++) {
      radG[i][j][k][mu] += radG_buf[i][j][k][mu];
      radG_buf[i][j][k][mu] = 0.;
    }
  }
}
void sync_radG() {
  memset(&radG_buf[0][0][0][0], 0,
      (N1 + 2 * NG) * (N2 + 2 * NG) * (N3 + 2 * NG) * NDIM * sizeof(double));

  MPI_Status status;

  MPI_Sendrecv(&radG[N1 + NG][NG - NGR][NG - NGR][0], 1, radG_face_subtype[1],
      neighbors[2][1][1], 0, &radG_buf[NG][NG - NGR][NG - NGR][0], 1,
      radG_face_subtype[1], neighbors[0][1][1], 0, comm, &status);
  MPI_Sendrecv(&radG[NG - NGR][NG - NGR][NG - NGR][0], 1, radG_face_subtype[1],
      neighbors[0][1][1], 1, &radG_buf[N1 + NG - NGR][NG - NGR][NG - NGR][0], 1,
      radG_face_subtype[1], neighbors[2][1][1], 1, comm, &status);
  reduce_radG_buf();

  MPI_Sendrecv(&radG[NG - NGR][N2 + NG][NG - NGR][0], 1, radG_face_subtype[2],
      neighbors[1][2][1], 2, &radG_buf[NG - NGR][NG][NG - NGR][0], 1,
      radG_face_subtype[2], neighbors[1][0][1], 2, comm, &status);
  MPI_Sendrecv(&radG[NG - NGR][NG - NGR][NG - NGR][0], 1, radG_face_subtype[2],
      neighbors[1][0][1], 3, &radG_buf[NG - NGR][N2 + NG - NGR][NG - NGR][0], 1,
      radG_face_subtype[2], neighbors[1][2][1], 3, comm, &status);
  reduce_radG_buf();

  MPI_Sendrecv(&radG[NG - NGR][NG - NGR][N3 + NG][0], 1, radG_face_subtype[3],
      neighbors[1][1][2], 4, &radG_buf[NG - NGR][NG - NGR][NG][0], 1,
      radG_face_subtype[3], neighbors[1][1][0], 4, comm, &status);
  MPI_Sendrecv(&radG[NG - NGR][NG - NGR][NG - NGR][0], 1, radG_face_subtype[3],
      neighbors[1][1][0], 5, &radG_buf[NG - NGR][NG - NGR][N3 + NG - NGR][0], 1,
      radG_face_subtype[3], neighbors[1][1][2], 5, comm, &status);
  reduce_radG_buf();
}

void reduce_Jrad_buf() {
  for (int n = 0; n < MAXNSCATT + 2; n++) {
    ZSLOOP(-NGR, N1 + NGR - 1, -NGR, N2 + NGR - 1, -NGR, N3 + NGR - 1) {
      Jrad[n][i][j][k] += Jrad_buf[n][i][j][k];
      Jrad_buf[n][i][j][k] = 0.;
    }
  }
}
void sync_Jrad() {
  memset(&Jrad_buf[0][0][0][0], 0, (MAXNSCATT + 2) * N123G * sizeof(double));

  MPI_Status status;

  MPI_Sendrecv(&Jrad[0][N1 + NG][NG - NGR][NG - NGR], 1, Jrad_face_subtype[1],
      neighbors[2][1][1], 0, &Jrad_buf[0][NG][NG - NGR][NG - NGR], 1,
      Jrad_face_subtype[1], neighbors[0][1][1], 0, comm, &status);
  MPI_Sendrecv(&Jrad[0][NG - NGR][NG - NGR][NG - NGR], 1, Jrad_face_subtype[1],
      neighbors[0][1][1], 1, &Jrad_buf[0][N1 + NG - NGR][NG - NGR][NG - NGR], 1,
      Jrad_face_subtype[1], neighbors[2][1][1], 1, comm, &status);
  reduce_Jrad_buf();

  MPI_Sendrecv(&Jrad[0][NG - NGR][N2 + NG][NG - NGR], 1, Jrad_face_subtype[2],
      neighbors[1][2][1], 2, &Jrad_buf[0][NG - NGR][NG][NG - NGR], 1,
      Jrad_face_subtype[2], neighbors[1][0][1], 2, comm, &status);
  MPI_Sendrecv(&Jrad[0][NG - NGR][NG - NGR][NG - NGR], 1, Jrad_face_subtype[2],
      neighbors[1][0][1], 3, &Jrad_buf[0][NG - NGR][N2 + NG - NGR][NG - NGR], 1,
      Jrad_face_subtype[2], neighbors[1][2][1], 3, comm, &status);
  reduce_Jrad_buf();

  MPI_Sendrecv(&Jrad[0][NG - NGR][NG - NGR][N3 + NG], 1, Jrad_face_subtype[3],
      neighbors[1][1][2], 4, &Jrad_buf[0][NG - NGR][NG - NGR][NG], 1,
      Jrad_face_subtype[3], neighbors[1][1][0], 4, comm, &status);
  MPI_Sendrecv(&Jrad[0][NG - NGR][NG - NGR][NG - NGR], 1, Jrad_face_subtype[3],
      neighbors[1][1][0], 5, &Jrad_buf[0][NG - NGR][NG - NGR][N3 + NG - NGR], 1,
      Jrad_face_subtype[3], neighbors[1][1][2], 5, comm, &status);
  reduce_Jrad_buf();

  for (int i = 0; i < NG; i++) {
    JLOOPALL KLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
  for (int i = N1 + NG; i < N1 + 2 * NG; i++) {
    JLOOPALL KLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
  for (int j = 0; j < NG; j++) {
    ILOOPALL KLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
  for (int j = N2 + NG; j < N2 + 2 * NG; j++) {
    ILOOPALL KLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
  for (int k = 0; k < NG; k++) {
    ILOOPALL JLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
  for (int k = N3 + NG; k < N3 + 2 * NG; k++) {
    ILOOPALL JLOOPALL JRADLOOP { Jrad[n][i][j][k] = 0.; }
  }
}

void reduce_radtype_buf(grid_radtype_type v) {
  TYPELOOP {
    ZSLOOP(-NGR, N1 + NGR - 1, -NGR, N2 + NGR - 1, -NGR, N3 + NGR - 1) {
      v[i][j][k][itp] += radtype_buf[i][j][k][itp];
      radtype_buf[i][j][k][itp] = 0.0;
    }
  }
}

void sync_radtype_vec(grid_radtype_type v) {
  memset(&radtype_buf[0][0][0][0], 0, RAD_NUM_TYPES * N123G * sizeof(double));

  MPI_Status status;

  MPI_Sendrecv(&v[N1 + NG][NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[1], neighbors[2][1][1], 0,
      &radtype_buf[NG][NG - NGR][NG - NGR][0], 1, grid_radtype_face_subtype[1],
      neighbors[0][1][1], 0, comm, &status);
  MPI_Sendrecv(&v[NG - NGR][NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[1], neighbors[0][1][1], 1,
      &radtype_buf[N1 + NG - NGR][NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[1], neighbors[2][1][1], 1, comm, &status);
  reduce_radtype_buf(v);

  MPI_Sendrecv(&v[NG - NGR][N2 + NG][NG - NGR][0], 1,
      grid_radtype_face_subtype[2], neighbors[1][2][1], 2,
      &radtype_buf[NG - NGR][NG][NG - NGR][0], 1, grid_radtype_face_subtype[2],
      neighbors[1][0][1], 2, comm, &status);
  MPI_Sendrecv(&v[NG - NGR][NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[2], neighbors[1][0][1], 3,
      &radtype_buf[NG - NGR][N2 + NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[2], neighbors[1][2][1], 3, comm, &status);
  reduce_radtype_buf(v);

  MPI_Sendrecv(&v[NG - NGR][NG - NGR][N3 + NG][0], 1,
      grid_radtype_face_subtype[3], neighbors[1][1][2], 4,
      &radtype_buf[NG - NGR][NG - NGR][NG][0], 1, grid_radtype_face_subtype[3],
      neighbors[1][1][0], 4, comm, &status);
  MPI_Sendrecv(&v[NG - NGR][NG - NGR][NG - NGR][0], 1,
      grid_radtype_face_subtype[3], neighbors[1][1][0], 5,
      &radtype_buf[NG - NGR][NG - NGR][N3 + NG - NGR][0], 1,
      grid_radtype_face_subtype[3], neighbors[1][1][2], 5, comm, &status);
  reduce_radtype_buf(v);

  for (int i = 0; i < NG; i++) {
    JLOOPALL KLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
  for (int i = N1 + NG; i < N1 + 2 * NG; i++) {
    JLOOPALL KLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
  for (int j = 0; j < NG; j++) {
    ILOOPALL KLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
  for (int j = N2 + NG; j < N2 + 2 * NG; j++) {
    ILOOPALL KLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
  for (int k = 0; k < NG; k++) {
    ILOOPALL JLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
  for (int k = N3 + NG; k < N3 + 2 * NG; k++) {
    ILOOPALL JLOOPALL TYPELOOP { v[i][j][k][itp] = 0.; }
  }
}

void sync_mpi_photons(
    struct of_photon **ph_mpi, grid_prim_type P, double t, double dt) {
  MPI_Status status;

  struct of_photon *ph_send_L[NDIM]  = {NULL, NULL, NULL, NULL};
  struct of_photon *ph_send_R[NDIM]  = {NULL, NULL, NULL, NULL};
  struct of_photon *ph_recv          = NULL;
  int               nph_send_L[NDIM] = {0, 0, 0, 0};
  int               nph_send_R[NDIM] = {0, 0, 0, 0};
  int               nph_recv_L[NDIM] = {0, 0, 0, 0};
  int               nph_recv_R[NDIM] = {0, 0, 0, 0};
  int               nph_recv         = 0;

  struct of_photon *ph_tmp_recv = NULL;
  struct of_photon *ph          = *ph_mpi;
  struct of_photon *ph_send_buf, *ph_recv_buf, *to_free;
  double            X[NDIM], Kcov[NDIM], Kcon[NDIM];
  while (ph != NULL) {
    get_X_K_interp(ph, t + dt, P, X, Kcov, Kcon);
    step_sent++;
    step_tot--;
    for (int i = 1; i < NDIM; i++) {
      if (X[i] < startx_proc[i]) {
        swap_ph(&ph, &ph_send_L[i]);
        nph_send_L[i]++;
        break;
      } else if (X[i] > stopx_proc[i]) {
        swap_ph(&ph, &ph_send_R[i]);
        nph_send_R[i]++;
        break;
      } else if (i == NDIM - 1) {
        step_sent--;
        step_tot++;
        ph = ph->next;
      }
    }
  }

  // Send right receive left
  int neighbor_L, neighbor_R;
  for (int dir = 1; dir < NDIM; dir++) {
    if (dir == 1) {
      neighbor_L = neighbors[0][1][1];
      neighbor_R = neighbors[2][1][1];
    } else if (dir == 2) {
      neighbor_L = neighbors[1][0][1];
      neighbor_R = neighbors[1][2][1];
    } else if (dir == 3) {
      neighbor_L = neighbors[1][1][0];
      neighbor_R = neighbors[1][1][2];
    }

    MPI_Sendrecv(&nph_send_R[dir], 1, MPI_INT, neighbor_R, 0, &nph_recv_L[dir],
        1, MPI_INT, neighbor_L, 0, comm, &status);

    ph_send_buf = malloc(nph_send_R[dir] * sizeof(struct of_photon));
    ph_recv_buf = malloc(nph_recv_L[dir] * sizeof(struct of_photon));

    for (int n = 0; n < nph_send_R[dir]; n++) {
      memcpy(&ph_send_buf[n], ph_send_R[dir], sizeof(struct of_photon));
      to_free        = ph_send_R[dir];
      ph_send_R[dir] = ph_send_R[dir]->next;
      free(to_free);
    }

    MPI_Sendrecv(ph_send_buf, nph_send_R[dir], mpi_photon_type, neighbor_R, 0,
        ph_recv_buf, nph_recv_L[dir], mpi_photon_type, neighbor_L, 0, comm,
        &status);
    free(ph_send_buf);
    for (int n = 0; n < nph_recv_L[dir]; n++) {
      ph = malloc(sizeof(struct of_photon));
      memcpy(ph, &(ph_recv_buf[n]), sizeof(struct of_photon));
      swap_ph(&ph, &ph_tmp_recv);
    }
    free(ph_recv_buf);

    // Send left receive right
    MPI_Sendrecv(&nph_send_L[dir], 1, MPI_INT, neighbor_L, 0, &nph_recv_R[dir],
        1, MPI_INT, neighbor_R, 0, comm, &status);

    ph_send_buf = malloc(nph_send_L[dir] * sizeof(struct of_photon));
    ph_recv_buf = malloc(nph_recv_R[dir] * sizeof(struct of_photon));

    for (int n = 0; n < nph_send_L[dir]; n++) {
      memcpy(&ph_send_buf[n], ph_send_L[dir], sizeof(struct of_photon));
      to_free        = ph_send_L[dir];
      ph_send_L[dir] = ph_send_L[dir]->next;
      free(to_free);
    }

    MPI_Sendrecv(ph_send_buf, nph_send_L[dir], mpi_photon_type, neighbor_L, 0,
        ph_recv_buf, nph_recv_R[dir], mpi_photon_type, neighbor_R, 0, comm,
        &status);
    free(ph_send_buf);
    for (int n = 0; n < nph_recv_R[dir]; n++) {
      ph = malloc(sizeof(struct of_photon));
      memcpy(ph, &(ph_recv_buf[n]), sizeof(struct of_photon));
      swap_ph(&ph, &ph_tmp_recv);
    }
    free(ph_recv_buf);

    // Send or keep received superphotons in ph_tmp_recv?
    while (ph_tmp_recv != NULL) {
      get_X_K_interp(ph_tmp_recv, t + dt, P, X, Kcov, Kcon);
      bound_rad_transport(X, ph_tmp_recv, 1);
      if (dir == NDIM - 1) {
        swap_ph(&ph_tmp_recv, &ph_recv);
        nph_recv++;
      } else {

        for (int i = dir + 1; i < NDIM; i++) {
          if (X[i] < startx_proc[i]) {
            swap_ph(&ph_tmp_recv, &ph_send_L[i]);
            nph_send_L[i]++;
            break;
          } else if (X[i] > stopx_proc[i]) {
            swap_ph(&ph_tmp_recv, &ph_send_R[i]);
            nph_send_R[i]++;
            break;
          } else if (i == NDIM - 1) {
            swap_ph(&ph_tmp_recv, &ph_recv);
            nph_recv++;
          }
        }
      }
    }
  }

  step_rcvd += nph_recv;
  step_tot += nph_recv;

  // Distribute received superphotons equitably between threads
  int thread_start  = (int)(get_rand() * nthreads);
  int ph_per_thread = nph_recv / nthreads;
  if (ph_per_thread * nthreads < nph_recv) {
    ph_per_thread++;
  }
  for (int n = thread_start; n < thread_start + nthreads; n++) {
    int index    = n % nthreads;
    int ph_added = 0;

    while (ph_added < ph_per_thread && ph_recv != NULL) {
      swap_ph(&ph_recv, &(photon_lists[index]));
      ph_added++;
    }
  }
}
#endif // RADIATION

int mpi_nprocs() { return numprocs; }

double mpi_max(double f) {
  double fmax;
  MPI_Allreduce(&f, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
  return fmax;
}

int mpi_max_int(int f) {
  int fmax;
  MPI_Allreduce(&f, &fmax, 1, MPI_INT, MPI_MAX, comm);
  return fmax;
}

double mpi_min(double f) {
  double fmin;
  MPI_Allreduce(&f, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);
  return fmin;
}

double mpi_reduce(double f) {
  double fsum;
  MPI_Allreduce(&f, &fsum, 1, MPI_DOUBLE, MPI_SUM, comm);
  return fsum;
}

int mpi_reduce_int(int f) {
  int fsum;
  MPI_Allreduce(&f, &fsum, 1, MPI_INT, MPI_SUM, comm);
  return fsum;
}

void mpi_dbl_allreduce_array(double *A, int size) {
  MPI_Allreduce(MPI_IN_PLACE, A, size, MPI_DOUBLE, MPI_SUM, comm);
}

int mpi_io_proc() { return (rank == 0 ? 1 : 0); }

void mpi_int_broadcast(int *val) { MPI_Bcast(val, 1, MPI_INT, 0, comm); }

void mpi_int_broadcast_array(int *val, int size) {
  MPI_Bcast(val, size, MPI_INT, 0, comm);
}

void mpi_int_broadcast_proc(int *val, int root) {
  MPI_Bcast(val, 1, MPI_INT, root, comm);
}

void mpi_dbl_broadcast(double *val) { MPI_Bcast(val, 1, MPI_DOUBLE, 0, comm); }

void mpi_allgather_int1(int *buffer, int element) {
  MPI_Allgather(&element, 1, MPI_INT, buffer, 1, MPI_INT, comm);
}

int mpi_accumulate_int(int my_val) {
  int *buffer = safe_malloc(mpi_nprocs() * sizeof(int));

  mpi_allgather_int1(buffer, my_val);

  int out = 0;
  for (int i = 0; i < mpi_myrank(); i++) {
    out += buffer[i];
  }

  free(buffer);
  return out;
}

double mpi_io_reduce(double val) {

  double local;
  MPI_Reduce(&val, &local, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
  return local;
}

double mpi_io_max(double val) {

  double local;
  MPI_Reduce(&val, &local, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
  return local;
}

int mpi_myrank() { return rank; }

void mpi_sync_output() {
  fflush(stderr);
  MPI_Barrier(comm);
}

void mpi_barrier() { MPI_Barrier(comm); }

int mpi_is_periodic(int dir) { return periodic[dir - 1]; }

#if RADIATION
static double nuLnu_local[NULNU_IDX0][NTH][NPHI][NU_BINS_SPEC];
void          mpi_reduce_nuLnu() {
  MPI_Reduce(nuLnu, nuLnu_local, (NULNU_IDX0)*NTH * NPHI * NU_BINS_SPEC,
      MPI_DOUBLE, MPI_SUM, 0, comm);
  memcpy(nuLnu, nuLnu_local,
      (NULNU_IDX0)*NTH * NPHI * NU_BINS_SPEC * sizeof(double));
}
#endif // RADIATION
