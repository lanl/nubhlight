/******************************************************************************
 *                                                                            *
 * IO.C                                                                       *
 *                                                                            *
 * HDF5 OUTPUT AND RESTART                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <hdf5.h>
#include <hdf5_hl.h>

#define MAX_GRID_DIM (5)

void    write_array(void *data, const char *name, hsize_t rank, hsize_t *fdims,
       hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart,
       hsize_t type);
void    read_array(void *data, const char *name, hsize_t rank, hsize_t *fdims,
       hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart,
       hsize_t type);
hsize_t product_hsize_t(hsize_t a[], int size);

// Some macro tricks to drastically reduce number of lines of code
#define TYPE_FLOAT H5T_NATIVE_FLOAT
#define TYPE_DBL H5T_NATIVE_DOUBLE
#define TYPE_INT H5T_NATIVE_INT
#define TYPE_STR H5T_NATIVE_CHAR

#define WRITE_HDR(x, type)                                                     \
  write_array((void *)&x, #x, 1, fdims_hdr, fstart_hdr, fcount_hdr, mdims_hdr, \
      mstart_hdr, type)
#define READ_HDR(x, type)                                                     \
  read_array((void *)&x, #x, 1, fdims_hdr, fstart_hdr, fcount_hdr, mdims_hdr, \
      mstart_hdr, type)

#define WRITE_ARRAY(x, rank, fdims, fstart, fcount, mdims, mstart, type) \
  write_array((void *)x, #x, rank, fdims, fstart, fcount, mdims, mstart, type)
#define WRITE_GRID(x, type)                                           \
  write_array((void *)x, #x, 3, fdims_grid, fstart_grid, fcount_grid, \
      mdims_grid, mstart_grid, type)
#define WRITE_GRID_NO_GHOSTS(x, type)                                 \
  write_array((void *)x, #x, 3, fdims_grid, fstart_grid, fcount_grid, \
      mdims_grid_noghost, mstart_grid_noghost, type)
#define WRITE_PRIM(x, type)                                           \
  write_array((void *)x, #x, 4, fdims_prim, fstart_prim, fcount_prim, \
      mdims_prim, mstart_prim, type)
#define WRITE_VEC(x, type)                                                    \
  write_array((void *)x, #x, 4, fdims_vec, fstart_vec, fcount_vec, mdims_vec, \
      mstart_vec, type)
#define WRITE_RADG(x, type)                                          \
  write_array((void *)x, #x, 4, fdims_radg, fstart_vec, fcount_radg, \
      mdims_radg, mstart_vec, type)
#define WRITE_RADTYPEVEC(x, type)                             \
  write_array((void *)x, #x, 4, fdims_radtypevec, fstart_vec, \
      fcount_radtypevec, mdims_radtypevec, mstart_vec, type)
#define WRITE_VEC_NO_GHOSTS(x, type)                               \
  write_array((void *)x, #x, 4, fdims_vec, fstart_vec, fcount_vec, \
      mdims_vec_noghost, mstart_vec_noghost, type)
#define WRITE_TENSOR(x, type)                                         \
  write_array((void *)x, #x, 5, fdims_tens, fstart_tens, fcount_tens, \
      mdims_tens, mstart_tens, type)
#define WRITE_TENSOR_NO_GHOSTS(x, type)                               \
  write_array((void *)x, #x, 5, fdims_tens, fstart_tens, fcount_tens, \
      mdims_tens_noghost, mstart_tens_noghost, type)
#define READ_ARRAY(x, rank, fdims, fstart, fcount, mdims, mstart, type) \
  read_array((void *)x, #x, rank, fdims, fstart, fcount, mdims, mstart, type)
#define READ_GRID(x, type)                                           \
  read_array((void *)x, #x, 3, fdims_grid, fstart_grid, fcount_grid, \
      mdims_grid, mstart_grid, type)
#define READ_PRIM(x, type)                                           \
  read_array((void *)x, #x, 4, fdims_prim, fstart_prim, fcount_prim, \
      mdims_prim, mstart_prim, type)
#define READ_VEC(x, type)                                                    \
  read_array((void *)x, #x, 4, fdims_vec, fstart_vec, fcount_vec, mdims_vec, \
      mstart_vec, type)
#define READ_RADTYPEVEC(x, type)                             \
  read_array((void *)x, #x, 4, fdims_radtypevec, fstart_vec, \
      fcount_radtypevec, mdims_radtypevec, mstart_vec, type)

hid_t   file_id, plist_id;
hid_t   filespace, memspace;
hid_t   filespace_hdr, memspace_hdr;
hsize_t hdr_dims[1] = {1}, zero = 0, one = 1;
hsize_t mem_start[MAX_GRID_DIM], file_grid_start[MAX_GRID_DIM];
hsize_t file_grid_count[MAX_GRID_DIM], file_grid_dims[MAX_GRID_DIM];
hsize_t file_hdr_start[1] = {0}, file_hdr_count[1], file_hdr_dims[1] = {1};
hsize_t hdr_rank = 1, grid_rank;

hsize_t fdims_grid[3] = {N1TOT, N2TOT, N3TOT};
hsize_t fstart_grid[3];
hsize_t fcount_grid[3] = {N1, N2, N3};
hsize_t mdims_grid[3]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
hsize_t mstart_grid[3] = {NG, NG, NG};

hsize_t mdims_grid_noghost[3]  = {N1, N2, N3};
hsize_t mstart_grid_noghost[3] = {0, 0, 0};

hsize_t fdims_prim[4] = {N1TOT, N2TOT, N3TOT, NVAR};
hsize_t fstart_prim[4];
hsize_t fcount_prim[4] = {N1, N2, N3, NVAR};
hsize_t mdims_prim[4]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NVAR};
hsize_t mstart_prim[4] = {NG, NG, NG, 0};

hsize_t fdims_vec[4] = {N1TOT, N2TOT, N3TOT, NDIM};
hsize_t fstart_vec[4];
hsize_t fcount_vec[4] = {N1, N2, N3, NDIM};
hsize_t mdims_vec[4]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NDIM};
hsize_t mstart_vec[4] = {NG, NG, NG, 0};

hsize_t fdims_radtypevec[4]  = {N1TOT, N2TOT, N3TOT, RAD_NUM_TYPES};
hsize_t fcount_radtypevec[4] = {N1, N2, N3, RAD_NUM_TYPES};
hsize_t mdims_radtypevec[4]  = {
    N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, RAD_NUM_TYPES};

hsize_t fdims_radg[4]  = {N1TOT, N2TOT, N3TOT, NDIM + NRADCOMP};
hsize_t fcount_radg[4] = {N1, N2, N3, NDIM + NRADCOMP};
hsize_t mdims_radg[4]  = {
    N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NDIM + NRADCOMP};

hsize_t mdims_vec_noghost[4]  = {N1, N2, N3, NDIM};
hsize_t mstart_vec_noghost[4] = {0, 0, 0, 0};

hsize_t fdims_tens[5] = {N1TOT, N2TOT, N3TOT, NDIM, NDIM};
hsize_t fstart_tens[5];
hsize_t fcount_tens[5] = {N1, N2, N3, NDIM, NDIM};
hsize_t mdims_tens[5]  = {N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG, NDIM, NDIM};
hsize_t mstart_tens[5] = {NG, NG, NG, 0, 0};

hsize_t mdims_tens_noghost[5]  = {N1, N2, N3, NDIM, NDIM};
hsize_t mstart_tens_noghost[5] = {0, 0, 0, 0, 0};

hsize_t fdims_hdr[1]  = {1};
hsize_t fstart_hdr[1] = {0};
hsize_t fcount_hdr[1] = {1};
hsize_t mdims_hdr[1]  = {1};
hsize_t mstart_hdr[1] = {0};

hid_t phfiletype, phmemtype, trackphfiletype, trackphmemtype;

const char *vnams_base[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
const char *vnams[NVAR];
static char version[STRLEN];
static int  dump_id = 0, restart_id = 0, restart_perm_id = 0, fdump_id = 0;
#if RADIATION && TRACERS
static int dumptrace_id = 0;
#endif
static grid_double_type divb;
#if OUTPUT_EOSVARS
static grid_double_type PRESS, ENT, TEMP, CS2;
#endif
#if RADIATION
static grid_double_type tau_cool;
static grid_double_type dtau_scatt, dtau_tot;
static grid_double_type denominator_scatt, denominator_tot;
#endif

void init_io() {
  strcpy(dumpdir, "dumps/");
  strcpy(restartdir, "restarts/");
  strcpy(xmfdir, "dumps/xmf/");
#if RADIATION && TRACERS
  strcpy(tracerdir, "dumps/tracers/");
#endif
  strcpy(version, VERSION);
  int len = strlen(outputdir);
  memmove(dumpdir + len, dumpdir, strlen(dumpdir) + 1);
  memmove(restartdir + len, restartdir, strlen(restartdir) + 1);
  memmove(xmfdir + len, xmfdir, strlen(xmfdir) + 1);
#if RADIATION && TRACERS
  memmove(tracerdir + len, tracerdir, strlen(tracerdir) + 1);
#endif
  for (int n = 0; n < len; ++n) {
    dumpdir[n]    = outputdir[n];
    restartdir[n] = outputdir[n];
    xmfdir[n]     = outputdir[n];
#if RADIATION && TRACERS
    tracerdir[n] = outputdir[n];
#endif
  }

  if (mpi_io_proc()) {
    char mkdircall[STRLEN];
    strcpy(mkdircall, "mkdir -p ");
    strcat(mkdircall, dumpdir);
    strcat(mkdircall, " ");
    strcat(mkdircall, restartdir);
    strcat(mkdircall, " ");
    strcat(mkdircall, xmfdir);
#if RADIATION && TRACERS
    strcat(mkdircall, " ");
    strcat(mkdircall, tracerdir);
#endif
    safe_system(mkdircall);
  }

  fstart_grid[0] = global_start[1];
  fstart_grid[1] = global_start[2];
  fstart_grid[2] = global_start[3];
  fstart_prim[0] = global_start[1];
  fstart_prim[1] = global_start[2];
  fstart_prim[2] = global_start[3];
  fstart_prim[3] = 0;
  fstart_vec[0]  = global_start[1];
  fstart_vec[1]  = global_start[2];
  fstart_vec[2]  = global_start[3];
  fstart_vec[3]  = 0;
  fstart_tens[0] = global_start[1];
  fstart_tens[1] = global_start[2];
  fstart_tens[2] = global_start[3];
  fstart_tens[3] = 0;
  fstart_tens[4] = 0;
  if (!mpi_io_proc())
    fcount_hdr[0] = 0;

  filespace_hdr     = H5Screate_simple(1, file_hdr_dims, NULL);
  file_hdr_count[0] = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace_hdr, H5S_SELECT_SET, file_hdr_start, NULL,
      file_hdr_count, NULL);
  memspace_hdr = H5Screate_simple(hdr_rank, file_hdr_count, NULL);

#if RADIATION
  // Create custom datatypes for arrays inside struct
  hsize_t array_dim[2]   = {NSUP, NDIM};
  hid_t   array_tid      = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, array_dim);
  hsize_t int_vec_dim[1] = {NDIM};
  hid_t   int_vec_tid    = H5Tarray_create(H5T_NATIVE_INT, 1, int_vec_dim);

  // Memory contiguous in file
  phfiletype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  int offset = 0;
  // Function that accepts phfiletype, name, offset, size
  H5Tinsert(phfiletype, "X", offset, array_tid);
  offset += NSUP * NDIM * sizeof(double);
  H5Tinsert(phfiletype, "Kcov", offset, array_tid);
  offset += NSUP * NDIM * sizeof(double);
  H5Tinsert(phfiletype, "Kcon", offset, array_tid);
  offset += NSUP * NDIM * sizeof(double);
  H5Tinsert(phfiletype, "w", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "KdotKprev", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "type", offset, H5T_NATIVE_INT);
  offset += sizeof(int);
  H5Tinsert(phfiletype, "nscatt", offset, H5T_NATIVE_INT);
  offset += sizeof(int);
  H5Tinsert(phfiletype, "origin", offset, int_vec_tid);
  offset += NDIM * sizeof(int);
  H5Tinsert(phfiletype, "t0", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(phfiletype, "is_tracked", offset, H5T_NATIVE_INT);
  offset += sizeof(int);
  H5Tinsert(phfiletype, "osc_count", offset, H5T_NATIVE_INT);

  // Use HOFFSET to account for struct padding in memory
  phmemtype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  H5Tinsert(phmemtype, "X", HOFFSET(struct of_photon, X), array_tid);
  H5Tinsert(phmemtype, "Kcov", HOFFSET(struct of_photon, Kcov), array_tid);
  H5Tinsert(phmemtype, "Kcon", HOFFSET(struct of_photon, Kcon), array_tid);
  H5Tinsert(phmemtype, "w", HOFFSET(struct of_photon, w), H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "KdotKprev", HOFFSET(struct of_photon, KdotKprev),
      H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "type", HOFFSET(struct of_photon, type), H5T_NATIVE_INT);
  H5Tinsert(
      phmemtype, "nscatt", HOFFSET(struct of_photon, nscatt), H5T_NATIVE_INT);
  H5Tinsert(
      phmemtype, "origin", HOFFSET(struct of_photon, origin), int_vec_tid);
  H5Tinsert(phmemtype, "t0", HOFFSET(struct of_photon, t0), H5T_NATIVE_DOUBLE);
  H5Tinsert(phmemtype, "is_tracked", HOFFSET(struct of_photon, is_tracked),
      H5T_NATIVE_INT);
  H5Tinsert(phmemtype, "osc_count",
      HOFFSET(struct of_photon, osc_count), H5T_NATIVE_INT);

  trackphfiletype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_photon));
  offset          = 0;
  H5Tinsert(trackphfiletype, "X1", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "X2", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "X3", offset, H5T_NATIVE_DOUBLE);
  offset += sizeof(double);
  H5Tinsert(trackphfiletype, "nscatt", offset, H5T_NATIVE_INT);

  trackphmemtype = H5Tcreate(H5T_COMPOUND, sizeof(struct of_track_photon));
  H5Tinsert(trackphmemtype, "X1", HOFFSET(struct of_track_photon, X1),
      H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "X2", HOFFSET(struct of_track_photon, X2),
      H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "X3", HOFFSET(struct of_track_photon, X3),
      H5T_NATIVE_DOUBLE);
  H5Tinsert(trackphmemtype, "nscatt", HOFFSET(struct of_track_photon, nscatt),
      H5T_NATIVE_INT);
#endif // RADIATION

  BASELOOP { vnams[ip] = vnams_base[ip]; }
#if NVAR_PASSIVE > 0
  PASSLOOP { vnams[ipass] = PASSNAME(ipass); }
#endif
#if ELECTRONS
  vnams[KEL]  = "KEL";
  vnams[KTOT] = "KTOT";
#endif // ELECTRONS
}

#if RADIATION
void track_ph() {
#if TRACK_PH
  static int track_id = 0;
  // Count number of tracked photons on this processor
  int num_tracked = 0;
#pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    while (ph != NULL) {
      if (ph->is_tracked == 1) {
#pragma omp atomic
        num_tracked++;
      }
      ph = ph->next;
    }
  }

  // Create buffer for output
  struct of_track_photon *io_track =
      safe_malloc((size_t)num_tracked * sizeof(struct of_track_photon));
  int    nio = 0;
  double X[NDIM], Kcov[NDIM], Kcon[NDIM];
  for (int n = 0; n < nthreads; n++) {
    struct of_photon *ph = photon_lists[n];
    while (ph != NULL) {
      if (ph->is_tracked == 1) {
        get_X_K_interp(ph, t, P, X, Kcov, Kcon);

        io_track[nio].X1     = X[1];
        io_track[nio].X2     = X[2];
        io_track[nio].X3     = X[3];
        io_track[nio].nscatt = ph->nscatt;
        nio++;
      }
      ph = ph->next;
    }
  }

  if (Nph_to_track > 0.) {

    char name[STRLEN];
    char fname[STRLEN];
    sprintf(fname, "trackph_%08d.h5", track_id);
    strcpy(name, dumpdir);
    strcat(name, fname);

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    // One dataset per processor; each processor must create all datasets
    hid_t *ph_dsets = safe_malloc((size_t)mpi_nprocs() * sizeof(hid_t));
    for (int n = 0; n < mpi_nprocs(); n++) {
      char dsetnam[STRLEN];
      sprintf(dsetnam, "trackph_%08d", n);
      int num_tracked_buf = num_tracked;
      mpi_int_broadcast_proc(&num_tracked_buf, n);
      hsize_t dims[1] = {num_tracked_buf};
      hid_t   space   = H5Screate_simple(1, dims, NULL);
      ph_dsets[n]     = H5Dcreate(file_id, dsetnam, trackphfiletype, space,
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Sclose(space);
    }
    if (num_tracked > 0)
      H5Dwrite(ph_dsets[mpi_myrank()], trackphmemtype, H5S_ALL, H5S_ALL,
          H5P_DEFAULT, io_track);

    // Close and release resources
    free(io_track);
    for (int n = 0; n < mpi_nprocs(); n++) {
      H5Dclose(ph_dsets[n]);
    }
    free(ph_dsets);
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    track_id++;
  }
#endif // TRACK_PH
}
#endif // RADIATION

void dump_grid() {
  char name[STRLEN], fname[STRLEN];
  sprintf(fname, "grid.h5");
  strcpy(name, dumpdir);
  strcat(name, fname);
  if (mpi_io_proc()) {
    fprintf(stdout, "WRITING GEOMETRY TO %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Could not create grid file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  {
    double *Xharm = safe_malloc(N1 * N2 * N3 * (size_t)NDIM * sizeof(double));
    int     n     = 0;
    ZLOOP {
      coord(i, j, k, CENT, &Xharm[n]);
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xharm, TYPE_DBL);
    free(Xharm);
  }

  { // If metric == Minkowski, then Xcart == Xharm
    double *Xcart = safe_malloc(N1 * N2 * N3 * (size_t)NDIM * sizeof(double));
    double  X[NDIM], Xcart_loc[NDIM];
    int     n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      cart_coord(X, Xcart_loc);
      /*
      #if METRIC == MKS
      { // poles should be at theta = 0 and pi
  // TODO: is this right? Maybe causes a weird offset.
  if (global_start[2] == 0 && j == NG) Xcart_loc[2] = 0;
  if (global_stop[2] == N2TOT && j == N2 + NG) Xcart_loc[2] = M_PI;
      }
      #endif // METRIC
      */
      for (int l = 0; l < NDIM; l++)
        Xcart[n + l] = Xcart_loc[l];
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xcart, TYPE_DBL);
    free(Xcart);
  }

  { // Face locations, in HARM and Cartesoan coordinates
#define RANK (4)
    // metadata
    hsize_t fdims[RANK]  = {N1TOT + 1, N2TOT + 1, N3TOT + 1, NDIM};
    hsize_t fstart[RANK] = {
        global_start[1], global_start[2], global_start[3], 0};
    hsize_t fcount[RANK] = {N1, N2, N3, NDIM};
    hsize_t mdims[RANK]  = {N1, N2, N3, NDIM};
    for (int d = 0; d < 3; d++) {
      if (global_stop[d + 1] == fdims[d] - 1) {
        fcount[d]++;
        mdims[d]++;
      }
    }
    hsize_t mstart[RANK] = {0, 0, 0, 0};
    hsize_t memsize_tot  = product_hsize_t(mdims, RANK);
    // Values
    double *XFharm = safe_malloc(memsize_tot * sizeof(double));
    double *XFcart = safe_malloc(memsize_tot * sizeof(double));
    double  X[NDIM], Xcart_loc[NDIM];
    int     n = 0;
    ZSLOOP(0, mdims[0] - 1, 0, mdims[1] - 1, 0, mdims[2] - 1) {
      coord(i, j, k, CORN, X);
      cart_coord(X, Xcart_loc);
      for (int l = 0; l < NDIM; l++) {
        XFharm[n + l] = X[l];
        XFcart[n + l] = Xcart_loc[l];
      }
      n += 4;
    }
    WRITE_ARRAY(XFharm, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    WRITE_ARRAY(XFcart, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    free(XFharm);
    free(XFcart);
#undef RANK
  }

#if METRIC == MKS
  {
    double *Xbl = safe_malloc((size_t)N1 * N2 * N3 * NDIM * sizeof(double));
    int     n   = 0;
    ZLOOP {
      double X[NDIM], r, th;
      coord(i, j, k, CENT, X);
      bl_coord(X, &r, &th);
      Xbl[n + 1] = r;
      Xbl[n + 2] = th;
      Xbl[n + 3] = X[3];
      n += 4;
    }
    WRITE_VEC_NO_GHOSTS(Xbl, TYPE_DBL);
    free(Xbl);
  }
#endif

  {
    double *gcov =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      DLOOP2 {
        gcov[n] = ggeom[i][j][CENT].gcov[mu][nu];
        n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(gcov, TYPE_DBL);
    free(gcov);
  }

  {
    double *gcon =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      DLOOP2 {
        gcon[n] = ggeom[i][j][CENT].gcon[mu][nu];
        n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(gcon, TYPE_DBL);
    free(gcon);
  }

  {
    double *gdet = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double));
    int     n    = 0;
    ZLOOP {
      gdet[n] = ggeom[i][j][CENT].g;
      n++;
    }
    WRITE_GRID_NO_GHOSTS(gdet, TYPE_DBL);
    free(gdet);
  }

  {
    double *alpha = safe_malloc((size_t)N1 * N2 * N3 * sizeof(double));
    int     n     = 0;
    ZLOOP {
      alpha[n] = ggeom[i][j][CENT].alpha;
      n++;
    }
    WRITE_GRID_NO_GHOSTS(alpha, TYPE_DBL);
    free(alpha);
  }

#if METRIC == MKS
  {
// connection coefficients
#define RANK (6)
    hsize_t fdims[RANK]  = {N1TOT, N2TOT, N3TOT, NDIM, NDIM, NDIM};
    hsize_t fstart[RANK] = {
        global_start[1], global_start[2], global_start[3], 0, 0, 0};
    hsize_t fcount[RANK] = {N1, N2, N3, NDIM, NDIM, NDIM};
    hsize_t mdims[RANK]  = {N1, N2, N3, NDIM, NDIM, NDIM};
    hsize_t mstart[RANK] = {0, 0, 0, 0, 0, 0};
    double *Gamma =
        safe_malloc(N1 * N2 * N3 * NDIM * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      DLOOP3 {
        Gamma[n] = conn[i][j][mu][nu][sigma];
        n++;
      }
    }
    WRITE_ARRAY(Gamma, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }
#endif // METRIC == MKS

#if METRIC == MKS
  {
    double  X[NDIM];
    double  Lambda_con_local[NDIM][NDIM];
    double  Lambda_cov_local[NDIM][NDIM];
    double *Lambda_h2bl_con =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    double *Lambda_h2bl_cov =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_harm_to_bl(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
        Lambda_h2bl_con[n] = Lambda_con_local[mu][nu];
        Lambda_h2bl_cov[n] = Lambda_cov_local[mu][nu];
        n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2bl_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2bl_cov, TYPE_DBL);
    free(Lambda_h2bl_con);
    free(Lambda_h2bl_cov);
  }

  {
    double  X[NDIM];
    double  Lambda_con_local[NDIM][NDIM];
    double  Lambda_cov_local[NDIM][NDIM];
    double *Lambda_bl2cart_con =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    double *Lambda_bl2cart_cov =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_bl_to_cart(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
        Lambda_bl2cart_con[n] = Lambda_con_local[mu][nu];
        Lambda_bl2cart_cov[n] = Lambda_cov_local[mu][nu];
        n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_bl2cart_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_bl2cart_cov, TYPE_DBL);
    free(Lambda_bl2cart_con);
    free(Lambda_bl2cart_cov);
  }
#endif // METRIC == MKS

  {
    double  X[NDIM];
    double  Lambda_con_local[NDIM][NDIM];
    double  Lambda_cov_local[NDIM][NDIM];
    double *Lambda_h2cart_con =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    double *Lambda_h2cart_cov =
        safe_malloc((size_t)N1 * N2 * N3 * NDIM * NDIM * sizeof(double));
    int n = 0;
    ZLOOP {
      coord(i, j, k, CENT, X);
      jac_harm_to_cart(X, Lambda_cov_local, Lambda_con_local);
      DLOOP2 {
        Lambda_h2cart_con[n] = Lambda_con_local[mu][nu];
        Lambda_h2cart_cov[n] = Lambda_cov_local[mu][nu];
        n++;
      }
    }
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2cart_con, TYPE_DBL);
    WRITE_TENSOR_NO_GHOSTS(Lambda_h2cart_cov, TYPE_DBL);
    free(Lambda_h2cart_con);
    free(Lambda_h2cart_cov);
  }

#if RADIATION
#if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS
  {
    double *local_angles_Xharm = safe_malloc(
        (NDIM - 1) * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 * sizeof(double));
#if METRIC == MKS
    double *local_angles_Xbl = safe_malloc(
        (NDIM - 1) * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 * sizeof(double));
#endif // MKS
    double *local_angles_Xcart = safe_malloc(
        (NDIM - 1) * LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2 * sizeof(double));
    int    n = 0;
    double X[NDIM], Xcart[NDIM], r, th;
    for (int i = 0; i < LOCAL_ANGLES_NX1; ++i) {
      X[1] = startx_rad[1] + (i + 0.5) * local_dx1_rad; // startx is a face
      for (int j = 0; j < LOCAL_ANGLES_NX2; ++j) {
        X[2] = startx_rad[2] + (j + 0.5) * local_dx2_rad;
        cart_coord(X, Xcart);

        local_angles_Xharm[n + 0] = 0;
        local_angles_Xharm[n + 1] = X[1];
        local_angles_Xharm[n + 2] = X[2];

#if METRIC == MKS
        bl_coord(X, &r, &th);
        local_angles_Xbl[n + 0] = 0;
        local_angles_Xbl[n + 1] = r;
        local_angles_Xbl[n + 2] = th;
#endif // MKS

        local_angles_Xcart[n + 0] = 0;
        local_angles_Xcart[n + 1] = Xcart[1];
        local_angles_Xcart[n + 2] = Xcart[2];

        n += (NDIM - 1);
      }
    }
#define RANK (3)
    hsize_t fdims[RANK]  = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2, NDIM - 1};
    hsize_t fstart[RANK] = {0, 0, 0};
    hsize_t fcount[RANK] = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2, NDIM - 1};
    hsize_t mdims[RANK]  = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2, NDIM - 1};
    hsize_t mstart[RANK] = {0, 0, 0};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
      fcount[1] = 0;
      fcount[2] = 0;
    }
    WRITE_ARRAY(local_angles_Xharm, RANK, fdims, fstart, fcount, mdims, mstart,
        TYPE_DBL);
    free(local_angles_Xharm);

    WRITE_ARRAY(local_angles_Xcart, RANK, fdims, fstart, fcount, mdims, mstart,
        TYPE_DBL);
    free(local_angles_Xcart);

#if METRIC == MKS
    WRITE_ARRAY(
        local_angles_Xbl, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    free(local_angles_Xbl);
#endif // MKS
#undef RANK
  }

  {
    double *local_angles_mu = safe_malloc(LOCAL_ANGLES_NMU * sizeof(double));
    for (int i = 0; i < LOCAL_ANGLES_NMU; ++i) {
      local_angles_mu[i] = -1 + (i + 0.5) * local_dx_costh;
    }
#define RANK (1)
    hsize_t fdims[RANK]  = {LOCAL_ANGLES_NMU};
    hsize_t fstart[RANK] = {0};
    hsize_t fcount[RANK] = {LOCAL_ANGLES_NMU};
    hsize_t mdims[RANK]  = {LOCAL_ANGLES_NMU};
    hsize_t mstart[RANK] = {0};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
    }
    WRITE_ARRAY(
        local_angles_mu, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    free(local_angles_mu);
  }
#endif // RADIATION
#endif // LOCAL_ANGULAR_DISTRIBUTIONS

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
}

void dump() {
  timer_start(TIMER_OUT);

  char name[STRLEN], fname[STRLEN];
  sprintf(fname, "dump_%08d.h5", dump_id);
  strcpy(name, dumpdir);
  strcat(name, fname);
  if (mpi_io_proc()) {
    fprintf(stdout, "DUMP %s\n", name);
  }
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Could not create dump file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  WRITE_HDR(version, TYPE_STR);
  WRITE_HDR(metric, TYPE_STR);
  WRITE_HDR(tracers, TYPE_INT);
  WRITE_HDR(eos, TYPE_STR);
  int gamma_fallback = GAMMA_FALLBACK;
  WRITE_HDR(gamma_fallback, TYPE_INT);

  WRITE_HDR(nulnutype, TYPE_STR);

#if METRIC == MKS
  int derefine_poles = DEREFINE_POLES;
  WRITE_HDR(derefine_poles, TYPE_INT);
#endif

  int full_dump = (dump_id % DTf == 0);
  WRITE_HDR(full_dump, TYPE_INT);
  int electrons = ELECTRONS;
  WRITE_HDR(electrons, TYPE_INT);
  int radiation = RADIATION;
  WRITE_HDR(radiation, TYPE_INT);
  int nvar = NVAR;
  WRITE_HDR(nvar, TYPE_INT);
  int nvar_passive = NVAR_PASSIVE;
  WRITE_HDR(nvar_passive, TYPE_INT);
  int output_eosvars = OUTPUT_EOSVARS;
  WRITE_HDR(output_eosvars, TYPE_INT);

#if DIAGNOSTICS_USE_RADTYPES
  int diagnostics_use_radtypes = DIAGNOSTICS_USE_RADTYPES;
#else
  int diagnostics_use_radtypes = 0;
#endif
  WRITE_HDR(diagnostics_use_radtypes, TYPE_INT);

  WRITE_HDR(t, TYPE_DBL);
  WRITE_HDR(tf, TYPE_DBL);
  WRITE_HDR(nstep, TYPE_INT);
  int N1tot = N1TOT;
  WRITE_HDR(N1tot, TYPE_INT);
  int N2tot = N2TOT;
  WRITE_HDR(N2tot, TYPE_INT);
  int N3tot = N3TOT;
  WRITE_HDR(N3tot, TYPE_INT);
  WRITE_HDR(startx[1], TYPE_DBL);
  WRITE_HDR(startx[2], TYPE_DBL);
  WRITE_HDR(startx[3], TYPE_DBL);
  WRITE_HDR(dx[1], TYPE_DBL);
  WRITE_HDR(dx[2], TYPE_DBL);
  WRITE_HDR(dx[3], TYPE_DBL);

#if (METRIC == MKS)
  WRITE_HDR(Rin, TYPE_DBL);
  WRITE_HDR(Rout, TYPE_DBL);
  WRITE_HDR(Rout_vis, TYPE_DBL);
  WRITE_HDR(Reh, TYPE_DBL);
  WRITE_HDR(Risco, TYPE_DBL);
  WRITE_HDR(hslope, TYPE_DBL);
  WRITE_HDR(a, TYPE_DBL);
  WRITE_HDR(poly_xt, TYPE_DBL);
  WRITE_HDR(poly_alpha, TYPE_DBL);
  WRITE_HDR(mks_smooth, TYPE_DBL);
#if NEED_UNITS
  WRITE_HDR(Mbh, TYPE_DBL);
  WRITE_HDR(mbh, TYPE_DBL);
#endif // NEED_UNITS
#if RADIATION
  WRITE_HDR(Rout_rad, TYPE_DBL);
#endif // RADIATION
#endif // METRIC

#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  WRITE_HDR(gam, TYPE_DBL);
#endif // GAMMA
#if EOS == EOS_TYPE_POLYTROPE
  WRITE_HDR(poly_K, TYPE_DBL);
  WRITE_HDR(poly_gam, TYPE_DBL);
#endif // POLYTROPE
#if EOS == EOS_TYPE_TABLE
  WRITE_HDR(eospath, TYPE_STR);
#endif // TABLE

#if ELECTRONS
  WRITE_HDR(game, TYPE_DBL);
  WRITE_HDR(gamp, TYPE_DBL);
#endif

#if NEED_UNITS
  WRITE_HDR(L_unit, TYPE_DBL);
  WRITE_HDR(M_unit, TYPE_DBL);
  WRITE_HDR(RHO_unit, TYPE_DBL);
  WRITE_HDR(T_unit, TYPE_DBL);
  WRITE_HDR(U_unit, TYPE_DBL);
  WRITE_HDR(B_unit, TYPE_DBL);
#if EOS == EOS_TYPE_TABLE
  WRITE_HDR(TEMP_unit, TYPE_DBL);
#endif // EOS_TYPE_TABLE
#if RADIATION
  WRITE_HDR(tp_over_te, TYPE_DBL);
  WRITE_HDR(Ne_unit, TYPE_DBL);
  WRITE_HDR(Thetae_unit, TYPE_DBL);
#endif // RADIATION
#endif

  WRITE_HDR(cour, TYPE_DBL);
  WRITE_HDR(DTd, TYPE_DBL);
  WRITE_HDR(DTl, TYPE_DBL);
  WRITE_HDR(DTr, TYPE_DBL);
  WRITE_HDR(DNr, TYPE_INT);
  WRITE_HDR(DTp, TYPE_INT);
  WRITE_HDR(DTf, TYPE_INT);
  WRITE_HDR(dump_cnt, TYPE_INT);
  WRITE_HDR(dump_id, TYPE_INT);
  WRITE_HDR(dt, TYPE_DBL);
  WRITE_HDR(failed, TYPE_INT);
#if RADIATION
  int maxnscatt = MAXNSCATT;
  WRITE_HDR(maxnscatt, TYPE_INT);
  int nubins = NU_BINS;
  WRITE_HDR(nubins, TYPE_INT);
  int nth = NTH;
  WRITE_HDR(nth, TYPE_INT);
  int nphi = NPHI;
  WRITE_HDR(nphi, TYPE_INT);
  int nubins_spec = NU_BINS_SPEC;
  WRITE_HDR(nubins_spec, TYPE_INT);
  WRITE_HDR(numin, TYPE_DBL);
  WRITE_HDR(numax, TYPE_DBL);

  int neutrino_oscillations = NEUTRINO_OSCILLATIONS;
  WRITE_HDR(neutrino_oscillations, TYPE_INT);
  int force_equipartition = FORCE_EQUIPARTITION;
  WRITE_HDR(force_equipartition, TYPE_INT);
  int local_angular_distributions = LOCAL_ANGULAR_DISTRIBUTIONS;
  WRITE_HDR(local_angular_distributions, TYPE_INT);
  int rz_histograms = RZ_HISTOGRAMS;
  WRITE_HDR(rz_histograms, TYPE_INT);
#endif

#if NVAR_PASSIVE > 0
  { // passive_type
#define RANK (1)
    hsize_t fdims[RANK]  = {NVAR_PASSIVE};
    hsize_t fstart[RANK] = {0};
    hsize_t fcount[RANK] = {NVAR_PASSIVE};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
    }
    hsize_t mdims[RANK]  = {NVAR_PASSIVE};
    hsize_t mstart[RANK] = {0};
    WRITE_ARRAY(
        passive_type, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
#undef RANK
  }
#endif

  WRITE_PRIM(P, TYPE_FLOAT);

  {
    hid_t prim_plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    hid_t prim_dset     = H5Dopen(file_id, "P", prim_plist_id);
    hid_t strtype       = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, H5T_VARIABLE);
    hsize_t str_dims[1] = {NVAR};
    hid_t   prim_space  = H5Screate_simple(1, str_dims, NULL);
    hid_t   str_attr    = H5Acreate(
        prim_dset, "vnams", strtype, prim_space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(str_attr, strtype, vnams);
    H5Aclose(str_attr);
    H5Sclose(prim_space);
    H5Dclose(prim_dset);
    H5Pclose(prim_plist_id);
  }

  if (full_dump) {
#pragma omp parallel for collapse(3)
    ZLOOP divb[i][j][k] = flux_ct_divb(i, j, k);
    WRITE_GRID(divb, TYPE_FLOAT);

#if OUTPUT_EOSVARS
    {
#pragma omp parallel for collapse(3)
      ZLOOP {
#if EOS == EOS_TYPE_TABLE
        EOS_SC_fill(P[i][j][k], extra[i][j][k]);
#endif // EOS_TYPE_TABLE

        PRESS[i][j][k] = EOS_pressure_rho0_u(
            P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
        ENT[i][j][k] =
            EOS_entropy_rho0_u(P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
        TEMP[i][j][k] =
            EOS_temperature(P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
        CS2[i][j][k] = EOS_sound_speed_rho0_u(
            P[i][j][k][RHO], P[i][j][k][UU], extra[i][j][k]);
      }
      WRITE_GRID(PRESS, TYPE_FLOAT);
      WRITE_GRID(ENT, TYPE_FLOAT);
      WRITE_GRID(TEMP, TYPE_FLOAT);
      WRITE_GRID(CS2, TYPE_FLOAT);
    }
#endif

#if ELECTRONS
    WRITE_GRID(Qvisc, TYPE_FLOAT);
#if RADIATION
    WRITE_GRID(Qcoul, TYPE_FLOAT);
#endif
#endif

    WRITE_GRID(fail_save, TYPE_INT);

    current_calc();
    WRITE_VEC(jcon, TYPE_FLOAT);

#if RADIATION
    WRITE_GRID(Nsph, TYPE_INT);

    WRITE_GRID(nph, TYPE_FLOAT);

    {
#define RANK (4)
      hsize_t fdims[RANK]  = {MAXNSCATT + 2, N1TOT, N2TOT, N3TOT};
      hsize_t fstart[RANK] = {
          0, global_start[1], global_start[2], global_start[3]};
      hsize_t fcount[RANK] = {MAXNSCATT + 2, N1, N2, N3};
      hsize_t mdims[RANK]  = {
          MAXNSCATT + 2, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
      hsize_t mstart[RANK] = {0, NG, NG, NG};
      WRITE_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_FLOAT);
#undef RANK
    }

    {
      ZLOOP {
        dtau_scatt[i][j][k]        = 0.0;
        dtau_tot[i][j][k]          = dtau_avg[0][i][j][k];
        denominator_scatt[i][j][k] = SMALL;
        denominator_tot[i][j][k]   = en_int_avg[0][i][j][k] + SMALL;
        SCATTLOOP {
          dtau_scatt[i][j][k] += dtau_avg[iscatt + 1][i][j][k];
          dtau_tot[i][j][k] += dtau_avg[iscatt + 1][i][j][k];

          denominator_scatt[i][j][k] += en_int_avg[iscatt + 1][i][j][k];
          denominator_tot[i][j][k] += en_int_avg[iscatt + 1][i][j][k];
        }
        denominator_scatt[i][j][k] *= DTd;
        denominator_tot[i][j][k] *= DTd;
        dtau_scatt[i][j][k] /= denominator_scatt[i][j][k];
        dtau_tot[i][j][k] /= denominator_tot[i][j][k];
      }
      WRITE_GRID(dtau_scatt, TYPE_DBL);
      WRITE_GRID(dtau_tot, TYPE_DBL);
    }

    {
      for (int iscatt = 0; iscatt < RAD_SCATT_TYPES + 1; iscatt++) {
        ZLOOP {
          dtau_avg[iscatt][i][j][k] /=
              ((SMALL + en_int_avg[iscatt][i][j][k]) * DTd);
        }
      }
#define RANK (4)
      hsize_t fdims[RANK]  = {RAD_SCATT_TYPES + 1, N1TOT, N2TOT, N3TOT};
      hsize_t fstart[RANK] = {
          0, global_start[1], global_start[2], global_start[3]};
      hsize_t fcount[RANK] = {RAD_SCATT_TYPES + 1, N1, N2, N3};
      hsize_t mdims[RANK]  = {
          RAD_SCATT_TYPES + 1, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
      hsize_t mstart[RANK] = {0, NG, NG, NG};
      WRITE_ARRAY(
          dtau_avg, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    }

    WRITE_GRID(Nem, TYPE_INT);
    WRITE_GRID(Nabs, TYPE_INT);
    WRITE_GRID(Nsc, TYPE_INT);

    set_cooling_time(tau_cool, P, extra);
    WRITE_GRID(tau_cool, TYPE_FLOAT);

    set_Rmunu();
    WRITE_TENSOR(Rmunu, TYPE_FLOAT);

    ZLOOP {
      TYPELOOP {
        Nem_phys[i][j][k][itp] /= DTd;
        Nabs_phys[i][j][k][itp] /= DTd;
      }
    }
    WRITE_RADTYPEVEC(Nem_phys, TYPE_DBL);
    WRITE_RADTYPEVEC(Nabs_phys, TYPE_DBL);

    ZLOOP {
      DLOOP1 { radG_int[i][j][k][mu] /= DTd; }
#if RADIATION == RADTYPE_NEUTRINOS
      radG_int[i][j][k][RADG_YE] /= DTd;
      radG_int[i][j][k][RADG_YE_EM] /= DTd;
#endif
    }

    WRITE_RADG(radG_int, TYPE_FLOAT);

    // nuLnu reduced over MPI, now just io proc outputs it
    {
#if METRIC == MINKOWSKI
      bin_all_superphotons();
#endif
      mpi_reduce_nuLnu();
#define RANK (4)
      hsize_t fdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
      hsize_t fstart[RANK] = {0, 0, 0, 0};
      hsize_t fcount[RANK] = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
      hsize_t mstart[RANK] = {0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
      }
      WRITE_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    }

    // local angle histograms
#if RADIATION == RADTYPE_NEUTRINOS && LOCAL_ANGULAR_DISTRIBUTIONS
    accumulate_local_angles();
    {
#define RANK (5)
      hsize_t fdims[RANK]  = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, RAD_NUM_TYPES, LOCAL_ANGLES_NMU};
      hsize_t fstart[RANK] = {0, 0, 0, 0, 0};
      hsize_t fcount[RANK] = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, RAD_NUM_TYPES, LOCAL_ANGLES_NMU};
      hsize_t mdims[RANK]  = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, RAD_NUM_TYPES, LOCAL_ANGLES_NMU};
      hsize_t mstart[RANK] = {0, 0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
        fcount[4] = 0;
      }
      WRITE_ARRAY(
          local_angles, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    }

    {
#define RANK (4)
      hsize_t fdims[RANK]  = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, LOCAL_ANGLES_NMU};
      hsize_t fstart[RANK] = {0, 0, 0, 0};
      hsize_t fcount[RANK] = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, LOCAL_ANGLES_NMU};
      hsize_t mdims[RANK]  = {LOCAL_NUM_BASES, LOCAL_ANGLES_NX1,
          LOCAL_ANGLES_NX2, LOCAL_ANGLES_NMU};
      hsize_t mstart[RANK] = {0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
      }
      WRITE_ARRAY(Gnu, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
      WRITE_ARRAY(
          local_Ns, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
      WRITE_ARRAY(
          local_wsqr, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    }

    {
#define RANK (4)
      hsize_t fdims[RANK] = {
          LOCAL_NUM_BASES, 2, LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t fstart[RANK] = {0, 0, 0, 0};
      hsize_t fcount[RANK] = {
          LOCAL_NUM_BASES, 2, LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t mdims[RANK] = {
          LOCAL_NUM_BASES, 2, LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t mstart[RANK] = {0, 0, 0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
        fcount[2] = 0;
        fcount[3] = 0;
      }
      WRITE_ARRAY(
          local_moments, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
    }

#if NEUTRINO_OSCILLATIONS
    {
      mpi_dbl_allreduce_array(
          (double *)local_osc_count, LOCAL_ANGLES_NX1 * LOCAL_ANGLES_NX2);
      for (int ix1 = 0; ix1 < LOCAL_ANGLES_NX1; ++ix1) {
        for (int ix2 = 0; ix2 < LOCAL_ANGLES_NX2; ++ix2) {
          local_osc_count[ix1][ix2] /= DTd;
        }
      }

#define RANK (2)
      hsize_t fdims[RANK]  = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t fstart[RANK] = {0, 0};
      hsize_t fcount[RANK] = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t mdims[RANK]  = {LOCAL_ANGLES_NX1, LOCAL_ANGLES_NX2};
      hsize_t mstart[RANK] = {0, 0};
      if (!mpi_io_proc()) {
        fcount[0] = 0;
        fcount[1] = 0;
      }
      WRITE_ARRAY(local_osc_count, RANK, fdims, fstart, fcount, mdims, mstart,
          TYPE_DBL);
#undef RANK
    }
#endif // NEUTRINO_OSCILLATIONS
#endif // LOCAL_ANGULAR_DISTRIBUTIONS

#if RZ_HISTOGRAMS
  generate_rz_histograms();
  {
#define RANK (1)
    hsize_t fdims[RANK]  = {RZ_HISTOGRAMS_N};
    hsize_t fstart[RANK] = {0};
    hsize_t fcount[RANK] = {RZ_HISTOGRAMS_N};
    hsize_t mdims[RANK]  = {RZ_HISTOGRAMS_N};
    hsize_t mstart[RANK] = {0};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
    }
    WRITE_ARRAY(rz_r_orig_hist, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    WRITE_ARRAY(rz_z_orig_hist, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#if NEUTRINO_OSCILLATIONS
    WRITE_ARRAY(osc_rz_r_orig_hist, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    WRITE_ARRAY(osc_rz_z_orig_hist, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#endif
#undef RANK
  }
#endif // RZ_HISTOGRAMS
#endif // RADIATION
  }

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  if (mpi_io_proc()) {
    write_xml_file(dump_id, t, vnams);
  }

  dump_id++;
  fdump_id++;

  timer_stop(TIMER_OUT);
}

void restart_write(int restart_type) {
  timer_start(TIMER_OUT);

  char name[STRLEN], fname[STRLEN];
  if (restart_type == RESTART_TEMP) {
    sprintf(fname, "restart_%08d.h5", restart_id);
  } else if (restart_type == RESTART_PERM) {
    sprintf(fname, "restart_perm_%08d.h5", restart_perm_id);
    restart_perm_id++;
  }
  strcpy(name, restartdir);
  strcat(name, fname);

  if (restart_type == RESTART_TEMP) {
    char lastname[STRLEN];
    sprintf(fname, "restart.last");
    strcpy(lastname, restartdir);
    strcat(lastname, fname);

    restart_id++;
    FILE *fp = fopen(lastname, "w");
    fprintf(fp, "%s\n", name);
    fclose(fp);
  }

  if (mpi_io_proc()) {
    fprintf(stdout, "RESTART %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    printf("Could not create restart file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  WRITE_HDR(version, TYPE_STR);
  WRITE_HDR(metric, TYPE_STR);
  WRITE_HDR(nulnutype, TYPE_STR);
  int electrons = ELECTRONS;
  WRITE_HDR(electrons, TYPE_INT);
  int radiation = RADIATION;
  WRITE_HDR(radiation, TYPE_INT);
  int nvar = NVAR;
  WRITE_HDR(nvar, TYPE_INT);
  int nvar_passive = NVAR_PASSIVE;
  WRITE_HDR(nvar_passive, TYPE_INT);
  int output_eosvars = OUTPUT_EOSVARS;
  WRITE_HDR(output_eosvars, TYPE_INT);

#if DIAGNOSTICS_USE_RADTYPES
  int diagnostics_use_radtypes = DIAGNOSTICS_USE_RADTYPES;
#else
  int diagnostics_use_radtypes = 0;
#endif
  WRITE_HDR(diagnostics_use_radtypes, TYPE_INT);

  WRITE_HDR(t, TYPE_DBL);
  WRITE_HDR(tf, TYPE_DBL);
  WRITE_HDR(nstep, TYPE_INT);
  WRITE_HDR(dump_id, TYPE_DBL);
  WRITE_HDR(fdump_id, TYPE_DBL);
#if RADIATION && TRACERS
  WRITE_HDR(dumptrace_id, TYPE_INT);
#endif
  WRITE_HDR(restart_id, TYPE_INT);
  WRITE_HDR(restart_perm_id, TYPE_INT);
  WRITE_HDR(dt, TYPE_DBL);
  WRITE_HDR(lim, TYPE_INT);
  WRITE_HDR(failed, TYPE_INT);

  WRITE_HDR(DTd, TYPE_DBL);
  WRITE_HDR(DTl, TYPE_DBL);
  WRITE_HDR(DTr, TYPE_DBL);
  WRITE_HDR(DNr, TYPE_INT);
  WRITE_HDR(DTp, TYPE_INT);
  WRITE_HDR(DTf, TYPE_INT);
  WRITE_HDR(tdump, TYPE_DBL);
  WRITE_HDR(trestart, TYPE_DBL);
  WRITE_HDR(tlog, TYPE_DBL);

#if METRIC == MINKOWSKI
  WRITE_HDR(x1Min, TYPE_DBL);
  WRITE_HDR(x1Max, TYPE_DBL);
  WRITE_HDR(x2Min, TYPE_DBL);
  WRITE_HDR(x2Max, TYPE_DBL);
  WRITE_HDR(x3Min, TYPE_DBL);
  WRITE_HDR(x3Max, TYPE_DBL);
#elif METRIC == MKS
  WRITE_HDR(a, TYPE_DBL);
  WRITE_HDR(Rin, TYPE_DBL);
  WRITE_HDR(Rout, TYPE_DBL);
  WRITE_HDR(Rout_vis, TYPE_DBL);
  WRITE_HDR(hslope, TYPE_DBL);
  WRITE_HDR(Reh, TYPE_DBL);
  WRITE_HDR(Risco, TYPE_DBL);
  WRITE_HDR(poly_xt, TYPE_DBL);
  WRITE_HDR(poly_alpha, TYPE_DBL);
  WRITE_HDR(mks_smooth, TYPE_DBL);
#if RADIATION
  WRITE_HDR(mbh, TYPE_DBL);
  WRITE_HDR(Mbh, TYPE_DBL);
#endif // RADIATION
#endif // METRIC

// EOS
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  WRITE_HDR(gam, TYPE_DBL);
#endif // GAMMA EOS
#if EOS == EOS_TYPE_POLYTROPE
  WRITE_HDR(poly_K, TYPE_DBL);
  WRITE_HDR(poly_gam, TYPE_DBL);
#elif EOS == EOS_TYPE_TABLE
  WRITE_HDR(eospath, TYPE_STR);
#endif // EOS

// Units
#if NEED_UNITS
  WRITE_HDR(L_unit, TYPE_DBL);
  WRITE_HDR(M_unit, TYPE_DBL);
  WRITE_HDR(RHO_unit, TYPE_DBL);
  WRITE_HDR(T_unit, TYPE_DBL);
  WRITE_HDR(U_unit, TYPE_DBL);
  WRITE_HDR(B_unit, TYPE_DBL);
#if EOS == EOS_TYPE_TABLE
  WRITE_HDR(TEMP_unit, TYPE_DBL);
#endif // EOS_TYPE_TABLE
#if RADIATION
  WRITE_HDR(tp_over_te, TYPE_DBL);
  WRITE_HDR(Ne_unit, TYPE_DBL);
  WRITE_HDR(Thetae_unit, TYPE_DBL);
#endif // RADIATION
#endif

  // Fluid
  WRITE_HDR(cour, TYPE_DBL);

#if ELECTRONS
  WRITE_HDR(game, TYPE_DBL);
  WRITE_HDR(gamp, TYPE_DBL);
  WRITE_HDR(fel0, TYPE_DBL);
  WRITE_HDR(tptemin, TYPE_DBL);
  WRITE_HDR(tptemax, TYPE_DBL);
#endif // ELECTRONS

#if RADIATION
  WRITE_HDR(numin, TYPE_DBL);
  WRITE_HDR(numax, TYPE_DBL);
  WRITE_HDR(tune_emiss, TYPE_DBL);
  WRITE_HDR(tune_scatt, TYPE_DBL);
  WRITE_HDR(t_tune_emiss, TYPE_DBL);
  WRITE_HDR(t_tune_scatt, TYPE_DBL);
  WRITE_HDR(dt_tune_emiss, TYPE_DBL);
  WRITE_HDR(dt_tune_scatt, TYPE_DBL);
  WRITE_HDR(made_tune_proc, TYPE_INT);
  WRITE_HDR(abs_tune_proc, TYPE_INT);
  WRITE_HDR(scatt_tune_proc, TYPE_INT);
  WRITE_HDR(thetae_max, TYPE_DBL);
  WRITE_HDR(sigma_max, TYPE_DBL);
  WRITE_HDR(kdotk_tol, TYPE_DBL);
  WRITE_HDR(Nph_to_track, TYPE_DBL);
#if METRIC == MKS
  WRITE_HDR(Rout_rad, TYPE_DBL);
#endif
#if RADIATION == RADTYPE_NEUTRINOS
  WRITE_HDR(lepton_tot, TYPE_DBL);
  WRITE_HDR(lepton_last, TYPE_DBL);
  WRITE_HDR(lepton_lost, TYPE_DBL);
#endif // RADIATION
#endif

#if NVAR_PASSIVE > 0
  // WRITE_ARRAY(passive_type, TYPE_INT);
  {
#define RANK (1)
    hsize_t fdims[RANK]  = {NVAR_PASSIVE};
    hsize_t fstart[RANK] = {0};
    hsize_t fcount[RANK] = {NVAR_PASSIVE};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
    }
    hsize_t mdims[RANK]  = {NVAR_PASSIVE};
    hsize_t mstart[RANK] = {0};
    WRITE_ARRAY(
        passive_type, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
#undef RANK
  }
#endif

  WRITE_PRIM(P, TYPE_DBL);

#if RADIATION
  WRITE_GRID(Nsph, TYPE_INT);
  WRITE_GRID(nph, TYPE_DBL);
  WRITE_GRID(Nem, TYPE_INT);
  WRITE_GRID(Nabs, TYPE_INT);
  WRITE_GRID(Nsc, TYPE_INT);

  WRITE_RADTYPEVEC(Nem_phys, TYPE_DBL);
  WRITE_RADTYPEVEC(Nabs_phys, TYPE_DBL);

  {
#define RANK (4)
    hsize_t fdims[RANK]  = {MAXNSCATT + 2, N1TOT, N2TOT, N3TOT};
    hsize_t fstart[RANK] = {
        0, global_start[1], global_start[2], global_start[3]};
    hsize_t fcount[RANK] = {MAXNSCATT + 2, N1, N2, N3};
    hsize_t mdims[RANK]  = {
        MAXNSCATT + 2, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
    hsize_t mstart[RANK] = {0, NG, NG, NG};
    WRITE_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }

  {
#define RANK (4)
    hsize_t fdims[RANK]  = {RAD_SCATT_TYPES + 1, N1TOT, N2TOT, N3TOT};
    hsize_t fstart[RANK] = {
        0, global_start[1], global_start[2], global_start[3]};
    hsize_t fcount[RANK] = {RAD_SCATT_TYPES + 1, N1, N2, N3};
    hsize_t mdims[RANK]  = {
        RAD_SCATT_TYPES + 1, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
    hsize_t mstart[RANK] = {0, NG, NG, NG};
    WRITE_ARRAY(dtau_avg, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    WRITE_ARRAY(
        en_int_avg, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }

  {
    mpi_reduce_nuLnu();
#define RANK (4)
    hsize_t fdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t fstart[RANK] = {0, 0, 0, 0};
    hsize_t fcount[RANK] = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t mdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t mstart[RANK] = {0, 0, 0, 0};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
      fcount[1] = 0;
      fcount[2] = 0;
      fcount[3] = 0;
    }
    WRITE_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_FLOAT);
#undef RANK
  }

  // Superphoton data
  // compute sizes
  size_t npart_local  = count_particles_local();
  size_t npart        = mpi_reduce_int(npart_local);
  size_t npart_offset = mpi_accumulate_int(npart_local);

  // allocate
  struct of_photon *superphotons;
  superphotons = safe_malloc((size_t)(npart_local) * sizeof(struct of_photon));

  // Copy superphotons into buffer
  int nph = 0;
  for (int n = 0; n < nthreads; n++) {
    struct of_photon *ph = photon_lists[n];
    while (ph != NULL) {
      copy_photon(ph, &(superphotons[nph]));
      ph = ph->next;
      nph++;
    }
  }

  // write number of particles per rank
  {
#define RANK (1)
    hsize_t NTOT_CPU         = N1CPU * N2CPU * N3CPU;
    hsize_t my_rank          = mpi_myrank();
    hsize_t fdims[RANK]      = {NTOT_CPU};
    hsize_t fstart[RANK]     = {my_rank};
    hsize_t fcount[RANK]     = {1};
    hsize_t mdims[RANK]      = {1};
    hsize_t mstart[RANK]     = {0};
    int *   particle_offsets = safe_malloc(1 * sizeof(int));
    int *   particle_counts  = safe_malloc(1 * sizeof(int));
    particle_offsets[0]      = npart_offset;
    particle_counts[0]       = npart_local;
    WRITE_ARRAY(
        particle_offsets, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
    WRITE_ARRAY(
        particle_counts, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
    free(particle_offsets);
    free(particle_counts);
#undef RANK
  }

  // write superphoton datasets
  {
#define RANK (1)
    hsize_t fdims[RANK]  = {npart};
    hsize_t fstart[RANK] = {npart_offset};
    hsize_t fcount[RANK] = {npart_local};
    hsize_t mdims[RANK]  = {npart_local};
    hsize_t mstart[RANK] = {0};
    WRITE_ARRAY(
        superphotons, RANK, fdims, fstart, fcount, mdims, mstart, phmemtype);
#undef RANK
  }

  // Close and release resources
  free(superphotons);
#endif // RADIATION

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  // Keep only two most recent restart files
  if (mpi_io_proc() && restart_id >= 3 && restart_type == RESTART_TEMP) {
    char nametodel[STRLEN];
    sprintf(fname, "restart_%08d.h5", restart_id - 3);
    strcpy(nametodel, restartdir);
    strcat(nametodel, fname);
    remove(nametodel);
  }

  timer_stop(TIMER_OUT);
}

void restart_read(char *fname) {
  timer_start(TIMER_OUT);

  if (mpi_io_proc()) {
    fprintf(stderr, "Restarting from %s\n\n", fname);
  }
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "ERROR %s is not a valid restart file\n", fname);
    exit(-1);
  }
  H5Pclose(plist_id);
  READ_HDR(t, TYPE_DBL);
  READ_HDR(tf, TYPE_DBL);
  READ_HDR(nstep, TYPE_INT);
  READ_HDR(dump_id, TYPE_DBL);
  READ_HDR(fdump_id, TYPE_DBL);
#if RADIATION && TRACERS
  READ_HDR(dumptrace_id, TYPE_INT);
#endif
  READ_HDR(restart_id, TYPE_INT);
  READ_HDR(restart_perm_id, TYPE_INT);
  READ_HDR(dt, TYPE_DBL);
  READ_HDR(lim, TYPE_INT);
  READ_HDR(failed, TYPE_INT);

  READ_HDR(DTd, TYPE_DBL);
  READ_HDR(DTl, TYPE_DBL);
  READ_HDR(DTr, TYPE_DBL);
  READ_HDR(DNr, TYPE_INT);
  READ_HDR(DTp, TYPE_INT);
  READ_HDR(DTf, TYPE_INT);
  READ_HDR(tdump, TYPE_DBL);
  READ_HDR(trestart, TYPE_DBL);
  READ_HDR(tlog, TYPE_DBL);

#if METRIC == MINKOWSKI
  READ_HDR(x1Min, TYPE_DBL);
  READ_HDR(x1Max, TYPE_DBL);
  READ_HDR(x2Min, TYPE_DBL);
  READ_HDR(x2Max, TYPE_DBL);
  READ_HDR(x3Min, TYPE_DBL);
  READ_HDR(x3Max, TYPE_DBL);
#elif METRIC == MKS
  READ_HDR(a, TYPE_DBL);
  READ_HDR(Rin, TYPE_DBL);
  READ_HDR(Rout, TYPE_DBL);
  READ_HDR(Rout_vis, TYPE_DBL);
  READ_HDR(hslope, TYPE_DBL);
  READ_HDR(Reh, TYPE_DBL);
  READ_HDR(Risco, TYPE_DBL);
  READ_HDR(poly_xt, TYPE_DBL);
  READ_HDR(poly_alpha, TYPE_DBL);
  READ_HDR(mks_smooth, TYPE_DBL);
#if RADIATION
  READ_HDR(mbh, TYPE_DBL);
  READ_HDR(Mbh, TYPE_DBL);
#endif // RADIATION
#endif // METRIC

// Units
#if NEED_UNITS
  READ_HDR(L_unit, TYPE_DBL);
  READ_HDR(T_unit, TYPE_DBL);
  READ_HDR(M_unit, TYPE_DBL);
  READ_HDR(RHO_unit, TYPE_DBL);
  READ_HDR(U_unit, TYPE_DBL);
  READ_HDR(B_unit, TYPE_DBL);
#if EOS == EOS_TYPE_TABLE
  READ_HDR(TEMP_unit, TYPE_DBL);
#endif // EOS_TYPE_TABLE

#if RADIATION
  READ_HDR(tp_over_te, TYPE_DBL);
  READ_HDR(Ne_unit, TYPE_DBL);
  READ_HDR(Thetae_unit, TYPE_DBL);
#endif // RADIATION
#endif // NEED_UNITS

  READ_HDR(cour, TYPE_DBL);

#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  READ_HDR(gam, TYPE_DBL);
#endif // EOS

#if EOS == EOS_TYPE_POLYTROPE
  READ_HDR(poly_K, TYPE_DBL);
  READ_HDR(poly_gam, TYPE_DBL);
#endif // EOS

#if ELECTRONS
  READ_HDR(game, TYPE_DBL);
  READ_HDR(gamp, TYPE_DBL);
  READ_HDR(fel0, TYPE_DBL);
  READ_HDR(tptemin, TYPE_DBL);
  READ_HDR(tptemax, TYPE_DBL);
#endif // ELECTRONS

#if RADIATION
  READ_HDR(numin, TYPE_DBL);
  READ_HDR(numax, TYPE_DBL);
  READ_HDR(tune_emiss, TYPE_DBL);
  READ_HDR(tune_scatt, TYPE_DBL);
  READ_HDR(t_tune_emiss, TYPE_DBL);
  READ_HDR(t_tune_scatt, TYPE_DBL);
  READ_HDR(dt_tune_emiss, TYPE_DBL);
  READ_HDR(dt_tune_scatt, TYPE_DBL);
  READ_HDR(made_tune_proc, TYPE_INT);
  READ_HDR(abs_tune_proc, TYPE_INT);
  READ_HDR(scatt_tune_proc, TYPE_INT);
  READ_HDR(thetae_max, TYPE_DBL);
  READ_HDR(sigma_max, TYPE_DBL);
  READ_HDR(kdotk_tol, TYPE_DBL);
  READ_HDR(Nph_to_track, TYPE_DBL);
#if METRIC == MKS
  READ_HDR(Rout_rad, TYPE_DBL);
#endif
#if RADIATION == RADTYPE_NEUTRINOS
  READ_HDR(lepton_tot, TYPE_DBL);
  READ_HDR(lepton_last, TYPE_DBL);
  READ_HDR(lepton_lost, TYPE_DBL);
#endif // RADIATION
#endif

#if NVAR_PASSIVE > 0
  {
#define RANK (1)
    hsize_t fdims[RANK]  = {NVAR_PASSIVE};
    hsize_t fstart[RANK] = {0};
    hsize_t fcount[RANK] = {NVAR_PASSIVE};
    hsize_t mdims[RANK]  = {NVAR_PASSIVE};
    hsize_t mstart[RANK] = {0};
    READ_ARRAY(
        passive_type, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
#undef RANK
  }
#endif

  READ_PRIM(P, TYPE_DBL);

#if RADIATION
  READ_GRID(Nsph, TYPE_INT);
  READ_GRID(nph, TYPE_DBL);
  READ_GRID(Nem, TYPE_INT);
  READ_GRID(Nabs, TYPE_INT);
  READ_GRID(Nsc, TYPE_INT);

  READ_RADTYPEVEC(Nem_phys, TYPE_DBL);
  READ_RADTYPEVEC(Nabs_phys, TYPE_DBL);

  {
#define RANK (4)
    hsize_t fdims[RANK]  = {MAXNSCATT + 2, N1TOT, N2TOT, N3TOT};
    hsize_t fstart[RANK] = {
        0, global_start[1], global_start[2], global_start[3]};
    hsize_t fcount[RANK] = {MAXNSCATT + 2, N1, N2, N3};
    hsize_t mdims[RANK]  = {
        MAXNSCATT + 2, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
    hsize_t mstart[RANK] = {0, NG, NG, NG};
    READ_ARRAY(Jrad, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }

  {
#define RANK (4)
    hsize_t fdims[RANK]  = {RAD_SCATT_TYPES + 1, N1TOT, N2TOT, N3TOT};
    hsize_t fstart[RANK] = {
        0, global_start[1], global_start[2], global_start[3]};
    hsize_t fcount[RANK] = {RAD_SCATT_TYPES + 1, N1, N2, N3};
    hsize_t mdims[RANK]  = {
        RAD_SCATT_TYPES + 1, N1 + 2 * NG, N2 + 2 * NG, N3 + 2 * NG};
    hsize_t mstart[RANK] = {0, NG, NG, NG};
    READ_ARRAY(dtau_avg, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
    READ_ARRAY(
        en_int_avg, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }

  {
#define RANK (4)
    hsize_t fdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t fstart[RANK] = {0, 0, 0, 0};
    hsize_t fcount[RANK] = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t mdims[RANK]  = {NULNU_IDX0, NTH, NPHI, NU_BINS_SPEC};
    hsize_t mstart[RANK] = {0, 0, 0, 0};
    if (!mpi_io_proc()) {
      fcount[0] = 0;
      fcount[1] = 0;
      fcount[2] = 0;
      fcount[3] = 0;
    }
    READ_ARRAY(nuLnu, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_DBL);
#undef RANK
  }

  // Superphoton datasets
  // read number of particles per rank
  int npart_offset, npart_local;
  {
#define RANK (1)
    hsize_t NTOT_CPU         = N1CPU * N2CPU * N3CPU;
    hsize_t my_rank          = mpi_myrank();
    hsize_t fdims[RANK]      = {NTOT_CPU};
    hsize_t fstart[RANK]     = {my_rank};
    hsize_t fcount[RANK]     = {1};
    hsize_t mdims[RANK]      = {1};
    hsize_t mstart[RANK]     = {0};
    int *   particle_offsets = safe_malloc(1 * sizeof(int));
    int *   particle_counts  = safe_malloc(1 * sizeof(int));
    READ_ARRAY(
        particle_offsets, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
    READ_ARRAY(
        particle_counts, RANK, fdims, fstart, fcount, mdims, mstart, TYPE_INT);
    npart_offset = particle_offsets[0];
    npart_local  = particle_counts[0];
    free(particle_offsets);
    free(particle_counts);
#undef RANK
  }
  int npart = mpi_reduce_int(npart_local);

  // read full dataset
  struct of_photon *superphotons;
  superphotons = safe_malloc((size_t)(npart_local) * sizeof(struct of_photon));
  {
#define RANK (1)
    hsize_t fdims[RANK]  = {npart};
    hsize_t fstart[RANK] = {npart_offset};
    hsize_t fcount[RANK] = {npart_local};
    hsize_t mdims[RANK]  = {npart_local};
    hsize_t mstart[RANK] = {0};
    READ_ARRAY(
        superphotons, RANK, fdims, fstart, fcount, mdims, mstart, phmemtype);
#undef RANK
  }

  // Divide array among openmp processes equitably
  int thread_start = (int)(get_rand() * nthreads);
  for (int ip = 0; ip < npart_local; ++ip) {
    int               index = (thread_start + ip) % nthreads;
    struct of_photon *ph    = safe_malloc(sizeof(struct of_photon));
    copy_photon(&superphotons[ip], ph);
    swap_ph(&ph, &photon_lists[index]);
  }

  free(superphotons);
#endif // RADIATION

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  timer_stop(TIMER_OUT);
}

// Use fluid restart data as initial conditions, usually for a GRRMHD simulation
void init_fluid_restart() {
  hsize_t file_grid_dims[4], file_start[4], file_count[4];
  hsize_t mem_grid_dims[4], mem_start[4];
  int     fdims[3] = {N1TOT, N2TOT, N3TOT};
  int     mdims[3] = {N1, N2, N3};

  // MAKE SURE N1TOT, N2TOT, N3TOT are right!

  if (mpi_io_proc()) {
    fprintf(stderr, "Initializing fluid from %s\n\n", init_from_grmhd);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(init_from_grmhd, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // Read some header data? gam? game? gamp?

  // Read fluid primitive data
  {
    char *name = "p";
    for (int d = 0; d < 3; d++)
      file_grid_dims[d] = fdims[d];
    file_grid_dims[3] = NVAR; // For vectors
    hid_t filespace   = H5Screate_simple(4, file_grid_dims, NULL);
    for (int d = 0; d < 3; d++) {
      file_start[d] = global_start[d + 1];
      file_count[d] = mdims[d];
    }
    file_start[3] = 0;
    file_count[3] = NVAR;
    H5Sselect_hyperslab(
        filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

    for (int d = 0; d < 3; d++) {
      mem_grid_dims[d] = mdims[d] + 2 * NG;
    }
    mem_grid_dims[3] = NVAR;
    hid_t memspace   = H5Screate_simple(4, mem_grid_dims, NULL);
    for (int d = 0; d < 3; d++)
      mem_start[d] = NG;
    mem_start[3] = 0;
    H5Sselect_hyperslab(
        memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

    plist_id      = H5Pcreate(H5P_DATASET_ACCESS);
    hid_t dset_id = H5Dopen(file_id, name, plist_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
        &P[0][0][0][0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);
}

int restart_init() {
  char fname[STRLEN], lastname[STRLEN];
  sprintf(fname, "restart.last");
  strcpy(lastname, restartdir);
  strcat(lastname, fname);

#if RADIATION
  {
    photon_lists = safe_malloc((size_t)nthreads * sizeof(struct of_photon *));
    photon_mpi_lists =
        safe_malloc((size_t)nthreads * sizeof(struct of_photon *));
#pragma omp parallel
    {
      photon_lists[omp_get_thread_num()]     = NULL;
      photon_mpi_lists[omp_get_thread_num()] = NULL;
    }
  }
#endif // RADIATION

  FILE *fp = fopen(lastname, "r");
  if (fp == NULL) {
    if (mpi_io_proc()) {
      fprintf(stdout, "No restart file\n\n");
    }
    return 0;
  }

  if (mpi_io_proc())
    fprintf(stdout, "Loading restart file\n\n");
  zero_arrays();

  safe_fscanf(fp, "%s\n", fname);
  restart_read(fname);

  ZSLOOP(-NG, N1 - 1 + NG, -NG, N2 - 1 + NG, -NG, N3 - 1 + NG) {
    PLOOP { Ph[i][j][k][ip] = P[i][j][k][ip]; }
  }

  set_grid();

  bound_prim(P);

#if NEED_UNITS
  {
    set_units();

#if EOS == EOS_TYPE_TABLE && POLYTROPE_FALLBACK
    { rho_poly_thresh = EOS_SC_get_min_rho(); }
#endif

#if RADIATION
    {
      ZLOOP {
        sim_vol += ggeom[i][j][CENT].g * dx[1] * dx[2] * dx[3] * L_unit *
                   L_unit * L_unit;
      }
      sim_vol = mpi_reduce(sim_vol);

      init_emissivity();

      set_weight(P, extra);

#if SCATTERING
      { init_all_hotcross(); }
#endif

      // Check boundary conditions and put photons where they belong
      if (mpi_myrank() == 0) {
        printf("Shuffling particles\n");
      }
      // number needed for particle to cross whole domain
      int max_comm_steps = 2 * (N1CPU + N2CPU + N3CPU);
      for (int comm_steps = 0; comm_steps < max_comm_steps; ++comm_steps) {
        bound_superphotons(P, t, 0);
      }
      printf("[%d] nph after balancing = %lu\n", mpi_myrank(),
          (unsigned long int)count_particles_local());
    }
#endif // RADIATION
  }
#endif // NEED_UNITS

  reset_dump_variables();

  return 1;
}

#if RADIATION && TRACERS
void dump_tracers() {
  timer_start(TIMER_OUT);

  char name[STRLEN], fname[STRLEN];

  sprintf(fname, "tracers_%08d.h5part", dumptrace_id);
  strcpy(name, tracerdir);
  strcat(name, fname);
  if (mpi_io_proc()) {
    fprintf(stdout, "TRACERS %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (file_id < 0) {
    fprintf(stderr, "Could not create dump file! Exiting...\n");
    exit(-1);
  }
  H5Pclose(plist_id);

  hid_t group_id =
      H5Gcreate(file_id, "/Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);

  double Time = t;
  WRITE_HDR(Time, TYPE_DBL);
  WRITE_HDR(nstep, TYPE_INT);
  WRITE_HDR(dumptrace_id, TYPE_INT);
#if NEED_UNITS
  {
    WRITE_HDR(L_unit, TYPE_DBL);
    WRITE_HDR(M_unit, TYPE_DBL);
    WRITE_HDR(RHO_unit, TYPE_DBL);
    WRITE_HDR(T_unit, TYPE_DBL);
    WRITE_HDR(U_unit, TYPE_DBL);
    WRITE_HDR(B_unit, TYPE_DBL);
#if EOS == EOS_TYPE_TABLE
    WRITE_HDR(TEMP_unit, TYPE_DBL);
#endif // EOS_TYPE_TABLE
  }
#endif

  size_t ntracers_local     = count_tracers_local();
  size_t ntracers           = mpi_reduce_int(ntracers_local);
  size_t ntracers_before_me = mpi_accumulate_int(ntracers_local);
  WRITE_HDR(ntracers, TYPE_INT);

  // arrays to fil
  int *   id     = safe_malloc(ntracers_local * sizeof(int));
  int *   it     = safe_malloc(ntracers_local * sizeof(int));
  int *   active = safe_malloc(ntracers_local * sizeof(int));
  double *time   = safe_malloc(ntracers_local * sizeof(double));
  double *mass   = safe_malloc(ntracers_local * sizeof(double));
  double *Xharm  = safe_malloc(3 * ntracers_local * sizeof(double));
  double *Xcart  = safe_malloc(3 * ntracers_local * sizeof(double));
  double *x      = safe_malloc(ntracers_local * sizeof(double));
  double *y      = safe_malloc(ntracers_local * sizeof(double));
  double *z      = safe_malloc(ntracers_local * sizeof(double));
  double *rho    = safe_malloc(ntracers_local * sizeof(double));
  double *uu     = safe_malloc(ntracers_local * sizeof(double));
  double *s      = safe_malloc(ntracers_local * sizeof(double));
  double *Press  = safe_malloc(ntracers_local * sizeof(double));
  double *T      = safe_malloc(ntracers_local * sizeof(double));

#if EOS == EOS_TYPE_TABLE
  double *Ye    = safe_malloc(ntracers_local * sizeof(double));
  double *Ye_em = safe_malloc(ntracers_local * sizeof(double));
#endif

  double *rate_emitted =
      safe_malloc(RAD_NUM_TYPES * ntracers_local * sizeof(double));
  double *rate_absorbed =
      safe_malloc(RAD_NUM_TYPES * ntracers_local * sizeof(double));

  double *ucon = safe_malloc(NDIM * ntracers_local * sizeof(double));
  double *ucov = safe_malloc(NDIM * ntracers_local * sizeof(double));
  double *bcon = safe_malloc(NDIM * ntracers_local * sizeof(double));
  double *bcov = safe_malloc(NDIM * ntracers_local * sizeof(double));

  // horrible and inefficient loop through tracer particles
  size_t          idx, n = 0;
  double          X_local[NDIM];
  double          Xcart_local[NDIM];
  double          Kcon[NDIM];
  double          Kcov[NDIM];
  double          P_interp[NVAR];
  double          extra[EOS_NUM_EXTRA];
  double          Nem_interp[RAD_NUM_TYPES];
  double          Nabs_interp[RAD_NUM_TYPES];
  struct of_geom  g;
  struct of_state q;
  for (int thread = 0; thread < nthreads; thread++) {
    struct of_photon *ph = photon_lists[thread];
    while (ph != NULL) {
      if (ph->type == TYPE_TRACER) {

        // IDs
        id[n]     = tracer_get_id(ph);
        active[n] = 1;

        // time
        it[n]   = nstep;
        time[n] = t;

        // mass
        mass[n] = tracer_get_mass(ph);

        // Interpolate prim to tracer position
        tracer_get_X_u(ph, t, X_local, Kcov, Kcon);
        lagrange_interp_prim_3d(P_interp, X_local, P);
        cart_coord(X_local, Xcart_local);
        lagrange_interp_grid_3d(
            Nem_interp, X_local, &(Nem_phys[0][0][0][0]), RAD_NUM_TYPES);
        lagrange_interp_grid_3d(
            Nabs_interp, X_local, &(Nabs_phys[0][0][0][0]), RAD_NUM_TYPES);

        // Set metric at tracer position
        set_metric(X_local, &g);
        get_state(P_interp, &g, &q);

        // Coordinates
        idx = n * 3;
        for (int m = 0; m < 3; m++) {
          Xharm[idx + m] = X_local[m + 1];
        }
        for (int m = 0; m < 3; m++) {
          Xcart[idx + m] = Xcart_local[m + 1];
        }
        x[n] = Xcart[1];
        y[n] = Xcart[2];
        z[n] = Xcart[3];

        // prims
#if EOS == EOS_TYPE_TABLE
        EOS_SC_fill(P_interp, extra);
#endif
        rho[n]   = P_interp[RHO];
        uu[n]    = P_interp[UU];
        s[n]     = EOS_entropy_rho0_u(P_interp[RHO], P_interp[UU], extra);
        Press[n] = EOS_pressure_rho0_u(P_interp[RHO], P_interp[UU], extra);
        T[n]     = EOS_temperature(P_interp[RHO], P_interp[UU], extra);
#if EOS == EOS_TYPE_TABLE
        Ye[n]    = extra[EOS_YE];
        Ye_em[n] = P_interp[YE_EM];
#endif

        // Rates
        idx = n * RAD_NUM_TYPES;
        TYPELOOP {
          rate_emitted[idx + itp]  = Nem_interp[itp];
          rate_absorbed[idx + itp] = Nabs_interp[itp];
        }

        // vectors
        idx = n * NDIM;
        DLOOP1 { ucon[idx + mu] = q.ucon[mu]; }
        DLOOP1 { ucov[idx + mu] = q.ucov[mu]; }
        DLOOP1 { bcon[idx + mu] = q.bcon[mu]; }
        DLOOP1 { bcov[idx + mu] = q.bcon[mu]; }

        // Iterate!
        n++;
      }
      ph = ph->next;
    }
  }
  if (n != ntracers_local) {
    fprintf(stderr,
        "[dump_tracers]: INDEXING ERROR!\n"
        "\tn                  = %lu\n"
        "\tntracers_local     = %lu\n"
        "\tntracers_before_me = %lu\n"
        "\tntracers_total     = %lu\n",
        (unsigned long int)n, (unsigned long int)ntracers_local,
        (unsigned long int)ntracers_before_me, (unsigned long int)ntracers);
    exit(1);
  }

  { // HDF5 output, scalars
#define RANK (1)
    hsize_t fdims[RANK]  = {ntracers};
    hsize_t fstart[RANK] = {ntracers_before_me};
    hsize_t fcount[RANK] = {ntracers_local};
    hsize_t mdims[RANK]  = {ntracers_local};
    hsize_t mstart[RANK] = {0};
    write_array((void *)id, "/Step#0/id", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_INT);
    write_array((void *)it, "/Step#0/it", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_INT);
    write_array((void *)active, "/Step#0/active", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_INT);
    write_array((void *)time, "/Step#0/time", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)x, "/Step#0/x", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)y, "/Step#0/y", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)z, "/Step#0/z", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)mass, "/Step#0/mass", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)rho, "/Step#0/rho", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)uu, "/Step#0/uu", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)s, "/Step#0/s", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)Press, "/Step#0/Press", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)T, "/Step#0/T", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
#if EOS == EOS_TYPE_TABLE
    write_array((void *)Ye, "/Step#0/Ye", RANK, fdims, fstart, fcount, mdims,
        mstart, TYPE_FLOAT);
    write_array((void *)Ye_em, "/Step#0/Ye_em", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
#endif
#undef RANK
  }

  { // HDF5 output 3-vectors
#define RANK (2)
    hsize_t fdims[RANK]  = {ntracers, 3};
    hsize_t fstart[RANK] = {ntracers_before_me, 0};
    hsize_t fcount[RANK] = {ntracers_local, 3};
    hsize_t mdims[RANK]  = {ntracers_local, 3};
    hsize_t mstart[RANK] = {0, 0};
    write_array((void *)Xharm, "/Step#0/Xharm", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)Xcart, "/Step#0/Xcart", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
#undef RANK
  }

  { // HDF5 output rates
#define RANK (2)
    hsize_t fdims[RANK]  = {ntracers, RAD_NUM_TYPES};
    hsize_t fstart[RANK] = {ntracers_before_me, 0};
    hsize_t fcount[RANK] = {ntracers_local, RAD_NUM_TYPES};
    hsize_t mdims[RANK]  = {ntracers_local, RAD_NUM_TYPES};
    hsize_t mstart[RANK] = {0, 0};
    write_array((void *)rate_emitted, "/Step#0/rate_emitted", RANK, fdims,
        fstart, fcount, mdims, mstart, TYPE_DBL);
    write_array((void *)rate_absorbed, "/Step#0/rate_absorbed", RANK, fdims,
        fstart, fcount, mdims, mstart, TYPE_DBL);
  }

  { // HDF5 output 4-vectors
#define RANK (2)
    hsize_t fdims[RANK]  = {ntracers, NDIM};
    hsize_t fstart[RANK] = {ntracers_before_me, 0};
    hsize_t fcount[RANK] = {ntracers_local, NDIM};
    hsize_t mdims[RANK]  = {ntracers_local, NDIM};
    hsize_t mstart[RANK] = {0, 0};
    write_array((void *)ucon, "/Step#0/ucon", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)ucov, "/Step#0/ucov", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
#if TRACERS_SAVE_BCON
    write_array((void *)bcon, "/Step#0/bcon", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
    write_array((void *)bcov, "/Step#0/bcov", RANK, fdims, fstart, fcount,
        mdims, mstart, TYPE_FLOAT);
#endif
#undef RANK
  }

  free(id);
  free(active);
  free(it);
  free(time);
  free(mass);
  free(Xharm);
  free(Xcart);
  free(x);
  free(y);
  free(z);
  free(rho);
  free(uu);
  free(T);
  free(s);
  free(Press);
#if EOS == EOS_TYPE_TABLE
  free(Ye);
#endif
  free(rate_emitted);
  free(rate_absorbed);
  free(ucon);
  free(ucov);
  free(bcon);
  free(bcov);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  dumptrace_id++;

  timer_stop(TIMER_OUT);
}
#endif

void write_array(void *data, const char *name, hsize_t rank, hsize_t *fdims,
    hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart,
    hsize_t type) {
  hid_t filespace, memspace, plist_id, dset_id;
  filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount, NULL);
  memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount, NULL);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  if (type == TYPE_FLOAT) {
    size_t ntot = 1;
    for (size_t n = 0; n < rank; n++) {
      ntot *= mdims[n];
    }
    float *buf = safe_malloc(ntot * sizeof(float));
    for (size_t n = 0; n < ntot; n++) {
      buf[n] = (float)(((double *)data)[n]);
    }
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
        plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, buf);
    free(buf);
  } else if (type == TYPE_DBL) {
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace,
        H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  } else if (type == TYPE_INT) {
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT,
        plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
  } else if (type == TYPE_STR) {
    hid_t string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, strlen(data));
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id  = H5Dcreate(file_id, name, string_type, filespace, H5P_DEFAULT,
        plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id, string_type, memspace, filespace, plist_id, data);
  } else {
    dset_id = H5Dcreate(
        file_id, name, type, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, type, memspace, filespace, plist_id, data);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
  H5Sclose(filespace);
}

void read_array(void *data, const char *name, hsize_t rank, hsize_t *fdims,
    hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart,
    hsize_t type) {
  hid_t filespace, memspace, plist_id, dset_id;
  filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount, NULL);
  memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount, NULL);
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, name, plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);

  if (type == TYPE_DBL) {
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
    if (rank == 1 && fdims[0] == 1) {
      mpi_dbl_broadcast(data);
    }
  } else if (type == TYPE_INT) {
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
    if (rank == 1 && fdims[0] == 1) {
      mpi_int_broadcast(data);
    }
  } else {
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id, type, memspace, filespace, plist_id, data);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(memspace);
  H5Sclose(filespace);
}

hsize_t product_hsize_t(hsize_t a[], int size) {
  hsize_t out = 1;
  for (int i = 0; i < size; i++) {
    out *= a[i];
  }
  return out;
}
