/******************************************************************************
 *                                                                            *
 * EOS_STELLAR_COLLAPSE.c                                                     *
 *                                                                            *
 * IMPLEMENTS ROUTINES FOR READING EOS TABLES PROVIDED ON STELLARCOLLAPSE.ORG *
 *                                                                            *
 ******************************************************************************/

/*
 * TODO: Make EOS reader use gamma law as fallback if we fall off the table
 *       by going to too low densities.
 */

#include "eos_stellar_collapse.h"

#if EOS == EOS_TYPE_TABLE

// HDF5
#include <hdf5.h>

// Init
// ----------------------------------------------------------------------
void EOS_SC_init(char *name) {
  hsize_t file_grid_dims[3], file_start[3], file_count[3];
  hsize_t mem_grid_dims[3], mem_start[3];

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(name, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  hsize_t one       = 1;
  file_grid_dims[0] = 1;
  hid_t filespace   = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0]     = 0;
  file_count[0]     = 1;
  H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace = H5Screate_simple(1, &one, NULL);

  plist_id      = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dopen(file_id, "pointsrho", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  int status =
      H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &Nrho);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "pointstemp", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &NT);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "pointsye", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &NYe);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "energy_shift", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, &e_shift);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  int tab_size    = NYe * NT * Nrho;
  int tab_size_2d = NYe * Nrho;

  tab_lrho = safe_malloc(Nrho * sizeof(double));
  tab_lT   = safe_malloc(NT * sizeof(double));
  tab_Ye   = safe_malloc(NYe * sizeof(double));

  tab_tmp        = safe_malloc(tab_size * sizeof(double));
  tab_lP         = safe_malloc(tab_size * sizeof(double));
  tab_ent        = safe_malloc(tab_size * sizeof(double));
  tab_cs2        = safe_malloc(tab_size * sizeof(double));
  tab_Xa         = safe_malloc(tab_size * sizeof(double));
  tab_Xh         = safe_malloc(tab_size * sizeof(double));
  tab_Xn         = safe_malloc(tab_size * sizeof(double));
  tab_Xp         = safe_malloc(tab_size * sizeof(double));
  tab_Abar       = safe_malloc(tab_size * sizeof(double));
  tab_Zbar       = safe_malloc(tab_size * sizeof(double));
  tab_le         = safe_malloc(tab_size * sizeof(double));
  tab_lwmrho     = safe_malloc(tab_size * sizeof(double));
  tab_hm1        = safe_malloc(tab_size * sizeof(double));
  tab_lhm1       = safe_malloc(tab_size * sizeof(double));
  tab_dpderho    = safe_malloc(tab_size * sizeof(double));
  tab_dpdrhoe    = safe_malloc(tab_size * sizeof(double));
  tab_poly_gamma = safe_malloc(tab_size * sizeof(double));
  tab_poly_K     = safe_malloc(tab_size * sizeof(double));

  tab_le_min_2d     = safe_malloc(tab_size_2d * sizeof(double));
  tab_le_max_2d     = safe_malloc(tab_size_2d * sizeof(double));
  tab_lP_min_2d     = safe_malloc(tab_size_2d * sizeof(double));
  tab_lP_max_2d     = safe_malloc(tab_size_2d * sizeof(double));
  tab_lwmrho_min_2d = safe_malloc(tab_size_2d * sizeof(double));
  tab_lwmrho_max_2d = safe_malloc(tab_size_2d * sizeof(double));
  // tab_hm1_min_1d    = safe_malloc(NYe*sizeof(double));

  hsize_t Narr      = Nrho;
  file_grid_dims[0] = Narr;
  filespace         = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0]     = 0;
  file_count[0]     = Narr;
  H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace = H5Screate_simple(1, &Narr, NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "logrho", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lrho);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  Narr              = NT;
  file_grid_dims[0] = Narr;
  filespace         = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0]     = 0;
  file_count[0]     = Narr;
  H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace = H5Screate_simple(1, &Narr, NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "logtemp", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lT);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  Narr              = NYe;
  file_grid_dims[0] = Narr;
  filespace         = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0]     = 0;
  file_count[0]     = Narr;
  H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace = H5Screate_simple(1, &Narr, NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "ye", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Ye);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  file_grid_dims[0] = NYe;
  file_grid_dims[1] = NT;
  file_grid_dims[2] = Nrho;
  file_start[0]     = 0;
  file_start[1]     = 0;
  file_start[2]     = 0;
  file_count[0]     = NYe;
  file_count[1]     = NT;
  file_count[2]     = Nrho;

  filespace = H5Screate_simple(3, file_grid_dims, NULL);
  H5Sselect_hyperslab(
      filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = NYe;
  mem_grid_dims[1] = NT;
  mem_grid_dims[2] = Nrho;
  memspace         = H5Screate_simple(3, mem_grid_dims, NULL);
  mem_start[0]     = 0;
  mem_start[1]     = 0;
  mem_start[2]     = 0;
  H5Sselect_hyperslab(
      memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "logpress", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lP);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "logenergy", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_le);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "entropy", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_ent);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "dpdrhoe", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dread(
      dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_dpdrhoe);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "dpderho", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dread(
      dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_dpderho);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Xa", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Xa);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Xh", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Xh);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Xn", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Xn);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Xp", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Xp);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Abar", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Abar);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "Zbar", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status =
      H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Zbar);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id  = H5Dopen(file_id, "gamma", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  status = H5Dread(
      dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_poly_gamma);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  H5Fclose(file_id);

  // Manage energy shift
  h_shift_cs = 0.0;
#if !USE_ENERGY_SHIFT
  {
#if OFFSET_H_CS && !OFFSET_U
    h_shift_cs = e_shift;
#endif
    e_shift = 0.0;
  }
#endif

// Filter noisy table data
#if DO_SMOOTH
  memcpy(tab_tmp, tab_dpdrhoe, tab_size * sizeof(double));
  median_filter(tab_tmp, tab_dpdrhoe);
  memcpy(tab_tmp, tab_dpderho, tab_size * sizeof(double));
  median_filter(tab_tmp, tab_dpderho);
#endif

  // mins and maxes of indep. variables
  tab_lrho_min = tab_lrho[0];
  tab_lrho_max = tab_lrho[Nrho - 1];
  tab_dlrho    = (tab_lrho_max - tab_lrho_min) / (Nrho - 1);
  tab_rho_min  = pow(10., tab_lrho_min);
  tab_rho_max  = pow(10., tab_lrho_max);

  tab_lT_min = tab_lT[0];
  tab_lT_max = tab_lT[NT - 1];
  tab_dlT    = (tab_lT_max - tab_lT_min) / (NT - 1);
  tab_T_min  = pow(10., tab_lT_min);
  tab_T_max  = pow(10., tab_lT_max);

  tab_Ye_min = tab_Ye[0];
  tab_Ye_max = tab_Ye[NYe - 1];
  tab_dYe    = (tab_Ye_max - tab_Ye_min) / (NYe - 1);

  // mins and maxes of dependent variables
  tab_lP_min      = find_min(tab_lP, tab_size);
  tab_lP_max      = find_max(tab_lP, tab_size);
  tab_ent_min     = find_min(tab_ent, tab_size);
  tab_ent_max     = find_max(tab_ent, tab_size);
  tab_le_min      = find_min(tab_le, tab_size);
  tab_le_max      = find_max(tab_le, tab_size);
  tab_Xa_min      = find_min(tab_Xa, tab_size);
  tab_Xa_max      = find_max(tab_Xa, tab_size);
  tab_Xh_min      = find_min(tab_Xh, tab_size);
  tab_Xh_max      = find_max(tab_Xh, tab_size);
  tab_Xn_min      = find_min(tab_Xn, tab_size);
  tab_Xn_max      = find_max(tab_Xn, tab_size);
  tab_Xp_min      = find_min(tab_Xp, tab_size);
  tab_Xp_max      = find_max(tab_Xp, tab_size);
  tab_Abar_min    = find_min(tab_Abar, tab_size);
  tab_Abar_max    = find_max(tab_Abar, tab_size);
  tab_Zbar_min    = find_min(tab_Zbar, tab_size);
  tab_Zbar_max    = find_max(tab_Zbar, tab_size);
  tab_dpdrhoe_min = find_min(tab_dpdrhoe, tab_size);
  tab_dpdrhoe_max = find_max(tab_dpdrhoe, tab_size);
  tab_dpderho_min = find_min(tab_dpderho, tab_size);
  tab_dpderho_max = find_max(tab_dpderho, tab_size);

  tab_e_min = le2e(tab_le_min);
  tab_e_max = le2e(tab_le_max);
  tab_P_min = pow(10., tab_lP_min);
  tab_P_max = pow(10., tab_lP_max);
  fill_max_min_2d(tab_le_max_2d, tab_le_min_2d, tab_le);
  fill_max_min_2d(tab_lP_max_2d, tab_lP_min_2d, tab_lP);

#if OFFSET_U
  { // Find u shift
    tab_u_min = INFINITY;
    double lrho, rho, le, e, u;
    for (int irho = 0; irho < Nrho; irho++) {
      lrho = tab_lrho[irho];
      rho  = pow(10., lrho);
      for (int iT = 0; iT < NT; iT++) {
        for (int iY = 0; iY < NYe; iY++) {
          le = tab_le[EOS_ELEM(irho, iT, iY)];
          e  = le2e(le);
          u  = rho * e;
          if (u < tab_u_min) {
            tab_u_min = u;
          }
        }
      }
    }
    if (tab_u_min <= 0) {
      u_shift = -(1.0 + 1e-10) * tab_u_min + SMALL;
      tab_u_min += u_shift;
    } else {
      u_shift = 0.0;
    }
  }
#else
  u_shift            = 0.0;
#endif

  { // make enthalpy table
    double *tab_wmrho = safe_malloc(tab_size * sizeof(double));
    // double *tab_hm1   = safe_malloc(tab_size*sizeof(double));
    double lrho, rho, le, lP, e, u, P, w, lw, h, lh;
    // first tabulate the real enthalpy
    for (int irho = 0; irho < Nrho; irho++) {
      lrho = tab_lrho[irho];
      rho  = pow(10., lrho);
      for (int iT = 0; iT < NT; iT++) {
        for (int iY = 0; iY < NYe; iY++) {
          le = tab_le[EOS_ELEM(irho, iT, iY)];
          lP = tab_lP[EOS_ELEM(irho, iT, iY)];
          e  = le2e(le);
          P  = pow(10., lP);
          // use log enthalpy minus rho.
          u = e2u(rho, e);
          w = u + P; // contains u_shift and is > 0
          // h = w/rho;
          h = e + P / rho; // does not contain u_shift
          // lh = h2lh(h);
          tab_wmrho[EOS_ELEM(irho, iT, iY)] = w;
          tab_hm1[EOS_ELEM(irho, iT, iY)]   = h;
          // tab_lhm1[EOS_ELEM(irho,iT,iY)]  = lh;
        }
      }
    }

    // Calculate min and max of wmrho
    tab_wmrho_min = find_min(tab_wmrho, tab_size);
    tab_wmrho_max = find_max(tab_wmrho, tab_size);
    tab_hm1_min   = find_min(tab_hm1, tab_size);
    tab_hm1_max   = find_max(tab_hm1, tab_size);
    // tab_lhm1_min  = find_min(tab_lhm1,  tab_size);
    // tab_lhm1_max  = find_max(tab_lhm1,  tab_size);
    // fill_min_1d(tab_hm1_min_1d, tab_hm1);

    // Next check to see if we need a min value
    // arbitrary. Just need to make enthalpy + shift > 0
    if (tab_wmrho_min <= 0) {
      w_shift = -1.01 * (tab_wmrho_min - SMALL);
    } else {
      w_shift = 0.0;
    }
    if (tab_hm1_min <= 0) {
      h_shift = -1.01 * (tab_hm1_min - SMALL);
    } else {
      h_shift = 0.0;
    }

    // Set the log enthalpy
    for (int irho = 0; irho < Nrho; irho++) {
      for (int iT = 0; iT < NT; iT++) {
        for (int iY = 0; iY < NYe; iY++) {
          w  = tab_wmrho[EOS_ELEM(irho, iT, iY)];
          h  = tab_hm1[EOS_ELEM(irho, iT, iY)];
          lw = w2lw(w);
          lh = h2lh(h);
          tab_lwmrho[EOS_ELEM(irho, iT, iY)] = lw;
          tab_lhm1[EOS_ELEM(irho, iT, iY)]   = lh;
        }
      }
    }
    free(tab_wmrho);
    // free(tab_hm1);
  }
  tab_lwmrho_min = find_min(tab_lwmrho, tab_size);
  tab_lwmrho_max = find_max(tab_lwmrho, tab_size);
  tab_lhm1_min   = find_min(tab_lhm1, tab_size);
  tab_lhm1_max   = find_max(tab_lhm1, tab_size);
  fill_max_min_2d(tab_lwmrho_max_2d, tab_lwmrho_min_2d, tab_lwmrho);

  { // make cs2 table
    int elem;
    for (int irho = 0; irho < Nrho; irho++) {
      double lrho = tab_lrho[irho];
      double rho  = pow(10., lrho);
      for (int iT = 0; iT < NT; iT++) {
        for (int iY = 0; iY < NYe; iY++) {
          elem           = EOS_ELEM(irho, iT, iY);
          double hm1     = tab_hm1[elem];
          double h       = hm1 + CL * CL - h_shift_cs;
          double lP      = tab_lP[elem];
          double P       = pow(10., lP);
          double dpdrhoe = tab_dpdrhoe[elem];
          double dpderho = tab_dpderho[elem];
          double cs2     = (dpdrhoe + (P / (rho * rho)) * dpderho) / h;
          // Tables can have weird kinks in them messing up
          // this analytically guaranteed relationship.
          // Fix this by hand.
          // problems only emerge for densities >= 10^{14} g/cm^3
          cs2 = MY_MAX(cs2, 10. * SMALL);
#if DO_SMOOTH
          if (cs2 > CL * CL && iT < NT - 2 && irho < Nrho - 3) { // DEBUG
            if (mpi_io_proc()) {
              fprintf(stderr,
                  "cs2 > 1!\n"
                  "\t[irho,iT,iY] = [%d, %d, %d]\n",
                  irho, iT, iY);
            }
          }
#endif
          cs2           = MY_MIN(cs2, 1.0 - 10. * SMALL);
          tab_cs2[elem] = cs2;
        }
      }
    }
  }
  tab_cs2_min = find_min(tab_cs2, tab_size);
  tab_cs2_max = find_max(tab_cs2, tab_size);

  { // make polytrope table
    int elem;
    for (int irho = 0; irho < Nrho; irho++) {
      double lrho = tab_lrho[irho];
      for (int iT = 0; iT < NT; iT++) {
        for (int iY = 0; iY < NYe; iY++) {
          elem             = EOS_ELEM(irho, iT, iY);
          double lP        = tab_lP[elem];
          double Gam       = tab_poly_gamma[elem];
          double lK        = lP - Gam * lrho;
          double K         = pow(10., lK);
          tab_poly_K[elem] = K;
        }
      }
    }
  }

  // sanity checks
  if (is_practically_nan(tab_lwmrho_min)) {
    fprintf(stderr, "[EOS_SC_init]: log enthalpy is nan.\n");
    exit(1);
  }
  if (status) {
    fprintf(stderr, "[EOS_SC_init]: HDF5 Returned an error. Status = %d.\n",
        status);
    exit(1);
  }
}

void do_ye_fixup(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]) {
  pv[YE] = catch_ye(pv[YE]);
}

void do_ye_em_fixup(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]) {
  pv[YE_EM] = catch_ye(pv[YE_EM]);
}

void do_atm_fixup(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]) {
  if (pv[ATM] < IS_ATM)
    pv[ATM] = IS_ATM;
  if (pv[ATM] > NOT_ATM)
    pv[ATM] = NOT_ATM;
}
// ----------------------------------------------------------------------

// Front-facing API
// ----------------------------------------------------------------------
void EOS_SC_fill(double *restrict p, double *restrict eos) {
  double lTguess, leosTemp;
  double lrho, e, le;
  double u   = p[UU];
  double rho = p[RHO];
  double ye  = p[YE];
  // double yedens = p[YE];
  double lT = eos[EOS_LT];

  // double ye     = yedens / (fabs(rho) + SMALL);

  // into CGS
  u *= U_unit;
  rho *= RHO_unit;

  // don't fall off the table
  lrho = catch_rho(rho);
  ye   = catch_ye(ye);
  e    = u2e(rho, u);
  // e = u / rho;
  le = catch_e(e);
  le = catch_var_2d(lrho, ye, le, tab_le_min_2d, tab_le_max_2d);
#if SC_MONOTONE_SAFE
  le = catch_var_2d_monotone(lrho, ye, le, tab_le);
#endif // SC_MONOTONE_SAFE
  // Get a good guess
  lTguess = catch_lT(lT);

// crash and die if something went wrong here
#if SC_DEBUG
  if (is_practically_nan(le)) {
    fprintf(stderr, "[EOS_SC_fill %d]: NAN detected!\n", mpi_io_proc());
    fprintf(stderr, "rho     = %.10e\n", rho);
    fprintf(stderr, "u       = %.10e\n", u);
    fprintf(stderr, "e       = %.10e\n", e);
    fprintf(stderr, "lrho    = %.10f\n", lrho);
    fprintf(stderr, "ye      = %.10f\n", ye);
    fprintf(stderr, "le      = %.10f\n", le);
    fprintf(stderr, "lTguess = %.10f\n", lTguess);
    exit(1);
  }
#endif // SC_DEBUG

  int status = find_lT(lrho, lTguess, ye, tab_le, le, &leosTemp);
  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[EOS_SC_fill]: Failed to root find table!\n"
        "\trho      = %e\n"
        "\tu        = %e\n"
        "\te        = %e\n"
        "\tye       = %g\n"
        "\tlrho     = %g\n"
        "\tle       = %g\n"
        "\tlTguess  = %g\n"
        "\tleosTemp = %g\n",
        rho / RHO_unit, u / U_unit, e, ye, lrho, le, lTguess, leosTemp);
    exit(1); // TODO: Handle this more gracefully
  }

  eos[EOS_LRHO] = lrho;
  eos[EOS_LT]   = leosTemp;
  eos[EOS_YE]   = ye;
  return;
}

double EOS_SC_pressure_rho0_u(double lrho, double lT, double ye) {
  const double lP = EOS_SC_interp(lrho, lT, ye, tab_lP);
  return pow(10., lP) / U_unit;
}

double EOS_SC_specific_enthalpy_rho0_u(double lrho, double lT, double ye) {
  const double hm1   = EOS_SC_interp(lrho, lT, ye, tab_hm1);
  const double h_cgs = hm1 + CL * CL;
  const double h     = h_cgs / (CL * CL);
  return h;
}

double EOS_SC_sound_speed(double lrho, double lT, double ye) {
  double cs2 = EOS_SC_interp(lrho, lT, ye, tab_cs2);
  return sqrt(cs2);
}

double EOS_SC_entropy(double lrho, double lT, double ye) {
  double ent = EOS_SC_interp(lrho, lT, ye, tab_ent);
  return ent;
}

double EOS_SC_gamma(double lrho, double lT, double ye) {
  const double cs2    = EOS_SC_interp(lrho, lT, ye, tab_cs2);
  const double lP     = EOS_SC_interp(lrho, lT, ye, tab_lP);
  const double lgamma = lrho - lP;
  double       gamma  = cs2 * pow(10., lgamma);
  if (gamma < 1.1)
    gamma = 1.1;
  if (gamma > 2.0)
    gamma = 2.0;
  return gamma;
}

double EOS_SC_temperature(double lT) {
  // temperature is in MeV to start, which is a fine code unit
  return pow(10., lT); // / TEMP_unit;
}

double EOS_SC_get_u_of_T(double rho, double T, double ye) {
  if (T < tab_T_min)
    T = tab_T_min;
  if (T > tab_T_max)
    T = tab_T_max;
  const double lrho = catch_rho(rho);
  const double lT   = catch_lT(log10(T));
  ye                = catch_ye(ye);
  const double le   = EOS_SC_interp(lrho, lT, ye, tab_le);
  const double e    = le2e(le);
  const double u    = e2u(rho, e);
  // const double u = e*rho;
  return u / U_unit;
}

double EOS_SC_pressure_rho0_w(double rho, double w, double ye, double *lTold) {
  double lTguess = *lTold;
  double leosTemp;
  lTguess = catch_lT(lTguess);

  // Subtract off rest energy
  w -= rho;

  // convert to CGS
  rho *= RHO_unit;
  w *= U_unit;

  // dont' fall off the table
  double lrho = catch_rho(rho);
  double lw   = catch_w(w);
  ye          = catch_ye(ye);

  // force w back onto the table for chosen ilrho and iY
  lw = catch_var_2d(lrho, ye, lw, tab_lwmrho_min_2d, tab_lwmrho_max_2d);
#if SC_MONOTONE_SAFE
  lw = catch_var_2d_monotone(lrho, ye, lw, tab_lwmrho);
#endif // SC_MONOTONE_SAFE

// crash and die if something went wrong here
#if SC_DEBUG
  if (is_practically_nan(lw)) {
    fprintf(
        stderr, "[EOS_SC_Pressure_rho0_w %d]: NAN detected!\n", mpi_io_proc());
    fprintf(stderr, "rho     = %.10e\n", rho);
    fprintf(stderr, "w       = %.10e\n", w);
    fprintf(stderr, "lrho    = %.10f\n", lrho);
    fprintf(stderr, "lw      = %.10f\n", lw);
    fprintf(stderr, "ye      = %.10f\n", ye);
    fprintf(stderr, "lTguess = %.10f\n", lTguess);
    exit(1);
  }
#endif // SC_DEBUG

  int status = find_lT(lrho, lTguess, ye, tab_lwmrho, lw, &leosTemp);

  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[EOS_SC_pressure_rho0_w %d]: "
        "Failed to root find table\n"
        "\trho      = %g\n"
        "\tlrho     = %g\n"
        "\tlTguess  = %g\n"
        "\tye       = %g\n"
        "\tleosTemp = %g\n"
        "\tlwmrho   = %g\n",
        mpi_io_proc(), rho, lrho, lTguess, ye, leosTemp, lw);

    fprintf(stderr, "tab_lwwmrho = \n");
    temp_map(lrho, ye, tab_lwmrho);
    exit(1); // TODO: Handle this more gracefully
  }

  const double log_press = EOS_SC_interp(lrho, leosTemp, ye, tab_lP);
  double       press     = pow(10., log_press);

  // back into code units
  press /= U_unit;

#if SC_DEBUG
  if (is_practically_nan(press)) {
    // TODO: handle this more gracefully.
    fprintf(stderr, "press from enthalpy = NaN.\n");
    exit(1);
  }
#endif // SC_DEBUG

  *lTold = leosTemp;
  return press;
}

double EOS_SC_u_press(double press, double rho, double ye, double *lTold) {
  double lTguess = *lTold;
  double leosTemp; //, lrho, lp;
  lTguess = catch_lT(lTguess);

  // Units to CGS
  press *= U_unit;
  rho *= RHO_unit;

  const double lrho = catch_rho(rho);
  double       lp   = catch_press(press);
  ye                = catch_ye(ye);
  lp                = catch_var_2d(lrho, ye, lp, tab_lP_min_2d, tab_lP_max_2d);

// crash and die if something went wrong here
#if SC_DEBUG
  if (is_practically_nan(lp)) {
    fprintf(
        stderr, "[EOS_SC_Pressure_rho0_w %d]: NAN detected!\n", mpi_io_proc());
    fprintf(stderr, "rho     = %.10e\n", rho);
    fprintf(stderr, "press   = %.10e\n", press);
    fprintf(stderr, "lrho    = %.10f\n", lrho);
    fprintf(stderr, "lP      = %.10f\n", lp);
    fprintf(stderr, "ye      = %.10f\n", ye);
    fprintf(stderr, "lTguess = %.10f\n", lTguess);
    exit(1);
  }
#endif // SC_DEBUG

  int status = find_lT(lrho, lTguess, ye, tab_lP, lp, &leosTemp);
  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[EOS_SC_u_press]: "
        "Failed to root find table %g %g %g %g   %g %d\n",
        lrho, lTguess, ye, leosTemp, rho, mpi_io_proc());
    exit(1); // TODO: Handle this more gracefully
  }
  const double le = EOS_SC_interp(lrho, leosTemp, ye, tab_le);
  const double e  = le2e(le);
  double       u  = e2u(rho, e); // e*rho;

  // to code units
  u /= U_unit;

#if SC_DEBUG
  if (is_practically_nan(u)) {
    // TODO: handle this more gracefully.
    fprintf(stderr, "u from press = NaN.\n");
    fprintf(stderr, "press = %f\n", press);
    fprintf(stderr, "rho = %f\n", rho);
    fprintf(stderr, "ye = %f\n", ye);
    exit(1);
  }
#endif // SC_DEBUG

  *lTold = leosTemp;
  return u;
}

void EOS_SC_mass_fractions(double Xi[NUM_MASS_FRACTIONS], const double *extra) {
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double Ye   = extra[EOS_YE];
  Xi[MF_XA]   = EOS_SC_interp(lrho, lT, Ye, tab_Xa);
  Xi[MF_XH]   = EOS_SC_interp(lrho, lT, Ye, tab_Xh);
  Xi[MF_XN]   = EOS_SC_interp(lrho, lT, Ye, tab_Xn);
  Xi[MF_XP]   = EOS_SC_interp(lrho, lT, Ye, tab_Xp);
}

void EOS_SC_avg_ions(double *Abar, double *Zbar, const double *extra) {
  double lrho = extra[EOS_LRHO];
  double lT   = extra[EOS_LT];
  double Ye   = extra[EOS_YE];
  *Abar       = EOS_SC_interp(lrho, lT, Ye, tab_Abar);
  *Zbar       = EOS_SC_interp(lrho, lT, Ye, tab_Zbar);
}

void EOS_SC_set_floors(double scale, double rho, double u, double ye,
    double bsq, double *rhoflr, double *uflr) {
  *rhoflr = EOS_SC_rho_floor(scale, bsq);
  *uflr   = EOS_SC_u_floor(scale, bsq, ye);
#if BOUND_UORHO && (!OFFSET_U || OLD_U_FLOORS)
  *rhoflr = MY_MAX(*rhoflr, u / UORHOMAX);
#endif
}

double EOS_SC_rho_floor(double scale, double bsq) {
  double rhoflr = RHOMIN * scale;
#if GAMMA_FALLBACK
  double rhominlimit = RHOMINLIMIT;
#else
  double rhominlimit = MY_MAX(RHOMINLIMIT, EOS_SC_get_min_rho());
#endif
  rhoflr = MY_MAX(rhoflr, rhominlimit);
  rhoflr = MY_MAX(rhoflr, bsq / BSQORHOMAX);
  return rhoflr;
}

double EOS_SC_get_min_lrho() { return tab_lrho_min; }

double EOS_SC_get_min_rho() {
  double delrho      = tab_lrho_max - tab_lrho_min;
  double lrho_min    = tab_lrho_min + 0.01 * delrho;
  double rho_min_cgs = pow(10., lrho_min);
  double rho_min     = rho_min_cgs / RHO_unit;
  return rho_min;
}

double EOS_SC_get_min_lT() { return tab_lT_min; }

double EOS_SC_get_minu(double rho, double ye, double scale) {
  double u;
#if OLD_U_FLOORS
  u = MY_MAX(u, UUMIN * scale);
#else
  double lrho        = catch_rho(rho);
  double lT          = tab_lT_min;
  ye                 = catch_ye(ye);
  double le          = EOS_SC_interp(lrho, lT, ye, tab_le);
  double e           = le2e(le);
  u                  = e2u(rho, e); // e*rho;
  u                  = u / U_unit;
#endif
  return u;
}

double EOS_SC_u_floor(double scale, double bsq, double ye) {
  // return -INFINITY;
  double rhoflr = EOS_SC_rho_floor(scale, bsq) * RHO_unit;
  double minu   = EOS_SC_get_minu(rhoflr, ye, scale);
  if ((bsq / BSQOUMAX) > fabs(minu)) { // Good idea?
    minu = bsq / BSQOUMAX;
  }
  return minu;
}

void EOS_SC_print_table_mins() {
  if (mpi_io_proc()) {
    double rho_min_cgs = pow(10., tab_lrho_min);
    double rho_min     = rho_min_cgs / RHO_unit;
    // tab_e_min*rho_min_cgs;
    double u_min_cgs = tab_u_min;
    double u_min     = u_min_cgs / U_unit;
    fprintf(stderr,
        "[EOS_SC_init]: Table mins:\n"
        "\trho      = %g\n"
        "\trho_code = %g\n"
        "\tu        = %g\n"
        "\tu_code   = %g\n",
        rho_min_cgs, rho_min, u_min_cgs, u_min);
    // DEBUG
    /*
    sleep(1);
    double rhoflr,uflr;
    printf("FLOORS:\nr\t\trho\t\tu\n-------------------\n");
    for (double r = 5.0; r <= 250.; r += 5.0) {
      double scale = pow(r,-2);
      EOS_SC_set_floors(scale,1.0,0.0,0.1,0.0,&rhoflr,&uflr);
      printf("%3g\t%10g\t%10g\n",r,rhoflr,uflr);
    }
    exit(1);
    */
  }
}

void EOS_SC_get_polytrope(
    double lrho, double lT, double ye, double *poly_K, double *poly_gamma) {
  lrho          = catch_lrho(lrho);
  lT            = catch_lT(lT);
  ye            = catch_ye(ye);
  double K      = EOS_SC_interp(lrho, lT, ye, tab_poly_K);
  double Gam    = EOS_SC_interp(lrho, lT, ye, tab_poly_gamma);
  double K_unit = U_unit / pow(RHO_unit, Gam);
  *poly_K       = K / K_unit;
  *poly_gamma   = Gam;
}
//----------------------------------------------------------------------

// Root-finding
// ----------------------------------------------------------------------
static int find_lT(const double lrho, double lTguess, const double ye,
    double *restrict tab, const double val, double *lT) {
  struct of_lT_params p;
  p.lrho     = lrho;
  p.ye       = ye;
  p.tab      = tab;
  int status = find_root(&lT_f, &p, val, lTguess, tab_lT_min,
      tab_lT_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lT);
  return status;
}

static int find_adiabat_0d(
    double lrho, double lTguess, double ye, double s, double *lT) {
  lrho    = catch_lrho(lrho);
  lTguess = catch_lT(lTguess);
  ye      = catch_ye(ye);
  s       = catch_s(s);

  struct of_lT_params p;
  p.lrho     = lrho;
  p.ye       = ye;
  p.tab      = tab_ent;
  int status = find_root(&lT_f, &p, s, lTguess, tab_lT_min,
      tab_lT_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lT);
  return status;
}

static void adiabat_lrho_min_max(
    double s, double ye, double *lrho_min, double *lrho_max) {
  s                 = catch_s(s);
  ye                = catch_ye(ye);
  *lrho_min         = catch_lrho(*lrho_min);
  *lrho_max         = catch_lrho(*lrho_max);
  double lrho_guess = 0.5 * (*lrho_min + *lrho_max);

  // lrho_min
  struct of_ent_f_params p;
  p.lT       = tab_lT_min;
  p.ye       = ye;
  int status = find_root(&ent_f_adiabat, &p, s, lrho_guess, tab_lrho_min,
      tab_lrho_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lrho_min);
  if (status != ROOT_SUCCESS || *lrho_min < tab_lrho_min ||
      is_practically_nan(*lrho_min)) {
    *lrho_min = tab_lrho_min + TABLE_TOL;
  }

  // lrho_max
  // struct of_ent_f_params p;
  p.lT   = tab_lT_max;
  p.ye   = ye;
  status = find_root(&ent_f_adiabat, &p, s, lrho_guess, tab_lrho_min,
      tab_lrho_max - 2. * TABLE_TOL / 100., TABLE_TOL, TABLE_FTOL, lrho_max);
  if (status != ROOT_SUCCESS || *lrho_max > tab_lrho_max ||
      is_practically_nan(*lrho_max)) {
    *lrho_max = tab_lrho_max - TABLE_TOL;
  }

  // final fix
  *lrho_min = catch_lrho(*lrho_min);
  *lrho_max = catch_lrho(*lrho_max);
}

double EOS_SC_hm1_min_adiabat(const struct of_adiabat *a) {
  return a->hm1_min / (CL * CL);
}

int EOS_SC_find_adiabat_1d(double s, double ye, struct of_adiabat *a) {
  s               = catch_s(s);
  ye              = catch_ye(ye);
  double lrho_min = tab_lrho_min;
  double lrho_max = tab_lrho_max;
  adiabat_lrho_min_max(s, ye, &lrho_min, &lrho_max);

  // Conservatively choose the max index for lrho
  int ilrho_min = find_index(lrho_min, tab_lrho, Nrho);
  int ilrho_max = find_index(lrho_max, tab_lrho, Nrho) - 1;
  if (ilrho_min < 0)
    ilrho_min = 0;
  if (ilrho_max > Nrho - 2)
    ilrho_max = Nrho - 2;
  lrho_min = tab_lrho[ilrho_min];
  lrho_max = tab_lrho[ilrho_max];

  double  slope     = tab_dlT / tab_dlrho;
  double *lT_of_rho = safe_malloc(Nrho * sizeof(double));

  double hm1_min = INFINITY;
  double hm1_max = -INFINITY;

  for (int i = ilrho_min; i < ilrho_max + 1; i++) {
    double lTguess = slope * tab_lrho[i];
    double lrho    = tab_lrho[i];
    int    status  = find_adiabat_0d(lrho, lTguess, ye, s, &(lT_of_rho[i]));
    double hm1     = EOS_SC_interp(lrho, lT_of_rho[i], ye, tab_hm1);
    if (hm1 < hm1_min)
      hm1_min = hm1;
    if (status != ROOT_SUCCESS)
      return ROOT_FAIL;
    if (hm1 > hm1_max)
      hm1_max = hm1;
  }

  a->s        = s;
  a->ye       = ye;
  a->lrho_min = lrho_min;
  a->lrho_max = lrho_max;
  a->imin     = ilrho_min;
  a->imax     = ilrho_max;
  a->hm1_min  = hm1_min;
  a->hm1_max  = hm1_max;
  a->lT       = lT_of_rho;

  return ROOT_SUCCESS;
}

void EOS_SC_print_adiabat(const struct of_adiabat *a) {
  double ye = a->ye;
  fprintf(stdout, "\n");
  fprintf(stdout, "Adiabat =\n");
  fprintf(stdout, "#----------\n");
  fprintf(stdout, "#lrho\tlT\thm1\ts\n");
  for (int i = a->imin; i < a->imax; i++) {
    double lrho = tab_lrho[i];
    double lT   = a->lT[i];
    double hm1  = EOS_SC_interp(lrho, lT, ye, tab_hm1);
    double s    = EOS_SC_interp(lrho, lT, ye, tab_ent);
    fprintf(stdout, "%e\t%e\t%e\t%e\n", lrho, lT, hm1, s);
  }
  fprintf(stdout, "#----------\n");
  fprintf(stdout, "\n");
}

void EOS_SC_adiabat_free(struct of_adiabat *a) { free(a->lT); }

void EOS_SC_isoentropy_hm1(double hm1, const struct of_adiabat *a,
    double *lrho_guess, double *rho, double *u) {
  hm1         = catch_hm1(hm1 * CL * CL);
  *lrho_guess = catch_lrho(*lrho_guess);
  double s    = catch_s(a->s);
  double ye   = catch_ye(a->ye);
  double lrho;

  // catch hm1 more carefully
  double hm1_min = a->hm1_min;
  if (hm1 < hm1_min)
    hm1 = hm1_min + TABLE_TOL / 100.;
  double lhm1 = h2lh(hm1);

  // root finding
  struct of_lT_adiabat_params p = {tab_lhm1, a};
  int status = find_root(&lT_f_adiabat, &p, lhm1, *lrho_guess, a->lrho_min,
      a->lrho_max, TABLE_TOL, TABLE_FTOL, &lrho);
  if (status != ROOT_SUCCESS) {
    fprintf(stderr,
        "[EOS_SC_isoentropy_hm1]: Failed to find root!\n"
        "\thm1        = %e\n"
        "\thm1_min    = %e\n"
        "\ts          = %e\n"
        "\tye         = %e\n"
        "\tlrho_guess = %e\n",
        hm1, hm1_min, s, ye, *lrho_guess);

    printf("ADIABAT MAP\n");
    printf("rho\thm1\n");
    for (int i = a->imin; i < a->imax; i++) {
      printf("%e\t%e\n", tab_lrho[i], lT_f_adiabat(tab_lrho[i], &p));
    }
    exit(1);
  }

  // and get what we care about
  double lT = interp_1d(
      lrho, a->lrho_min, a->lrho_max, a->imin, a->imax, tab_lrho, a->lT);
  double le = EOS_SC_interp(lrho, lT, ye, tab_le);
  // if (le > tab_le_max) le = tab_le_max - TABLE_TOL;
  //  double s_tab = EOS_SC_interp(lrho,lT,ye,tab_ent);
  double e       = le2e(le);
  double rho_cgs = pow(10., lrho);
  double u_cgs   = e2u(rho_cgs, e); // rho_cgs*e;
  // and to cgs
  *rho = rho_cgs / RHO_unit;
  *u   = u_cgs / U_unit;

  // save initial guesses
  *lrho_guess = lrho;

  return;
}
// ----------------------------------------------------------------------

// Interpolation
// ----------------------------------------------------------------------
static double EOS_SC_interp(
    const double lrho, const double lT, const double Ye, double *restrict tab) {
  // We don't want our indices to fall off the table,
  // but we want to extrapolate appropriately
  double lT2   = catch_lT(lT);
  double Ye2   = catch_ye(Ye);
  double lrho2 = catch_lrho(lrho);

  // indices
  const int irho = (lrho2 - tab_lrho_min) / tab_dlrho;
  const int iT   = (lT2 - tab_lT_min) / tab_dlT;
  const int iY   = (Ye2 - tab_Ye_min) / tab_dYe;

  // assumes evenly spaced table
  const double delrho = (lrho - tab_lrho[irho]) / tab_dlrho;
  const double delT   = (lT - tab_lT[iT]) / tab_dlT;
  const double delY   = (Ye - tab_Ye[iY]) / tab_dYe;

  // one-sided, trilinear interpolation
  return ((1 - delrho) * (1 - delT) * (1 - delY) * tab[EOS_ELEM(irho, iT, iY)] +
          delrho * (1 - delT) * (1 - delY) * tab[EOS_ELEM(irho + 1, iT, iY)] +
          (1 - delrho) * delT * (1 - delY) * tab[EOS_ELEM(irho, iT + 1, iY)] +
          delrho * delT * (1 - delY) * tab[EOS_ELEM(irho + 1, iT + 1, iY)] +
          (1 - delrho) * (1 - delT) * delY * tab[EOS_ELEM(irho, iT, iY + 1)] +
          delrho * (1 - delT) * delY * tab[EOS_ELEM(irho + 1, iT, iY + 1)] +
          (1 - delrho) * delT * delY * tab[EOS_ELEM(irho, iT + 1, iY + 1)] +
          delrho * delT * delY * tab[EOS_ELEM(irho + 1, iT + 1, iY + 1)]);
}

static double interp_2d(
    const double lrho, const double Ye, const double *tab_2d) {
  // indices
  const int irho = (lrho - tab_lrho_min) / tab_dlrho;
  const int iY   = (Ye - tab_Ye_min) / tab_dYe;

  // assumes evenly spaced table
  const double delrho = (lrho - tab_lrho[irho]) / tab_dlrho;
  const double delY   = (Ye - tab_Ye[iY]) / tab_dYe;

  // one-sided, bilinear interpolation
  return ((1 - delrho) * (1 - delY) * tab_2d[MMA_ELEM(irho, iY)] +
          delrho * (1 - delY) * tab_2d[MMA_ELEM(irho + 1, iY)] +
          (1 - delrho) * delY * tab_2d[MMA_ELEM(irho, iY + 1)] +
          delrho * delY * tab_2d[MMA_ELEM(irho + 1, iY + 1)]);
}
// ----------------------------------------------------------------------

// Error functions
// ----------------------------------------------------------------------
static double lT_f(const double lT, const void *params) {
  struct of_lT_params *p = (struct of_lT_params *)params;
  double *restrict tab   = p->tab;
  const double     lrho  = p->lrho;
  const double     ye    = p->ye;
  return EOS_SC_interp(lrho, lT, ye, tab);
}

static double lT_f_adiabat(const double lrho, const void *params) {
  struct of_lT_adiabat_params *p = (struct of_lT_adiabat_params *)params;
  const struct of_adiabat *    a = p->a;
  double *restrict tab           = p->tab;

  const double ye = a->ye;
  const double lT = interp_1d(
      lrho, a->lrho_min, a->lrho_max, a->imin, a->imax, tab_lrho, a->lT);
#if SC_DEBUG
  if (isnan(lrho) || isnan(lT)) {
    fprintf(stderr,
        "[lT_f_adiabat]: NaN detected!\n"
        "\tlrho      = %f\n"
        "\tlT        = %f\n"
        "\tlrho_min  = %f\n"
        "\tlrho_max = %f\n"
        "\timin      = %d\n"
        "\timax      = %d\n"
        "\n",
        lrho, lT, a->lrho_min, a->lrho_max, a->imin, a->imax);
    exit(1);
  }
#endif // SC_DEBUG
  return EOS_SC_interp(lrho, lT, ye, tab);
}

static double ent_f_adiabat(double lrho, const void *params) {
  struct of_ent_f_params *p  = (struct of_ent_f_params *)params;
  double                  lT = catch_lT(p->lT);
  double                  ye = catch_ye(p->ye);
  lrho                       = catch_lrho(lrho);
  double ent                 = EOS_SC_interp(lrho, lT, ye, tab_ent);
  return ent;
}

void EOS_SC_get_bounds(struct of_tablebounds *b) {
  b->Nrho = Nrho;
  b->NT   = NT;
  b->NYe  = NYe;

  b->lrho_min = tab_lrho_min;
  b->lrho_max = tab_lrho_max;
  b->dlrho    = tab_dlrho;

  b->lT_min = tab_lT_min;
  b->lT_max = tab_lT_max;
  b->dlT    = tab_dlT;

  b->Ye_min = tab_Ye_min;
  b->Ye_max = tab_Ye_max;
  b->dYe    = tab_dYe;
}
// ----------------------------------------------------------------------

// Utilities
// ----------------------------------------------------------------------

#if DO_SMOOTH
/*
 * Fills buffer with data from tab contained in
 * box of wide width 2*width+1
 * and centered on (iYe,iT,irho)
 */
void fill_mdn_buffer(double buffer[], int width, int iY, int iT, int irho,
    double *restrict tab) {
  // fill buffer
  int i = 0;
  for (int iWy = -width; iWy < width + 1; iWy++) {
    for (int iWt = -width; iWt < width + 1; iWt++) {
      for (int iWr = -width; iWr < width + 1; iWr++) {
        int elem  = EOS_ELEM(irho + iWr, iT + iWt, iY + iWy);
        buffer[i] = tab[elem];
        i++;
      }
    }
  }
}

/*
 * Uses a median filter to smooth out tab_in and copy it to tab_out
 * Filter has width of MF_W points
 * filter is applied only if data is EPSSMOOTH away from average
 */
static void median_filter(double *restrict tab_in, double *restrict tab_out) {
  double buffer[MF_S];
  // filter, overwriting as needed
  for (int iY = MF_W; iY < NYe - MF_W; iY++) {
    for (int iT = MF_W; iT < NT - MF_W; iT++) {
      for (int irho = MF_W; irho < Nrho - MF_W; irho++) {

        // copy tab_in to tab_out
        int elem      = EOS_ELEM(irho, iT, iY);
        tab_out[elem] = tab_in[elem];

        // overwrite tab_out as needed
        fill_mdn_buffer(buffer, MF_W, iY, iT, irho, tab_in);
        double pnt = tab_in[elem];
        double avg = find_median(buffer, MF_S);
        int    bad = fabs(avg - pnt) / fabs(avg) > EPSSMOOTH;
        if (bad)
          tab_out[elem] = avg;
      }
    }
  }
}
#endif

/*
static void fill_min_1d(double* tab_min_1d, double* tab)
{
  int elem;
  for (int iY = 0; iY < NYe; iY++) {
    double min = INFINITY;
    for (int iT = 0; iT < NT; iT++) {
      for (int irho = 0; irho < Nrho; irho++) {
  elem = EOS_ELEM(irho,iT,iY);
  if (tab[elem] < min) min = tab[elem];
      }
    }
    tab_min_1d[iY] = min;
  }
}
*/

static void fill_max_min_2d(
    double *tab_max_2d, double *tab_min_2d, double *tab) {
  int elem, elem2d;
  for (int irho = 0; irho < Nrho; irho++) {
    for (int iY = 0; iY < NYe; iY++) {
      elem2d     = MMA_ELEM(irho, iY);
      double max = -INFINITY;
      double min = INFINITY;
      for (int iT = 0; iT < NT; iT++) {
        elem = EOS_ELEM(irho, iT, iY);
        if (tab[elem] > max)
          max = tab[elem];
        if (tab[elem] < min)
          min = tab[elem];
      }
      tab_max_2d[elem2d] = max;
      tab_min_2d[elem2d] = min;
    }
  }
}

static double catch_var_2d(const double lrho, const double Ye, const double var,
    const double *tab_min_2d, const double *tab_max_2d) {
#if SC_DEBUG
  {
    if (isnan(lrho) || isnan(Ye) || isnan(var)) {
      fprintf(stderr,
          "[EOS_SC::catch_var_2d]: NaN detected!\n"
          "\tlrho = %e\n"
          "\tYe   = %e\n"
          "\tvar  = %e\n",
          lrho, Ye, var);
    }
  }
#endif // SC_DEBUG

  const double varmin = interp_2d(lrho, Ye, tab_min_2d);
  if (var <= varmin)
    return varmin + TABLE_TOL / 100.;
  const double varmax = interp_2d(lrho, Ye, tab_max_2d);
  if (var >= varmax)
    return varmax - TABLE_TOL / 100.;
  return var;
}

static double catch_var_2d_monotone(const double lrho, const double Ye,
    const double var, double *restrict tab) {
  const double varmin = EOS_SC_interp(lrho, tab_lT_min, Ye, tab);
  if (var <= varmin)
    return varmin + TABLE_TOL / 100.;
  const double varmax = EOS_SC_interp(lrho, tab_lT_max, Ye, tab);
  if (var >= varmax)
    return varmax - TABLE_TOL / 100.;
  return var;
}

static void temp_map(double lrho, double Ye, const double *tab) {
  // Don't fall off table
  if (lrho < tab_lrho_min)
    lrho = tab_lrho_min;
  if (lrho >= tab_lrho_max)
    lrho = tab_lrho_max - TABLE_TOL / 100.;
  if (Ye < tab_Ye_min)
    Ye = tab_Ye_min;
  if (Ye >= tab_Ye_max)
    Ye = tab_Ye_max - TABLE_TOL / 100.;
  // Indices
  const int irho = (lrho - tab_lrho_min) / tab_dlrho;
  const int iY   = (Ye - tab_Ye_min) / tab_dYe;
  for (int iT = 0; iT < NT; iT++) {
    fprintf(stderr, "%d\t%.10f\t%.10f\n", iT, tab_lT[iT],
        tab[EOS_ELEM(irho, iT, iY)]);
  }
}

static double u2e(const double rho, const double u) {
  return (u - u_shift) / rho;
}

static double e2u(const double rho, const double e) {
  return e * rho + u_shift;
}

static double le2e(const double le) { return pow(10., le) - e_shift; }

static double e2le(const double e) { return log10(e + e_shift); }

/*
static double lh2h(const double lh)
{
  return pow(10.,lh) - e_shift;
}
*/

static double h2lh(const double h) { return log10(h + h_shift); }

/*
static double lw2w(const double lw)
{
  return pow(10., lw) - w_shift;
}
*/

static double w2lw(const double w) { return log10(w + w_shift); }

static double catch_rho(const double rho) {
  // dont' fall off the table
  if (rho < tab_rho_min) {
    return tab_lrho_min;
  } else if (rho >= tab_rho_max) {
    return tab_lrho_max - TABLE_TOL / 100.;
  } else {
    return log10(rho);
  }
}

static double catch_lrho(const double lrho) {
  if (lrho < tab_lrho_min) {
    return tab_lrho_min;
  } else if (lrho > tab_lrho_max) {
    return tab_lrho_max - TABLE_TOL / 100.;
  } else {
    return lrho;
  }
}

static double catch_e(const double e) {
  if (e < tab_e_min) {
    return tab_le_min;
  } else if (e >= tab_e_max) {
    return tab_le_max - TABLE_TOL / 100.;
  } else {
    return e2le(e);
  }
}

static double catch_press(const double press) {
  if (press < tab_P_min) {
    return tab_lP_min;
  } else if (press >= tab_P_max) {
    return tab_lP_max - TABLE_TOL / 100.;
  } else {
    return log10(press);
  }
}

static double catch_w(const double w) {
  if (w < tab_wmrho_min) {
    return tab_lwmrho_min;
  } else if (w >= tab_wmrho_max) {
    return tab_lwmrho_max - TABLE_TOL / 100.;
  } else {
    return w2lw(w);
  }
}

/*
static double catch_h(const double h)
{
  if ( h < tab_hm1_min ) {
    return tab_hm1_min;
  } else if ( h >= tab_hm1_max ) {
    return tab_hm1_max - TABLE_TOL/100.;
  } else {
    return h;
  }
}
*/

static double catch_ye(const double ye) {
  if (ye < tab_Ye_min) {
    return tab_Ye_min;
  } else if (ye >= tab_Ye_max) {
    return tab_Ye_max - TABLE_TOL / 100.;
  } else {
    return ye;
  }
}

/*
static double catch_temp(const double temp)
{
  if ( temp < tab_T_min ) {
    return tab_lT_min;
  } else if ( temp >= tab_T_max ) {
    return tab_lT_max - TABLE_TOL/100.;
  } else {
    return log10(temp);
  }
}
*/

static double catch_lT(const double lT) {
  if (lT < tab_lT_min) {
    return tab_lT_min;
  } else if (lT >= tab_lT_max) {
    return tab_lT_max - TABLE_TOL / 100.;
    // } else if (isnan(lT)) {
    //   return 0.5*(tab_lT_min + tab_lT_max);
  } else {
    return lT;
  }
}

static double catch_s(const double s) {
  if (s < tab_ent_min) {
    return tab_ent_min;
  } else if (s >= tab_ent_max) {
    return tab_ent_max - TABLE_TOL / 100.;
  } else {
    return s;
  }
}

static double catch_hm1(const double hm1) {
  if (hm1 < tab_hm1_min) {
    return tab_hm1_min;
  } else if (hm1 >= tab_hm1_max) {
    return tab_hm1_max - TABLE_TOL / 100.;
  } else {
    return hm1;
  }
}
// ----------------------------------------------------------------------
#endif // EOS == EOS_TYPE_TABLE
