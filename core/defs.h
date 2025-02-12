/******************************************************************************
 *                                                                            *
 * DEFS.H                                                                     *
 *                                                                            *
 * GLOBAL VARIABLE DEFINITIONS                                                *
 *                                                                            *
 ******************************************************************************/

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
grid_prim_type       P;
grid_prim_type       F1;
grid_prim_type       F2;
grid_prim_type       F3;
grid_prim_type       Ph;
grid_int_type        pflag;
grid_int_type        fail_save;
grid_int_type        fixup_required;
grid_prim_type       Psave;
grid_fourvector_type jcon;
grid_eosvar_type     extra;
#if RADIATION
grid_radg_type     radG;     // Radiation four-force
grid_radg_type     radG_int; // ...integrated
grid_radg_type     radG_buf; // ...buffer for communication
grid_tensor_type   Rmunu;
grid_int_type      Nsph;
grid_double_type   nph;
#if RZ_HISTOGRAMS
rz_hist_type rz_r_orig_hist, rz_z_orig_hist;
#if NEUTRINO_OSCILLATIONS
rz_hist_type osc_rz_r_orig_hist, osc_rz_z_orig_hist;
#endif // NEUTRINO_OSCILLATIONS
#endif // RZ_HISTOGRAMS

struct of_photon **photon_lists;
struct of_photon **photon_mpi_lists;

double nuLnu[NULNU_IDX0][NTH][NPHI][NU_BINS_SPEC];

#if RADIATION == RADTYPE_NEUTRINOS
double rad_type_counts[RAD_NUM_TYPES];
double lepton_tot, lepton_last, dlepton_tot, dlepton_perc;
double lepton_gas, lepton_rad;
double lepton_lost, lepton_lost_step, lepton_lost_local;
#endif

double Jrad[MAXNSCATT + 2][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
double Jrad_buf[MAXNSCATT + 2][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
double dtau_avg[RAD_SCATT_TYPES + 1][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
double en_int_avg[RAD_SCATT_TYPES + 1][N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
grid_int_type     Nem, Nabs, Nsc;
grid_radtype_type Nem_phys, Nabs_phys, radtype_buf;

grid_int_type    Nsuper;
grid_double_type Esuper;
grid_prim_type   psupersave;

#if LOCAL_ANGULAR_DISTRIBUTIONS
grid_local_angles_type local_angles;
double local_dx1_rad, local_dx2_rad, local_dx_costh;
#if RAD_NUM_TYPES >= 4
grid_Gnu_type Gnu, local_Ns, local_wsqr;
grid_local_moment_type local_moments;
grid_local_basis_idx_type local_b_osc;
grid_local_count_type local_osc_count;
#endif //  RAD_NUM_TYPES
#endif // LOCAL_ANGULAR_DISTRIBUTIONS
#endif // RADIATION

#if ELECTRONS
grid_double_type Qvisc, Qcoul;
#endif // ELECTRONS
#if NVAR_PASSIVE > 0
int                passive_type[NVAR_PASSIVE];
char               passive_name[NVAR_PASSIVE][STRLEN];
passive_init_ftype do_passive_fixup[NVAR_PASSIVE];
#endif

double         conn[N1 + 2 * NG][N2 + 2 * NG][NDIM][NDIM][NDIM];
struct of_geom ggeom[N1 + 2 * NG][N2 + 2 * NG][NPG];
#if RADIATION
double                 dt_light[N1 + 2 * NG][N2 + 2 * NG], dt_light_min;
struct of_microphysics m_grd[N1 + 2 * NG][N2 + 2 * NG][N3 + 2 * NG];
grid_fourvector_type   Ucon_grd, Ucov_grd, Bcon_grd, Bcov_grd;
#endif

/*******************************************************************************
    GLOBAL VARIABLES
*******************************************************************************/
char outputdir[STRLEN], dumpdir[STRLEN], restartdir[STRLEN];
char xmfdir[STRLEN];
#if RADIATION && TRACERS
char tracerdir[STRLEN];
#endif
char init_from_grmhd[STRLEN];
char metric[STRLEN], reconstruction[STRLEN];
char eos[STRLEN], nulnutype[STRLEN];
int  tracers;

#if EOS == EOS_TYPE_TABLE
char eospath[STRLEN];
#endif

// opacity table paths
#if RADIATION
#if (RADIATION == RADTYPE_NEUTRINOS)
#if BURROWS_OPACITIES
char opac_param_file[STRLEN];
char opac_file[STRLEN];
#endif // BURROWS
#if HDF5_OPACITIES
char opac_file[STRLEN];
#endif // hdf5
#endif // neutrinos
#endif // radiation

#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
double gam;
#endif
#if EOS == EOS_TYPE_POLYTROPE
double poly_K, poly_gam;
#endif
#if POLYTROPE_FALLBACK
double rho_poly_thresh;
#endif

double a;
double M_unit;
double Reh;
double Risco;

#if NEED_UNITS
double mbh, Mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
#endif

#if EOS == EOS_TYPE_TABLE
double TEMP_unit;
#endif

#if RADIATION
double Ne_unit, Thetae_unit, kphys_to_num;
double tp_over_te, thetae_max, sigma_max, kdotk_tol;
#endif

double Rin, Rout, Rout_vis, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth;
#if RADIATION
double Rout_rad;
double nph_per_proc;
double tune_emiss, t_tune_emiss, dt_tune_emiss;
double tune_scatt, t_tune_scatt, dt_tune_scatt;
int    made_tune_proc, abs_tune_proc, scatt_tune_proc;
double numin, numax;
double kappa;
double startx_rad[NDIM], stopx_rad[NDIM];
double wgtC;
int    step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
int    tracer_tot, tracer_tot_all, ntcr_per_proc = 0;
int    step_sent, step_rcvd, step_fail;
int    step_made_all, step_abs_all, step_scatt_all, step_lost_all;
int    step_rec_all, step_sent_all, step_rcvd_all, step_tot_all, step_fail_all;
int    step_tot_max, step_tot_min;
double load_imbalance;
double Nph_to_track;
double sim_vol;
#if FLATEMISS
double cnu_flat;
#endif // MULTISCATT_TEST
#if MULTISCATT_TEST
double ms_theta_nu0, ms_delta0;
#endif // MULTISCATT_TEST
#if RZ_HISTOGRAMS
double rz_rmax, rz_zmax, delta_rcyl, delta_z;
#endif // RZ_HISTOGRAMS
#endif // RADIATION
#if ELECTRONS
double tptemin, tptemax;
#endif
double cour;
#if RADIATION
double cour_cool;
#endif
double dV, dx[NDIM], startx[NDIM], stopx[NDIM], startx_proc[NDIM],
    stopx_proc[NDIM];
double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
double dt, dtsave;
double t, tf;
double rcurr, hcurr;
int    istart, istop, jstart, jstop;
int    nstep;
int    is_restart;

double DTd;
double DTl;
double DTr;
int    DNr;
int    DTp;
int    DTf;
double DTw;
int    dump_cnt;
int    rdump_cnt;
double tdump, trestart, tlog;
double root_fcount[FCOUNT_NBINS];

int failed;
int lim;

double mdot = 0., mdot_eh = 0.;
double edot = 0., edot_eh = 0.;
double ldot = 0., ldot_eh = 0.;
int    icurr, jcurr, kcurr;

int nthreads;

#if ELECTRONS
double game, gamp;
double fel0;
#endif

int global_start[NDIM];
int global_stop[NDIM];
