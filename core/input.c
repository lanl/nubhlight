/******************************************************************************
 *                                                                            *
 * INPUT.C                                                                    *
 *                                                                            *
 * READ IN PARAMETER FILES                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <ctype.h>

struct param {
  char *key;
  void *data;
};

#define MAXPARAMS (1000)
static struct param table[MAXPARAMS];
static int          nparam    = 0;
static int          nparamset = 0;

void set_param(char *key, void *data) {
  table[nparam].key  = key;
  table[nparam].data = data;
  // printf("Will require param %d: %s\n",nparam, key); // debug
  nparam++;
}

int get_param(char *key, void **data) {
  int n = 0;
  while (strncmp(key, table[n].key, strlen(key)) != 0) {
    n++;
    if (n >= nparam) {
      *data = NULL;
      return 0;
    }
  }
  *data = table[n].data;

  return 1;
}

// char metric[STRLEN], reconstuction[STRLEN], eos[STRLEN];
void set_core_params() {
#if RECONSTRUCTION == LINEAR
  sprintf(reconstruction, "linear");
#elif RECONSTRUCTION == PPM
  sprintf(reconstruction, "ppm");
#elif RECONSTRUCTION == WENO
  sprintf(reconstruction, "weno");
#elif RECONSTRUCTION == MP5
  sprintf(reconstruction, "mp5");
#endif

  set_param("tf", &tf);
  set_param("dt", &dt);
#if METRIC == MINKOWSKI
  sprintf(metric, "MINKOWSKI");
  sprintf(nulnutype, "kstatistics");
  set_param("x1Min", &x1Min);
  set_param("x1Max", &x1Max);
  set_param("x2Min", &x2Min);
  set_param("x2Max", &x2Max);
  set_param("x3Min", &x3Min);
  set_param("x3Max", &x3Max);
#elif METRIC == MKS
  sprintf(metric, "MKS");
  sprintf(nulnutype, "camera");
  set_param("a", &a);
  set_param("hslope", &hslope);
  set_param("poly_xt", &poly_xt);
  set_param("poly_alpha", &poly_alpha);
  set_param("mks_smooth", &mks_smooth);
  if (N2 < NG)
    hslope = 1.;
  set_param("Rout", &Rout);
  set_param("Rout_vis", &Rout_vis);
#if RADIATION
  set_param("Rout_rad", &Rout_rad);
#endif // RADIATION
#endif

#if RADIATION && TRACERS
  tracers = 1;
#else
  tracers = 0;
#endif

  set_param("cour", &cour);
#if RADIATION
  set_param("cour_cool", &cour_cool);
#endif

// EOS
#if EOS == EOS_TYPE_GAMMA || GAMMA_FALLBACK
  set_param("gam", &gam);
#endif
#if EOS == EOS_TYPE_GAMMA
  sprintf(eos, "GAMMA");
#elif EOS == EOS_TYPE_POLYTROPE
  sprintf(eos, "POLYTROPE");
  set_param("poly_K", &poly_K);
  set_param("poly_gam", &poly_gam);
#elif EOS == EOS_TYPE_TABLE
  sprintf(eos, "TABLE");
  set_param("eospath", &eospath);
#endif // EOS

#if NEED_UNITS
#if METRIC == MINKOWSKI
  set_param("L_unit", &L_unit);
  set_param("M_unit", &M_unit);
#elif METRIC == MKS
  set_param("mbh", &mbh);
  set_param("M_unit", &M_unit);
#endif // METRIC
#endif // NEED_UNITS

#if ELECTRONS
  set_param("game", &game);
  set_param("gamp", &gamp);
  set_param("fel0", &fel0);
  set_param("tptemin", &tptemin);
  set_param("tptemax", &tptemax);
#endif

#if RADIATION
#if !ELECTRONS
  set_param("tp_over_te", &tp_over_te);
#endif // !ELECTRONS
  set_param("nph_per_proc", &nph_per_proc);
  set_param("numin", &numin);
  set_param("numax", &numax);
  set_param("tune_emiss", &tune_emiss);
  set_param("tune_scatt", &tune_scatt);
  set_param("t0_tune_emiss", &t_tune_emiss);
  set_param("t0_tune_scatt", &t_tune_scatt);
  set_param("thetae_max", &thetae_max);
  set_param("sigma_max", &sigma_max);
  set_param("kdotk_tol", &kdotk_tol);
  set_param("Nph_to_track", &Nph_to_track);
#if FLATEMISS
  set_param("cnu_flat", &cnu_flat);
#endif // FLATEMISS
#if MULTISCATT_TEST
  set_param("ms_theta_nu0", &ms_theta_nu0);
  set_param("ms_delta0", &ms_delta0);
#endif // MULTISCATT_TEST
#if RZ_HISTOGRAMS
  set_param("rz_rmax", &rz_rmax);
  set_param("rz_zmax", &rz_zmax);
#endif // RZ_HISTOGRAMS
#if (RADIATION == RADTYPE_NEUTRINOS)
#if BURROWS_OPACITIES
  set_param("opac_param_file", &opac_param_file);
  set_param("opac_file", &opac_file);
#endif // BURROWS_OPACITIES
#if HDF5_OPACITIES
  set_param("opac_file", &opac_file);
#endif // HDF5_OPACITIES
#endif // RADTYPE_NEUTRINOS
#endif // RADIATION

  set_param("init_from_grmhd", &init_from_grmhd);

  set_param("DTd", &DTd);
  set_param("DTl", &DTl);
  set_param("DTr", &DTr);
  set_param("DNr", &DNr);
  set_param("DTp", &DTp);
  set_param("DTf", &DTf);
  set_param("outputdir", &outputdir);
}

void read_params(char *pfname) {
  void *ptr;

  FILE *fp = fopen(pfname, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open parameter file: %s\n", pfname);
    exit(-1);
  }

  char line[STRLEN];
  while (fgets(line, STRLEN, fp)) {
    // Ignore comments, newlines, and leading whitespace
    if (line[0] == '#' || line[0] == '\n' || isspace(line[0]))
      continue;

    // Is key in dictionary, and is variable empty?
    char test[STRLEN], key[STRLEN];
    test[0] = '\0';
    sscanf(line, "%*s %s %*s %s", key, test);
    char *word = test;
    while (isspace(*word)) {
      word++;
    }
    if (word[0] == '\0') {
      continue;
    }

    // Read in parameter depending on datatype
    char type[6];
    strncpy(type, line, 5);
    type[5] = 0;
    if (get_param(key, &ptr)) {
      if (!strncmp(type, "[int]", 5)) {
        int buf;
        sscanf(line, "%*s %s %*s %d", key, &buf);
        *((int *)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[dbl]", 5)) {
        double buf;
        sscanf(line, "%*s %s %*s %lf", key, &buf);
        *((double *)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[str]", 5)) {
        char buf[STRLEN];
        sscanf(line, "%*s %s %*s %s", key, buf);
        strcpy((char *)ptr, buf);
        nparamset++;
      }
    }
  }

#if (METRIC == MKS) && NEED_UNITS
  Mbh = mbh * MSUN;
#endif

  if (nparamset != nparam && mpi_io_proc()) {
    fprintf(stderr, "Set %i parameters, needed %i!\n", nparamset, nparam);
    fprintf(stderr, "Params are:\n");
    for (int i = 0; i < nparam; i++) {
      fprintf(stderr, "\t%d\t%s\n", i, table[i].key);
    }
    exit(-1);
  }

  fclose(fp);

  if (mpi_io_proc())
    fprintf(stdout, "Parameter file read\n\n");
}

void init_params(char *pfname) {
  set_core_params();
  set_problem_params();
  read_params(pfname);
  // TODO: Not the best place for this?
#if RZ_HISTOGRAMS
  delta_rcyl = rz_rmax / RZ_HISTOGRAMS_N;
  delta_z = rz_zmax / RZ_HISTOGRAMS_N;
#endif // RZ_HISTOGRAMS
}
