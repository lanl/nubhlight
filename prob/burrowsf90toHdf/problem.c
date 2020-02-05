/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR TEST OF TABULATED EOS READER                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static char           fname[STRLEN];
struct of_tablebounds b;

void set_problem_params() {
  set_param("outfile", &fname);
  set_param("lrho_min", &(b.lrho_min));
  set_param("lrho_max", &(b.lrho_max));
  set_param("numrho_out", &(b.Nrho));
  set_param("lT_min", &(b.lT_min));
  set_param("lT_max", &(b.lT_max));
  set_param("numT_out", &(b.NT));
  set_param("Ye_min", &(b.Ye_min));
  set_param("Ye_max", &(b.Ye_max));
  set_param("numYe_out", &(b.NYe));
}

void init_prob() {

#if EOS != EOS_TYPE_TABLE
  fprintf(stderr, "This problem must be run with a tabulated EOS!\n");
  exit(1);
#endif

#if RADIATION != RADTYPE_NEUTRINOS
  fprintf(stderr, "This problem must be run with with neutrinos active.\n");
  exit(1);
#endif

#ifndef BURROWS_OPACITIES
  fprintf(stderr, "This problem must be run with Burrows opacities.\n");
  exit(1);
#endif

#if !BURROWS_OPACITIES
  fprintf(stderr, "This problem must be run with Burrows opacities.\n");
  exit(1);
#endif

#if NVAR_PASSIVE < 2
  fprintf(stderr,
      "This problem must be run with more than %d passive variables!\n"
      "There are %d passive vars.\n",
      2, NVAR_PASSIVE);
  exit(1);
#endif

  opac_emis_to_hdf(fname, &b);
}
