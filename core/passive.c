/******************************************************************************
 *                                                                            *
 * PASSIVE.C                                                                  *
 *                                                                            *
 * PASSIVE SCALARS                                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if NVAR_PASSIVE > 0

void fixup_passive(
    int i, int j, int k, double pv[NVAR], double pv_prefloor[NVAR]) {
  PASSLOOP { // advected numbers should not be changed by fixup
    if (PASSTYPE(ipass) == PASSTYPE_NUMBER) {
      pv[ipass] = pv_prefloor[ipass];
    }
    if (do_passive_fixup[PASSELEM(ipass)] != NULL) {
      do_passive_fixup[PASSELEM(ipass)](i, j, k, pv, pv_prefloor);
    }
  }
}

void __attribute__((weak)) init_passives() {
  PASSLOOP { do_passive_fixup[PASSELEM(ipass)] = NULL; }

#if EOS == EOS_TYPE_TABLE
  do_passive_fixup[PASSELEM(YE)]    = do_ye_fixup;
  do_passive_fixup[PASSELEM(YE_EM)] = do_ye_em_fixup;
#if METRIC == MKS
  do_passive_fixup[PASSELEM(ATM)] = do_atm_fixup;
#endif
#endif
}

void __attribute__((weak)) name_passives() {
  char name[STRLEN];
  PASSLOOP {
    sprintf(name, "var%d", PASSELEM(ipass));
    strcpy(PASSNAME(ipass), name);
  }
#if EOS == EOS_TYPE_TABLE
  // The first passive scalar is always Ye, the electron/proton fraction
  strcpy(PASSNAME(YE), "Ye");
  strcpy(PASSNAME(YE_EM), "Ye_em");
#if METRIC == MKS
  strcpy(PASSNAME(ATM), "ATM");
#endif // METRIC
#endif // EOS_TYPE_TABLE
}

#endif
