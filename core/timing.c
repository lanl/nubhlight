/******************************************************************************
 *                                                                            *
 * TIMING.C                                                                   *
 *                                                                            *
 * PERFORMANCE TIMING AND REPORTING                                           *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double timers[NUM_TIMERS];
static double times[NUM_TIMERS];
static int    nstep_start;

void time_set(int n, double val) { times[n] = val; }

double time_read(int n) { return times[n]; }

void time_init() {
  for (int n = 0; n < NUM_TIMERS; n++) {
    times[n] = 0.;
  }
  nstep_start = nstep;
}

void        timer_start(int timerCode) {
#pragma omp barrier

#pragma omp single
  {
    mpi_barrier();
    timers[timerCode] = omp_get_wtime();
  }
}

void        timer_stop(int timerCode) {
#pragma omp barrier

#pragma omp single
  {
    mpi_barrier();
    times[timerCode] += (omp_get_wtime() - timers[timerCode]);
  }
}

void timers_reset() {
  for (int n = 0; n < NUM_TIMERS; n++) {
    times[n] = 0.;
  }
}

double get_time_per_step(int timerCode) {
  int dnstep = nstep - nstep_start;
  return times[timerCode] / dnstep;
}

// Report a running average of timing data
void report_performance() {
  if (mpi_myrank() == 0) {
    int dnstep = nstep - nstep_start;
    fprintf(stdout, "\n********** PERFORMANCE **********\n");
    fprintf(stdout, "   UPDATE:   %8.4g s (%.4g %%)\n",
        times[TIMER_UPDATE] / dnstep,
        100. * times[TIMER_UPDATE] / times[TIMER_ALL]);
    fprintf(stdout, "   FLUXCALC: %8.4g s (%.4g %%)\n",
        times[TIMER_FLUXCALC] / dnstep,
        100. * times[TIMER_FLUXCALC] / times[TIMER_ALL]);
    fprintf(stdout, "   FIXUP:    %8.4g s (%.4g %%)\n",
        times[TIMER_FIXUP] / dnstep,
        100. * times[TIMER_FIXUP] / times[TIMER_ALL]);
    fprintf(stdout, "   BOUND:    %8.4g s (%.4g %%)\n",
        times[TIMER_BOUND] / dnstep,
        100. * times[TIMER_BOUND] / times[TIMER_ALL]);
    fprintf(stdout, "   DIAG:     %8.4g s (%.4g %%)\n",
        times[TIMER_DIAG] / dnstep,
        100. * times[TIMER_DIAG] / times[TIMER_ALL]);
    fprintf(stdout, "   OUTPUT:   %8.4g s (%.4g %%)\n",
        times[TIMER_OUT] / dnstep, 100. * times[TIMER_OUT] / times[TIMER_ALL]);
#if ELECTRONS
    fprintf(stdout, "   ELECTRON: %8.4g s (%.4g %%)\n",
        times[TIMER_ELECTRON] / dnstep,
        100. * times[TIMER_ELECTRON] / times[TIMER_ALL]);
#endif
#if RADIATION
    fprintf(stdout, "   MAKE:     %8.4g s (%.4g %%)\n",
        times[TIMER_MAKE] / dnstep,
        100. * times[TIMER_MAKE] / times[TIMER_ALL]);
    fprintf(stdout, "   PUSH:     %8.4g s (%.4g %%)\n",
        times[TIMER_PUSH] / dnstep,
        100. * times[TIMER_PUSH] / times[TIMER_ALL]);
    fprintf(stdout, "   INTERACT: %8.4g s (%.4g %%)\n",
        times[TIMER_INTERACT] / dnstep,
        100. * times[TIMER_INTERACT] / times[TIMER_ALL]);
#if NEUTRINO_OSCILLATIONS || LOCAL_ANGULAR_DISTRIBUTIONS
    fprintf(stdout, "   OSCILL:   %8.4g s (%.4g %%)\n",
        times[TIMER_OSCILLATIONS] / dnstep,
        100. * times[TIMER_OSCILLATIONS] / times[TIMER_ALL]);
#endif
    fprintf(stdout, "   MICRPHYS: %8.4g s (%.4g %%)\n",
        times[TIMER_MICRO] / dnstep,
        100. * times[TIMER_MICRO] / times[TIMER_ALL]);
#endif
    fprintf(stdout, "   ALL:      %8.4g s\n", times[TIMER_ALL] / dnstep);
    fprintf(stdout, "   ZONE CYCLES PER\n");
    fprintf(stdout, "     CORE-SECOND: %e\n",
        N1 * N2 * N3 / (times[TIMER_ALL] * nthreads / dnstep));
    fprintf(stdout, "     NODE-SECOND: %e\n",
        N1 * N2 * N3 / (times[TIMER_ALL] / dnstep));
    fprintf(stdout, "*********************************\n\n");
  }
}
