/******************************************************************************
 *                                                                            *
 * ROOT_FINDING.c                                                             *
 *                                                                            *
 * DEFINES HAND-TUNED ROUTIENS FOR ROOT_FINDING FOR PROBLEMS LIKE EOS TABLES  *
 *                                                                            *
 ******************************************************************************/

/* Implementation based on gsl root finder API */

#include "decs.h"

#define ROOT_DEBUG (0)
#define ROOT_VERBOSE (0)
#define ROOT_NAN_OK (0)
#define SECANT_NITER_MAX (10)

// PROTOTYPES
// ----------------------------------------------------------------------
static int root_secant(double (*f)(const double, const void *),
    const void *params, const double ytarget, const double xguess,
    const double xmin, const double xmax, const double xtol, const double ytol,
    double *xroot);
static int root_bisect(double (*f)(const double, const void *),
    const void *params, const double ytarget, const double xguess,
    const double xmin, const double xmax, const double xtol, const double ytol,
    double *xroot);
// ----------------------------------------------------------------------

// Implementation
// ----------------------------------------------------------------------
int find_root(double (*f)(const double, const void *), const void *params,
    const double ytarget, double xguess, const double xmin, const double xmax,
    const double xtol, const double ytol, double *xroot) {
  int status;

  // first check if we're at the max or min values
  const double fmax   = (*f)(xmax, params);
  const double errmax = fabs(fmax - ytarget) / (fabs(fmax) + ytol);
  if (errmax < ytol) {
    *xroot = xmax;
    root_fcount[0]++;
    return ROOT_SUCCESS;
  }
  const double fmin   = (*f)(xmin, params);
  const double errmin = fabs(fmin - ytarget) / (fabs(fmin) + ytol);
  if (errmin < ytol) {
    *xroot = xmin;
    root_fcount[0]++;
    return ROOT_SUCCESS;
  }

  if (xguess >= xmax)
    xguess = xmax - xtol;
  if (xguess < xmin)
    xguess = xmin;

  // Next try Secant
  status =
      root_secant(f, params, ytarget, xguess, xmin, xmax, xtol, ytol, xroot);
  if (status == ROOT_SUCCESS)
    return ROOT_SUCCESS;

#if ROOT_DEBUG
  if (isnan(*xroot)) {
    fprintf(stderr, "xroot is nan after secant\n");
  }
#endif

// Secant failed. Try bisection.
#if ROOT_VERBOSE
  fprintf(stderr,
      "\n\nRoot finding. Secant failed. Trying bisection.\n"
      "\txguess  = %.10g\n"
      "\tytarget = %.10g\n"
      "\txmin    = %.10g\n"
      "\txmax    = %.10g\n",
      xguess, ytarget, xmin, xmax);
#endif
  status =
      root_bisect(f, params, ytarget, xguess, xmin, xmax, xtol, ytol, xroot);

  // Check for something horrible happening
  if (isnan(*xroot) || isinf(*xroot)) {
#if ROOT_DEBUG
    fprintf(stderr, "xroot is nan after bisection\n");
#endif
    return ROOT_FAIL;
  }
  if (*xroot < xmin)
    return ROOT_FAIL;
  if (*xroot > xmax)
    return ROOT_FAIL;

  return status;
}

static int root_secant(double (*f)(const double, const void *),
    const void *params, const double ytarget, const double xguess,
    const double xmin, const double xmax, const double xtol, const double ytol,
    double *xroot) {
  double dx;
  double x_last, y, yp, ym, dyNum, dyDen, dy;

  double       x    = xguess;
  unsigned int iter = 0;
  do {
    x_last = x;
    dx     = fabs(1.e-7 * x) + xtol;
    y      = (*f)(x, params) - ytarget;
    yp     = (*f)(x + dx, params);
    ym     = (*f)(x - dx, params);
    dyNum  = yp - ym;
    dyDen  = (2. * dx);
    dy     = dyNum / dyDen;
    x -= y / dy;
    iter++;
    if (isnan(x) || isinf(x)) {
// can't recover from this
#if ROOT_DEBUG
      fprintf(stderr,
          "\n\n[root_secant]: NAN or out-of-bounds detected!\n"
          "\txguess  = %.10e\n"
          "\tytarget = %.10e\n"
          "\tx       = %.10e\n"
          "\tx_last  = %.10e\n"
          "\tx_min   = %.10e\n"
          "\tx_max   = %.10e\n"
          "\ty       = %.10e\n"
          "\tdx      = %.10e\n"
          "\typ      = %.10e\n"
          "\tym      = %.10e\n"
          "\tdyNum   = %.10e\n"
          "\tdyDen   = %.10e\n"
          "\tdy      = %.10e\n"
          "\titer    = %d\n"
          "\tsign x  = %d\n",
          xguess, ytarget, x, x_last, xmin, xmax, y, dx, yp, ym, dyNum, dyDen,
          dy, iter, (int)MY_SIGN(x));
#endif
#if ROOT_NAN_OK
      if (isinf(x)) {
        if (x < xmin)
          x = xmin;
        if (x > xmax)
          x = xmax;
      } else {
        root_fcount[FCOUNT_MORE]++;
        return ROOT_FAIL;
      }
#else
      root_fcount[FCOUNT_MORE]++;
      return ROOT_FAIL;
#endif
      if (x < xmin)
        x = xmin;
      if (x > xmax)
        x = xmax;
    }
  } while (
      iter < SECANT_NITER_MAX && fabs(x - x_last) / (fabs(x) + xtol) > xtol);

  if (iter < FCOUNT_NBINS)
    root_fcount[iter]++;
  else
    root_fcount[FCOUNT_MORE]++;

  *xroot = x;

  y                       = (*f)(x, params);
  const double frac_error = fabs(y - ytarget) / (fabs(y) + ytol);
#if ROOT_DEBUG
  if (frac_error > ytol) {
    fprintf(stderr,
        "\n\n[root_secant]: Failed via too large yerror.\n"
        "\tfractional error = %.10e\n"
        "\tx                = %.10e\n"
        "\ty                = %.10e\n"
        "\tytarget          = %.10e\n"
        "\typ               = %.10e\n"
        "\tym               = %.10e\n"
        "\tdy               = %.10e\n"
        "\tdx               = %.10e\n"
        "\titer             = %d\n",
        frac_error, x, y, ytarget, yp, ym, dy, dx, iter);
  }
  if (fabs(x - x_last) > xtol) {
    fprintf(stderr,
        "\n\n[root_secant]: failed via dx too big.\n"
        "\tfractional error = %.10e\n"
        "\tx                = %.10e\n"
        "\tx_last           = %.10e\n"
        "\tdx               = %.10e\n"
        "\ty                = %.10e\n"
        "\tytarget          = %.10e\n"
        "\typ               = %.10e\n"
        "\tym               = %.10e\n"
        "\tdy               = %.10e\n"
        "\tdx               = %.10e\n"
        "\titer             = %d\n",
        frac_error, x, x_last, fabs(x - x_last), y, ytarget, yp, ym, dy, dx,
        iter);
  }
#endif

  const int secant_failed =
      ((fabs(x - x_last) > xtol && fabs(frac_error) > ytol) || isnan(x) ||
          isinf(x));
  return secant_failed ? ROOT_FAIL : ROOT_SUCCESS;
}

static int root_bisect(double (*f)(const double, const void *),
    const void *params, const double ytarget, const double xguess,
    const double xmin, const double xmax, const double xtol, const double ytol,
    double *xroot) {
  double xl, xr, fl, fr, dx;

  double grow = 0.01;
  double x    = xguess;
  if (fabs(x) < xtol)
    x += 2. * xtol;
  do { // Try to find reasonable region for bisection
    dx = fabs(grow * x);
    xl = x - dx;
    xr = x + dx;
    fl = (*f)(xl, params) - ytarget;
    fr = (*f)(xr, params) - ytarget;
    grow *= 1.1;
  } while (fl * fr > 0 && xl >= xmin && xr <= xmax);

  // force back onto the bisection region
  if (xr > xmax) {
    xr = xmax;
    fr = (*f)(xr, params) - ytarget;
  }
  if (xl < xmin) {
    xl = xmin;
    fl = (*f)(xl, params) - ytarget;
  }

  // if they have the same sign, change that.
  // if we can't fix it, fail.
  if (fl * fr > 0) {
    xl = xmin;
    fl = (*f)(xl, params) - ytarget;
    if (fl * fr > 0) {
      xr = xmax;
      fr = (*f)(xr, params) - ytarget;
      if (fl * fr > 0) {
#if ROOT_DEBUG
        double il = (*f)(xl, params);
        double ir = (*f)(xr, params);
        fprintf(stderr,
            "\n\n[root_bisect]: fl*fr > 0!\n"
            "\txguess  = %.10e\n"
            "\tytarget = %.10e\n"
            "\txl      = %.10e\n"
            "\txr      = %.10e\n"
            "\tfl      = %.10e\n"
            "\tfr      = %.10e\n"
            "\til      = %.10e\n"
            "\tir      = %.10e\n",
            xguess, ytarget, xl, xr, fl, fr, il, ir);
        int    nx = 300;
        double dx = (xmax - xmin) / (nx - 1);
        fprintf(stderr, "Area map:\nx\ty\n");
        for (int i = 0; i < nx; i++) {
          fprintf(stderr, "%.4f\t%.4e\n", x + i * dx, (*f)(x + i * dx, params));
        }
#endif
        return ROOT_FAIL;
      }
    }
  }

  do { // bisection algorithm
    double xm = 0.5 * (xl + xr);
    double fm = (*f)(xm, params) - ytarget;
    if (fl * fm <= 0) {
      xr = xm;
      fr = fm;
    } else {
      xl = xm;
      fl = fm;
    }
  } while (xr - xl > xtol);

  *xroot = 0.5 * (xl + xr);

  if (isnan(*xroot)) {
#if ROOT_DEBUG
    double il = (*f)(xl, params);
    double ir = (*f)(xr, params);
    fprintf(stderr,
        "\n\n[root_bisect]: NAN DETECTED!\n"
        "\txguess  = %.10e\n"
        "\tytarget = %.10e\n"
        "\txl      = %.10e\n"
        "\txr      = %.10e\n"
        "\tdx      = %.10e\n"
        "\tgrow    = %.10e\n"
        "\txtol    = %.10e\n"
        "\tfl      = %.10e\n"
        "\tfr      = %.10e\n"
        "\til      = %.10e\n"
        "\tir      = %.10e\n"
        "\txmin    = %.10e\n"
        "\txmax    = %.10e\n",
        xguess, ytarget, xl, xr, dx, grow, xtol, fl, fr, il, ir, xmin, xmax);
#endif
    return ROOT_FAIL;
  }

  return ROOT_SUCCESS;
}

void initialize_root_fcounts() {
  for (int i = 0; i < FCOUNT_NBINS; i++)
    root_fcount[i] = 0.0;
}

void print_root_fcounts() {
  double fcount_tot = 0.0;
  double global_fcount[FCOUNT_NBINS];
  double fcount_percs[FCOUNT_NBINS];

  for (int i = 0; i < FCOUNT_NBINS; i++) {
    global_fcount[i] = mpi_reduce(root_fcount[i]);
    fcount_tot += global_fcount[i];
  }
  for (int i = 0; i < FCOUNT_NBINS; i++) {
    fcount_percs[i] = (100. * global_fcount[i]) / fcount_tot;
  }
  if (mpi_io_proc()) {
    fprintf(stdout, "\n********** ROOT FINDING *********\n");
    fprintf(stdout, "   ITERATIONS          PERCENTAGE\n");
    for (int i = 0; i < FCOUNT_NBINS; i++) {
      if (i == FCOUNT_NBINS - 1) {
        fprintf(stdout, "         more          %.2e %%\n", fcount_percs[i]);
      } else {
        fprintf(stdout, "          %3d          %.2e %%\n", i, fcount_percs[i]);
      }
    }
    fprintf(stdout, "*********************************\n\n");
  }
  // reset counts
  initialize_root_fcounts();
}
// ----------------------------------------------------------------------
