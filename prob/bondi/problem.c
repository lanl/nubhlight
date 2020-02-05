/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR BONDI INFLOW                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Local functions
void coord_transform(double *Pr, int i, int j);
double lfish_calc(double rmax) ;

void set_problem_params()
{
}

// Rootfinding for analytic Bondi solution
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
struct params {
  double r, mdot, n, C1, C2;
};
double rfunc(double T, void *params)
{
  struct params *p = (struct params *) params;
  //double r = p->r;
  double r = MY_MAX(2.5, p->r); // Solution breaks inside event horizon
  double mdot = p->mdot;
  double n = p->n;
  double C1 = p->C1;
  double C2 = p->C2;

  double resid = (pow(1. + (1. + n)*T,2)*(1. - 2.*mdot/r + pow(C1,2)/
         (pow(r,4)*pow(T,2*n))) - C2);

  //printf("r = %e T = %e resid = %e\n", r, T, resid);
  return resid;
  //return (pow(1. + (1. + n)*T,2)*(1. - 2.*mdot/r + pow(C1,2)/
  //       (pow(r,4)*pow(T,2*n))) - C2);
}

double C1, C2, T_bondi[N1+2*NG][N2+2*NG], n;

// Adapted from M. Chandra
double get_Tfunc(double T, double r)
{
  return pow(1.+(1.+n)*T,2.)*(1.-2./r+pow(C1/r/r/pow(T,n),2.))-C2;
}

double get_T(double r)
{
  double rtol = 1.e-12;
  double ftol = 1.e-14;
  double Tmin = 0.6*(sqrt(C2) - 1.)/(n + 1);
  double Tmax = pow(C1*sqrt(2./r/r/r),1./n);
  double f0, f1, fh;
  double T0, T1, Th;
  T0 = 0.6*Tmin;
  f0 = get_Tfunc(T0, r);
  T1 = Tmax;
  f1 = get_Tfunc(T1, r);

  if (f0*f1 > 0.) {
    printf("Failed solving for T at r = %e C1 = %e C2 = %e\n", r, C1, C2);
    exit(-1);
  }

  Th = (f1*T0 - f0*T1)/(f1 - f0);
  fh = get_Tfunc(Th, r);
  double epsT = rtol*(Tmin + Tmax);
  while (fabs(Th - T0) > epsT && fabs(Th - T1) > epsT && fabs(fh) > ftol) {
    if (fh*f0 < 0.) {
      T0 = Th;
      f0 = fh;
    } else {
      T1 = Th;
      f1 = fh;
    }

    Th = (f1*T0 - f0*T1)/(f1 - f0);
    fh = get_Tfunc(Th, r);
  }

  return Th;
}

void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM];
  DLOOP2 trans[mu][nu] = 0.;
  DLOOP1 trans[mu][mu] = 1.;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a);
  trans[3][1] = a/(r*r - 2.*r + a*a);

  DLOOP1 ucon_ks[mu] = 0.;
  DLOOP2 ucon_ks[mu] += trans[mu][nu]*ucon_bl[nu];
}

void fourvel_to_prim(double ucon[NDIM], double prim[NVAR],
  struct of_geom *geom)
{
  double alpha, beta[NDIM], gamma;
  alpha = 1.0/sqrt(-geom->gcon[0][0]);
  beta[1] = alpha*alpha*geom->gcon[0][1];
  beta[2] = alpha*alpha*geom->gcon[0][2];
  beta[3] = alpha*alpha*geom->gcon[0][3];
  gamma = ucon[0]*alpha;

  prim[U1] = ucon[1] + beta[1]*gamma/alpha;
  prim[U2] = ucon[2] + beta[2]*gamma/alpha;
  prim[U3] = ucon[3] + beta[3]*gamma/alpha;
}

void set_ut(double ucon[NDIM], struct of_geom *geom)
{
  double AA, BB, CC;

  AA = geom->gcov[0][0];
  BB = 2.*(geom->gcov[0][1]*ucon[1] +
           geom->gcov[0][2]*ucon[2] +
           geom->gcov[0][3]*ucon[3]);
  CC = 1. + geom->gcov[1][1]*ucon[1]*ucon[1] +
            geom->gcov[2][2]*ucon[2]*ucon[2] +
            geom->gcov[3][3]*ucon[3]*ucon[3] +
       2. *(geom->gcov[1][2]*ucon[1]*ucon[2] +
            geom->gcov[1][3]*ucon[1]*ucon[3] +
            geom->gcov[2][3]*ucon[2]*ucon[3]);

  double discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
}

void get_prim_bondi(int i, int j, int k, double P[NVAR])
{
  double r, th, X[NDIM];
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  while (r < Reh) {
    i++;
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
  }

  double T = T_bondi[i][j];
  double ur = -C1/(pow(T,n)*pow(r,2));
  double rho = pow(T,n);
  double u = rho*T/(gam - 1.);
  double ucon_bl[NDIM], ucon_ks[NDIM], ucon_mks[NDIM];
  struct of_geom geom_bl, *geom_mks;

  blgset(i, j, &geom_bl);

  DLOOP1 {
    ucon_bl[mu] = 0.;
    ucon_ks[mu] = 0.;
    ucon_mks[mu] = 0.;
  }
  ucon_bl[1] = ur;

  set_ut(ucon_bl, &geom_bl);
  bl_to_ks(X, ucon_bl, ucon_ks);

  double dxdX[NDIM][NDIM], dXdx[NDIM][NDIM];
  set_dxdX(X, dxdX);
  invert(&dxdX[0][0], &dXdx[0][0]);
  DLOOP2 {
    ucon_mks[mu] += dXdx[mu][nu]*ucon_ks[nu];
  }

  geom_mks = get_geometry(i, j, k, CENT);
  fourvel_to_prim(ucon_mks, P, geom_mks);
  //printf("[%i %i %i] r = %e ucons = %e %e %e vr = %e\n", 
  //  i, j, k, r, ucon_bl[1], ucon_ks[1], ucon_mks[1], P[U1]);

  P[RHO] = rho;
  P[UU] = u;
  P[B1] = 0.;
  P[B2] = 0.;
  P[B3] = 0.;

  /*if (j == N2/2){
    printf("\nu[1]: %e %e %e\n", ucon_bl[1], ucon_ks[1], ucon_mks[1]);
  printf("T rho ur u = %e %e %e %e\n", T, rho, ur, u);
  printf("P[] = %e %e %e %e\n", P[0], P[1], P[2], P[3]);
  }*/

  /*double gcov[NDIM][NDIM];
  bl_gcov_func(r, th, gcov);
  double ut = -sqrt((-1. - gcov[1][1]*ur*ur)/gcov[0][0]);
  double ucon[NDIM] = {ut, ur, 0, 0};



//  if ( j == NG)
//  printf("ucon[] = %e %e %e %e\n", ucon[0], ucon[1], ucon[2], ucon[3]);

  double trans[NDIM][NDIM];
  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  double tmp[NDIM];
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
    tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  set_dxdX(X, trans);
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
    tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  struct of_geom *geom = get_geometry(i, j, k, CENT);
  double alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  double gamma = ucon[0] * alpha;

  double beta[NDIM];
  beta[1] = alpha * alpha * geom->gcon[0][1];
  beta[2] = alpha * alpha * geom->gcon[0][2];
  beta[3] = alpha * alpha * geom->gcon[0][3];

  P[RHO] = rho;
  P[UU] = u;
  P[U1] = ucon[1] + beta[1] * gamma / alpha;
  P[U2] = 0.;
  P[U3] = 0.;
  P[B1] = 0.;
  P[B2] = 0.;
  P[B3] = 0.;
  if (r < 5.) P[U1] = 0.;

  //P[U1] *= 0.;*/
}

/*void get_bondi_prim(double r, double P[NDIM])
{
  double T, rho, ur, ut;
  if (r > 2. + 1.e-8) {
    T = get_T(r);
    rho = pow(T, n);
    ur = -C1/r/r/pow(T,n);
    ut = (2.*ur/r+sqrt(ur*ur+1.-2./r))/(1.-2./r);
  } else {
    T = 1.;
    rho = RHOMIN;
    ur = 0.;
    ut = 1.;
  }

  P[RHO] = rho;
  P[UU] = n*rho*T;
  P[U1] = (ur + 2.*ut/(r + 2.))/r;
  P[U2] = 0.;
  P[U3] = 0.;
  P[B1] = 0.;
  P[B2] = 0.;
  P[B3] = 0.;
}*/

void init_prob()
{
  double r, th, X[NDIM];

  double mdot = 1.;
  double rs = 8.;
  n = 1./(gam - 1.);

  // Solution constants
  double uc = sqrt(mdot/(2.*rs));
  double Vc = -sqrt(pow(uc,2)/(1. - 3.*pow(uc,2)));
  double Tc = -n*pow(Vc,2)/((n + 1.)*(n*pow(Vc,2) - 1.));
  C1 = uc*pow(rs,2)*pow(Tc,n);
  C2 = pow(1. + (1. + n)*Tc,2)*(1. - 2.*mdot/rs + pow(C1,2)/
       (pow(rs,4)*pow(Tc,2*n)));

  printf("a = %e Reh = %e\n", a, Reh);

  printf("mdot = %e\n", mdot);
  printf("rs   = %e\n", rs);
  printf("n    = %e\n", n);
  printf("uc   = %e\n", uc);
  printf("Vc   = %e\n", Vc);
  printf("Tc   = %e\n", Tc);
  printf("C1   = %e\n", C1);
  printf("C2   = %e\n", C2);

  ZSLOOP(-NG,N1+NG-1,-NG,N2+NG-1,-NG,N3+NG-1) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    T_bondi[i][j] = get_T(r);
  }

  ZLOOP {
    get_prim_bondi(i, j, k, P[i][j][k]);
  } // ZSLOOP
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(double *Pr, int ii, int jj)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct of_geom *geom, blgeom;

  coord(ii, jj, 0, CENT, X);
  bl_coord(X, &r, &th);
  blgset(ii, jj, &blgeom);

  ucon[1] = Pr[U1];
  ucon[2] = Pr[U2];
  ucon[3] = Pr[U3];

  AA = blgeom.gcov[0][0];
  BB = 2. * (blgeom.gcov[0][1] * ucon[1] +
       blgeom.gcov[0][2] * ucon[2] +
       blgeom.gcov[0][3] * ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1] * ucon[1] * ucon[1] +
      blgeom.gcov[2][2] * ucon[2] * ucon[2] +
      blgeom.gcov[3][3] * ucon[3] * ucon[3] +
      2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
      blgeom.gcov[1][3] * ucon[1] * ucon[3] +
      blgeom.gcov[2][3] * ucon[2] * ucon[3]);

  discr = BB * BB - 4. * AA * CC;
  ucon[0] = (-BB - sqrt(discr)) / (2. * AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  set_dxdX(X, trans);
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  //ucon[1] *= (1. / (r - R0));
  //ucon[1] /= dr_dx(X[1]);
  //ucon[2] *=
  //    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
  //ucon[2] /= dth_dx(X[2]);

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  geom = get_geometry(ii, jj, 0, CENT) ;
  alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  gamma = ucon[0] * alpha;

  beta[1] = alpha * alpha * geom->gcon[0][1];
  beta[2] = alpha * alpha * geom->gcon[0][2];
  beta[3] = alpha * alpha * geom->gcon[0][3];

  Pr[U1] = ucon[1] + beta[1] * gamma / alpha;
  Pr[U2] = ucon[2] + beta[2] * gamma / alpha;
  Pr[U3] = ucon[3] + beta[3] * gamma / alpha;
}

double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}

void bound_gas_prob_x1r(int i, int j, int k, grid_prim_type P)
{
  get_prim_bondi(i, j, k, P[i][j][k]);
}

