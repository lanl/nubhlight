/******************************************************************************
 *                                                                            *
 * RECONSTRUCTION.C                                                           *
 *                                                                            *
 * RECONSTRUCTION ALGORITHMS                                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Performs the slope-limiting for the numerical flux calculation
double slope_lim(double y1, double y2, double y3) {
  double Dqm, Dqp, Dqc, s;

  // Limiter choice now hard-coded
  lim = MC;

  // Woodward, or monotonized central, slope limiter
  if (lim == MC) {
    Dqm = 2. * (y2 - y1);
    Dqp = 2. * (y3 - y2);
    Dqc = 0.5 * (y3 - y1);
    s   = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else {
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
        return (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
        return (Dqp);
      else
        return (Dqc);
    }
  }

  // van Leer slope limiter
  else if (lim == VANL) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s   = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else
      return (2. * s / (Dqm + Dqp));
  }

  // Minmod slope limiter (crude but robust)
  else if (lim == MINM) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s   = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else if (fabs(Dqm) < fabs(Dqp))
      return Dqm;
    else
      return Dqp;
  }

  fprintf(stderr, "unknown slope limiter\n");
  exit(10);

  return (0.);
}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout) {
  double Dqm, Dqp, Dqc, s;

  Dqm = 2. * (x2 - x1);
  Dqp = 2. * (x3 - x2);
  Dqc = 0.5 * (x3 - x1);

  s = Dqm * Dqp;

  if (s <= 0.)
    s = 0.;
  else {
    if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      s = Dqm;
    else if (fabs(Dqp) < fabs(Dqc))
      s = Dqp;
    else
      s = Dqc;
  }

  // Reconstruct left, right
  *lout = x2 - 0.5 * s;
  *rout = x2 + 0.5 * s;
}

// Colella & Sekora 2007 PPM extension to reduce clipping at extrema
void para_cs(double x1, double x2, double x3, double x4, double x5,
    double *lout, double *rout) {}

/*
#define PARA2LIM 0
void para2(double x1, double x2, double x3, double x4, double x5,
  double *lout, double *rout)
{
  int mm ;
  double y[5];
  y[0] = x1;
  y[1] = x2;
  y[2] = x3;
  y[3] = x4;
  y[4] = x5;
  double dq0[5];
  double dq[5];
  // PROBLEM! dq0 is uninitialized!
  for (int i = 0; i < 5; i++) dq[i] = dq0[i];
  double Dqc,l,r,qa, qd, qe;

#if PARA2LIM == 0 || PARA2LIM == 1
  double s, Dqm, Dqp, aDqm, aDqp , aDqc;
#endif

  // shifted dq
  //dq=dq0+2;

  // CW1.7
  for(mm=-1+2 ; mm<=1+2 ; mm++) {
    Dqc = 0.5 *(y[mm+1]-y[mm-1]);
#if PARA2LIM == 0 || PARA2LIM == 1
    Dqm = 2.0 *(y[mm]-y[mm-1]);
    Dqp = 2.0 *(y[mm+1]-y[mm]);
    aDqm = fabs(Dqm);
    aDqp = fabs(Dqp);
    aDqc = fabs(Dqc);
    s = Dqm*Dqp;
#endif

#if(PARA2LIM == 0)  // VANL
    double Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    double aDqvanl=fabs(Dqvanl);
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=MY_MIN(MY_MIN(aDqc,aDqvanl),MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
#elif(PARA2LIM == 1) // PMC
    if (s <=0.) dq[mm]=0.;       //CW1.8
    else dq[mm]=MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
#elif(PARA2LIM == 2) // MC
    dq[mm] =Dqc;
#endif
  }
  // CW1.6

  l=0.5*(y[0+2]+y[-1+2])-(dq[0+2]-dq[-1+2])/6.0;
  r=0.5*(y[1+2]+y[0+2])-(dq[1+2]-dq[0+2])/6.0;


  //  l=max(min(y[0],y[-1]),l);
  //  l=min(max(y[0],y[-1]),l);
  //  r=max(min(y[0],y[1]),r);
  //  r=min(max(y[0],y[1]),r);


  qa=(r-y[0+2])*(y[0+2]-l);
  qd=(r-l);
  qe=6.0*(y[0+2]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[0+2];
    r=y[0+2];
  }

  else if (qd*(qd-qe)<0.0) l=3.0*y[0+2]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[0+2]-2.0*l;


  *lout=l;   //a_L,j
  *rout=r;
  // *dw=r-l;                      //CW1.5
  // *w6=6.0*(y[0]-0.5*(l+r));
}
*/

// Parabolic interpolation (Colella & Woodward 1984). Implemented by F. Foucart
/*void para(double y0, double y1, double y2, double y3, double y4, double *lout,
double *rout)
{
  // Slopes
  double d0 = 2.*(y1-y0);
  double d1 = 2.*(y2-y1);
  double d2 = 2.*(y3-y2);
  double d3 = 2.*(y4-y3);
  double D1 = 0.5*(y2-y0);
  double D2 = 0.5*(y3-y1);
  double D3 = 0.5*(y4-y2);

  double condZeroSlope1 = (d1*d0<=0.);
  double sign1 = MY_SIGN(D1);//(D1>0.)*2.-1.;
  double DQ1 =
(1-condZeroSlope1)*sign1*MY_MIN(fabs(D1),MY_MIN(fabs(d0),fabs(d1))); double
condZeroSlope2 = (d2*d1<=0.); double sign2 = MY_SIGN(D2);//(D2>0.)*2.-1.; double
DQ2 = (1.-condZeroSlope2)*sign2*MY_MIN(fabs(D2),MY_MIN(fabs(d1),fabs(d2)));
  double condZeroSlope3 = (d3*d2<=0.);
  double sign3 = MY_SIGN(D3);//(D3>0.)*2.-1.;
  double DQ3 =
(1.-condZeroSlope3)*sign3*MY_MIN(fabs(D3),MY_MIN(fabs(d2),fabs(d3)));

  // Fiducial reconstruction
  double leftV = 0.5*(y2+y1)-1./6.*(DQ2-DQ1);
  double rightV = 0.5*(y3+y2)-1./6.*(DQ3-DQ2);

  // Corrections
  double corr1 = ((rightV-y2)*(y2-leftV)<=0.);
  double qd = rightV-leftV;
  double qe = 6.*(y2-0.5*(rightV+leftV));
  double corr2 = (qd*(qd-qe)<0.);
  double corr3 = (qd*(qd+qe)<0.);
  leftV = leftV*(1.-corr1)+corr1*y2;
  rightV = rightV*(1.-corr1)+corr1*y2;

  *lout = leftV*(1.-corr2)+corr2*(3.*y2-2.*rightV);
  *rout =
rightV*corr2+(1.-corr2)*rightV*(1.-corr3)+(1.-corr2)*corr3*(3.*y2-2.*leftV);
}*/

// Parabolic interpolation (see Colella & Woodward 1984; CW)
// Implemented by Xiaoyue Guan
/*void para(double x1, double x2, double x3, double x4, double x5, double *lout,
double *rout)
{
  double y[5], dq[5];
  double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  y[0] = x1;
  y[1] = x2;
  y[2] = x3;
  y[3] = x4;
  y[4] = x5;

  // CW 1.7
  for(int i = 1; i < 4; i++) {
    Dqm = 2. *(y[i  ] - y[i-1]);
    Dqp = 2. *(y[i+1] - y[i  ]);
    Dqc = 0.5*(y[i+1] - y[i-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    // CW 1.8
    if (s <= 0.) {
      dq[i] = 0.;
    } else {
      dq[i] = MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
    }
  }

  // CW 1.6
  l = 0.5*(y[2] + y[1]) - (dq[2] - dq[1])/6.0;
  r = 0.5*(y[3] + y[2]) - (dq[3] - dq[2])/6.0;

  qa = (r - y[2])*(y[2] - l);
  qd = (r - l);
  qe = 6.0*(y[2] - 0.5*(l + r));

  if (qa <= 0.) {
    l = y[2];
    r = y[2]; // r = y[3]? Mignone+ 2005
  }

  // CW 1984 limiter (1st order at extrema)
  if (qd*(qd - qe) < 0.0) {
    l = 3.0*y[2] - 2.0*r;
  } else if (qd*(qd + qe) < 0.0) {
    r = 3.0*y[2] - 2.0*l;
  }

  //if ((l-r)*(y[2] - (l+r)/2.) < -pow(l-r,2)/6.) {
  //  l = 3.0*y[2] - 2.0*r;
  //}

  //if ((l-r)*(y[2] - (l+r)/2.) >  pow(l-r,2)/6.) {
  //  r = 3.0*y[2] - 2.0*l;
  //}

  *lout = l;
  *rout = r;
}*/

// Colella & Sekora 2008 (CS08) extremum-preserving PPM
void para(double x1, double x2, double x3, double x4, double x5, double *lout,
    double *rout) {
  double C = 1.25;
  int    i = 2;
  double a[5], amin, amax, D2a, D2aL, D2aR, D2alim,
      s; //, D2ma, D2ca, D2pa;//alphamax;
  a[0] = x1;
  a[1] = x2;
  a[2] = x3;
  a[3] = x4;
  a[4] = x5;

  // CS08 Eqn. 16
  double aiph = 7. / 12. * (a[i] + a[i + 1]) - 1. / 12. * (a[i - 1] + a[i + 2]);
  double aimh = 7. / 12. * (a[i - 1] + a[i]) - 1. / 12. * (a[i - 2] + a[i + 1]);

  // Does aiph lie between a[i] and a[i+1]? CS08 Eqn. 13
  amin = MY_MIN(a[i], a[i + 1]);
  amax = MY_MAX(a[i], a[i + 1]);
  if (aiph <= amin || aiph >= amax) {
    // Extremum-preserving limiter
    D2a  = 3. * (a[i] - 2. * aiph + a[i + 1]);
    D2aL = a[i - 1] - 2. * a[i] + a[i + 1];
    D2aR = a[i] - 2. * a[i + 1] + a[i + 2];
    if (D2a * D2aL > 0. && D2a * D2aR > 0.) {
      s      = MY_SIGN(D2a);
      D2alim = s * MY_MIN(MY_MIN(C * fabs(D2aL), C * fabs(D2aR)), fabs(D2a));
    } else {
      D2alim = 0.;
    }
    aiph = 0.5 * (a[i] + a[i + 1]) - 1. / 3. * D2alim;
  }

  // Does aimh lie between a[i-1] and a[i]? CS08 Eqn. 13
  amin = MY_MIN(a[i - 1], a[i]);
  amax = MY_MAX(a[i - 1], a[i]);
  if (aimh < amin || aimh > amax) {
    // Extremum-preserving limiter
    D2a  = 3. * (a[i - 1] - 2. * aimh + a[i]);
    D2aL = a[i - 2] - 2. * a[i - 1] + a[i];
    D2aR = a[i - 1] - 2. * a[i] + a[i + 1];
    if (D2a * D2aL > 0. && D2a * D2aR > 0.) {
      s      = MY_SIGN(D2a);
      D2alim = s * MY_MIN(MY_MIN(C * fabs(D2aL), C * fabs(D2aR)), fabs(D2a));
    } else {
      D2alim = 0.;
    }
    aimh = 0.5 * (a[i - 1] + a[i]) - 1. / 3. * D2alim;
  }

  // Does a[i+1/2] lie between a[i] and a[i+1]? CS08 Eqn. 13
  /*amin = MY_MIN(a[i], a[i+1]);
  amax = MY_MAX(a[i], a[i+1]);
  if (aiph < amin || aiph > amax) {
    // Extremum-preserving limiter
    D2a  = 3.*(a[i] - 2.*aiph + a[i+1]);
    D2aL = a[i-1] - 2.*a[i] + a[i+1];
    D2aR = a[i] - 2.*a[i+1] + a[i+2];
    if (D2a*D2aL > 0. && D2a*D2aR > 0.) {
      s = MY_SIGN(D2a);
      D2alim = s*MY_MIN(MY_MIN(C*fabs(D2aL), C*fabs(D2aR)), fabs(D2a));
    } else {
      D2alim = 0.;
    }
    aiph = 0.5*(a[i] + a[i+1]) - 1./3.*D2alim;
  }*/

  /*if ((aiph - a[i])*(aimh - a[i]) > 0.) {
    aiph = a[i];
    aimh = a[i];
*/

  if ((a[i] - aimh) * (aiph - a[i]) <= 0. ||
      (a[i + 1] - a[i]) * (a[i] - a[i - 1]) <= 0.) {
    // aiph = a[i];
    // aimh = a[i];

    // THIS BIT OF CODE KILLS EVERYTHING
    D2a         = 6. * (aimh - 2. * a[i] + aiph);
    double D2aL = a[i - 2] - 2. * a[i - 1] + a[i];
    double D2aC = a[i - 1] - 2. * a[i] + a[i + 1];
    double D2aR = a[i] - 2. * a[i + 1] + a[i + 2];
    // fprintf(stdout, "%e %e %e %e\n", D2a, D2aL, D2aC, D2aR);
    double D2alim = 0.;
    if (MY_SIGN(D2a) == MY_SIGN(D2aL) && MY_SIGN(D2a) == MY_SIGN(D2aC) &&
        MY_SIGN(D2a) == MY_SIGN(D2aR)) {
      s      = MY_SIGN(D2a);
      D2alim = s * MY_MIN(MY_MIN(fabs(D2a), C * fabs(D2aL)),
                       MY_MIN(C * fabs(D2aC), C * fabs(D2aR)));
    }

    aiph = a[i] + (aiph - a[i]) * D2alim / D2a;
    aimh = a[i] + (aimh - a[i]) * D2alim / D2a;
  } else {
    if (fabs(aiph - a[i]) >= 2. * fabs(aimh - a[i])) {
      aiph = a[i] - 2. * (aimh - a[i]);
    }
    if (fabs(aimh - a[i]) >= 2. * fabs(aiph - a[i])) {
      aimh = a[i] - 2. * (aiph - a[i]);
    }
  }

  *lout = aimh;
  *rout = aiph;
}

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka
void weno(double x1, double x2, double x3, double x4, double x5, double *lout,
    double *rout) {
  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] = (3. / 8.) * x1 - (5. / 4.) * x2 + (15. / 8.) * x3;
  vr[1] = (-1. / 8.) * x2 + (3. / 4.) * x3 + (3. / 8.) * x4;
  vr[2] = (3. / 8.) * x3 + (3. / 4.) * x4 - (1. / 8.) * x5;

  vl[0] = (3. / 8.) * x5 - (5. / 4.) * x4 + (15. / 8.) * x3;
  vl[1] = (-1. / 8.) * x4 + (3. / 4.) * x3 + (3. / 8.) * x2;
  vl[2] = (3. / 8.) * x3 + (3. / 4.) * x2 - (1. / 8.) * x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13. / 12.) * pow(x1 - 2. * x2 + x3, 2) +
            (1. / 4.) * pow(x1 - 4. * x2 + 3. * x3, 2);
  beta[1] =
      (13. / 12.) * pow(x2 - 2. * x3 + x4, 2) + (1. / 4.) * pow(x4 - x2, 2);
  beta[2] = (13. / 12.) * pow(x3 - 2. * x4 + x5, 2) +
            (1. / 4.) * pow(x5 - 4. * x4 + 3. * x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;
  eps = 1.e-26;

  den = eps + beta[0];
  den *= den;
  wtr[0] = (1. / 16.) / den;
  den    = eps + beta[1];
  den *= den;
  wtr[1] = (5. / 8.) / den;
  den    = eps + beta[2];
  den *= den;
  wtr[2] = (5. / 16.) / den;
  Wr     = wtr[0] + wtr[1] + wtr[2];
  wr[0]  = wtr[0] / Wr;
  wr[1]  = wtr[1] / Wr;
  wr[2]  = wtr[2] / Wr;

  den = eps + beta[2];
  den *= den;
  wtl[0] = (1. / 16.) / den;
  den    = eps + beta[1];
  den *= den;
  wtl[1] = (5. / 8.) / den;
  den    = eps + beta[0];
  den *= den;
  wtl[2] = (5. / 16.) / den;
  Wl     = wtl[0] + wtl[1] + wtl[2];
  wl[0]  = wtl[0] / Wl;
  wl[1]  = wtl[1] / Wl;
  wl[2]  = wtl[2] / Wl;

  *lout = vl[0] * wl[0] + vl[1] * wl[1] + vl[2] * wl[2];
  *rout = vr[0] * wr[0] + vr[1] * wr[1] + vr[2] * wr[2];
}

// MP5 reconstruction from PLUTO
// Imported by Mani Chandra
#define MINMOD(a, b) ((a) * (b) > 0.0 ? (fabs(a) < fabs(b) ? (a) : (b)) : 0.0)
double Median(double a, double b, double c) {
  return (a + MINMOD(b - a, c - a));
}
double mp5_subcalc(
    double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2) {
  double        f, d2, d2p, d2m;
  double        dMMm, dMMp;
  double        scrh1, scrh2, Fmin, Fmax;
  double        fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f = 2.0 * Fjm2 - 13.0 * Fjm1 + 47.0 * Fj + 27.0 * Fjp1 - 3.0 * Fjp2;
  f /= 60.0;

  fMP = Fj + MINMOD(Fjp1 - Fj, alpha * (Fj - Fjm1));

  if ((f - Fj) * (f - fMP) <= epsm)
    return f;

  d2m = Fjm2 + Fj - 2.0 * Fjm1; // Eqn. 2.19
  d2  = Fjm1 + Fjp1 - 2.0 * Fj;
  d2p = Fj + Fjp2 - 2.0 * Fjp1; // Eqn. 2.19

  scrh1 = MINMOD(4.0 * d2 - d2p, 4.0 * d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1, scrh2); // Eqn. 2.27

  scrh1 = MINMOD(4.0 * d2m - d2, 4.0 * d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1, scrh2); // Eqn. 2.27

  fUL = Fj + alpha * (Fj - Fjm1);                   // Eqn. 2.8
  fAV = 0.5 * (Fj + Fjp1);                          // Eqn. 2.16
  fMD = fAV - 0.5 * dMMp;                           // Eqn. 2.28
  fLC = 0.5 * (3.0 * Fj - Fjm1) + 4.0 / 3.0 * dMMm; // Eqn. 2.29

  scrh1 = fmin(Fj, Fjp1);
  scrh1 = fmin(scrh1, fMD);
  scrh2 = fmin(Fj, fUL);
  scrh2 = fmin(scrh2, fLC);
  Fmin  = fmax(scrh1, scrh2); // Eqn. (2.24a)

  scrh1 = fmax(Fj, Fjp1);
  scrh1 = fmax(scrh1, fMD);
  scrh2 = fmax(Fj, fUL);
  scrh2 = fmax(scrh2, fLC);
  Fmax  = fmin(scrh1, scrh2); // Eqn. 2.24b

  f = Median(f, Fmin, Fmax); // Eqn. 2.26
  return f;
}

void mp5(double x1, double x2, double x3, double x4, double x5, double *lout,
    double *rout) {
  *rout = mp5_subcalc(x1, x2, x3, x4, x5);
  *lout = mp5_subcalc(x5, x4, x3, x2, x1);
}
#undef MINMOD

void reconstruct_lr_lin(double Ptmp[NMAX + 2 * NG][NVAR], int N,
    double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
  double dqtmp[NMAX + 2 * NG][NVAR];

  ISLOOP(-1, N) {
    PLOOP { dqtmp[i][ip] = 0.; }
  }

  // Calculate slopes
  ISLOOP(-1, N) {
    PLOOP {
      dqtmp[i][ip] = slope_lim(Ptmp[i - 1][ip], Ptmp[i][ip], Ptmp[i + 1][ip]);
    }
  }

  // Reconstruct left
  // ISLOOP(0, N) {
  ISLOOP(-1, N) {
    PLOOP { P_l[i][ip] = Ptmp[i][ip] - 0.5 * dqtmp[i][ip]; }
  }

  // Reconstruct right
  // ISLOOP(-1, N - 1) {
  ISLOOP(-1, N) {
    PLOOP { P_r[i][ip] = Ptmp[i][ip] + 0.5 * dqtmp[i][ip]; }
  }
}

void reconstruct_lr_par(double Ptmp[NMAX + 2 * NG][NVAR], int N,
    double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
  ISLOOP(-1, N) {
    PLOOP {
      para(Ptmp[i - 2][ip], Ptmp[i - 1][ip], Ptmp[i][ip], Ptmp[i + 1][ip],
          Ptmp[i + 2][ip], &P_l[i][ip], &P_r[i][ip]);
      /*para2(Ptmp[i-2][ip],
           Ptmp[i-1][ip],
           Ptmp[i][ip],
           Ptmp[i+1][ip],
           Ptmp[i+2][ip],
           &P_l[i][ip],
           &P_r[i][ip]);*/
    }
  }
}

// Extremum-preserving PPM -- implementation from athena++
double dph[NMAX + 2 * NG];
double d2qc[NMAX + 2 * NG];
double qminus[NMAX + 2 * NG];
double qplus[NMAX + 2 * NG];
double dqf_minus[NMAX + 2 * NG], dqf_plus[NMAX + 2 * NG], d2qf[NMAX + 2 * NG];
void   reconstruct_lr_ppm(double q[NMAX + 2 * NG][NVAR], int N,
      double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
  double C = 1.25;
  double qa, qb, qc, qd, qe, rho;

  PLOOP {
    // Reconstruct interface averages
    ISLOOP(-2, N + 1) {
      dph[i] =
          (7. * (q[i - 1][ip] + q[i][ip]) - (q[i + 1][ip] + q[i - 2][ip])) /
          12.;
      d2qc[i] = q[i - 1][ip] - 2. * q[i][ip] + q[i + 1][ip];
    }

    // Limit interpolated interface states
    ISLOOP(-1, N + 1) {
      double qmin = MY_MIN(q[i - 1][ip], q[i][ip]);
      double qmax = MY_MAX(q[i - 1][ip], q[i][ip]);
      /*if (dph[i] < qmin) dph[i] = qmin;
      if (dph[i] > qmax) dph[i] = qmax;
      continue;*/
      if (dph[i] < qmin || dph[i] > qmax) {
        qa = 3.0 * (q[i - 1][ip] - 2. * dph[i] + q[i][ip]);
        qb = d2qc[i - 1];
        qc = d2qc[i];
        qd = 0.;
        if (MY_SIGN(qa) == MY_SIGN(qb) && MY_SIGN(qa) == MY_SIGN(qc)) {
          qd = MY_SIGN(qa) *
               MY_MIN(C * fabs(qb), MY_MIN(C * fabs(qc), fabs(qa)));
        }
        dph[i] = 0.5 * (q[i - 1][ip] + q[i][ip]) - qd / 6.;
      }
      /*continue;

      qa = dph[i] - q[i-1][ip];
      qb = q[i][ip] - dph[i];
      if (qa*qb < 0.) { // Local extrema at i-1/2
        qa = 3.0*(q[i-1][ip] - 2.*dph[i] + q[i][ip]);
        qb = d2qc[i-1];
        qc = d2qc[i];
        qd = 0.;
        if (MY_SIGN(qa) == MY_SIGN(qb) && MY_SIGN(qa) == MY_SIGN(qc)) {
          qd = MY_SIGN(qa)*MY_MIN(C*fabs(qb), MY_MIN(C*fabs(qc), fabs(qa)));
        }
        dph[i] = 0.5*(q[i-1][ip] + q[i][ip]) - qd/6.;
      }*/
    }

    ISLOOP(-1, N) {
      qminus[i] = dph[i];
      qplus[i]  = dph[i + 1];
    }

    ISLOOP(-1, N) {
      dqf_minus[i] = q[i][ip] - qminus[i];
      dqf_plus[i]  = qplus[i] - q[i][ip];
      d2qf[i]      = 6. * (dph[i] - 2. * q[i][ip] + dph[i + 1]);
    }

    // CS08 limiters
    ISLOOP(-1, N) {
      qa = dqf_minus[i] * dqf_plus[i];
      qb = (q[i + 1][ip] - q[i][ip]) * (q[i][ip] - q[i - 1][ip]);

      if (qa <= 0. || qb <= 0.) {
        qa = d2qc[i - 1];
        qb = d2qc[i];
        qc = d2qc[i + 1];
        qd = d2qf[i];
        if (MY_SIGN(qa) == MY_SIGN(qb) && MY_SIGN(qa) == MY_SIGN(qc) &&
            MY_SIGN(qa) == MY_SIGN(qd)) {
          qe = MY_SIGN(qd) * MY_MIN(MY_MIN(C * fabs(qa), C * fabs(qb)),
                                 MY_MIN(C * fabs(qc), fabs(qd)));
        } else {
          qe = 0.;
        }

        // Is 2nd derivative close to roundoff error?
        qa = MY_MAX(fabs(q[i - 1][ip]), fabs(q[i - 2][ip]));
        qb = MY_MAX(
            fabs(q[i][ip]), MY_MAX(fabs(q[i + 1][ip]), fabs(q[i + 2][ip])));
        if (fabs(qd) <= (1.e-12) * MY_MAX(qa, qb)) {
          rho = 0.;
        } else {
          rho = qe / qd;
        }
        if (rho <= (1. - (1.e-12))) {
          qminus[i] = q[i][ip] - rho * dqf_minus[i];
          qplus[i]  = q[i][ip] + rho * dqf_plus[i];
        }
        // qminus[i] = q[i][ip];
        // qplus[i] = q[i][ip];
      } else {
        if (fabs(dqf_minus[i]) >= 2. * fabs(dqf_plus[i])) {
          qminus[i] = q[i][ip] - 2. * dqf_plus[i];
          // qminus[i] = q[i][ip];
        }
        if (fabs(dqf_plus[i]) >= 2. * fabs(dqf_minus[i])) {
          qplus[i] = q[i][ip] + 2. * dqf_minus[i];
          // qplus[i] = q[i][ip];
        }
      }
    }

    ISLOOP(-1, N) {
      P_l[i][ip] = qminus[i];
      P_r[i][ip] = qplus[i];
    }
  }
}

void reconstruct_lr_weno(double Ptmp[NMAX + 2 * NG][NVAR], int N,
    double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
  ISLOOP(-1, N) {
    PLOOP {
      weno(Ptmp[i - 2][ip], Ptmp[i - 1][ip], Ptmp[i][ip], Ptmp[i + 1][ip],
          Ptmp[i + 2][ip], &P_l[i][ip], &P_r[i][ip]);
    }
  }
}

void reconstruct_lr_mp5(double Ptmp[NMAX + 2 * NG][NVAR], int N,
    double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
  ISLOOP(-1, N) {
    PLOOP {
      mp5(Ptmp[i - 2][ip], Ptmp[i - 1][ip], Ptmp[i][ip], Ptmp[i + 1][ip],
          Ptmp[i + 2][ip], &P_l[i][ip], &P_r[i][ip]);
    }
  }
}

void reconstruct(double Ptmp[NMAX + 2 * NG][NVAR], int N,
    double P_l[NMAX + 2 * NG][NVAR], double P_r[NMAX + 2 * NG][NVAR]) {
#if RECONSTRUCTION == LINEAR
  reconstruct_lr_lin(Ptmp, N, P_l, P_r);
#elif RECONSTRUCTION == PPM
  fprintf(stderr, "ERROR PPM currently broken!\n");
  exit(-1);
  // reconstruct_lr_par(Ptmp, N, P_l, P_r);
  reconstruct_lr_ppm(Ptmp, N, P_l, P_r);
#elif RECONSTRUCTION == WENO
  reconstruct_lr_weno(Ptmp, N, P_l, P_r);
#elif RECONSTRUCTION == MP5
  reconstruct_lr_mp5(Ptmp, N, P_l, P_r);
#else
  fprintf(stderr, "Reconstruction algorithm not specified!\n");
  exit(-1);
#endif
}
