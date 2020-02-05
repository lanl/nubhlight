/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR TEST OF TABULATED EOS READER                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <stdio.h>
#define EXTRA_VARS  (14)
#define NRHO        (N1TOT)
#define NU          (N2TOT)
#define NYE         (N3TOT)
#define VARSTART    (YE_EM)
#define PENTH       (VARSTART+1)
#define PTEMP       (VARSTART+2)
#define TEMP_TRUE   (VARSTART+3)
#define TEMP_ERR    (VARSTART+4)
#define PRESS_TRUE  (VARSTART+5)
#define PRESS_U     (VARSTART+6)
#define PRESS_U_ERR (VARSTART+7)
#define PRESS_W     (VARSTART+8)
#define PRESS_W_ERR (VARSTART+9)
#define PCS         (VARSTART+10)
#define CS_TRUE     (VARSTART+11)
#define CS_ERR      (VARSTART+12)
#define UPRESS      (VARSTART+13)
#define UPRESS_ERR  (VARSTART+14)

#define MAX_PERCENT_ERROR     (1e-7)
// much higher because of interpolation errors
// log(cs2) is not linear at large temperatures
#define MAX_PERCENT_ERROR_CS2 (1e-3) 

static double lrhomin, lrhomax, lumin, lumax, yemin, yemax;
static char name[STRLEN], fname[STRLEN];
static FILE *errfile;
void set_problem_params()
{
  set_param("lrhomin", &lrhomin);
  set_param("lrhomax", &lrhomax);
  set_param("lumin", &lumin);
  set_param("lumax", &lumax);
  set_param("yemin", &yemin);
  set_param("yemax", &yemax);
}

void init_prob()
{

  sprintf(fname, "table_errors.txt");
  strcpy(name, dumpdir);
  strcat(name, fname);

  #if EOS != EOS_TYPE_TABLE
  fprintf(stderr,"This problem must be run with a tabulated EOS!\n");
  exit(1);
  #endif

  #if NVAR_PASSIVE != EXTRA_VARS + VARSTART - PASSIVE_START + 1
  fprintf(stderr,
    "This problem must be run with exactly %d passive variables!\n"
    "There are %d passive vars.\n",
    EXTRA_VARS + VARSTART - PASSIVE_START + 1,NVAR_PASSIVE);
  exit(1);
  #endif

  double X[NDIM];
  double U[NVAR];
  static double extra1[EOS_NUM_EXTRA];
  static double extra2[EOS_NUM_EXTRA];
  static double extra3[EOS_NUM_EXTRA];
  struct of_geom  *geom;
  struct of_state state;

  const double v0    = 0.0;
  //const double gamma = 0.0;
  const double b0    = 0.0;  

  ZLOOP { // fill arrays
    coord(i, j, k, CENT, X);
    geom = get_geometry(i, j, k, CENT);

    double lrho     = X[1];
    double rho      = pow(10., lrho);
    P[i][j][k][RHO] = rho;

    double lu       = X[2];
    double u        = pow(10., lu);
    P[i][j][k][UU]  = u;

    P[i][j][k][U1]  = v0;
    P[i][j][k][U2]  = v0;
    P[i][j][k][U3]  = v0;
    P[i][j][k][B1]  = b0;
    P[i][j][k][B2]  = b0;
    P[i][j][k][B3]  = b0;

    double ye         = X[3];
    P[i][j][k][YE]    = ye;
    P[i][j][k][YE_EM] = YE_EM;

    double enth       = EOS_Gamma_enthalpy_rho0_u(rho,u);
    P[i][j][k][PENTH] = enth;

    double temp_true  = EOS_Gamma_temp(rho,u);
    EOS_SC_fill(P[i][j][k], extra1);\
    // make temperature unitless
    double temp = EOS_temperature(rho,u,extra1)*TEMP_unit/(MP*CL*CL);
    P[i][j][k][PTEMP]       = temp;
    P[i][j][k][TEMP_TRUE]   = temp_true;
    P[i][j][k][TEMP_ERR]    = temp - temp_true;
    double eps = u/rho;
    double eps_cgs = eps*CL*CL;
    if ( fabs(P[i][j][k][TEMP_ERR]/temp_true) > MAX_PERCENT_ERROR ) {
      fprintf(stderr,
	      "ERROR IN TEMPERATURE:\n"
	      "TEMP_unit     = %e\n"
	      "\tT           = %e\n"
	      "\tT_TRUE      = %e\n"
	      "\tPercent_err = %e\n"
	      "\tT_CGS       = %e\n"
	      "\tlT_cgs      = %f\n"
	      "\tlT_true_cgs = %f\n"
	      "\tRHO         = %e\n"
	      "\tU           = %e\n"
	      "\tEPS         = %e\n"
	      "\te_cgs       = %e\n"
	      "\tle_cgs      = %f\n",
	      TEMP_unit,
	      P[i][j][k][PTEMP],
	      temp_true,
	      100.*(P[i][j][k][PTEMP]-temp_true)/temp_true,
	      temp*TEMP_unit,
	      log10(temp*TEMP_unit),
	      log10(temp_true*TEMP_unit),
	      rho,u,eps,eps_cgs,
	      log10(eps_cgs));
      exit(1);
    }

    double press_true = EOS_Gamma_pressure_rho0_u(rho,u);
    double press = EOS_pressure_rho0_u(rho,u,extra1);
    P[i][j][k][PRESS_TRUE]  = press_true;
    P[i][j][k][PRESS_U]     = press;
    P[i][j][k][PRESS_U_ERR] = press - press_true;
    if ( fabs(P[i][j][k][PRESS_U_ERR]/press_true) > MAX_PERCENT_ERROR ) {
      fprintf(stderr,
	      "ERROR IN PRESSURE:\n"
	      "\tP             = %e\n"
	      "\tP_TRUE        = %e\n"
	      "\tPercent Error = %e\n"
	      "\tT_TRUE        = %e\n"
	      "\tRHO           = %e\n"
	      "\tU             = %e\n"
	      "\tYE            = %e\n",
	      P[i][j][k][PRESS_U],
	      press_true,
	      100*fabs(P[i][j][k][PRESS_U_ERR]/press_true),
	      temp_true,
	      rho,u,ye);
      exit(1);
    }
    
    get_state(P[i][j][k], geom, &state);
    primtoflux(P[i][j][k], &state, 0, 0, geom, U);
    extra2[EOS_YE] = U[YE];
    
    P[i][j][k][PRESS_W] = EOS_SC_pressure_rho0_w(rho, enth, ye, &extra2[EOS_LT]);
    P[i][j][k][PRESS_W_ERR] = P[i][j][k][PRESS_W] - press_true;
    
    double cs_true = EOS_Gamma_sound_speed_rho0_u(rho,u);
    double ef_true = EOS_Gamma_enthalpy_rho0_u(rho,u);
    double cs = EOS_sound_speed_rho0_u(rho,u,extra1);
    P[i][j][k][PCS]            = cs;
    P[i][j][k][CS_ERR]         = cs - cs_true;
    P[i][j][k][CS_TRUE]        = cs_true;
    if (fabs(P[i][j][k][CS_ERR]/cs_true) > MAX_PERCENT_ERROR_CS2) {
      fprintf(stderr,
	      "ERROR IN CS:\n"
	      "U_UNIT           = %e\n"
	      "\tCS             = %e\n"
	      "\tCS_TRUE        = %e\n"
	      "\tCS_ERR         = %e\n"
	      "\tCS_PERCENT_ERR = %e\n"
	      "\tT              = %e\n"
	      "\tT_TRUE         = %e\n"
	      "\tT_CGS          = %e\n"
	      "\tlT_cgs         = %f\n"
	      "\tlT_true_cgs    = %f\n"
	      "\tRHO            = %e\n"
	      "\trho_cgs        = %e\n"
	      "\tlrho_cgs       = %e\n"
	      "\tU              = %e\n"
	      "\tU_cgs          = %e\n"
	      "\tEF             = %e\n"
	      "\tEF_cgs         = %e\n"
	      "\tP              = %e\n"
	      "\tP_cgs          = %e\n"
	      "\tEPS            = %e\n"
	      "\te_cgs          = %e\n"
	      "\tle_cgs         = %f\n",
	      U_unit,
	      P[i][j][k][PCS],
	      cs_true,
	      P[i][j][k][CS_ERR],
	      100*P[i][j][k][CS_ERR]/cs_true,
	      P[i][j][k][PTEMP],
	      temp_true,
	      temp*TEMP_unit,
	      log10(temp*TEMP_unit),
	      log10(temp_true*TEMP_unit),
	      rho,
	      rho*RHO_unit,
	      log10(rho*RHO_unit),
	      u,
	      u*U_unit,
	      ef_true,
	      ef_true*U_unit,
	      press_true,
	      press_true*U_unit,
	      eps,
	      eps_cgs,
	      log10(eps_cgs));      
      exit(1);
    }
    
    extra3[EOS_YE] = ye;
    P[i][j][k][UPRESS] = EOS_u_press(press_true, rho, extra3);
    P[i][j][k][UPRESS_ERR] = P[i][j][k][UPRESS] - u;
  }

  double lerr_pu = 0.;
  double err_pu  = 0.;
  double lerr_pw = 0.;
  double err_pw  = 0.;
  double lerr_cs = 0.;
  double err_cs  = 0.;
  double lerr_up = 0.;
  double err_up  = 0.;
  ZLOOP {
    if ( fabs(P[i][j][k][PRESS_U_ERR]) > lerr_pu ) {
      lerr_pu = fabs(P[i][j][k][PRESS_U_ERR]);
    }
    lerr_pw = MY_MAX(fabs(P[i][j][k][PRESS_W_ERR]),lerr_pw);
    lerr_cs = MY_MAX(fabs(P[i][j][k][CS_ERR]),lerr_cs);
    lerr_up = MY_MAX(fabs(P[i][j][k][UPRESS_ERR]),lerr_up);
  }
  MPI_Allreduce(&lerr_pu, &err_pu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&lerr_pw, &err_pw, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&lerr_cs, &err_cs, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&lerr_up, &err_up, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (mpi_io_proc()) {
    fprintf(stdout,
	    "Max Errors are:\n"
	    "\tPress(rho,u) = %e\n"
	    "\tPress(rho,w) = %e\n"
	    "\tCS           = %e\n"
	    "\tU(rho,P)     = %e\n",
	    err_pu, err_pw, err_cs, err_up);

    fprintf(stdout,"Saving to file: %s\n",name);
    errfile = fopen(name, "w");
    fprintf(errfile,"%e\n%e\n%e\n%e\n",
	    err_pu, err_pw, err_cs, err_up);
    fclose(errfile);
  }
}
