#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <omp.h>

#include "include_global.h"
#include "background.h"
#include "general_purpose.h"

/********************************************************************/
/*       DIFFERENTIAL EQUATIONS                                     */
/********************************************************************/
double Xbprime(double DB, double DC, double DN, double XB, double XC, double XN, double fnu, double Bfunc, double Afunc)
{
  double DCB = (OB0/(OB0+OC0))*DB + (OC0/(OB0+OC0))*DC;
  return -1.0*(-Bfunc*(fnu*DN + (1.-fnu)*DCB) - Afunc*XB); // !!!!
}

double Xcprime(double DB, double DC, double DN, double XB, double XC, double XN, double fnu, double Bfunc, double Afunc)
{
  double DCB = (OB0/(OB0+OC0))*DB + (OC0/(OB0+OC0))*DC;
  return -1.0*(-Bfunc*(fnu*DN + (1.-fnu)*DCB) - Afunc*XC);
}

double Xnprime(double DB, double DC, double DN, double XB, double XC, double XN, double fnu, double k, double Bfunc, double Afunc,double KJ2)
{
  double DCB = (OB0/(OB0+OC0))*DB + (OC0/(OB0+OC0))*DC;

  if (M_nu==0.0) return -1.0*(-Bfunc*((1.-fnu)*DCB) -Afunc*XN);
  return -1.0*(-Bfunc*((fnu - (k*k)/KJ2)*DN + (1.-fnu)*DCB) -Afunc*XN);
}

double Dbprime(double XB)
{
  return -1.0*XB;   // !!!!
}

double Dcprime(double XC)
{
  return -1.0*XC;
}

double Dnprime(double XN)
{
  return -1.0*XN;
}

/********************************************************************/
/*   RUNGE-KUTTA 4 MODULE                                           */
/*                                                                  */
/* Needs to be passed the (integer) total number of wavenumbers that*/
/* are required and a pointer to the vector of the wavenumbers      */
/* themselves. It initializes the integration to the values at z0   */
/* and procedes through nstep iterations up to z1. It writes the    */
/* values of k,Dc,Dn,fc,fn,Dm,fm at z=z1 in an output file.         */
/********************************************************************/
void RK (int k_num, double *k,
         double *BETAB, double *BETANU,
         double *FB, double *FC, double *FN,
         double **Delta_b, double **Delta_c, double **Delta_n, double **Delta_m,
         double **growth_b, double **growth_c, double **growth_n, double **growth_m)
{
  int j;
  size_t start = time(NULL);
  size_t stop;

  double Db_k_z[k_num][output_number];
  double Dc_k_z[k_num][output_number];
  double Dn_k_z[k_num][output_number];
  double Xb_k_z[k_num][output_number];
  double Xc_k_z[k_num][output_number];
  double Xn_k_z[k_num][output_number];

  printf("Begin computation! \n");

  #pragma omp parallel for shared(k_num,k,BETAB,BETANU,FB,FC,FN,OM0,OB0,OC0,OX0,OG0,OR0,h,M_nu,tau_reio,ns,As,kmax,N_nu,Neff,wrong_nu,z_final,z_initial,z_output,output_number,mode,alpha,Dc_k_z,Dn_k_z,Xc_k_z,Xn_k_z) schedule(dynamic)
  for (j = 0; j < k_num; j++)
  {
    double step;
    int index_out;
    int rk_index_i,rk_index_j;
    int n_rk_step = 0;

    double z0 = z_final;
    double z1 = z_initial;

    double K[4];
    double L[4];
    double M[4];
    double N[4];
    double O[4];
    double P[4];
    double K_increment,L_increment,M_increment,N_increment,O_increment,P_increment;
    double db_current,dc_current,dn_current,xb_current,xc_current,xn_current;

    double xcurrent;
    double a_rk;
    double E2_rk;
    double OR_rk;
    double OX_rk;
    double ON_E2_rk;
    double A_rk;
    double B_rk;
    double FNU;
    double kj2_rk;
    double rk_1_3w;
    double rk_1_plus_3w;
    double OCB_rk;
    double ON_rk;

    double db_0,dc_0,dn_0,xb_0,xc_0,xn_0;
    double db_1,dc_1,dn_1,xb_1,xc_1,xn_1;
    double x_0,x_1;

    //Fehlberg-Runge-Kutta with Dormand-Prince coefficients
    //(See Numerical Recipes page 912)

    double c_coeff[4] = {0.,1./2.,1./2.,1.};
    double b_coeff[4] = {1./6.,1./3.,1./3.,1./6.};
    double a_coeff[4][3] = { {0.,0.,0.},
                             {1./2.,0.,0.},
                             {0.,1./2.,0.},
                             {0.,0.,1.}      };

    int count100 = 0;

    int n_step_per_period = 50;
    double omega2,period,tmp_step;

    count100++;
    if (count100 == 50)
    {
      printf("    %.0lf%% completed\n",((double)j/k_num)*100.);
      count100=0;
    }
    // printf("%i from thread %i\n",j,omp_get_thread_num());

    // boundary conditions
    db_1  = BETAB[j]*1.;
    xb_1  = -FB[j]*db_1;
    dc_1  = 1.;
    xc_1  = -FC[j]*dc_1;
    dn_1  = BETANU[j]*1.;
    xn_1  = -FN[j]*dn_1;

    // Compute optimal stepsize for any given k
    if (M_nu==0.0)
    {
      step = (log(1./(1.+z0))-log(1./(1.+z1)))/(20000-1);
    }
    else
    {
  		a_rk = 1./(1.+z_initial);
  		ON_E2_rk = ONE2(a_rk);
  		E2_rk = E2(a_rk,ON_E2_rk);
  		OCB_rk = OCB(a_rk,E2_rk);
  		ON_rk = ON(ON_E2_rk,E2_rk);
  		OR_rk = OR(a_rk,E2_rk);
  		OX_rk = OX(a_rk,E2_rk);
  		rk_1_3w = func_1_3w(a_rk);
  		B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
  		kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
  		FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
  		if (M_nu==0.) kj2_rk=1.;
      omega2 = B_rk*(FNU - (k[j]*k[j])/kj2_rk);

      if (omega2 <= 0.0) step = (log(1./(1.+z0))-log(1./(1.+z1)))/(20000-1);
      else
      {
      	period = (2.0*M_PI)/sqrt(omega2);
      	step = period/n_step_per_period;
      	tmp_step = (log(1./(1.+z0))-log(1./(1.+z1)))/(20000-1);
      	if (fabs(step) > fabs(tmp_step)) step = tmp_step;
      }
    }

    // Begin of actual RK computation
    x_0 = log(1./(1.+z1));
    x_1 = log(1./(1.+z1));
    index_out=output_number-1;

    while(exp(x_1) < (1./(1.+z0)))
    {
      x_0 = x_1;
      db_0 = db_1;
      dc_0 = dc_1;
      dn_0 = dn_1;
      xb_0 = xb_1;
      xc_0 = xc_1;
      xn_0 = xn_1;

      for (rk_index_i=0; rk_index_i<4; rk_index_i++)
      {
        xcurrent = x_0 + c_coeff[rk_index_i]*step;
        a_rk = exp(xcurrent);
        // if (j==10) printf("j = %i, a = %lf\n",j,a_rk);
        ON_E2_rk = ONE2(a_rk);
        E2_rk = E2(a_rk,ON_E2_rk);
        OCB_rk = OCB(a_rk,E2_rk);
        ON_rk = ON(ON_E2_rk,E2_rk);
        OR_rk = OR(a_rk,E2_rk);
        OX_rk = OX(a_rk,E2_rk);
        rk_1_3w = func_1_3w(a_rk);
        rk_1_plus_3w = 2.0 - rk_1_3w;
        A_rk = A_func(a_rk,OR_rk,OCB_rk,OX_rk,rk_1_plus_3w,ON_rk,E2_rk);
        B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
        kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
        FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
        if (M_nu==0.) kj2_rk=1.;

        db_current = db_0;
        dc_current = dc_0;
        dn_current = dn_0;
        xb_current = xb_0;
        xc_current = xc_0;
        xn_current = xn_0;

        for (rk_index_j=0; rk_index_j<rk_index_i; rk_index_j++)
        {
          db_current += a_coeff[rk_index_i][rk_index_j]*K[rk_index_j];
          dc_current += a_coeff[rk_index_i][rk_index_j]*L[rk_index_j];
          dn_current += a_coeff[rk_index_i][rk_index_j]*M[rk_index_j];
          xb_current += a_coeff[rk_index_i][rk_index_j]*N[rk_index_j];
          xc_current += a_coeff[rk_index_i][rk_index_j]*O[rk_index_j];
          xn_current += a_coeff[rk_index_i][rk_index_j]*P[rk_index_j];
        }

        K[rk_index_i] = (step) * Dbprime(xb_current);
        L[rk_index_i] = (step) * Dcprime(xc_current);
        M[rk_index_i] = (step) * Dnprime(xn_current);
        N[rk_index_i] = (step) * Xbprime(db_current,dc_current,dn_current,xb_current,xc_current,xn_current,FNU,B_rk,A_rk);
        O[rk_index_i] = (step) * Xcprime(db_current,dc_current,dn_current,xb_current,xc_current,xn_current,FNU,B_rk,A_rk);
        P[rk_index_i] = (step) * Xnprime(db_current,dc_current,dn_current,xb_current,xc_current,xn_current,FNU,k[j],B_rk,A_rk,kj2_rk);
      }

      K_increment=0.;
      L_increment=0.;
      M_increment=0.;
      N_increment=0.;
      O_increment=0.;
      P_increment=0.;

      for (rk_index_i=0; rk_index_i<4; rk_index_i++)
      {
        K_increment += b_coeff[rk_index_i]*K[rk_index_i];
        L_increment += b_coeff[rk_index_i]*L[rk_index_i];
        M_increment += b_coeff[rk_index_i]*M[rk_index_i];
        N_increment += b_coeff[rk_index_i]*N[rk_index_i];
        O_increment += b_coeff[rk_index_i]*O[rk_index_i];
        P_increment += b_coeff[rk_index_i]*P[rk_index_i];
      }

      db_1 = db_0 + K_increment ;
      dc_1 = dc_0 + L_increment ;
      dn_1 = dn_0 + M_increment ;
      xb_1 = xb_0 + N_increment ;
      xc_1 = xc_0 + O_increment ;
      xn_1 = xn_0 + P_increment ;

      x_1 = x_0 + step;

      if (exp(x_1)>1./(1.+z_output[index_out]))
      {
        Db_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             db_0,db_1);
        Dc_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             dc_0,dc_1);
        Dn_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             dn_0,dn_1);
        Xb_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             xb_0,xb_1);
        Xc_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             xc_0,xc_1);
        Xn_k_z[j][index_out] =
          lin_interp_between(z_output[index_out],
                             1./exp(x_0)-1.0,1./exp(x_1)-1.0,
                             xn_0,xn_1);
        index_out--;
      }
      n_rk_step++;
    }
  }

  double delta_b,delta_c,delta_n,delta_m,f_m,f_b,f_c,f_n,x_b,x_c,x_n;
  int index_out;

  double a_rk;
  double E2_rk;
  double ON_E2_rk;
  double FNU;
  // double rk_1_3w;
  double OCB_rk;
  double ON_rk;

  double DCB,XCB,normalization;

  for (j=0; j<k_num; j++)
  {
    a_rk=1./(1.+z_output[0]);
    ON_E2_rk = ONE2(a_rk);
    E2_rk = E2(a_rk,ON_E2_rk);
    OCB_rk = OCB(a_rk,E2_rk);
    ON_rk = ON(ON_E2_rk,E2_rk);
    // rk_1_3w = func_1_3w(a_rk);
    // if (strcmp(boltzmann_code,"camb")==0) FNU = ON_rk/(OCB_rk+ON_rk);
    // else FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
    FNU = ON_rk/(OCB_rk+ON_rk);

    DCB = (OB0/(OB0+OC0))*Db_k_z[j][0] + (OC0/(OB0+OC0))*Dc_k_z[j][0];
    normalization = 1./((1.-FNU)*DCB + FNU*Dn_k_z[j][0]);

    for(index_out=0; index_out<output_number; index_out++)
    {
      delta_b = Db_k_z[j][index_out]*normalization;
      delta_c = Dc_k_z[j][index_out]*normalization;
      delta_n = Dn_k_z[j][index_out]*normalization;
      x_b = Xb_k_z[j][index_out];
      x_c = Xc_k_z[j][index_out];
      x_n = Xn_k_z[j][index_out];

      a_rk=1./(1.+z_output[index_out]);
      ON_E2_rk = ONE2(a_rk);
      E2_rk = E2(a_rk,ON_E2_rk);
      OCB_rk = OCB(a_rk,E2_rk);
      ON_rk = ON(ON_E2_rk,E2_rk);
      // rk_1_3w = func_1_3w(a_rk);
      // if (strcmp(boltzmann_code,"camb")==0) FNU = ON_rk/(OCB_rk+ON_rk);
      // else FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
      FNU = ON_rk/(OCB_rk+ON_rk);

      DCB = (OB0/(OB0+OC0))*delta_b + (OC0/(OB0+OC0))*delta_c;
      XCB = (OB0/(OB0+OC0))*x_b + (OC0/(OB0+OC0))*x_c;

      delta_m = FNU*delta_n+(1.-FNU)*DCB;
      f_b = -x_b/Db_k_z[j][index_out];
      f_c = -x_c/Dc_k_z[j][index_out];
      f_n = -x_n/Dn_k_z[j][index_out];
      f_m = -((1.-FNU)*XCB+FNU*x_n)/((1.-FNU)*((OB0/(OB0+OC0))*Db_k_z[j][index_out]+(OC0/(OB0+OC0))*Dc_k_z[j][index_out])+FNU*Dn_k_z[j][index_out]);

      Delta_b[index_out][j] = delta_b;
      Delta_c[index_out][j] = delta_c;
      Delta_n[index_out][j] = delta_n;
      Delta_m[index_out][j] = delta_m;
      growth_b[index_out][j] = f_b;
      growth_c[index_out][j] = f_c;
      growth_n[index_out][j] = f_n;
      growth_m[index_out][j] = f_m;
    }
  }

  stop = time(NULL);
  printf("   %.0lf%% completed, end of computation! It took %i s.\n",100.,(int)(stop-start));
}
