#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"
extern double lin_interp_between(double x, double x0, double x1, double y0, double y1);
extern double E2(double A, double ON_CURRENT);
extern double func_1_3w(double A);
extern double OCB(double A, double E2);
extern double ONE2(double A);
extern double ON(double ON_e2, double E2);
extern double OR(double A, double E2);
extern double OL(double A, double E2);
extern double det(double **a_in,int n);

/********************************************************************/
/*       DIFFERENTIAL EQUATIONS                                     */
/********************************************************************/
double Xcprime(double DC, double DN, double XC, double XN, double fnu, double Bfunc, double Afunc)
{
  return -Bfunc*(fnu*DN + (1.-fnu)*DC) - Afunc*XC;
}

double Xnprime(double DC, double DN, double XC, double XN, double fnu, double k, double Bfunc, double Afunc,double KJ2)
{
  if (M_nu==0.0) return -Bfunc*((1.-fnu)*DC) -Afunc*XN;
  return -Bfunc*((fnu - (k*k)/KJ2)*DN + (1.-fnu)*DC) -Afunc*XN;
}

double Dcprime(double XC)
{
  return XC;
}

double Dnprime(double XN)
{
  return XN;
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
void RK (int k_num, double *k, double *BETA, double *FC, double *FN)
{
  int i,j,ic_set,index_out;
  int rk_index_i,rk_index_j;
  int n_rk_step[4] = {0,0,0,0};

  double K[4];
  double L[4];
  double M[4];
  double N[4];
  double K_increment,L_increment,M_increment,N_increment;
  double dc_current,dn_current,xc_current,xn_current;

  double delta_c,delta_n,delta_m,f_m,f_c,f_n,x_c,x_n;
  double detM;
  double coeff[4];
  double t;

  double xcurrent;
  double a_rk;
  double E2_rk;
  double OM_rk;
  double OR_rk;
  double OL_rk;
  double ON_E2_rk;
  double A_rk;
  double B_rk;
  double FNU;
  double kj2_rk;
  double rk_1_3w;
  double rk_1_plus_3w;
  double OCB_rk;
  double ON_rk;
  double ff,gg,y;

  size_t start = time(NULL);
  size_t stop;

  double dc_0[4],dn_0[4],xc_0[4],xn_0[4];
  double dc_1[4],dn_1[4],xc_1[4],xn_1[4];
  double x_0[4],x_1[4];

  double Dc_k_z[k_num][output_number][4];
  double Dn_k_z[k_num][output_number][4];
  double Xc_k_z[k_num][output_number][4];
  double Xn_k_z[k_num][output_number][4];

  char out_files[200];
  FILE *o[output_number];

  for (index_out=0; index_out<output_number; index_out++)
  {
    sprintf(out_files,"%s_znum%i.txt",outputfile,index_out);
    o[index_out] = fopen(out_files,"w");
  }

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

  printf("Begin computation! \n");

  for (j = 0; j < k_num; j++)
  {
    count100++;
    if (count100 == 50)
    {
      printf("    %.0lf%% completed\n",((double)j/k_num)*100.);
      count100=0;
    }

    dc_1[0]  = 1.;
    xc_1[0]  = 0.;
    dn_1[0]  = 1.;
    xn_1[0]  = 0.;

    dc_1[1]  = 0.;
    xc_1[1]  = 1.;
    dn_1[1]  = 0.;
    xn_1[1]  = 1.;

    dc_1[2]  = 0.;
    xc_1[2]  = 1.;
    dn_1[2]  = 1.;
    xn_1[2]  = 0.;

    dc_1[3]  = 1.;
    xc_1[3]  = 0.;
    dn_1[3]  = 0.;
    xn_1[3]  = -1.;

    // Compute optimal stepsize for any given k
    if (M_nu==0.0)
    {
      step = (log(1.+z1)-log(1.+z0))/(20000-1);
    }
    else
    {
  		a_rk = 1./(1.+z_beta_in);
  		ON_E2_rk = ONE2(a_rk);
  		E2_rk = E2(a_rk,ON_E2_rk);
  		OCB_rk = OCB(a_rk,E2_rk);
  		ON_rk = ON(ON_E2_rk,E2_rk);
  		OR_rk = OR(a_rk,E2_rk);
  		OL_rk = OL(a_rk,E2_rk);
  		rk_1_3w = func_1_3w(a_rk);
  		B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
  		kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
  		FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
  		if (M_nu==0.) kj2_rk=1.;
      omega2 = B_rk*(FNU - (k[j]*k[j])/kj2_rk);

      if (omega2 <= 0.0) step = (log(1.+z1)-log(1.+z0))/(20000-1);
      else
      {
      	period = (2.0*M_PI)/sqrt(omega2);
      	step = period/n_step_per_period;
      	tmp_step = (log(1.+z1)-log(1.+z0))/(20000-1);
      	if (step > tmp_step) step = tmp_step;
      }
    }

    // Begin of actual RK computation
    for (ic_set=0; ic_set<4; ic_set++)
    {
      x_0[ic_set] = log(1.+z0);
      x_1[ic_set] = log(1.+z0);
      index_out=0;

      while(exp(x_1[ic_set])-1. < z1)
      {
        x_0[ic_set] = x_1[ic_set];
        dc_0[ic_set] = dc_1[ic_set];
        dn_0[ic_set] = dn_1[ic_set];
        xc_0[ic_set] = xc_1[ic_set];
        xn_0[ic_set] = xn_1[ic_set];

        for (rk_index_i=0; rk_index_i<4; rk_index_i++)
        {
          xcurrent = x_0[ic_set] + c_coeff[rk_index_i]*step;
          a_rk = 1./exp(xcurrent);
          ON_E2_rk = ONE2(a_rk);
          E2_rk = E2(a_rk,ON_E2_rk);
          OCB_rk = OCB(a_rk,E2_rk);
          ON_rk = ON(ON_E2_rk,E2_rk);
          OR_rk = OR(a_rk,E2_rk);
          OL_rk = OL(a_rk,E2_rk);
          rk_1_3w = func_1_3w(a_rk);
          rk_1_plus_3w = 2.0 - rk_1_3w;
          A_rk = 0.5*(OCB_rk+2.*OR_rk-2.*OL_rk+rk_1_plus_3w*ON_rk-2.);
          B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
          kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
          FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
          if (M_nu==0.) kj2_rk=1.;

          dc_current = dc_0[ic_set];
          dn_current = dn_0[ic_set];
          xc_current = xc_0[ic_set];
          xn_current = xn_0[ic_set];

          for (rk_index_j=0; rk_index_j<rk_index_i; rk_index_j++)
          {
            dc_current += a_coeff[rk_index_i][rk_index_j]*K[rk_index_j];
            dn_current += a_coeff[rk_index_i][rk_index_j]*L[rk_index_j];
            xc_current += a_coeff[rk_index_i][rk_index_j]*M[rk_index_j];
            xn_current += a_coeff[rk_index_i][rk_index_j]*N[rk_index_j];
          }

          K[rk_index_i] = (step) * Dcprime(xc_current);
          L[rk_index_i] = (step) * Dnprime(xn_current);
          M[rk_index_i] = (step) * Xcprime(dc_current,dn_current,xc_current,xn_current,FNU,B_rk,A_rk);
          N[rk_index_i] = (step) * Xnprime(dc_current,dn_current,xc_current,xn_current,FNU,k[j],B_rk,A_rk,kj2_rk);
        }

        K_increment=0.;
        L_increment=0.;
        M_increment=0.;
        N_increment=0.;

        for (rk_index_i=0; rk_index_i<4; rk_index_i++)
        {
          K_increment += b_coeff[rk_index_i]*K[rk_index_i];
          L_increment += b_coeff[rk_index_i]*L[rk_index_i];
          M_increment += b_coeff[rk_index_i]*M[rk_index_i];
          N_increment += b_coeff[rk_index_i]*N[rk_index_i];
        }

        dc_1[ic_set] = dc_0[ic_set] + K_increment ;
        dn_1[ic_set] = dn_0[ic_set] + L_increment ;
        xc_1[ic_set] = xc_0[ic_set] + M_increment ;
        xn_1[ic_set] = xn_0[ic_set] + N_increment ;

        x_1[ic_set] = x_0[ic_set] + step;

        if (exp(x_1[ic_set])-1.>z_output[index_out])
        {
          Dc_k_z[j][index_out][ic_set] =
            lin_interp_between(z_output[index_out],
                               exp(x_0[ic_set])-1.0,exp(x_1[ic_set])-1,
                               dc_0[ic_set],dc_1[ic_set]);
          Dn_k_z[j][index_out][ic_set] =
            lin_interp_between(z_output[index_out],
                               exp(x_0[ic_set])-1.0,exp(x_1[ic_set])-1,
                               dn_0[ic_set],dn_1[ic_set]);
          Xc_k_z[j][index_out][ic_set] =
            lin_interp_between(z_output[index_out],
                               exp(x_0[ic_set])-1.0,exp(x_1[ic_set])-1,
                               xc_0[ic_set],xc_1[ic_set]);
          Xn_k_z[j][index_out][ic_set] =
            lin_interp_between(z_output[index_out],
                               exp(x_0[ic_set])-1.0,exp(x_1[ic_set])-1,
                               xn_0[ic_set],xn_1[ic_set]);
          index_out++;
        }
        n_rk_step[ic_set]++;
      }
    }
  }

  for (j=0; j<k_num; j++)
  {
    a_rk=1.0;
    ON_E2_rk = ONE2(a_rk);
    E2_rk = E2(a_rk,ON_E2_rk);
    OCB_rk = OCB(a_rk,E2_rk);
    ON_rk = ON(ON_E2_rk,E2_rk);
    rk_1_3w = func_1_3w(a_rk);
    FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);

    matrix_M[0][0] = 1.;
    matrix_M[0][1] = 0.;
    matrix_M[0][2] = FNU;
    matrix_M[0][3] = 1.-FNU;

    for (ic_set=0; ic_set<4; ic_set++)
    {
      bin_beta_in = output_number-1;   //FIX THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      bin_fc_in = output_number-1;
      bin_fn_in = output_number-1;
      matrix_M[1][ic_set] =
        Xc_k_z[j][bin_fc_in][ic_set] + FC[j]*Dc_k_z[j][bin_fc_in][ic_set];
      matrix_M[2][ic_set] =
        BETA[j]*Dc_k_z[j][bin_beta_in][ic_set] - Dn_k_z[j][bin_beta_in][ic_set];
      matrix_M[3][ic_set] =
        Xn_k_z[j][bin_fn_in][ic_set] + FN[j]*Dn_k_z[j][bin_fn_in][ic_set];
    }

    detM = det(matrix_M,4);

    int jj;
    for(ic_set=0; ic_set<4; ic_set++)
    {
      for(i=0; i<4; i++)
      {
        for(jj=0; jj<4; jj++)
        {
          if (jj==ic_set)
          {
            if (i==0) matrix_M2[i][jj] = alpha;
            else matrix_M2[i][jj]=0.;
          }
          else matrix_M2[i][jj] = matrix_M[i][jj];
        }
      }
      coeff[ic_set]=det(matrix_M2,4)/detM;
    }

    for(index_out=0; index_out<output_number; index_out++)
    {
      delta_c=delta_n=x_c=x_n=0.;
      for (ic_set=0; ic_set<4; ic_set++)
      {
        delta_c += coeff[ic_set]*Dc_k_z[j][index_out][ic_set];
        delta_n += coeff[ic_set]*Dn_k_z[j][index_out][ic_set];
        x_c += coeff[ic_set]*Xc_k_z[j][index_out][ic_set];
        x_n += coeff[ic_set]*Xn_k_z[j][index_out][ic_set];
      }

      a_rk=1./(1.+z_output[index_out]);
      ON_E2_rk = ONE2(a_rk);
      E2_rk = E2(a_rk,ON_E2_rk);
      OCB_rk = OCB(a_rk,E2_rk);
      ON_rk = ON(ON_E2_rk,E2_rk);
      rk_1_3w = func_1_3w(a_rk);
      FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);

      delta_m = FNU*delta_n+(1.-FNU)*delta_c;
      f_c = -x_c/delta_c;
      f_n = -x_n/delta_n;
      f_m = -(FNU*x_n+(1.-FNU)*x_c)/delta_m;
      delta_c = delta_c;
      delta_n = delta_n;

      fprintf(o[index_out],
        "%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
        k[j],delta_c,delta_n,delta_m,f_c,f_n,f_m);
    }
  }

  stop = time(NULL);
  printf("   %.0lf%% completed, end of computation! It took %i s.\n",100.,(int)(stop-start));

  for (index_out=0; index_out<output_number; index_out++) fclose(o[index_out]);
}
