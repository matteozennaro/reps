#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

/*INITIAL SETTINGS *************************************************/

  char input_file[200];
  char outputfile[200];
  char class_folder[500];
  char ic_mode[200];
  char normfile[200];
  char do_rescaled_ps;
  char neutrino_pressure;
  char fixed_neutrino_frac;
  char only_backg_rel;
  char print_hubble;
  char generate_boundary;
  char fixed_f_in_bc;

  int    nstep;
  double step,z0,z1;

  double z_beta_in,z_fc_in,z_fn_in;

  double *z_output;

  int output_number;
  int bin_beta_in,bin_fc_in,bin_fn_in;
  int bin_out;
  int mode = 0;

  // int nreallocs = 0;
  double **matrix_M;
  double **matrix_M2;
  // double **dc;
  // double **dn;
  // double **xc;
  // double **xn;

  double H0 = 100.;
  double OM0,OB0,OC0,OL0,OG0,OR0,h,M_nu,tau_reio,ns,As,kmax,N_nu,Neff;

  double Tcmb_0 = 2.7255;      //K
  double Kb = 8.617342e-5;     //eV K^-1
  double Gamma_nu = 0.71611;   //Tnu0/Tcmb0

  double *ytab;
  double *FFtab;
  double *GGtab;
  double ytab_max;
  double ytab_min;
  double F_inf,F_0;
  double G_inf,G_0;
  double ytabstep;
  double ps_norm_at_z;
  int ntab;

  //Hard coded boundary condition: the amplitude of the total matter
  //growth factor at redshift z0 (which is usually z0=0) is fixed to
  //be 1.
  double alpha = 1.;

/********************************************************************/
/*    DECLARATION OF FUNCTIONS                                      */
/********************************************************************/
void fscanf_error(int n);
double det(double **a_in,int n);
double **allocate_matrix(int rows, int columns);
void deallocate_matrix(double **m, int rows, int columns);
int find_z_bin(double Z, double *vec, int n_elems);
void read_parameter_file(char parfile[]);
void RK (int k_num, double *k, double *BETA, double *FC, double *FN);
void rescale_ps(int knum, double *k);
void create_class_ini_file (char dir_chain[]);
void read_D(char filename[],int n, double *dc, double *dn, double *dm);
void read_GG_FF_tabs();
double FF(double Y);
double GG(double Y);
void print_hubble_table();
void reallocate_matrix(double **m,int rows,int newcols);
void print_kfs();

/********************************************************************/
/*    BACKGROUND                                                    */
/********************************************************************/
double func_1_3w(double A)
{
  if (neutrino_pressure == 'F')
    return 1.;
  else
  {
    if (M_nu==0.) return 1.;
    double y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
    return y*(GG(y)/FF(y));
  }
}

double OCB(double A, double E2)
{
  return ((OC0+OB0)/(A*A*A))/E2;
}

double ONE2(double A)
{
  if (neutrino_pressure == 'T')
  {
    if (M_nu==0.) return 0.;
    double Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
    double pi4 = M_PI*M_PI*M_PI*M_PI;
    double a4 = A*A*A*A;
    double y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);  //if(A==1.) printf("%e %e \n",y,FF(y));   //FILE *yyy = fopen("yvalues.dat","a"); fprintf(yyy,"%e\n",y); fclose(yyy);
    return ((15.*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
  }
  else
  return (M_nu/(93.14*h*h))/(A*A*A);
}

double ON(double ON_e2, double E2)
{
  return ON_e2/E2;
}

double OR(double A, double E2)
{
    return (OR0/(A*A*A*A))/E2;
}

double OL(double A, double E2)
{
  return (OL0)/E2;
}

#include "boundary_conditions_module/BC.h"

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

/******************************************************************************/
/*        MAIN                                                                */
/******************************************************************************/

int main(int argc, char *argv[])
{
  read_GG_FF_tabs();

  read_parameter_file(argv[1]);

  if (generate_boundary=='T')
  generate_boundary_conditions();

  step  = (log(1.+z1)-log(1.+z0))/(nstep-1.);

  int lines = 0;
  char test = 'a';
  int j;
  int fscanfcheck;

  FILE *fin = fopen(input_file,"r");
    if (fin == NULL) {printf("Problem loading file...\n"); exit(1);}
  while (fscanf(fin,"%c",&test)!=EOF)
    if (test == '\n') lines++;

  double *k_requested = malloc(lines*sizeof(double));
  double *fc_init = malloc(lines*sizeof(double));
  double *fn_init = malloc(lines*sizeof(double));
  double *beta = malloc(lines*sizeof(double));
    if (k_requested == NULL) {printf("Bad memory allocation!\n"); exit(1);}
    if (fc_init == NULL)     {printf("Bad memory allocation!\n"); exit(1);}
    if (fn_init == NULL)     {printf("Bad memory allocation!\n"); exit(1);}
    if (beta == NULL)        {printf("Bad memory allocation!\n"); exit(1);}

  fseek(fin,0,SEEK_SET);

  for (j = 0; j < lines; j++)
  {
    fscanfcheck=fscanf(fin,"%lf %lf %lf %lf",&k_requested[j],&beta[j],&fc_init[j],&fn_init[j]);
    if (fscanfcheck!=4) fscanf_error(4);
  }
  fclose(fin);

  //Allocating the two 4x4 matrixes
  matrix_M = allocate_matrix(4,4);
  matrix_M2 = allocate_matrix(4,4);

  // dc = allocate_matrix(4,1e6+nreallocs*5e5);
  // dn = allocate_matrix(4,1e6+nreallocs*5e5);
  // xc = allocate_matrix(4,1e6+nreallocs*5e5);
  // xn = allocate_matrix(4,1e6+nreallocs*5e5);

  RK(lines,k_requested,beta,fc_init,fn_init);

  if(do_rescaled_ps=='T')
    rescale_ps(lines,k_requested);
  else
    printf("\nRescaling was not requested. Skip.\n");

  if(print_hubble=='T')
    print_hubble_table();
  else
    printf("\nPrinting the table of H values was not requested. Skip.\n");

  //print_hubble_table();
  print_kfs();
/*Freeing allocated memory*************************************************/
  free(k_requested);
  free(beta);
  free(fc_init);
  free(fn_init);
  free(z_output);
  free(ytab);
  free(GGtab);
  free(FFtab);

  deallocate_matrix(matrix_M,4,4);
  deallocate_matrix(matrix_M2,4,4);

  // deallocate_matrix(dc,4,1e6+nreallocs*1e5);
  // deallocate_matrix(dn,4,1e6+nreallocs*1e5);
  // deallocate_matrix(xc,4,1e6+nreallocs*1e5);
  // deallocate_matrix(xn,4,1e6+nreallocs*1e5);

/*Returning null value*****************************************************/
  return 0;
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
  double FNU = (M_nu/(93.14*h*h))/OM0;

  int i,j,ic_set,index_out;
  int rk_index_i,rk_index_j;
  int n_rk_step[4] = {0,0,0,0};

  double K[4];
  double L[4];
  double M[4];
  double N[4];
  double K_increment,L_increment,M_increment,N_increment;
  double dc_current,dn_current,xc_current,xn_current;
  double DeltaK,DeltaL,DeltaM,DeltaN;
  double scaleK,scaleL,scaleM,scaleN;
  double err,err_i[4];
  double abs_tol = 1.e-4;
  double rel_tol = 1.e-4;
  double safety_factor = 0.95;
  double **x;

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

  //if (x==NULL) {printf("Bad memory allocation!\n"); exit(1);}

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
		a_rk = 1./(1.+z_beta_in);
		ON_E2_rk = ONE2(a_rk);
		E2_rk = OR0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ON_E2_rk;
		OCB_rk = OCB(a_rk,E2_rk);
		if (only_backg_rel == 'T') ON_rk=((M_nu/(93.14*h*h))/(a_rk*a_rk*a_rk))/E2_rk;
		else ON_rk = ON(ON_E2_rk,E2_rk);
		OR_rk = OR(a_rk,E2_rk);
		OL_rk = OL(a_rk,E2_rk);
		rk_1_3w = func_1_3w(a_rk);
		if (only_backg_rel == 'T') B_rk = -3.*(OCB_rk+ON_rk)/2.;
		else B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
		kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
		if (fixed_neutrino_frac=='F') FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
		if (only_backg_rel == 'T') FNU = ON_rk/(OCB_rk+rk_1_3w*ON_rk);
		if (M_nu==0.) kj2_rk=1.;
    omega2 = B_rk*(FNU - (k[j]*k[j])/kj2_rk);

    if (omega2 <= 0.0) step = (log(1.+z1)-log(1.+z0))/(20000-1);
    else
    {
    	period = (2.0*M_PI)/sqrt(omega2);
    	step = period/n_step_per_period;
    	tmp_step = (log(1.+z1)-log(1.+z0))/(20000-1);
    	if (step > tmp_step) step = tmp_step;
    	if (M_nu==0.0) step=tmp_step;
    }

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
          E2_rk = OR0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ON_E2_rk;
          OCB_rk = OCB(a_rk,E2_rk);
          if (only_backg_rel == 'T') ON_rk=((M_nu/(93.14*h*h))/(a_rk*a_rk*a_rk))/E2_rk;
          else ON_rk = ON(ON_E2_rk,E2_rk);
          OR_rk = OR(a_rk,E2_rk);
          OL_rk = OL(a_rk,E2_rk);
          rk_1_3w = func_1_3w(a_rk);
          rk_1_plus_3w = 2.0 - rk_1_3w;
          A_rk = 0.5*(OCB_rk+2.*OR_rk-2.*OL_rk+rk_1_plus_3w*ON_rk-2.);
          if (only_backg_rel == 'T') B_rk = -3.*(OCB_rk+ON_rk)/2.;
          else B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
          kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
          if (fixed_neutrino_frac=='F') FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
          if (only_backg_rel == 'T') FNU = ON_rk/(OCB_rk+rk_1_3w*ON_rk);
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
          Dc_k_z[j][index_out][ic_set] = dc_1[ic_set] + ((dc_1[ic_set]-dc_0[ic_set])/(exp(x_1[ic_set])-exp(x_0[ic_set])))*(z_output[index_out]-(exp(x_1[ic_set])-1.));
          Dn_k_z[j][index_out][ic_set] = dn_1[ic_set] + ((dn_1[ic_set]-dn_0[ic_set])/(exp(x_1[ic_set])-exp(x_0[ic_set])))*(z_output[index_out]-(exp(x_1[ic_set])-1.));
          Xc_k_z[j][index_out][ic_set] = xc_1[ic_set] + ((xc_1[ic_set]-xc_0[ic_set])/(exp(x_1[ic_set])-exp(x_0[ic_set])))*(z_output[index_out]-(exp(x_1[ic_set])-1.));
          Xn_k_z[j][index_out][ic_set] = xn_1[ic_set] + ((xn_1[ic_set]-xn_0[ic_set])/(exp(x_1[ic_set])-exp(x_0[ic_set])))*(z_output[index_out]-(exp(x_1[ic_set])-1.));
          //if (j==15) printf("%e   %e    %e\n",exp(x_1[ic_set])-1.,z_output[index_out],exp(x_0[ic_set])-1.);
          index_out++;
        }
        n_rk_step[ic_set]++;
      }
    }
  }

  for (j=0; j<k_num; j++)
  {
      //a_rk=1./(1.+z_output[index_out]);
      a_rk=1.0;
      ON_E2_rk = ONE2(a_rk);
      E2_rk = OR0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ON_E2_rk;
      OCB_rk = OCB(a_rk,E2_rk);
      if (only_backg_rel == 'T') ON_rk=((M_nu/(93.14*h*h))/(a_rk*a_rk*a_rk))/E2_rk;
      else ON_rk = ON(ON_E2_rk,E2_rk);
      rk_1_3w = func_1_3w(a_rk);
      FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
      if (only_backg_rel == 'T') FNU = ON_rk/(OCB_rk+rk_1_3w*ON_rk);

    matrix_M[0][0] = 1.;
    matrix_M[0][1] = 0.;
    matrix_M[0][2] = FNU;
    matrix_M[0][3] = 1.-FNU;

    for (ic_set=0; ic_set<4; ic_set++)
    {
      bin_beta_in = output_number-1;   //FIX THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      bin_fc_in = output_number-1;
      bin_fn_in = output_number-1;
      matrix_M[1][ic_set] = Xc_k_z[j][bin_fc_in][ic_set] + FC[j]*Dc_k_z[j][bin_fc_in][ic_set];
      matrix_M[2][ic_set] = BETA[j]*Dc_k_z[j][bin_beta_in][ic_set] - Dn_k_z[j][bin_beta_in][ic_set];
      matrix_M[3][ic_set] = Xn_k_z[j][bin_fn_in][ic_set] + FN[j]*Dn_k_z[j][bin_fn_in][ic_set];
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
      if (fixed_neutrino_frac=='F')
      {
        a_rk=1./(1.+z_output[index_out]);
        ON_E2_rk = ONE2(a_rk);
        E2_rk = OR0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ON_E2_rk;
        OCB_rk = OCB(a_rk,E2_rk);
        if (only_backg_rel == 'T') ON_rk=((M_nu/(93.14*h*h))/(a_rk*a_rk*a_rk))/E2_rk;
        else ON_rk = ON(ON_E2_rk,E2_rk);
        rk_1_3w = func_1_3w(a_rk);
        FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
        if (only_backg_rel == 'T') FNU = ON_rk/(OCB_rk+rk_1_3w*ON_rk);
      }
      delta_m = FNU*delta_n+(1.-FNU)*delta_c;
      f_c = -x_c/delta_c;
      f_n = -x_n/delta_n;
      f_m = -(FNU*x_n+(1.-FNU)*x_c)/delta_m;
      delta_c = delta_c;
      delta_n = delta_n;

      fprintf(o[index_out],"%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",k[j],delta_c,delta_n,delta_m,f_c,f_n,f_m);
    }
  }

  stop = time(NULL);
  printf("   %.0lf%% completed, end of computation! It took %i s.\n",100.,(int)(stop-start));

  for (index_out=0; index_out<output_number; index_out++) fclose(o[index_out]);
}

double det(double **a_in,int n)
{
  int i,j,j1,j2;
  double determinant = 0;
  double **m = NULL;

  if (n < 1)
  {
    printf("You tried to compute the determinant of \n");
    printf("a matrix of negative order?\n");
    exit(1);
  }
  else if (n == 1)
  {
     determinant = a_in[0][0];
  }
  else if (n == 2)
  {
     determinant = a_in[0][0] * a_in[1][1] - a_in[1][0] * a_in[0][1];
  }
  else
  {
     determinant = 0;
     for (j1=0; j1<n; j1++)
     {
        m = malloc((n-1)*sizeof(double *));
        for (i=0; i<n-1; i++)
           m[i] = malloc((n-1)*sizeof(double));
        for (i=1; i<n; i++)
        {
           j2 = 0;
           for (j=0;j<n;j++)
           {
              if (j == j1)
                 continue;
              m[i-1][j2] = a_in[i][j];
              j2++;
           }
        }
        determinant += pow(-1.,2.+j1) * a_in[0][j1] * det(m,n-1);
        for (i=0; i<n-1; i++)
           free(m[i]);
        free(m);
     }
  }
  return determinant;
}

int find_z_bin(double Z, double *vec, int n_elems)
{
  int i = 0;
  if (Z < z0 || Z > z1) exit(1);
  while (log(1.+Z) > vec[i]) {i++; if(i>n_elems) {printf("Search exceded array dimension\n"); exit(1);} }
  return i;
}

double ** allocate_matrix(int rows, int columns)
{
  int i;
  double **m;
  m = malloc(rows*sizeof(double *));
  if (m==NULL)
  {
    printf("Bad memory allocation!\n");
    exit(1);
  }
  for(i=0; i<rows; i++)
  {
    m[i] = malloc(columns*sizeof(double));
    if (m==NULL)
    {
      printf("Bad memory allocation!\n");
      exit(1);
    }
  }
  return m;
}

void deallocate_matrix(double **m, int rows, int columns)
{
  int i;
  for(i=0; i<rows; i++)
  {
    free(m[i]);
  }
  free(m);
}

void read_parameter_file(char parfile[])
{
  int i,j;
  int fscanfcheck=0;
  int Nread = 32;
  char do_rescaled_ps_str[100];
  char neutrino_pressure_str[100];
  char print_hubble_str[100];
  char fixed_neutrino_frac_str[100];
  char only_backg_rel_str[100];
  char generate_boundary_str[100];
  char fixed_f_in_bc_str[100];
  char param_names[Nread][100];

  sprintf(param_names[0],"input_file");
  sprintf(param_names[1],"outputfile");
  sprintf(param_names[2],"nstep");
  sprintf(param_names[3],"z0");
  sprintf(param_names[4],"z1");
  sprintf(param_names[5],"z_beta_in");
  sprintf(param_names[6],"z_fc_in");
  sprintf(param_names[7],"z_fn_in");
  sprintf(param_names[8],"output_number");
  sprintf(param_names[9],"z_output");
  sprintf(param_names[10],"h");
  sprintf(param_names[11],"OB0");
  sprintf(param_names[12],"OC0");
  sprintf(param_names[13],"OG0");
  sprintf(param_names[14],"M_nu");
  sprintf(param_names[15],"do_rescaled_ps");
  sprintf(param_names[16],"As");
  sprintf(param_names[17],"ns");
  sprintf(param_names[18],"tau_reio");
  sprintf(param_names[19],"kmax");
  sprintf(param_names[20],"neutrino_pressure");
  sprintf(param_names[21],"N_nu");
  sprintf(param_names[22],"Neff");
  sprintf(param_names[23],"print_hubble");
  sprintf(param_names[24],"fixed_neutrino_frac");
  sprintf(param_names[25],"only_backg_rel");
  sprintf(param_names[26],"class_folder");
  sprintf(param_names[27],"generate_boundary");
  sprintf(param_names[28],"ps_norm_at_z");
  sprintf(param_names[29],"fixed_f_in_bc");
  sprintf(param_names[30],"ic_mode");
  sprintf(param_names[31],"normfile");
  int check[Nread]; for (i=0; i<Nread; i++) check[i]=0;
  int check_output_number = 0;
  char test_str[100];
  char c;

  FILE *input = fopen(parfile,"r");
  if (input==NULL)
  {
    printf("You need to specify a parameter file containing the following (%i) entries:\n",Nread);
    for(i=0; i<Nread; i++)
    printf("   %s\n",param_names[i]);
    exit(1);
  }

  while (!feof(input))
  {
    c = getc(input);
    if (isspace(c)==0)
    {
      if (c!='#')
      {
        ungetc(c,input);
        fscanfcheck=fscanf(input,"%s", test_str);
        for (i = 0; i < Nread; i++)
        {
          if (strcmp(test_str,param_names[i]) == 0)
          {
            if (i==0)
            {
              fscanfcheck=fscanf(input,"%s",input_file);
              check[i] = 1;
            }

            if (i==1)
            {
              fscanfcheck=fscanf(input,"%s",outputfile);
              check[i] = 1;
            }

            if (i==2)
            {
              fscanfcheck=fscanf(input,"%i",&nstep);
              check[i] = 1;
            }

            if (i==3)
            {
              fscanfcheck=fscanf(input,"%lf",&z0);
              check[i] = 1;
            }

            if (i==4)
            {
              fscanfcheck=fscanf(input,"%lf",&z1);
              check[i] = 1;
            }

            if (i==5)
            {
              fscanfcheck=fscanf(input,"%lf",&z_beta_in);
              check[i] = 1;
            }

            if (i==6)
            {
              fscanfcheck=fscanf(input,"%lf",&z_fc_in);
              check[i] = 1;
            }

            if (i==7)
            {
              fscanfcheck=fscanf(input,"%lf",&z_fn_in);
              check[i] = 1;
            }

            if (i==8)
            {
              fscanfcheck=fscanf(input,"%i",&output_number);
              if (output_number<=0)
              {
                printf("You should specify at least 1 output redshift.\n");
                exit(1);
              }
              else
              {
                check_output_number=1;
                check[i] = 1;
              }
            }

            if (i==9)
            {
              if (check_output_number==0)
              {
                printf("I need you specify the number of required outputs\n");
                printf("before you tell me the output redshifts.\n");
                exit(1);
              }
              else
              {
                z_output = malloc(output_number*sizeof(double));
                if (z_output==NULL)
                { printf("Bad memory allocation...\n"); exit(1); }
                for (j=0; j<output_number; j++)
                  fscanfcheck=fscanf(input,"%lf",&z_output[j]);
                check[i] = 1;
              }
            }

            if (i==10)
            {
              fscanfcheck=fscanf(input,"%lf",&h);
              check[i] = 1;
            }

            if (i==11)
            {
              fscanfcheck=fscanf(input,"%lf",&OB0);
              check[i] = 1;
            }

            if (i==12)
            {
             fscanfcheck=fscanf(input,"%lf",&OC0);
             check[i] = 1;
            }

            if (i==13)
            {
              fscanfcheck=fscanf(input,"%lf",&OG0);
              OG0/=(h*h);
              check[i] = 1;
            }

            if (i==14)
            {
              fscanfcheck=fscanf(input,"%lf",&M_nu);
              check[i] = 1;
            }

            if (i==15)
            {
              fscanfcheck=fscanf(input,"%s",do_rescaled_ps_str);
              if(strcmp(do_rescaled_ps_str,"T")==0||strcmp(do_rescaled_ps_str,"True")==0||strcmp(do_rescaled_ps_str,"t")==0||strcmp(do_rescaled_ps_str,"true")==0)
                do_rescaled_ps='T';
              else if(strcmp(do_rescaled_ps_str,"F")==0||strcmp(do_rescaled_ps_str,"False")==0||strcmp(do_rescaled_ps_str,"f")==0||strcmp(do_rescaled_ps_str,"false")==0)
                do_rescaled_ps='F';
              else
              {
                printf("\nYou wrote \'do_rescaled_ps %s\'.\n",do_rescaled_ps_str);
                printf("Please, set \'do_rescaled_ps\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==16)
            {
              fscanfcheck=fscanf(input,"%lf",&As);
              check[i] = 1;
            }

            if (i==17)
            {
              fscanfcheck=fscanf(input,"%lf",&ns);
              check[i] = 1;
            }

            if (i==18)
            {
              fscanfcheck=fscanf(input,"%lf",&tau_reio);
              check[i] = 1;
            }

            if (i==19)
            {
              fscanfcheck=fscanf(input,"%lf",&kmax);
              check[i] = 1;
            }

            if (i==20)
            {
              fscanfcheck=fscanf(input,"%s",neutrino_pressure_str);
              if(strcmp(neutrino_pressure_str,"T")==0||strcmp(neutrino_pressure_str,"True")==0||strcmp(neutrino_pressure_str,"t")==0||strcmp(neutrino_pressure_str,"true")==0)
                neutrino_pressure='T';
              else if(strcmp(neutrino_pressure_str,"F")==0||strcmp(neutrino_pressure_str,"False")==0||strcmp(neutrino_pressure_str,"f")==0||strcmp(neutrino_pressure_str,"false")==0)
                neutrino_pressure='F';
              else
              {
                printf("\nYou wrote \'neutrino_pressure %s\'.\n",neutrino_pressure_str);
                printf("Please, set \'neutrino_pressure\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==21)
            {
              fscanfcheck=fscanf(input,"%lf",&N_nu);
              check[i] = 1;
            }

            if (i==22)
            {
              fscanfcheck=fscanf(input,"%lf",&Neff);
              check[i] = 1;
            }

            if (i==23)
            {
              fscanfcheck=fscanf(input,"%s",print_hubble_str);
              if(strcmp(print_hubble_str,"T")==0||strcmp(print_hubble_str,"True")==0||strcmp(print_hubble_str,"t")==0||strcmp(print_hubble_str,"true")==0)
                print_hubble='T';
              else if(strcmp(print_hubble_str,"F")==0||strcmp(print_hubble_str,"False")==0||strcmp(print_hubble_str,"f")==0||strcmp(print_hubble_str,"false")==0)
                print_hubble='F';
              else
              {
                printf("\nYou wrote \'print_hubble %s\'.\n",print_hubble_str);
                printf("Please, set \'print_hubble\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==24)
            {
              fscanfcheck=fscanf(input,"%s",fixed_neutrino_frac_str);
              if(strcmp(fixed_neutrino_frac_str,"T")==0||strcmp(fixed_neutrino_frac_str,"True")==0||strcmp(fixed_neutrino_frac_str,"t")==0||strcmp(fixed_neutrino_frac_str,"true")==0)
                fixed_neutrino_frac='T';
              else if(strcmp(fixed_neutrino_frac_str,"F")==0||strcmp(fixed_neutrino_frac_str,"False")==0||strcmp(fixed_neutrino_frac_str,"f")==0||strcmp(fixed_neutrino_frac_str,"false")==0)
                fixed_neutrino_frac='F';
              else
              {
                printf("\nYou wrote \'fixed_neutrino_frac %s\'.\n",fixed_neutrino_frac_str);
                printf("Please, set \'fixed_neutrino_frac\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==25)
            {
              fscanfcheck=fscanf(input,"%s",only_backg_rel_str);
              if(strcmp(only_backg_rel_str,"T")==0||strcmp(only_backg_rel_str,"True")==0||strcmp(only_backg_rel_str,"t")==0||strcmp(only_backg_rel_str,"true")==0)
                only_backg_rel='T';
              else if(strcmp(only_backg_rel_str,"F")==0||strcmp(only_backg_rel_str,"False")==0||strcmp(only_backg_rel_str,"f")==0||strcmp(only_backg_rel_str,"false")==0)
                only_backg_rel='F';
              else
              {
                printf("\nYou wrote \'only_backg_rel %s\'.\n",only_backg_rel_str);
                printf("Please, set \'only_backg_rel\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              if (only_backg_rel=='T'&&neutrino_pressure=='F')
              {
                printf("Incompatible parameters: since you want \'only_backg_rel=T\'\n\'neutrino_pressure\' will be set to T\n");
                neutrino_pressure='T';
              }
              check[i] = 1;
            }

            if (i==26)
            {
              fscanfcheck=fscanf(input,"%s",class_folder);
              check[i] = 1;
            }

            if (i==27)
            {
              fscanfcheck=fscanf(input,"%s",generate_boundary_str);
              if(strcmp(generate_boundary_str,"T")==0||strcmp(generate_boundary_str,"True")==0||strcmp(generate_boundary_str,"t")==0||strcmp(generate_boundary_str,"true")==0)
                generate_boundary='T';
              else if(strcmp(generate_boundary_str,"F")==0||strcmp(generate_boundary_str,"False")==0||strcmp(generate_boundary_str,"f")==0||strcmp(generate_boundary_str,"false")==0)
                generate_boundary='F';
              else
              {
                printf("\nYou wrote \'generate_boundary %s\'.\n",generate_boundary_str);
                printf("Please, set \'generate_boundary\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }
            if (i==28)
            {
              fscanfcheck=fscanf(input,"%lf",&ps_norm_at_z);
              if (ps_norm_at_z!=0. && ps_norm_at_z!=99.)
              {
                printf("Error! Normalization can only be at z=0 or 99.\n");
                printf("Please, set \'ps_norm_at_z\' to a legal value.\n");
                exit(1);
              }
              check[i] = 1;
            }
            if (i==29)
            {
              fscanfcheck=fscanf(input,"%s",fixed_f_in_bc_str);
              if(strcmp(fixed_f_in_bc_str,"T")==0||strcmp(fixed_f_in_bc_str,"True")==0||strcmp(fixed_f_in_bc_str,"t")==0||strcmp(fixed_f_in_bc_str,"true")==0)
                fixed_f_in_bc='T';
              else if(strcmp(fixed_f_in_bc_str,"F")==0||strcmp(fixed_f_in_bc_str,"False")==0||strcmp(fixed_f_in_bc_str,"f")==0||strcmp(fixed_f_in_bc_str,"false")==0)
                fixed_f_in_bc='F';
              else
              {
                printf("\nYou wrote \'fixed_f_in_bc  %s\'.\n",generate_boundary_str);
                printf("Please, set \'fixed_f_in_bc\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }
            if (i==30)
            {
              fscanfcheck=fscanf(input,"%s",ic_mode);
              check[i] = 1;
            }
            if (i==31)
            {
              fscanfcheck=fscanf(input,"%s",normfile);
              check[i] = 1;
            }
          }
        }
      }
      else while(c!='\n') c = getc(input);
    }
  }
  fclose(input);

  for(i=0; i<Nread; i++)
  {
    if (check[i]==0)
    {
      printf("You didn't specify a value for %s.\n",param_names[i]);
      exit(1);
    }
  }

  OR0 = (Neff*(7./8.)*pow(4./11.,4./3.)+1.)*(OG0);

  printf("\nRead parameters:\n");
  printf("  %s = %s\n",param_names[0],input_file);
  printf("  %s = %s\n",param_names[1],outputfile);
  printf("  %s = %i\t\t\t",param_names[2],nstep);
  printf("  %s = %lf\n",param_names[3],z0);
  printf("  %s = %lf\t\t",param_names[4],z1);
  printf("  %s = %lf\n",param_names[5],z_beta_in);
  printf("  %s = %lf\t\t",param_names[6],z_fc_in);
  printf("  %s = %lf\n",param_names[7],z_fn_in);
  printf("  %s = %i\n",param_names[8],output_number);
  printf("  %s =",param_names[9]);
  for (i=0; i<output_number; i++)
    printf(" %.2lf",z_output[i]); printf("\n");
  printf("  %s = %lf\t\t\t",param_names[10],h);
  printf("  %s = %lf\n",param_names[11],OB0);
  printf("  %s = %lf\t\t",param_names[12],OC0);
  printf("  %s = %e\n","OR0",OR0);
  printf("  %s = %lf\t\t",param_names[14],M_nu);
  printf("  %s = %e\n",param_names[16],As);
  printf("  %s = %lf\t\t\t",param_names[17],ns);
  printf("  %s = %lf\n",param_names[18],tau_reio);
  printf("  %s = %c\n",param_names[15],do_rescaled_ps);
  printf("  %s = %c\n",param_names[20],neutrino_pressure);
  printf("  %s = %c\n",param_names[24],fixed_neutrino_frac);
  printf("  %s = %c\n",param_names[25],only_backg_rel);

  OM0 = OB0+OC0+ONE2(1.);
  OL0 = 1.-OM0-OR0;

/*  OM0=0.3175;*/
/*  OL0=0.6825;*/

  printf("\nDerived parameters: OM0 = %lf, OL0 = %lf\n\n",OM0,OL0);
}

void fscanf_error(int n)
{
  printf("Error reading file. %i values were expected.\n",n);
  exit(1);
}

void rescale_ps(int knum, double *k)
{
  printf("\nRescaling of the PS requested.\n");

  char dir_chain[1000];
  if (getcwd(dir_chain,sizeof(dir_chain))==NULL)
  {
    printf("\nError retrieving current directory path.\n");
    exit(1);
  }

  create_class_ini_file(dir_chain);

  printf("\nCalling camb and generating the P(K) and T(k) at the \n"
         "requested output redshifts.\n");
  char command[200];
  sprintf(command,"%scamb %s/class_tabs/power.ini",class_folder,dir_chain);
  FILE *cterm=popen(command,"w");
  pclose(cterm);

  // power spectra at z=0,99 used for normalization
  mode = 3;
  create_class_ini_file(dir_chain);

  sprintf(command,"%scamb %s/class_tabs/power_norm.ini",class_folder,dir_chain);
  cterm=popen(command,"w");
  pclose(cterm);
  mode = 0;

  if (ps_norm_at_z == 0.)
  {
    char test='a';
    int n=0;
    sprintf(command,"%s/class_tabs/power_norm_z1_pk.dat",dir_chain);
    FILE *spectrumz0 = fopen(command,"r");
    if (spectrumz0 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    while (fscanf(spectrumz0,"%c",&test)!=EOF) if(test=='\n') n++;
    if ((n-1)!=knum)
    {
      printf("\nThe number of requested ks doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    double *P = malloc(knum*sizeof(double));
    if (P==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dcb = malloc(knum*sizeof(double));
    if (Dcb==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dnu = malloc(knum*sizeof(double));
    if (Dnu==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dm = malloc(knum*sizeof(double));
    if (Dm==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    fseek(spectrumz0,0,SEEK_SET);
    int i;
    int fscanfcheck=0;
    char *dummy;
    double val,valk,Tc,Tb,Tn,Tm;
    double kpivot=0.05;

    char buf[5000];
    dummy=fgets (buf, sizeof(buf), spectrumz0);

    //FILE *ccc = fopen("testps.txt","w");

    if (N_nu != 0)
    {
      for(i=0; i<knum; i++)
      {
        fscanfcheck=fscanf(spectrumz0,"%lf %lf",&valk,&P[i]);
        if (fscanfcheck!=2) fscanf_error(2);
	//fprintf(ccc,"%e\t%e\n",k[i],P[i]);
      }
    }
    else
    {
      for(i=0; i<knum; i++)
      {
        fscanfcheck=fscanf(spectrumz0,"%lf %lf",&valk,&P[i]);
        if (fscanfcheck!=2) fscanf_error(2);
      }
    }
    fclose(spectrumz0);
    //fclose(ccc);

    int index_out=0;
    char outfile[200];
    char Dfile[200];
    printf("\n");
    for(index_out=0; index_out<output_number; index_out++)
    {
      sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/class_tabs/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Dcb[i]*Dcb[i]*P[i]);
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/class_tabs/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Dnu[i]*Dnu[i]*P[i]);
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/class_tabs/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledm = fopen(outfile,"w");
      // double acurr=1./(z_output[index_out]+1.);
      // double FNU = (ONE2(acurr))/((OC0+OB0)/pow(acurr,3)+ONE2(acurr));
      // double p_cb_z,p_n_z;
      for (i=0; i<knum; i++)
      {
        // p_cb_z=Dcb[i]*Dcb[i]*P[i];
        // p_n_z=Dnu[i]*Dnu[i]*P[i];
        fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Dm[i]*Dm[i]*P[i]);
        // fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],((1.-FNU)*(1.-FNU))*p_cb_z+(FNU*FNU)*p_n_z+FNU*(1.-FNU)*sqrt(p_cb_z*p_n_z));
      }
      fclose(outrescaledm);
    }

    free(P);
    free(Dcb);
    free(Dnu);
    free(Dm);
  }
  else
  {
    char test='a';
    int n=0;
    // if (M_nu==0.6) sprintf(command,"%s/class_tabs/Pcb_Mnu_0.6_planck_correct_rescaled_znum17.txt",dir_chain);
    // else if (M_nu==0.0) sprintf(command,"%s/class_tabs/Pcb_Mnu_0_planck_correct_rescaled_znum17.txt",dir_chain);
    // else sprintf(command,"%s/class_tabs/Pcb_Mnu_0.2_planck_correct_rescaled_znum17.txt",dir_chain);
    sprintf(command,"%s/class_tabs/Pcb_%s.txt",dir_chain,normfile);
    FILE *spectrumz99 = fopen(command,"r");
    if (spectrumz99 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    while (fscanf(spectrumz99,"%c",&test)!=EOF) if(test=='\n') n++;
    if (n!=knum)
    {
      printf("\nThe number of ks requested doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    double *P = malloc(knum*sizeof(double));
    if (P==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}

    sprintf(command,"%s/class_tabs/power_norm_z1_pk.dat",dir_chain);
    FILE *spectrumz0 = fopen(command,"r");
    if (spectrumz0 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    char *dummy;
    double val,valk,Tc,Tb,Tn,Tm;
    double kpivot=0.05;
    char buf[5000];
    int i;
    int fscanfcheck=0;
    dummy=fgets (buf, sizeof(buf), spectrumz0);

    if (N_nu != 0)
    {
      for(i=0; i<knum; i++)
      {
        fscanfcheck=fscanf(spectrumz0,"%lf %lf",&valk,&P[i]);
        if (fscanfcheck!=2) fscanf_error(2);
      }
    }
    else
    {
      for(i=0; i<knum; i++)
      {
        fscanfcheck=fscanf(spectrumz0,"%lf %lf",&valk,&P[i]);
        if (fscanfcheck!=2) fscanf_error(2);
      }
    }
    fclose(spectrumz0);

    double *Pc = malloc(knum*sizeof(double));
    if (Pc==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Pn = malloc(knum*sizeof(double));
    if (Pn==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Pm = malloc(knum*sizeof(double));
    if (Pm==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dcb = malloc(knum*sizeof(double));
    if (Dcb==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dnu = malloc(knum*sizeof(double));
    if (Dnu==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *Dm = malloc(knum*sizeof(double));
    if (Dm==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *D99cb = malloc(knum*sizeof(double));
    if (D99cb==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *D99nu = malloc(knum*sizeof(double));
    if (D99nu==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    double *D99m = malloc(knum*sizeof(double));
    if (D99m==NULL)
    {printf("Bad memory allocation.\n"); exit(1);}
    fseek(spectrumz99,0,SEEK_SET);

    for(i=0; i<knum; i++)
    {
      fscanfcheck=fscanf(spectrumz99,"%lf %lf",&valk,&Pc[i]);
      if (fscanfcheck!=2) fscanf_error(2);
    }
    fclose(spectrumz99);

    // if (M_nu==0.6) sprintf(command,"%s/class_tabs/Pnu_Mnu_0.6_planck_correct_rescaled_znum17.txt",dir_chain);
    // else if (M_nu==0.0) sprintf(command,"%s/class_tabs/Pcb_Mnu_0_planck_correct_rescaled_znum17.txt",dir_chain);
    // else sprintf(command,"%s/class_tabs/Pnu_Mnu_0.2_planck_correct_rescaled_znum17.txt",dir_chain);
    sprintf(command,"%s/class_tabs/Pnu_%s.txt",dir_chain,normfile);
    spectrumz99 = fopen(command,"r");
    if (spectrumz99 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    for(i=0; i<knum; i++)
    {
      fscanfcheck=fscanf(spectrumz99,"%lf %lf",&valk,&Pn[i]);
      if (fscanfcheck!=2) fscanf_error(2);
    }
    fclose(spectrumz99);

    sprintf(command,"%s/class_tabs/Pm_%s.txt",dir_chain,normfile);
    spectrumz99 = fopen(command,"r");
    if (spectrumz99 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    for(i=0; i<knum; i++)
    {
      fscanfcheck=fscanf(spectrumz99,"%lf %lf",&valk,&Pm[i]);
      if (fscanfcheck!=2) fscanf_error(2);
    }
    fclose(spectrumz99);

    int index_out=0;
    char outfile[200];
    char Dfile[200];
    char D99file[200];

    //sprintf(D99file,"Mnu_0.6_planck_correct_znum17.txt");
    sprintf(D99file,"%s_znum%i.txt",outputfile,output_number-1);
    read_D(D99file,knum,D99cb,D99nu,D99m);

    printf("\n");
    for(index_out=0; index_out<output_number; index_out++)
    {
      sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/class_tabs/Pcb_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pcb_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        //if (index_out==output_number-1) printf("%lf\n",(Dcb[i]*Dcb[i])/(D99cb[i]*D99cb[i]));
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Pc[i]*((Dcb[i]*Dcb[i])/(D99cb[i]*D99cb[i])));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/class_tabs/Pnu_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pnu_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Pn[i]*((Dnu[i]*Dnu[i])/(D99nu[i]*D99nu[i])));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/class_tabs/Pm_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pm_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledm = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Pm[i]*((Dm[i]*Dm[i])/(D99m[i]*D99m[i])));
      }
      fclose(outrescaledm);

      // norm z0

      sprintf(outfile,"%s/class_tabs/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        //if (index_out==output_number-1) printf("%lf\n",(Dcb[i]*Dcb[i])/(D99cb[i]*D99cb[i]));
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],P[i]*(Dcb[i]*Dcb[i]));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/class_tabs/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],P[i]*(Dnu[i]*Dnu[i]));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/class_tabs/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file class_tabs/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescaledm = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],P[i]*(Dm[i]*Dm[i]));
      }
      fclose(outrescaledm);
    }

    free(P);
    free(Pc);
    free(Pn);
    free(Pm);
    free(Dcb);
    free(Dnu);
    free(Dm);
    free(D99cb);
    free(D99nu);
    free(D99m);
  }
}

void read_D(char filename[],int n, double *dc, double *dn, double *dm)
{
  FILE *f = fopen(filename,"r");
  if (f == NULL) {printf("Problem opening file %s!\n",filename); exit(1);}

  int i; double val;
  int fscanfcheck=0;
  for (i=0; i<n; i++)
  {
    fscanfcheck=fscanf(f,"%lf %lf %lf %lf %lf %lf %lf",&val,&dc[i],&dn[i],&dm[i],&val,&val,&val);
    if (fscanfcheck!=7) fscanf_error(7);
  }

  fclose(f);
}

void create_class_ini_file (char dir_chain[])
{
  char fileopen[200];

  if (mode==0)
  {
    printf("Creating camb input_file\n   %s/class_tabs/power.ini\n",dir_chain);
    sprintf(fileopen,"%s/class_tabs/power.ini",dir_chain);
  }
  if (mode==1 || mode==2)
  {
    printf("Creating camb input_file\n   %s/boundary_conditions_module/tabs/power.ini\n",dir_chain);
    sprintf(fileopen,"%s/boundary_conditions_module/tabs/power.ini",dir_chain);
  }
  if (mode==3)
  {
    sprintf(fileopen,"%s/class_tabs/power_norm.ini",dir_chain);
  }

  double Nur,N_eigenstates,N_degeneracies;
  if (M_nu==0.)
  {
    Nur = 3.046;
    N_eigenstates = 0.0;
    N_degeneracies = 0.0;
  }
  else
  {
    Nur = 0.046;
    N_eigenstates = 1.0;
    N_degeneracies = 3.0;
  }

  int i;

FILE *out = fopen(fileopen,"w");
if (out == NULL)
{
  if (mode==0 || mode==3)
  {
    printf("Does the folder class_tabs exist? Creating...\n");
    FILE *foldercreation = popen("mkdir class_tabs","w");
    pclose(foldercreation);
    out = fopen(fileopen,"w");
    if (out==NULL)
    {
      printf("Error creating params file...\n");
      exit(1);
    }
  }
  else
  {
    printf("Does the folder boundary_conditions_module/tabs exist? Creating...\n");
    FILE *foldercreation = popen("mkdir boundary_conditions_module/tabs","w");
    pclose(foldercreation);
    out = fopen(fileopen,"w");
    if (out==NULL)
    {
      printf("Error creating params file...\n");
      exit(1);
    }
  }
}

int znum=0;

if (mode==0) {znum=output_number; fprintf(out,"output_root = %s/class_tabs/power\n",dir_chain);}
if (mode==1) {znum=50; fprintf(out,"output_root = %s/boundary_conditions_module/tabs/power\n",dir_chain);}
if (mode==2) {znum=1; fprintf(out,"output_root = %s/boundary_conditions_module/tabs/power_zin\n",dir_chain);}
if (mode==3) {znum=2; fprintf(out,"output_root = %s/class_tabs/power_norm\n",dir_chain);}

fprintf(out,
"get_scalar_cls = F \n"
"get_vector_cls = F\n"
"get_tensor_cls = F\n"
"get_transfer = T\n"
"do_lensing= F\n"
"do_nonlinear = 0\n"
"l_max_scalar	   = 2000\n"
"k_eta_max_scalar  = 4000\n"
"l_max_tensor	   = 1500\n"
"k_eta_max_tensor  = 3000\n"
"												     \n"
"#Main cosmological parameters, neutrino masses are assumed degenerate\n"
"# If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k\n"
"use_physical	= F\n"
"												     \n"
"omk= 0.0\n"
"hubble= %lf\n"
"w= -1\n"
"cs2_lam= 1\n"
"omega_baryon= %.10lf\n"
"omega_cdm= %.10lf\n"
"omega_lambda= %.10lf\n"
"omega_neutrino = %.10lf\n"
"												     \n"
"temp_cmb= 2.7255\n"
"helium_fraction= 0.24\n"
"massless_neutrinos = %lf\n"
"massive_neutrinos  = %i\n"
"												     \n"
"												     \n"
"nu_mass_eigenstates = %i\n"
"nu_mass_degeneracies = %i\n"
"nu_mass_fractions = 1\n"
"												     \n"
"initial_power_num = 1\n"
"pivot_scalar = 0.05\n"
"pivot_tensor = 0.05\n"
"scalar_amp(1) = %e\n"
"scalar_spectral_index(1) = %.8lf\n"
"scalar_nrun(1) = 0\n"
"tensor_spectral_index(1)  = 0\n"
"initial_ratio(1)= 1\n"
"												     \n"
"reionization= T\n"
"												     \n"
"re_use_optical_depth = T\n"
"re_optical_depth= 0.0925\n"
"re_redshift= 11\n"
"re_delta_redshift= 0.5\n"
"re_ionization_frac= -1\n"
"												     \n"
"RECFAST_fudge = 1.14\n"
"RECFAST_fudge_He = 0.86\n"
"RECFAST_Heswitch = 6\n"
"RECFAST_Hswitch  = T\n"
"												     \n"
"initial_condition= 1\n"
"initial_vector = -1 0 0 0 0\n"
"												     \n"
"vector_mode = 0\n"
"												     \n"
"COBE_normalize = F\n"
"CMB_outputscale = 7.42835025e12 \n"
"												     \n"
"transfer_power_var = 7\n"
"transfer_high_precision = T\n"
"transfer_kmax = 50\n"
"transfer_k_per_logint= 10\n"
"												     \n"
"transfer_num_redshifts = %i\n"
"												     \n"
"transfer_interp_matterpower = F\n"
"												     \n"
"scalar_output_file = scalCls.dat\n"
"vector_output_file = vecCls.dat\n"
"tensor_output_file = tensCls.dat\n"
"total_output_file  = totCls.dat\n"
"lensed_output_file = lensedCls.dat\n"
"lensed_total_output_file  =lensedtotCls.dat\n"
"lens_potential_output_file = lenspotentialCls.dat\n"
"FITS_filename      = scalCls.fits\n"
"												     \n"
"feedback_level = 0\n"
"												     \n"
"lensing_method = 1\n"
"accurate_BB = F\n"
"												     \n"
"massive_nu_approx = 1\n"
"												     \n"
"accurate_polarization   = T\n"
"accurate_reionization   = T\n"
"do_tensor_neutrinos	 = F\n"
"												     \n"
"#Whether to turn off small-scale late time radiation hierarchies (save time,v. accurate)\n"
"do_late_rad_truncation = F\n"
"high_accuracy_default = F\n"
"use_spline_template = F\n"
"												     \n"
"number_of_threads= 0\n"
"accuracy_boost= 2\n"
"l_accuracy_boost= 2\n"
"l_sample_boost= 2\n",
h*100.0,
OB0,
OC0,
1.0-OB0-OC0-ONE2(1.0),
ONE2(1.),
Nur,
(int)round(N_nu),
(int)round(N_eigenstates),
(int)round(N_degeneracies),
As,
ns,
znum);

if (mode==0) {
for (i=0; i< output_number; i++)
{
	fprintf(out,
	"transfer_redshift(%i) = %lf \n"
	"transfer_filename(%i) = z%i_tk.dat \n"
	"transfer_matterpower(%i) = z%i_pk.dat \n",
	i+1,
	z_output[output_number-1-i],
	i+1,
	output_number-i,
	i+1,
	output_number-i);
} }

if (mode==1) {
for (i=0; i<50; i++)
{
	fprintf(out,
	"transfer_redshift(%i) = %lf \n"
	"transfer_filename(%i) = z%i_tk.dat \n"
	"transfer_matterpower(%i) = z%i_pk.dat \n",
	i+1,
	bc_zz[49-i],
	i+1,
	50-i,
	i+1,
	50-i);
} }

if (mode==2) {
fprintf(out,
"transfer_redshift(%i) = %lf \n"
"transfer_filename(%i) = tk.dat \n"
"transfer_matterpower(%i) = pk.dat \n",
1,z_fc_in,1,1,1,1);
}

if (mode==3) {
fprintf(out,
"transfer_redshift(%i) = %lf \n"
"transfer_filename(%i) = z%i_tk.dat \n"
"transfer_matterpower(%i) = z%i_pk.dat \n"
"transfer_redshift(%i) = %lf \n"
"transfer_filename(%i) = z%i_tk.dat \n"
"transfer_matterpower(%i) = z%i_pk.dat \n",
	1,z_fc_in,1,2,1,2,
2,0.0,2,1,2,1);
}

fclose(out);
}

void read_GG_FF_tabs()
{
  FILE *gft = fopen("tabulated_functions/FF_GG_func_tab.dat","r");
  if (gft==NULL)
  {
    char dir[1000];
    if (getcwd(dir,sizeof(dir))==NULL)
    {
      printf("\nError retrieving current directory path.\n");
      exit(1);
    }
    printf("Ooops! A folder named \'tabulated_functions\' was expected to be found\n"
    "here %s\n"
    "and it was expected to contain 1 file named \'FF_GG_func_tab.dat\' but \n"
    "something went wrong.\n",dir);
    exit(1);
  }

  int i; int l=0;
  int fscanfcheck=0;
  char test='a';
  while(fscanf(gft,"%c",&test)!=EOF) if (test=='\n') l++;
  fseek(gft,0,SEEK_SET);

  ytab = malloc(l*sizeof(double));
  GGtab = malloc(l*sizeof(double));
  FFtab = malloc(l*sizeof(double));

  for (i=0; i<l; i++)
  {
    fscanfcheck=fscanf(gft,"%lf %lf %lf",&ytab[i],&FFtab[i],&GGtab[i]);
    if (fscanfcheck!=3) fscanf_error(3);
  }
  fclose(gft);
  ytab_min=ytab[0];
  ytab_max=ytab[l-1];  //printf("ytabminmax= %e %e\n",ytab_min,ytab_max);
  ntab=l;
  ytabstep = log(ytab[1])-log(ytab[0]);

  double x=0.;
  double xmin=0.;
  double xmax=20.;
  double xstep=(xmax-xmin)/(100000.-1.);

  F_inf=0.;
  for(x=xmin; x<=xmax; x+=xstep)
  {
    F_inf += (x*x/(1.+exp(x)))*xstep;
  }

  G_inf=F_inf;

  F_0=0.;
  for(x=xmin; x<=xmax; x+=xstep)
  {
    F_0 += (x*x*x/(1.+exp(x)))*xstep;
  }

  G_0=0.;

  double Fa1=0;
  xmin=0.; xmax=1000.; xstep=(xmax-xmin)/(1000000.-1.);
  double y = (M_nu/(Kb*Gamma_nu*N_nu*Tcmb_0));
  for(x=xmin; x<=xmax; x+=xstep)
  {
    Fa1 += ((x*x*sqrt(x*x+y*y))/(1.+exp(x)))*xstep;
  }
}

double FF(double Y)
{
  int i=0;
  double m,q;
  if (Y >= ytab_max)
  {
    return Y*F_inf;
  }
  else if (Y < ytab_min)
  {
    return F_0;
  }
  else
  {
    i = 0;
    i = floor((log(Y)-log(ytab_min))/(ytabstep));
    m = (FFtab[i]-FFtab[i+1])/(ytab[i]-ytab[i+1]);
    q = FFtab[i] - m*ytab[i];
    return (m*Y+q);
  }
}

double GG(double Y)
{
  int i=0;
  double m,q;
  if (Y >= ytab_max)
  {
    return G_inf;
  }
  else if (Y < ytab_min)
  {
    i=0;
    m=(GGtab[i]-0.)/(ytab[i]-0.);
    q = 0.;
    return (m*Y+q);
  }
  else
  {
    i=floor((log(Y)-log(ytab_min))/(ytabstep));
    m = (GGtab[i]-GGtab[i+1])/(ytab[i]-ytab[i+1]);
    q = GGtab[i] - m*ytab[i];
    return (m*Y+q);
  }
}

void print_hubble_table()
{
  char openf[200];
  sprintf(openf,"%s_hubble.txt",outputfile);

  printf("\nHubble table written to %s\n",openf);

  FILE *hub = fopen(openf,"w");

  int z_nstep=1000;
  double z_step=(double)((log(1.0+z1)-log(1.0+z0))/(z_nstep-1));

  double zz,E2;

  int i;
  for (i=0; i<z_nstep; i++)
  {
    zz = exp(log(1.+z0)+i*z_step);
    E2 = OR0*(zz*zz*zz*zz)+(OC0+OB0)*(zz*zz*zz)+OL0+ONE2(1./zz);

    fprintf(hub,"%.10e\t%.10e\n",zz-1.,h*H0*sqrt(E2));
  }
  printf("Compare ONU0 = %.10lf    %.10lf\n",ONE2(1.),(M_nu/(93.14*h*h)));

}

void reallocate_matrix(double **m,int rows,int newcols)
{
  int i;
  for (i=0; i<rows; i++)
  {
    if (realloc(m[i],newcols)==NULL);
    {
      printf("Some vector reallocation has been attempted\n");
      printf("but was rejected due to lack of memory.\n");
      printf("Please, try again with less accuracy.\n");
      exit(1);
    }
  }
}


void print_kfs()
{
  char openf[200];
  sprintf(openf,"k_fs.txt");
  FILE *fs = fopen(openf,"w");

  double zz,E2;

  int fs_Nstep = 1e2;
  double fs_step = (log(1.0+z1)-log(1.0+z0))/(fs_Nstep-1);

  double xcurrent,a_rk,ON_E2_rk,E2_rk,OCB_rk,ON_rk,OR_rk,OL_rk,rk_1_3w,
         A_rk,B_rk,kj2_rk,FNU;

  double rk_1_plus_3w;

  int i;

  zz=99.0+1.0;
  a_rk=1.0/zz;
  double fixed_onue2_z99 = ((M_nu/(93.14*h*h))*(zz*zz*zz));
  double E2z99 = OG0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ONE2(a_rk);

  printf("\n\n\t%.10lf  %.10lf  %.10lf\n\n",
          ONE2(1.0),
          func_1_3w(1.0)*ONE2(1.0),
          M_nu/(93.14*h*h));
  printf("\t%.10lf  %.10lf  %.10lf\n\n",
          ONE2(0.01),
          func_1_3w(0.01)*ONE2(0.01),
          fixed_onue2_z99);
  printf("\t%.10lf  %.10lf  %.10lf\n\n",
          ONE2(0.01)/E2z99,
          (func_1_3w(0.01)*ONE2(0.01))/E2z99,
          fixed_onue2_z99/E2z99);

  for (i=0; i<fs_Nstep; i++)
  {
    zz = exp(log(1.+z0)+i*fs_step);

    xcurrent = zz;
    a_rk = 1./xcurrent;
    ON_E2_rk = ONE2(a_rk);
    //E2 = OR0*(zz*zz*zz*zz)+(OC0+OB0)*(zz*zz*zz)+OL0+ONE2(1./zz);
    E2_rk = OG0/(a_rk*a_rk*a_rk*a_rk)+(OC0+OB0)/(a_rk*a_rk*a_rk)+OL0+ON_E2_rk;
    OCB_rk = OCB(a_rk,E2_rk);
    if (only_backg_rel == 'T') ON_rk=((M_nu/(93.14*h*h))/(a_rk*a_rk*a_rk))/E2_rk;
    else ON_rk = ON(ON_E2_rk,E2_rk);
    OR_rk = OR(a_rk,E2_rk);
    OL_rk = OL(a_rk,E2_rk);
    rk_1_3w = func_1_3w(a_rk);
    rk_1_plus_3w = 2.0 - rk_1_3w;
    A_rk = 0.5*(OCB_rk+2.*OR_rk-2.*OL_rk+rk_1_plus_3w*ON_rk-2.);
    if (only_backg_rel == 'T') B_rk = -3.*(OCB_rk+ON_rk)/2.;
    else B_rk = -3.*(OCB_rk+rk_1_3w*ON_rk)/2.;
    kj2_rk = -(E2_rk/(1.34423*1.34423))*(a_rk*a_rk*a_rk*a_rk)*B_rk*(M_nu*M_nu/9.);
    if (fixed_neutrino_frac=='F') FNU = rk_1_3w*ON_rk/(OCB_rk+rk_1_3w*ON_rk);
    if (only_backg_rel == 'T') FNU = ON_rk/(OCB_rk+ON_rk);
    if (M_nu==0.) kj2_rk=1.;

    //printf("%lf %lf %e %e\n",zz-1.0,OCB(1./zz,E2_rk),E2_rk,E2);
    fprintf(fs,"%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n",
                zz-1.,sqrt(kj2_rk),OR_rk,OCB_rk,rk_1_3w,ON_rk,A_rk,B_rk,FNU);
  }


  fclose(fs);
}
