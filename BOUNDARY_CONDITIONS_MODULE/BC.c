#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

/******************************************************************************/
/*    INITIAL SETTINGS - GLOBAL VARIABLES                                     */
/******************************************************************************/
#include "global_variables.h"

/******************************************************************************/
/*    DECLARATION OF FUNCTIONS                                                */
/******************************************************************************/
#include "boltzmann_solver.h"
#include "general_purpose.h"
#include "neutrino_distribution_function.h"
#include "read_ini_file.h"
#include "background.h"

void which_k (double *wn, int n_k, char filename[]);
void lsq(double *x,double **y,int znum,int kindex,double *m,double *q);
void read_ith_pk(double z, int n_k, double *PB, double *PC, double *PN, char psname[], char tname[]);
void num_deriv(int lines, double *FB, double *FC, double *FN, double *Pbminus, double *Pbplus, double *Pcminus, double *Pcplus, double *Pnminus, double *Pnplus, double zminus, double zplus);
void print_help();
int wrong_ic = 0;

/******************************************************************************/
/*        MAIN                                                                */
/******************************************************************************/
int main(int argc, char *argv[])
{
  if (argc!=3)
  {
    if (argc==2)
    {
      if (strcmp(argv[1], "-h")==0 || strcmp(argv[1],"--help")==0 ||
          strcmp(argv[1],"--h")==0 || strcmp(argv[1], "-help")==0 )
      {
        print_help();
        exit(-1);
      }
      printf("  Please, choose the kind of boundary conditions\n");
      printf("  required, between: \n");
      printf("    0] correct -numeric- from boltzmann code\n");
      printf("    1] fcb = Omega_m ^ 0.55\n");
      printf("    2] fnu = Omega_m ^ 0.55\n");
      printf("    3] fcb and fnu = Omega_m ^ 0.55\n");
      printf("\n  Your choice: ");
      scanf("%i",&wrong_ic); printf("\n");
      if (wrong_ic != 0 &&
          wrong_ic != 1 &&
          wrong_ic != 2 &&
          wrong_ic != 3)
      {
        printf("Error! Illegal value!\n");
        exit(-1);
      }
    }
    else
    {
        printf("Error! You should run ./BC --help and read\n");
        printf("       the help page to have info on the usage of\n");
        printf("       this code.\n");
        exit(-1);
    }
  }
  else
  {
    wrong_ic = atoi(argv[2]);
    if (wrong_ic != 0 &&
        wrong_ic != 1 &&
        wrong_ic != 2 &&
        wrong_ic != 3)
    {
      printf("Error! Illegal value!\n");
      exit(-1);
    }
  }

  read_GG_FF_tabs();
  read_parameter_file(argv[1]);

  printf("Generating the boundary conditions...\n");

  int bc_nstep=50;
  double zmin = z_initial-2.;
  double zmax = z_initial+2.;
  double xmin = log(1.+zmin);
  double xmax = log(1.+zmax);
  double bc_zz[50];
  double bc_step = (xmax-xmin)/(50.-1.);

  int i;
  for (i=0; i<50; i++)
    bc_zz[i] = exp(xmin + i*bc_step) - 1.;

  char dir_chain[1000];
  if (getcwd(dir_chain,sizeof(dir_chain))==NULL)
  {
    frame("Error retrieving current directory path.\n");
    exit(-1);
  }

  mode=1;
  create_boltzmann_ini_file(dir_chain);
  mode=0;

  printf("Calling %s and creating a tab of power spectra \n",boltzmann_code);
  printf("distributed around z=%lf...\n",z_initial);

  char command[200];
  sprintf(command,"%s%s %s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini > boltzmann.log",
                  boltzmann_folder,boltzmann_code,dir_chain);
  system(command);

  char powfile[200];
  char tfile[200];

  sprintf(powfile,"BOUNDARY_CONDITIONS_MODULE/tabs/power_z1_pk.dat");

  int knum = count_lines(powfile);
  int header_lines = count_header_lines(powfile);
  knum -= header_lines;

  double **Pb,**Pc,**Pn,*k;
  Pb = allocate_matrix(bc_nstep,knum);
  Pc = allocate_matrix(bc_nstep,knum);
  Pn = allocate_matrix(bc_nstep,knum);
  k = allocate_double_vec(knum);
  which_k(k,knum,powfile);

  for (i=0; i<bc_nstep; i++)
  {
    sprintf(powfile,"BOUNDARY_CONDITIONS_MODULE/tabs/power_z%i_pk.dat",i+1);
    sprintf(tfile,"BOUNDARY_CONDITIONS_MODULE/tabs/power_z%i_tk.dat",i+1);
    read_ith_pk(bc_zz[i],knum,Pb[i],Pc[i],Pn[i],powfile,tfile);
  }

  double z_der[bc_nstep-1];
  double **fb, **fc, **fn;

  fb = allocate_matrix(bc_nstep-1,knum);
  fc = allocate_matrix(bc_nstep-1,knum);
  fn = allocate_matrix(bc_nstep-1,knum);

  for (i=0; i<(bc_nstep-1); i++)
  {
    z_der[i] = 0.5*(bc_zz[i]+bc_zz[i+1]);
    num_deriv(knum,fb[i],fc[i],fn[i],Pb[i],Pb[i+1],Pc[i],Pc[i+1],Pn[i],Pn[i+1],bc_zz[i],bc_zz[i+1]);
  }

  double m_b[knum];
  double q_b[knum];
  double m_c[knum];
  double q_c[knum];
  double m_nu[knum];
  double q_nu[knum];

  for (i=0; i<knum; i++)
  {
    lsq(z_der,fb,bc_nstep-1,i,m_b,q_b);
    lsq(z_der,fc,bc_nstep-1,i,m_c,q_c);
    lsq(z_der,fn,bc_nstep-1,i,m_nu,q_nu);
  }

  mode=2;
  create_boltzmann_ini_file(dir_chain);
  mode=0;

  printf("Calling %s for computing beta at z=%lf...\n",boltzmann_code,z_initial);

  sprintf(command,"%s%s %s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini > boltzmann.log",
                  boltzmann_folder,boltzmann_code,dir_chain);
  system(command);

  double PPb[knum],PPc[knum],PPn[knum];
  read_ith_pk(99., knum, PPb, PPc, PPn, "./BOUNDARY_CONDITIONS_MODULE/tabs/power_zin_pk.dat", "./BOUNDARY_CONDITIONS_MODULE/tabs/power_zin_tk.dat");

  double beta_b[knum];
  double beta_n[knum];

  for (i=0; i<knum; i++) beta_b[i]=sqrt(PPb[i]/PPc[i]);
  for (i=0; i<knum; i++) beta_n[i]=sqrt(PPn[i]/PPc[i]);

  double A, OM99;
  char outfinal_file[300];
  sprintf(outfinal_file,"%s",input_file);
  FILE *outfinal = fopen(outfinal_file,"w");
  if (outfinal==NULL)
  {
    printf("Error creating file %s\n",outfinal_file);
    exit(-1);
  }
  if (wrong_ic==0)
  {
    for (i=0; i<knum; i++)
    {
      fprintf(outfinal,"%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
      k[i],beta_b[i],beta_n[i],m_b[i]*z_initial+q_b[i],m_c[i]*z_initial+q_c[i],m_nu[i]*z_initial+q_nu[i]);
    }
  }
  else if (wrong_ic==1)
  {
    A = 1./(1.+z_output[output_number-1]);
    OM99 = pow((OM0/(A*A*A))/E2(A,ONE2(A)),0.55);
    for (i=0; i<knum; i++)
    {
      fprintf(outfinal,"%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
      k[i],beta_b[i],beta_n[i],OM99,OM99,m_nu[i]*z_initial+q_nu[i]);
    }
  }
  else if (wrong_ic==2)
  {
    A = 1./(1.+z_output[output_number-1]);
    OM99 = pow((OM0/(A*A*A))/E2(A,ONE2(A)),0.55);
    for (i=0; i<knum; i++)
    {
      fprintf(outfinal,"%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
      k[i],beta_b[i],beta_n[i],m_b[i]*z_initial+q_b[i],m_c[i]*z_initial+q_c[i],OM99);
    }
  }
  else if (wrong_ic==3)
  {
    A = 1./(1.+z_output[output_number-1]);
    OM99 = pow((OM0/(A*A*A))/E2(A,ONE2(A)),0.55);
    for (i=0; i<knum; i++)
    {
      fprintf(outfinal,"%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
      k[i],beta_b[i],beta_n[i],OM99,OM99,OM99);
    }
  }

  fclose(outfinal);

  printf("Boundary conditions written in file %s\n",input_file);
  system("rm -rf BOUNDARY_CONDITIONS_MODULE/tabs");

  deallocate_matrix(Pb,bc_nstep,knum);
  deallocate_matrix(Pc,bc_nstep,knum);
  deallocate_matrix(Pn,bc_nstep,knum);
  free(k);
  deallocate_matrix(fb,bc_nstep-1,knum);
  deallocate_matrix(fc,bc_nstep-1,knum);
  deallocate_matrix(fn,bc_nstep-1,knum);

  return 0;
}

void which_k (double *wn, int n_k, char filename[])
{
  int i; double val;
  int fscanfcheck=0;
  char buf[5000];
  char *dummy;
  int header_lines = count_header_lines(filename);
  int ncol = count_number_of_columns(filename,header_lines);
  if (ncol!=2)
  {
    char error[2000];
    sprintf(error,"Error! The power spectrum file is expected to contain\n"
           "two columns: k [h Mpc^{-1}] and P(k) [h^{-3} Mpc^3]\n"
           "while %i columns were detected\n",ncol);
    frame(error);
    exit(-1);
  }
  int hdl=0;

  FILE *f = fopen(filename,"r");
  while(hdl<header_lines)
  {
    dummy=fgets (buf, sizeof(buf), f);
    if (dummy==NULL) exit(-1);
    hdl++;
  }

  for (i=0; i<n_k; i++)
  {
    fscanfcheck=fscanf(f,"%lf %lf",&wn[i],&val);
    if (fscanfcheck!=2) fscanf_error(2);
  }
  fclose(f);
}

void lsq(double *x,double **y,int znum,int kindex,double *m,double *q)
{
  double Sum_x = 0;
  double Sum_y = 0;
  double Sum_x2 = 0;
  double Sum_x_y = 0;

  int i;
  for (i=0; i < znum; i++)
  {
    Sum_x += x[i];
    Sum_y += y[i][kindex];
    Sum_x_y += x[i]*y[i][kindex];
    Sum_x2 += pow(x[i],2);
  }

  int N = znum;

  m[kindex] = (N*Sum_x_y-Sum_x*Sum_y)/(N*Sum_x2-pow(Sum_x,2));

  q[kindex] = (Sum_y*Sum_x2-Sum_x*Sum_x_y)/(N*Sum_x2-pow(Sum_x,2));
}

void num_deriv(int lines, double *FB, double *FC, double *FN, double *Pbminus, double *Pbplus, double *Pcminus, double *Pcplus, double *Pnminus, double *Pnplus, double zminus, double zplus)
{
  int i;
  for (i=0; i < lines; i++)
  {
    FB[i] = (log(sqrt(Pbplus[i]/Pbminus[i]))/log((1.+zminus)/(1.+zplus)));
    FC[i] = (log(sqrt(Pcplus[i]/Pcminus[i]))/log((1.+zminus)/(1.+zplus)));
    if (N_nu!=0) FN[i] = (log(sqrt(Pnplus[i]/Pnminus[i]))/log((1.+zminus)/(1.+zplus))); else FN[i] = 0.;
  }
}

void read_ith_pk(double z, int n_k, double *PB, double *PC, double *PN, char psname[], char tname[])
{
  int fscanfcheck=0;
  int ncol = count_number_of_columns(tname,count_header_lines(tname));

  FILE *ft = fopen(tname,"r");
  if (ft == NULL)
  {
    char error[1000];
    sprintf(error,"Problem opening file %s!\n",tname);
    frame(error);
    exit(-1);
  }

  int i; double val;
  double Tc[n_k];
  double Tb[n_k];
  double Tn[n_k];
  double k[n_k];

  char buf[5000];
  char *dummy;

  dummy = fgets(buf, sizeof(buf), ft);
  if (dummy==NULL) exit(-1);

  if (strcmp(boltzmann_code,"camb")==0)
  {
    if (ncol!=13)
    {
      char error[1000];
      sprintf(error,"Error! In the transfer function file there are\n"
             "%i columns, while 13 were expected.\n"
             "Please, check this before continuing.\n",ncol);
      frame(error);
      exit(-1);
    }
    for(i = 0; i < n_k; i++)
    {
      fscanfcheck=fscanf(ft,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      &k[i],&Tc[i],&Tb[i],&val,&val,&Tn[i],&val,&val,&val,&val,&val,&val,&val);
      if (fscanfcheck!=13) fscanf_error(13);
    }

    for (i=0; i < n_k; i++)
    {
      PB[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI*k[i]*(h*h*h*h))*(Tb[i]*Tb[i]);
      PC[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI*k[i]*(h*h*h*h))*(Tc[i]*Tc[i]);
      PN[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI*k[i]*(h*h*h*h))*(Tn[i]*Tn[i]);
    }
  }
  else
  {
    if (N_nu==0)
    {
      if (ncol!=6)
      {
        char error[1000];
        sprintf(error,"Error! In the transfer function file there are\n"
               "%i columns, while 6 were expected.\n"
               "Please, check this before continuing.\n",ncol);
        frame(error);
        exit(-1);
      }
      for(i = 0; i < n_k; i++)
      {
        fscanfcheck=fscanf(ft,"%lf %lf %lf %lf %lf %lf",&k[i],&val,&Tb[i],&Tc[i],&val,&val);
        Tn[i] = 0.0;
      }
    }
    else
    {
      if (ncol!=9)
      {
        char error[1000];
        sprintf(error,"Error! In the transfer function file there are\n"
               "%i columns, while 9 were expected.\n"
               "Please, check this before continuing.\n",ncol);
        frame(error);
        exit(-1);
      }
      for(i = 0; i < n_k; i++)
      {
        fscanfcheck=fscanf(ft,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&k[i],&val,&Tb[i],&Tc[i],&val,&Tn[i],&val,&val,&val);
      }
    }
    for (i=0; i < n_k; i++)
    {
      PB[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI/(k[i]*k[i]*k[i]))*(Tb[i]*Tb[i]);
      PC[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI/(k[i]*k[i]*k[i]))*(Tc[i]*Tc[i]);
      PN[i] = As*pow((k[i]*h)/(kpivot),ns-1.)*(2.0*M_PI*M_PI/(k[i]*k[i]*k[i]))*(Tn[i]*Tn[i]);
    }
  }
  fclose(ft);
}

void print_help()
{
  printf("\n\tUsage:\n");
  printf("\t./BC param_file wrong_ic\n");
  printf("\n\tparam_file = the same parameter file needed for rkps\n");
  printf("\twrong_ic = a number chosen among:\n");
  printf("\t  0] correct -numeric- from boltzmann code\n");
  printf("\t  1] fcb = Omega_m ^ 0.55\n");
  printf("\t  2] fnu = Omega_m ^ 0.55\n");
  printf("\t  3] fcb and fnu = Omega_m ^ 0.55\n");
  printf("\n");
}
