#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_global.h"
#include "boltzmann_solver.h"
#include "general_purpose.h"

void read_power_spectrum(char file[],int knum,double *k,double *P)
{
  int tot_nl = count_lines(file);
  int headerlines = count_header_lines(file);
  int ncol = count_number_of_columns(file,headerlines);
  if (ncol!=2)
  {
    char error[2000];
    sprintf(error,"Error! The power spectrum file is expected to contain\n"
           "two columns: k [h Mpc^{-1}] and P(k) [h^{-3} Mpc^3]\n"
           "while %i columns were detected\n",ncol);
    frame(error);
    exit(-1);
  }

  if (tot_nl-headerlines!=knum)
  {
    frame("Error! The numebr of wavenumbers in your P(k) file\n"
          "doesn't match the expected one\n");
    exit(-1);
  }

  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[1000];
    int current_knum=0;
    int i;
    for (i=0;i<headerlines;i++) if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
    for (i=headerlines;i<tot_nl;i++)
    {
      if (current_knum==knum) exit(-1);
      if(fscanf(f,"%lf %lf",&k[current_knum],&P[current_knum])!=2) fscanf_error(2);
      current_knum++;
    }
    fclose(f);
  }
  else
  {
    char error[1000];
    sprintf(error,"Error! The file %s doesn't exist\n",file);
    frame(error);
    exit(-1);
  }
}

void print_power_spectrum(char file[],int knum,double *k,double *Delta,double *Pmz0)
{
  FILE *f = fopen(file,"w");
  if (f!=NULL)
  {
    int i;
    for (i=0; i<knum; i++)
    {
      fprintf(f,"%.10e\t%.10e\n",k[i],(Delta[i]*Delta[i])*Pmz0[i]);
    }
    fclose(f);
    if (verb>1) printf("Written file %s\n",file);
  }
  else
  {
    char error[1000];
    sprintf(error,"Error! Creating the file %s\n"
           "raised a fatal exception\n",file);
    frame(error);
    exit(-1);
  }
}

void print_growth_rate(char file[],int knum,double *k,double *fk)
{
  FILE *f = fopen(file,"w");
  if (f!=NULL)
  {
    int i;
    for (i=0; i<knum; i++)
    {
      fprintf(f,"%.10e\t%.10e\n",k[i],fk[i]);
    }
    fclose(f);
    if (verb>1) printf("Written file %s\n",file);
  }
  else
  {
    char error[1000];
    sprintf(error,"Error! Creating the file %s\n"
           "raised a fatal exception\n",file);
    frame(error);
    exit(-1);
  }
}

void retrieve_Pm_z0(int knum,double *k, double *P)
{
  char dir_chain[1000];
  if (getcwd(dir_chain,sizeof(dir_chain))==NULL)
  {
    frame("Error retrieving current directory path.\n");
    exit(-1);
  }

  create_boltzmann_ini_file(dir_chain);
  printf("\nCalling %s and generating the P(K) and T(k) at the \n"
         "requested output redshifts.\n",boltzmann_code);
  char command[1000];
  sprintf(command,"%s%s %s/PK_TABS/power.ini > boltzmann.log",
        boltzmann_folder,boltzmann_code,dir_chain);
  system(command);

  // power spectra at z=0,z_initial used for normalization
  mode = 3;
  create_boltzmann_ini_file(dir_chain);
  sprintf(command,"%s%s %s/PK_TABS/power_norm.ini > boltzmann.log",
          boltzmann_folder,boltzmann_code,dir_chain);
  system(command);
  mode = 0;

  char filename[500];
  sprintf(filename,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
  read_power_spectrum(filename,knum,k,P);
}

void write_output(int knum, double *k,
             double **Delta_b,double **Delta_c,double **Delta_n,double **Delta_m,
             double **growth_b,double **growth_c,double **growth_n,double **growth_m)
{
  double ktmp[knum];
  double Pmz0[knum];
  retrieve_Pm_z0(knum,ktmp,Pmz0);

  char currentfile[1000];

  // int i;
  // for (i=0;i<knum;i++) printf("%e   %e\n",k[i],Delta_c[output_number-1][i]);

  int index_out;
  for (index_out=0; index_out<output_number; index_out++)
  {
    sprintf(currentfile,"%s_Pb_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
    print_power_spectrum(currentfile,knum,k,Delta_b[index_out],Pmz0);

    sprintf(currentfile,"%s_Pc_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
    print_power_spectrum(currentfile,knum,k,Delta_c[index_out],Pmz0);

    sprintf(currentfile,"%s_Pn_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
    print_power_spectrum(currentfile,knum,k,Delta_n[index_out],Pmz0);

    sprintf(currentfile,"%s_Pm_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
    print_power_spectrum(currentfile,knum,k,Delta_m[index_out],Pmz0);

    sprintf(currentfile,"%s_fb_z%.4lf.txt",outputfile,z_output[index_out]);
    print_growth_rate(currentfile,knum,k,growth_b[index_out]);

    sprintf(currentfile,"%s_fc_z%.4lf.txt",outputfile,z_output[index_out]);
    print_growth_rate(currentfile,knum,k,growth_c[index_out]);

    sprintf(currentfile,"%s_fn_z%.4lf.txt",outputfile,z_output[index_out]);
    print_growth_rate(currentfile,knum,k,growth_n[index_out]);

    sprintf(currentfile,"%s_fm_z%.4lf.txt",outputfile,z_output[index_out]);
    print_growth_rate(currentfile,knum,k,growth_m[index_out]);
  }
}
