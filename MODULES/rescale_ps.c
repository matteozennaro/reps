#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"
extern void create_boltzmann_ini_file (char dir_chain[]);
extern void fscanf_error(int n);
extern double *allocate_double_vec(int n_elems);
extern int count_lines(char file[]);

void read_D(char filename[],int n, double *db, double *dc, double *dn, double *dm)
{
  FILE *f = fopen(filename,"r");
  if (f == NULL) {printf("Problem opening file %s!\n",filename); exit(1);}

  int i; double val;
  int fscanfcheck=0;
  for (i=0; i<n; i++)
  {
    fscanfcheck=fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&val,&db[i],&dc[i],&dn[i],&dm[i],&val,&val,&val,&val);
    if (fscanfcheck!=9) fscanf_error(9);
  }
  fclose(f);
}

void rescale_camb_ps(int knum, double *k)
{
  printf("\nRescaling of the PS requested.\n");

  char dir_chain[1000];
  if (getcwd(dir_chain,sizeof(dir_chain))==NULL)
  {
    printf("\nError retrieving current directory path.\n");
    exit(1);
  }

  create_boltzmann_ini_file(dir_chain);

  printf("\nCalling camb and generating the P(K) and T(k) at the \n"
         "requested output redshifts.\n");
  char command[200];
  sprintf(command,"%scamb %s/PK_TABS/power.ini > boltzmann.log",boltzmann_folder,dir_chain);
  system(command);

  // power spectra at z=0,99 used for normalization
  mode = 3;
  create_boltzmann_ini_file(dir_chain);

  sprintf(command,"%scamb %s/PK_TABS/power_norm.ini > boltzmann.log",boltzmann_folder,dir_chain);
  system(command);
  mode = 0;

  int n=0;
  sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
  n = count_lines(command);
  if ((n-1)!=knum)
  {
    printf("\nThe number of requested ks doesn\'t match the one in the PS file!\n");
    exit(1);
  }
  double *P = allocate_double_vec(knum);
  double *Dc = allocate_double_vec(knum);
  double *Db = allocate_double_vec(knum);
  double *Dn = allocate_double_vec(knum);
  double *Dm = allocate_double_vec(knum);

  int i;
  int fscanfcheck=0;
  char *dummy;
  double valk;

  char buf[5000];
  FILE *spectrumz0 = fopen(command,"r");
  if (spectrumz0 == NULL)
  {
    printf("\nError opening file %s\n",command);
    exit(1);
  }
  dummy=fgets (buf, sizeof(buf), spectrumz0);
  if (dummy==NULL) exit(-1);

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

  int index_out=0;
  char outfile[200];
  char Dfile[200];
  printf("\n");
  for(index_out=0; index_out<output_number; index_out++)
  {
    sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
    read_D(Dfile,knum,Db,Dc,Dn,Dm);

    sprintf(outfile,"%s/PK_TABS/Pb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledb = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledb,"%.12e\t%.12e\n",k[i],Db[i]*Db[i]*P[i]);
    }
    fclose(outrescaledb);

    sprintf(outfile,"%s/PK_TABS/Pc_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pc_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledc = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledc,"%.12e\t%.12e\n",k[i],Dc[i]*Dc[i]*P[i]);
    }
    fclose(outrescaledc);

    sprintf(outfile,"%s/PK_TABS/Pn_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pn_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescalednu = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Dn[i]*Dn[i]*P[i]);
    }
    fclose(outrescalednu);

    sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledm = fopen(outfile,"w");

    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Dm[i]*Dm[i]*P[i]);
    }
    fclose(outrescaledm);
  }

  free(P);
  free(Db);
  free(Dc);
  free(Dn);
  free(Dm);
}

void rescale_class_ps(int knum, double *k)
{
  printf("\nRescaling of the PS requested.\n");

  char dir_chain[1000];
  if (getcwd(dir_chain,sizeof(dir_chain))==NULL)
  {
    printf("\nError retrieving current directory path.\n");
    exit(1);
  }

  create_boltzmann_ini_file(dir_chain);

  printf("\nCalling class and generating the P(K) and T(k) at the \n"
         "requested output redshifts.\n");
  char command[200];
  sprintf(command,"%sclass %s/PK_TABS/power.ini > boltzmann.log",boltzmann_folder,dir_chain);
  system(command);

  // power spectra at z=0,99 used for normalization
  mode = 3;
  create_boltzmann_ini_file(dir_chain);

  sprintf(command,"%sclass %s/PK_TABS/power_norm.ini > boltzmann.log",boltzmann_folder,dir_chain);
  system(command);
  mode = 0;

  int n=0;
  sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
  n = count_lines(command);
  if ((n)!=knum)
  {
    printf("\nThe number of requested ks doesn\'t match the one in the PS file!\n");
    exit(1);
  }
  double *P = allocate_double_vec(knum);
  double *Dc = allocate_double_vec(knum);
  double *Db = allocate_double_vec(knum);
  double *Dn = allocate_double_vec(knum);
  double *Dm = allocate_double_vec(knum);

  int i;
  int fscanfcheck=0;
  double valk;

  FILE *spectrumz0 = fopen(command,"r");
  if (spectrumz0 == NULL)
  {
    printf("\nError opening file %s\n",command);
    exit(1);
  }

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

  int index_out=0;
  char outfile[200];
  char Dfile[200];
  printf("\n");
  for(index_out=0; index_out<output_number; index_out++)
  {
    sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
    read_D(Dfile,knum,Db,Dc,Dn,Dm);

    sprintf(outfile,"%s/PK_TABS/Pb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledb = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledb,"%.12e\t%.12e\n",k[i],Db[i]*Db[i]*P[i]);
    }
    fclose(outrescaledb);

    sprintf(outfile,"%s/PK_TABS/Pc_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pc_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledc = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledc,"%.12e\t%.12e\n",k[i],Dc[i]*Dc[i]*P[i]);
    }
    fclose(outrescaledc);

    sprintf(outfile,"%s/PK_TABS/Pn_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pn_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescalednu = fopen(outfile,"w");
    for (i=0; i<knum; i++)
    {
      fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Dn[i]*Dn[i]*P[i]);
    }
    fclose(outrescalednu);

    sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
    printf("Written file PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
    FILE*outrescaledm = fopen(outfile,"w");

    for (i=0; i<knum; i++)
    {
      fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Dm[i]*Dm[i]*P[i]);
    }
    fclose(outrescaledm);
  }

  free(P);
  free(Db);
  free(Dc);
  free(Dn);
  free(Dm);
}
