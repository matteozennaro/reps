#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"
extern void create_boltzmann_ini_file (char dir_chain[]);
extern void fscanf_error(int n);
extern double *allocate_double_vec(int n_elems);

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

  if (ps_norm_at_z == 0.)
  {
    char test='a';
    int n=0;
    sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
    n = count_lines(command);
    if ((n-1)!=knum)
    {
      printf("\nThe number of requested ks doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    double *P = allocate_double_vec(knum);
    double *Dcb = allocate_double_vec(knum);
    double *Dnu = allocate_double_vec(knum);
    double *Dm = allocate_double_vec(knum);

    int i;
    int fscanfcheck=0;
    char *dummy;
    double val,valk,Tc,Tb,Tn,Tm;
    double kpivot=0.05;

    char buf[5000];
    FILE *spectrumz0 = fopen(command,"r");
    if (spectrumz0 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
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

    int index_out=0;
    char outfile[200];
    char Dfile[200];
    printf("\n");
    for(index_out=0; index_out<output_number; index_out++)
    {
      sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Dcb[i]*Dcb[i]*P[i]);
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Dnu[i]*Dnu[i]*P[i]);
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
    free(Dcb);
    free(Dnu);
    free(Dm);
  }
  else
  {
    int n=0;
    sprintf(command,"%s/PK_TABS/Pcb_%s.txt",dir_chain,normfile);
    n = count_lines(command);
    if (n!=knum)
    {
      printf("\nThe number of ks requested doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    FILE *spectrumz99 = fopen(command,"r");
    if (spectrumz99 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    double *P = allocate_double_vec(knum);

    sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
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

    double *Pc = allocate_double_vec(knum);
    double *Pn = allocate_double_vec(knum);
    double *Pm = allocate_double_vec(knum);
    double *Dcb = allocate_double_vec(knum);
    double *Dnu = allocate_double_vec(knum);
    double *Dm = allocate_double_vec(knum);
    double *D99cb = allocate_double_vec(knum);
    double *D99nu = allocate_double_vec(knum);
    double *D99m = allocate_double_vec(knum);

    for(i=0; i<knum; i++)
    {
      fscanfcheck=fscanf(spectrumz99,"%lf %lf",&valk,&Pc[i]);
      if (fscanfcheck!=2) fscanf_error(2);
    }
    fclose(spectrumz99);

    sprintf(command,"%s/PK_TABS/Pnu_%s.txt",dir_chain,normfile);
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

    sprintf(command,"%s/PK_TABS/Pm_%s.txt",dir_chain,normfile);
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

    sprintf(D99file,"%s_znum%i.txt",outputfile,output_number-1);
    read_D(D99file,knum,D99cb,D99nu,D99m);

    printf("\n");
    for(index_out=0; index_out<output_number; index_out++)
    {
      sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Pc[i]*((Dcb[i]*Dcb[i])/(D99cb[i]*D99cb[i])));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Pn[i]*((Dnu[i]*Dnu[i])/(D99nu[i]*D99nu[i])));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pm_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledm = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Pm[i]*((Dm[i]*Dm[i])/(D99m[i]*D99m[i])));
      }
      fclose(outrescaledm);

      // norm z0

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],P[i]*(Dcb[i]*Dcb[i]));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],P[i]*(Dnu[i]*Dnu[i]));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
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

  if (ps_norm_at_z == 0.)
  {
    char test='a';
    int n=0;
    sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
    n = count_lines(command);
    if ((n)!=knum)
    {
      printf("\nThe number of requested ks doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    double *P = allocate_double_vec(knum);
    double *Dcb = allocate_double_vec(knum);
    double *Dnu = allocate_double_vec(knum);
    double *Dm = allocate_double_vec(knum);

    int i;
    int fscanfcheck=0;
    char *dummy;
    double val,valk,Tc,Tb,Tn,Tm;
    double kpivot=0.05;

    char buf[5000];
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
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Dcb[i]*Dcb[i]*P[i]);
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Dnu[i]*Dnu[i]*P[i]);
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
    free(Dcb);
    free(Dnu);
    free(Dm);
  }
  else
  {
    int n=0;
    sprintf(command,"%s/PK_TABS/Pcb_%s.txt",dir_chain,normfile);
    n = count_lines(command);
    if (n!=knum)
    {
      printf("\nThe number of ks requested doesn\'t match the one in the PS file!\n");
      exit(1);
    }
    FILE *spectrumz99 = fopen(command,"r");
    if (spectrumz99 == NULL)
    {
      printf("\nError opening file %s\n",command);
      exit(1);
    }
    double *P = allocate_double_vec(knum);

    sprintf(command,"%s/PK_TABS/power_norm_z1_pk.dat",dir_chain);
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

    double *Pc = allocate_double_vec(knum);
    double *Pn = allocate_double_vec(knum);
    double *Pm = allocate_double_vec(knum);
    double *Dcb = allocate_double_vec(knum);
    double *Dnu = allocate_double_vec(knum);
    double *Dm = allocate_double_vec(knum);
    double *D99cb = allocate_double_vec(knum);
    double *D99nu = allocate_double_vec(knum);
    double *D99m = allocate_double_vec(knum);

    for(i=0; i<knum; i++)
    {
      fscanfcheck=fscanf(spectrumz99,"%lf %lf",&valk,&Pc[i]);
      if (fscanfcheck!=2) fscanf_error(2);
    }
    fclose(spectrumz99);

    sprintf(command,"%s/PK_TABS/Pnu_%s.txt",dir_chain,normfile);
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

    sprintf(command,"%s/PK_TABS/Pm_%s.txt",dir_chain,normfile);
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

    sprintf(D99file,"%s_znum%i.txt",outputfile,output_number-1);
    read_D(D99file,knum,D99cb,D99nu,D99m);

    printf("\n");
    for(index_out=0; index_out<output_number; index_out++)
    {
      sprintf(Dfile,"%s_znum%i.txt",outputfile,index_out);
      read_D(Dfile,knum,Dcb,Dnu,Dm);

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],Pc[i]*((Dcb[i]*Dcb[i])/(D99cb[i]*D99cb[i])));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],Pn[i]*((Dnu[i]*Dnu[i])/(D99nu[i]*D99nu[i])));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm99_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pm_%s_rescaled_norm99_znum%i.txt\n",outputfile,index_out);
      FILE*outrescaledm = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledm,"%.12e\t%.12e\n",k[i],Pm[i]*((Dm[i]*Dm[i])/(D99m[i]*D99m[i])));
      }
      fclose(outrescaledm);

      // norm z0

      sprintf(outfile,"%s/PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pcb_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescaledcb = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescaledcb,"%.12e\t%.12e\n",k[i],P[i]*(Dcb[i]*Dcb[i]));
      }
      fclose(outrescaledcb);

      sprintf(outfile,"%s/PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pnu_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
      outrescalednu = fopen(outfile,"w");
      for (i=0; i<knum; i++)
      {
        fprintf(outrescalednu,"%.12e\t%.12e\n",k[i],P[i]*(Dnu[i]*Dnu[i]));
      }
      fclose(outrescalednu);

      sprintf(outfile,"%s/PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt",dir_chain,outputfile,index_out);
      printf("Written file PK_TABS/Pm_%s_rescaled_norm00_znum%i.txt\n",outputfile,index_out);
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
