#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"

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
