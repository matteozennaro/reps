#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_global.h"

void frame(char str[])
{
  int add_newline=0;
  int i=0;
  while (str[i]!='\0') i++; if(str[i-1]!='\n') add_newline=1;
  printf("---------------------------------------------------------------------\n");
  printf("%s",str);
  if (add_newline==1) printf("\n");
  printf("---------------------------------------------------------------------\n");
}

double * allocate_double_vec(int n_elems)
{
  double *vec = (double*)malloc(n_elems*sizeof(double));
  if (vec==NULL)
  {
    frame("Bad memory allocation!\n");
    exit(-1);
  }
  else return vec;
}

double ** allocate_matrix(int rows, int columns)
{
  int i;
  double **m;
  m = malloc(rows*sizeof(double *));
  if (m==NULL)
  {
    frame("Bad memory allocation!\n");
    exit(-1);
  }
  for(i=0; i<rows; i++)
  {
    m[i] = malloc(columns*sizeof(double));
    if (m==NULL)
    {
      frame("Bad memory allocation!\n");
      exit(-1);
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

void reallocate_matrix(double **m,int rows,int newcols)
{
  int i;
  for (i=0; i<rows; i++)
  {
    if (realloc(m[i],newcols)==NULL);
    {
      frame("Some vector reallocation has been attempted\n"
            "but was rejected due to lack of memory.\n"
            "Please, try again with less accuracy.\n");
      exit(-1);
    }
  }
}

int find_z_bin(double Z, double *vec, int n_elems)
{
  int i = 0;
  if (Z < z_final || Z > z_initial) exit(-1);
  while (log(1.+Z) > vec[i])
  {
    i++;
    if(i>n_elems)
    {
      frame("Search exceded array dimension\n");
      exit(-1);
    }
  }
  return i;
}

double det(double **a_in,int n)
{
  int i,j,j1,j2;
  double determinant = 0;
  double **m = NULL;

  if (n < 1)
  {
    frame("You tried to compute the determinant of \n"
          "a matrix of negative order?\n");
    exit(-1);
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

void sort_double_vec(int n_elems, double * vec)
{
  double tmp;
  int i,j;
  for (i=0;i<n_elems;i++)
  {
    for (j=i;j<n_elems;j++)
    {
      if (vec[j] < vec[i])
      {
        tmp = vec[i];
        vec[i] = vec[j];
        vec[j] = tmp;
      }
    }
  }
}

int count_lines(char file[])
{
  int l=0;
  char test;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(fscanf(f,"%c",&test)!=EOF) if (test=='\n') l++;
    fclose(f);
  }
  else
  {
    char error[1000];
    sprintf(error,"Problem loading file %s!\n",file);
    frame(error);
    exit(-1);
  }
  return l;
}

int count_header_lines(char file[])
{
  int hdl=0;
  char buf[2000];
  char * token;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
      token = strtok(buf," \t");
      if (token[0]=='#' || token[0]=='/' || token[0]=='!') hdl++;
      else break;
    }
    fclose(f);
  }
  else
  {
    char error[1000];
    sprintf(error,"Problem loading file %s!\n",file);
    frame(error);
    exit(-1);
  }
  return hdl;
}

int count_number_of_columns(char file[], int number_of_header_lines)
{
  int ncol=0;
  int line=0;
  char buf[2000];
  char * token;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(line<number_of_header_lines)
    {
      if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
      line++;
    }
    if(fgets(buf,sizeof(buf),f)==NULL) exit(-1);
    token = strtok(buf," \t");
    if(token!=NULL)
    {
      while (token!=NULL)
      {
        ncol++;
        token = strtok(NULL," \t");
      }
    }
    fclose(f);
  }
  else
  {
    char error[1000];
    sprintf(error,"Problem loading file %s!\n",file);
    frame(error);
    exit(-1);
  }
  return ncol;
}

void fscanf_error(int n)
{
  char error[1000];
  sprintf(error,"Error reading file. %i values were expected.\n",n);
  frame(error);
  exit(-1);
}

double lin_interp_between(double x, double x0, double x1, double y0, double y1)
{
  return y1 + ((y1-y0)/(x1-x0))*(x-x1);
}

double lin_interp_at(double x0, double *x, double *y, int n)
{
  int i=0;

  if (x0 > x[n-1])
  {
    i=n-1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else if (x0 < x[0])
  {
    i=1;
    return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
  else
  {
    while (x0 > x[i]) i++;
    if (i==0) return y[i];
    else return ((y[i]-y[i-1])/(x[i]-x[i-1])*x0 + y[i] - (y[i]-y[i-1])/(x[i]-x[i-1])*x[i]);
  }
}
