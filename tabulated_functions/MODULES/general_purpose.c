#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"

double * allocate_double_vec(int n_elems)
{
  double *vec = (double*)malloc(n_elems*sizeof(double));
  if (vec==NULL)
  {
    printf("Bad memory allocation!\n");
    exit(1);
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

int find_z_bin(double Z, double *vec, int n_elems)
{
  int i = 0;
  if (Z < z0 || Z > z1) exit(1);
  while (log(1.+Z) > vec[i]) {i++; if(i>n_elems) {printf("Search exceded array dimension\n"); exit(1);} }
  return i;
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

int count_lines(char file[])
{
  int l=0;
  char test;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    while(fscanf(f,"%c",&test)!=EOF) if (test=='\n') l++;
  }
  else
  {
    printf("Problem loading file %s!\n",file);
    exit(1);
  }
  return l;
}

void fscanf_error(int n)
{
  printf("Error reading file. %i values were expected.\n",n);
  exit(1);
}

double lin_interp_between(double x, double x0, double x1, double y0, double y1)
{
  return y1 + ((y1-y0)/(x1-x0))*(x-x1);
}
