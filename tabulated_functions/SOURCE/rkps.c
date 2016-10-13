#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
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
#include "include.h"

/******************************************************************************/
/*        MAIN                                                                */
/******************************************************************************/
int main(int argc, char *argv[])
{
  read_GG_FF_tabs();

  read_parameter_file(argv[1]);

  int lines = 0;
  char test = 'a';
  int j;
  int fscanfcheck;

  lines = count_lines(input_file);

  double *k_requested = allocate_double_vec(lines);
  double *fc_init = allocate_double_vec(lines);
  double *fn_init = allocate_double_vec(lines);
  double *beta = allocate_double_vec(lines);

  FILE *fin = fopen(input_file,"r");
  if (fin!=NULL)
  {
    for (j = 0; j < lines; j++)
    {
      fscanfcheck=fscanf(fin,"%lf %lf %lf %lf",
        &k_requested[j],&beta[j],&fc_init[j],&fn_init[j]);
      if (fscanfcheck!=4) fscanf_error(4);
    }
    fclose(fin);
  }
  else
  {
    printf("Error opening file %s\n",input_file);
    exit(1);
  }

  //Allocating the two 4x4 matrixes
  matrix_M = allocate_matrix(4,4);
  matrix_M2 = allocate_matrix(4,4);

  RK(lines,k_requested,beta,fc_init,fn_init);

  if(do_rescaled_ps=='T')
  {
    if (strcmp(boltzmann_code,"camb")==0)
      rescale_camb_ps(lines,k_requested);
    else if (strcmp(boltzmann_code,"class")==0)
      rescale_class_ps(lines,k_requested);
    else
    {
      printf("Illegal Boltzmann code selected.\n");
      printf("The computation of the grwoth rates and factors\n");
      printf("has been performed, but the rescaling will be skipped.\n");
    }
  }
  else
    printf("\nRescaling was not requested. Skip.\n");

  if(print_hubble=='T')
    print_hubble_table();
  else
    printf("\nPrinting the table of H values was not requested. Skip.\n");

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

/*Returning null value*****************************************************/
  return 0;
}
