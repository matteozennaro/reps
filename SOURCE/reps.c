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
#include "neutrino_distribution_function.h"
#include "read_ini_file.h"
#include "general_purpose.h"
#include "RungeKutta_solver.h"
#include "write_output.h"
#include "background.h"

/******************************************************************************/
/*        MAIN                                                                */
/******************************************************************************/
int main(int argc, char *argv[])
{
  // printf("Welcome in\n");
  // printf("                               REPS\n");
  // printf("REscaled Power Spectra for initial conditions with massive neutrinos\n\n");

  frame("Welcome in\n"
        "                               REPS\n"
        "REscaled Power Spectra for initial conditions with massive neutrinos\n\n");

  read_GG_FF_tabs();

  read_parameter_file(argv[1]);

  int lines = 0;
  int j;
  int fscanfcheck;

  lines = count_lines(input_file);

  double *k_requested = allocate_double_vec(lines);
  double *fb_init = allocate_double_vec(lines);
  double *fc_init = allocate_double_vec(lines);
  double *fn_init = allocate_double_vec(lines);
  double *beta_b = allocate_double_vec(lines);
  double *beta_n = allocate_double_vec(lines);

  FILE *fin = fopen(input_file,"r");
  if (fin!=NULL)
  {
    for (j = 0; j < lines; j++)
    {
      fscanfcheck=fscanf(fin,"%lf %lf %lf %lf %lf %lf",
        &k_requested[j],&beta_b[j],&beta_n[j],&fb_init[j],&fc_init[j],&fn_init[j]);
      if (fscanfcheck!=6) fscanf_error(6);
    }
    fclose(fin);
  }
  else
  {
    printf("Error opening file %s\n",input_file);
    exit(1);
  }

  double **Delta_b = allocate_matrix(output_number,lines);
  double **Delta_c = allocate_matrix(output_number,lines);
  double **Delta_n = allocate_matrix(output_number,lines);
  double **Delta_m = allocate_matrix(output_number,lines);
  double **growth_b = allocate_matrix(output_number,lines);
  double **growth_c = allocate_matrix(output_number,lines);
  double **growth_n = allocate_matrix(output_number,lines);
  double **growth_m = allocate_matrix(output_number,lines);

  RK(lines,k_requested,beta_b,beta_n,fb_init,fc_init,fn_init,
     Delta_b,Delta_c,Delta_n,Delta_m,
     growth_b,growth_c,growth_n,growth_m);

  write_output(lines,k_requested,
               Delta_b,Delta_c,Delta_n,Delta_m,
               growth_b,growth_c,growth_n,growth_m);

  if(print_hubble=='T')
    print_hubble_table();
  else
    printf("\nPrinting the table of H values was not requested. Skip.\n");

/*Freeing allocated memory*************************************************/
  deallocate_matrix(Delta_b,output_number,lines);
  deallocate_matrix(Delta_c,output_number,lines);
  deallocate_matrix(Delta_n,output_number,lines);
  deallocate_matrix(Delta_m,output_number,lines);
  deallocate_matrix(growth_b,output_number,lines);
  deallocate_matrix(growth_c,output_number,lines);
  deallocate_matrix(growth_n,output_number,lines);
  deallocate_matrix(growth_m,output_number,lines);

  free(k_requested);
  free(beta_b);
  free(beta_n);
  free(fb_init);
  free(fc_init);
  free(fn_init);
  free(z_output);
  free(ytab);
  free(GGtab);
  free(FFtab);

/* Exit ******************************************************************/
  return 0;
}
