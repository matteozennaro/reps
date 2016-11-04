#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <ctype.h>

#include "include_global.h"
#include "background.h"
#include "general_purpose.h"

void read_err(char paramname[])
{
  char error[1000];
  sprintf(error,"Error! Parameter %s hasn't been specified\n",paramname);
  frame(error);
  exit(-1);
}

int read_double_from_file(char file[], char paramname[], double *val)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[1000];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) return reading_success;
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              *val = atof(token);
              reading_success = 1;
              if (verb>0)
              printf("  %s = %lf\n",paramname,*val);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    exit(-1);
  }
  return reading_success;
}

int read_int_from_file(char file[], char paramname[], int *val)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[1000];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) return reading_success;
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              *val = atoi(token);
              reading_success = 1;
              if (verb>0)
              printf("  %s = %i\n",paramname,*val);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    exit(-1);
  }
  return reading_success;
}

int read_bool_from_file(char file[], char paramname[], char *truth)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[1000];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) return reading_success;
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              if (token[0]=='T' || token[0]=='t')
                *truth = 'T';
              else
                *truth = 'F';
              if (verb>0)
              printf("  %s = %c\n",paramname,*truth);
              reading_success = 1;
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    exit(-1);
  }
  return reading_success;
}

int read_string_from_file(char file[], char paramname[], char *str)
{
  int reading_success=0;
  FILE *f = fopen(file,"r");
  if (f!=NULL)
  {
    char buf[1000];
    char *token;

    while(!feof(f))
    {
      if(fgets(buf,sizeof(buf),f)==NULL) return reading_success;
      token = strtok(buf," \t=\n");
      if (token!=NULL)
      {
        if(token[0]!='#')
        {
          while (token!=NULL)
          {
            if (strcmp(token,paramname)==0)
            {
              token = strtok(NULL," \t=\n");
              sprintf(str,"%s",token);
              reading_success = 1;
              if (verb>0)
              printf("  %s = %s\n",paramname,str);
              break;
            }
            else
            {
              token = strtok(NULL," \t=\n");
            }
          }
        }
      }
      if (reading_success==1) break;
    }
    fclose(f);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Unable to open param file %s\n"
          "Make sure you have specified a valid param file\n",file);
    frame(error);
    exit(-1);
  }
  return reading_success;
}

int read_output_redshifts(char file[])
{
  int reading_success = 0;
  int index_out;
  double *z_tmp;
  if(read_int_from_file(file,"output_number",&output_number)!=1)
  {
    printf("warning:: Auto-selecting output redshifts\n");
    output_number=2;
    z_tmp = allocate_double_vec(2);
    z_tmp[0]=z_final;
    z_tmp[1]= z_initial;
  }
  else
  {
    z_tmp = allocate_double_vec(output_number);

    FILE *f = fopen(file,"r");
    if (f!=NULL)
    {
      char buf[1000];
      char *token;

      while(!feof(f))
      {
        if(fgets(buf,sizeof(buf),f)==NULL) return reading_success;
        token = strtok(buf," \t=\n");
        if (token!=NULL)
        {
          if(token[0]!='#')
          {
            while (token!=NULL)
            {
              if (strcmp(token,"z_output")==0)
              {
                for (index_out=0; index_out<output_number; index_out++)
                {
                  token = strtok(NULL," \t=\n");
                  z_tmp[index_out] = atof(token);
                }
                reading_success = 1;
                break;
              }
              else
              {
                token = strtok(NULL," \t=\n");
              }
            }
          }
        }
        if (reading_success==1) break;
      }
      fclose(f);
    }
    else
    {
      char error[2000];
      sprintf(error,"Error! Unable to open param file %s\n"
            "Make sure you have specified a valid param file\n",file);
      frame(error);
      exit(-1);
    }
  }

  sort_double_vec(output_number,z_tmp);

  int add_z_initial = 0;
  int add_z_final = 0;

  if (z_tmp[output_number-1]!=z_initial)
  {
    output_number++;
    add_z_initial=1;
  }

  if (z_tmp[0]!=z_final)
  {
    output_number++;
    add_z_final=1;
  }

  z_output = allocate_double_vec(output_number);
  for (index_out=0;index_out<(output_number-add_z_initial-add_z_final);index_out++)
  {
    z_output[index_out+add_z_final] = z_tmp[index_out];
  }
  if (add_z_final==1) z_output[0] = z_final;
  if (add_z_initial==1) z_output[output_number-1] = z_initial;

  if (verb>0)
  {
    printf("  Updated output_number = %i\n",output_number);
    printf("  Redshifts:\n  ");
    for(index_out=0;index_out<output_number;index_out++)
    printf("%.4lf ",z_output[index_out]);
    printf("\n");
  }

  free(z_tmp);
  return reading_success;
}

void read_parameter_file(char parfile[])
{
  if(read_int_from_file(parfile,"verb",&verb)!=1) verb=1;

  if(read_string_from_file(parfile,"input_file",input_file)!=1) read_err("input_file");
  if(read_string_from_file(parfile,"outputfile",outputfile)!=1) read_err("outputfile");
  if(read_string_from_file(parfile,"output_format",output_format)!=1) read_err("output_format");
  if(read_string_from_file(parfile,"boltzmann_code",boltzmann_code)!=1) read_err("boltzmann_code");
  if(read_string_from_file(parfile,"boltzmann_folder",boltzmann_folder)!=1) read_err("boltzmann_folder");

  if(read_double_from_file(parfile,"z_initial",&z_initial)!=1) read_err("z_initial");
  if(read_double_from_file(parfile,"z_final",&z_final)!=1) read_err("z_final");
  if(read_double_from_file(parfile,"h",&h)!=1) read_err("h");
  if(read_double_from_file(parfile,"OB0",&OB0)!=1) read_err("OB0");
  if(read_double_from_file(parfile,"OC0",&OC0)!=1) read_err("OC0");
  if(read_double_from_file(parfile,"OG0",&OG0)!=1) read_err("OG0");
  if(read_double_from_file(parfile,"M_nu",&M_nu)!=1) read_err("M_nu");
  if(read_double_from_file(parfile,"As",&As)!=1) read_err("As");
  if(read_double_from_file(parfile,"ns",&ns)!=1) read_err("ns");
  if(read_double_from_file(parfile,"tau_reio",&tau_reio)!=1) read_err("tau_reio");
  if(read_double_from_file(parfile,"kmax",&kmax)!=1) read_err("kmax");
  if(read_double_from_file(parfile,"N_nu",&N_nu)!=1) read_err("N_nu");
  if(read_double_from_file(parfile,"Neff",&Neff)!=1) read_err("Neff");
  if(read_double_from_file(parfile,"w0",&w0)!=1) read_err("w0");
  if(read_double_from_file(parfile,"wa",&wa)!=1) read_err("wa");
  if(read_double_from_file(parfile,"kpivot",&kpivot)!=1) kpivot=0.05;

  if(read_int_from_file(parfile,"k_per_logint_camb",&k_per_logint_camb)!=1)
  k_per_logint_camb = 10;

  if(read_int_from_file(parfile,"wrong_nu",&wrong_nu)!=1) read_err("wrong_nu");

  if(read_bool_from_file(parfile,"print_hubble",&print_hubble)!=1) read_err("print_hubble");

  if (read_bool_from_file(parfile,"compute_Pk_0",&compute_Pk_0)==1)
  {
    if(compute_Pk_0=='F')
    if (read_string_from_file(parfile,"file_Pk_0_in",file_Pk_0_in)!=1)
    read_err("file_Pk_0_in");
  }
  else compute_Pk_0='T';

  read_output_redshifts(parfile);

  OG0/=(h*h);

  OR0 = (Neff*(7./8.)*pow(4./11.,4./3.)+1.)*(OG0);

  OM0 = OB0+OC0+set_ON0();
  OX0 = 1.-OM0-OR0;

  printf("\nDerived parameters: OM0 = %lf, OX0 = %lf\n\n",OM0,OX0);
}
