#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <ctype.h>

#include "include_extern.h"
extern double set_ON0();
extern int count_lines(char file[]);

void read_parameter_file(char parfile[])
{
  int i,j;
  int fscanfcheck=0;
  int Nread = 24;
  char do_rescaled_ps_str[100];
  char print_hubble_str[100];
  char param_names[Nread][100];

  sprintf(param_names[0],"input_file");
  sprintf(param_names[1],"outputfile");
  sprintf(param_names[2],"z_initial");
  sprintf(param_names[3],"z_final");
  sprintf(param_names[4],"output_number");
  sprintf(param_names[5],"z_output");
  sprintf(param_names[6],"h");
  sprintf(param_names[7],"OB0");
  sprintf(param_names[8],"OC0");
  sprintf(param_names[9],"OG0");
  sprintf(param_names[10],"M_nu");
  sprintf(param_names[11],"do_rescaled_ps");
  sprintf(param_names[12],"As");
  sprintf(param_names[13],"ns");
  sprintf(param_names[14],"tau_reio");
  sprintf(param_names[15],"kmax");
  sprintf(param_names[16],"N_nu");
  sprintf(param_names[17],"Neff");
  sprintf(param_names[18],"print_hubble");
  sprintf(param_names[19],"boltzmann_folder");
  sprintf(param_names[20],"wrong_nu");
  sprintf(param_names[21],"boltzmann_code");
  sprintf(param_names[22],"w0");
  sprintf(param_names[23],"wa");
  int check[Nread]; for (i=0; i<Nread; i++) check[i]=0;
  int check_output_number = 0;
  char test_str[100];
  char c;

  FILE *input = fopen(parfile,"r");
  if (input==NULL)
  {
    printf("You need to specify a parameter file containing the following (%i) entries:\n",Nread);
    for(i=0; i<Nread; i++)
    printf("   %s\n",param_names[i]);
    exit(1);
  }

  while (!feof(input))
  {
    c = getc(input);
    if (isspace(c)==0)
    {
      if (c!='#')
      {
        ungetc(c,input);
        fscanfcheck=fscanf(input,"%s", test_str);
        for (i = 0; i < Nread; i++)
        {
          if (strcmp(test_str,param_names[i]) == 0)
          {
            if (i==0)
            {
              fscanfcheck=fscanf(input,"%s",input_file);
              check[i] = 1;
            }

            if (i==1)
            {
              fscanfcheck=fscanf(input,"%s",outputfile);
              check[i] = 1;
            }

            if (i==2)
            {
              fscanfcheck=fscanf(input,"%lf",&z_initial);
              check[i] = 1;
            }

            if (i==3)
            {
              fscanfcheck=fscanf(input,"%lf",&z_final);
              check[i] = 1;
            }

            if (i==4)
            {
              fscanfcheck=fscanf(input,"%i",&output_number);
              if (output_number<=0)
              {
                printf("You should specify at least 1 output redshift.\n");
                exit(1);
              }
              else
              {
                check_output_number=1;
                check[i] = 1;
              }
            }

            if (i==5)
            {
              if (check_output_number==0)
              {
                printf("I need you specify the number of required outputs\n");
                printf("before you tell me the output redshifts.\n");
                exit(1);
              }
              else
              {
                z_output = malloc(output_number*sizeof(double));
                if (z_output==NULL)
                { printf("Bad memory allocation...\n"); exit(1); }
                for (j=0; j<output_number; j++)
                  fscanfcheck=fscanf(input,"%lf",&z_output[j]);
                check[i] = 1;
              }
            }

            if (i==6)
            {
              fscanfcheck=fscanf(input,"%lf",&h);
              check[i] = 1;
            }

            if (i==7)
            {
              fscanfcheck=fscanf(input,"%lf",&OB0);
              check[i] = 1;
            }

            if (i==8)
            {
             fscanfcheck=fscanf(input,"%lf",&OC0);
             check[i] = 1;
            }

            if (i==9)
            {
              fscanfcheck=fscanf(input,"%lf",&OG0);
              OG0/=(h*h);
              check[i] = 1;
            }

            if (i==10)
            {
              fscanfcheck=fscanf(input,"%lf",&M_nu);
              check[i] = 1;
            }

            if (i==11)
            {
              fscanfcheck=fscanf(input,"%s",do_rescaled_ps_str);
              if(strcmp(do_rescaled_ps_str,"T")==0||strcmp(do_rescaled_ps_str,"True")==0||strcmp(do_rescaled_ps_str,"t")==0||strcmp(do_rescaled_ps_str,"true")==0)
                do_rescaled_ps='T';
              else if(strcmp(do_rescaled_ps_str,"F")==0||strcmp(do_rescaled_ps_str,"False")==0||strcmp(do_rescaled_ps_str,"f")==0||strcmp(do_rescaled_ps_str,"false")==0)
                do_rescaled_ps='F';
              else
              {
                printf("\nYou wrote \'do_rescaled_ps %s\'.\n",do_rescaled_ps_str);
                printf("Please, set \'do_rescaled_ps\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==12)
            {
              fscanfcheck=fscanf(input,"%lf",&As);
              check[i] = 1;
            }

            if (i==13)
            {
              fscanfcheck=fscanf(input,"%lf",&ns);
              check[i] = 1;
            }

            if (i==14)
            {
              fscanfcheck=fscanf(input,"%lf",&tau_reio);
              check[i] = 1;
            }

            if (i==15)
            {
              fscanfcheck=fscanf(input,"%lf",&kmax);
              check[i] = 1;
            }

            if (i==16)
            {
              fscanfcheck=fscanf(input,"%lf",&N_nu);
              check[i] = 1;
            }

            if (i==17)
            {
              fscanfcheck=fscanf(input,"%lf",&Neff);
              check[i] = 1;
            }

            if (i==18)
            {
              fscanfcheck=fscanf(input,"%s",print_hubble_str);
              if(strcmp(print_hubble_str,"T")==0||strcmp(print_hubble_str,"True")==0||strcmp(print_hubble_str,"t")==0||strcmp(print_hubble_str,"true")==0)
                print_hubble='T';
              else if(strcmp(print_hubble_str,"F")==0||strcmp(print_hubble_str,"False")==0||strcmp(print_hubble_str,"f")==0||strcmp(print_hubble_str,"false")==0)
                print_hubble='F';
              else
              {
                printf("\nYou wrote \'print_hubble %s\'.\n",print_hubble_str);
                printf("Please, set \'print_hubble\' to either\n");
                printf("T,True,t,true or F,False,f,false. \n\n");
                exit(1);
              }
              check[i] = 1;
            }

            if (i==19)
            {
              fscanfcheck=fscanf(input,"%s",boltzmann_folder);
              check[i] = 1;
            }

            if (i==20)
            {
              fscanfcheck=fscanf(input,"%i",&wrong_nu);
              if (wrong_nu!=0 && wrong_nu!=1 && wrong_nu!=2)
              {
                printf("Invalid value of parameter wrong_nu.\n");
                printf("Please, choose among:\n");
                printf("  0] Correct solution\n");
                printf("  1] case05 - rel neutrinos only in background\n");
                printf("  2] case45 - no rel neutrinos\n");
              }
              check[i] = 1;
            }

            if (i==21)
            {
              fscanfcheck=fscanf(input,"%s",boltzmann_code);
              check[i] = 1;
            }

            if (i==22)
            {
              fscanfcheck=fscanf(input,"%lf",&w0);
              check[i] = 1;
            }

            if (i==23)
            {
              fscanfcheck=fscanf(input,"%lf",&wa);
              check[i] = 1;
            }
          }
        }
      }
      else while(c!='\n') c = getc(input);
    }
  }
  fclose(input);

  for(i=0; i<Nread; i++)
  {
    if (check[i]==0)
    {
      printf("You didn't specify a value for %s.\n",param_names[i]);
      exit(1);
    }
  }
  if (fscanfcheck==0) exit(-1);

  //if (wrong_nu==2) Neff = 3.046;
  OR0 = (Neff*(7./8.)*pow(4./11.,4./3.)+1.)*(OG0);

  OM0 = OB0+OC0+set_ON0();
  OX0 = 1.-OM0-OR0;

  printf("\nDerived parameters: OM0 = %lf, OX0 = %lf\n\n",OM0,OX0);
}
