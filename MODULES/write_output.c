#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_global.h"
#include "background.h"
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

void print_camblike_transfer(char file[],int knum,double *k,double *Db,
                        double *Dc,double *Dn,double *Dm,double *fb,
                        double *fc,double *fn,double *fm,double *P0,
                        double z)
{
  // 1. Create power spectra for each species
  double Pb[knum];
  double Pc[knum];
  double Pn[knum];
  double Pm[knum];
  int i;
  for (i=0;i<knum;i++)
  {
    Pb[i] = Db[i]*Db[i]*P0[i];
    Pc[i] = Dc[i]*Dc[i]*P0[i];
    Pn[i] = Dn[i]*Dn[i]*P0[i];
    Pm[i] = Dm[i]*Dm[i]*P0[i];
  }
  // 2. Compute camb normalization and use it to compute transfers
  double Tb[knum];
  double Tc[knum];
  double Tn[knum];
  double Tm[knum];
  double Tcb;
  double Norm;
  for (i=0;i<knum;i++)
  {
    Norm = As*pow((k[i]*h)/kpivot,ns-1.0)*(2.0*M_PI*M_PI)*(h*h*h*h*k[i]);
    Tb[i] = sqrt(Pb[i]/Norm);
    Tc[i] = sqrt(Pc[i]/Norm);
    Tn[i] = sqrt(Pn[i]/Norm);
    Tm[i] = sqrt(Pm[i]/Norm);
  }

  // 3. Compute baryon and cdm velocities
  double Vb[knum];
  double Vc[knum];
  double a = 1./(1.+z);
  double H = 100.0*sqrt(E2(a,ONE2(a)));
  double DDb,DDc;
  for (i=0;i<knum;i++)
  {
    // INCORRECT!
    Norm = sqrt(As*pow((k[i]*h)/kpivot,ns-1.0)*(2.0*M_PI*M_PI)*pow(h,4)*k[i]);
    DDb = sqrt(Pb[i]);
    DDc = sqrt(Pc[i]);
    // DDb = Tb[i];
    // DDc = Tc[i];
    // sembra funzionare a z=99
    Vb[i] = pow(2.*M_PI,2)/pow(h,5./2.)* ((a*H*fb[i]*DDb /k[i]) * k[i]/H) /Norm ;
    Vc[i] = pow(2.*M_PI,2)/pow(h,5./2.)* ((a*H*fc[i]*DDc /k[i]) * k[i]/H) /Norm ;

    Vb[i] = (a*H*DDb*fb[i]) * ((k[i])/H) /pow(h*k[i],2);
    Vc[i] = (a*H*DDc*fc[i]) * ((k[i])/H) /pow(h*k[i],2);

    Vb[i] = fb[i]*Tb[i];
    Vc[i] = fc[i]*Tc[i];
  }

  // IMPORTANT !!!!
  // The following is relative to CAMB - version May 2016
  // Hopefully, these things won't change too much among different versions
  // of CAMB, but SOMETIMES they do.
  // Please check what your code really wants as an input (it might even be
  // differnt from the format of the latest CAMB realease) and, in case,
  // modify this part accordingly (then 'make clean' and 'make' again).
  // I hope this is general enough, so it should be simple to adjust to
  // everyone's needs.

  // Specify the total number of columns you want your file to have
  // Irrelevant columns will be set to 0
  int Number_of_columns_in_transfer_file = 13;

  // Now specify in which column you want each output
  // (Columns are numbered starting from 0)
  int k_col = 0;
  int cdm_transfer_col = 1;
  int baryon_transfer_col = 2;
  int neutrino_transfer_col = 5;
  int matter_transfer_col = 6;
  int cdm_baryon_transfer_col = 7;
  int cdm_vel_col = 10;
  int baryon_vel_col = 11;
  int rel_vel_col = 12;

  double **OutputTable = allocate_matrix(knum,Number_of_columns_in_transfer_file);
  int j;
  for (i=0;i<Number_of_columns_in_transfer_file;i++)
  {
    if (i==k_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = k[j];
      }
    }
    else if(i==cdm_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tc[j];
      }
    }
    else if(i==baryon_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tb[j];
      }
    }
    else if(i==neutrino_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tn[j];
      }
    }
    else if(i==matter_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tm[j];
      }
    }
    else if(i==cdm_baryon_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        Tcb = (OB0/(OB0+OC0))*Tb[j] + (OC0/(OB0+OC0))*Tc[j];
        OutputTable[j][i] = Tcb;
      }
    }
    else if(i==cdm_vel_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Vc[j];
      }
    }
    else if(i==baryon_vel_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Vb[j];
      }
    }
    else if(i==rel_vel_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Vb[j]-Vc[i];
      }
    }
    else
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = 0.0;
      }
    }
  }

  FILE * f = fopen(file,"w");
  if(f!=NULL)
  {
    for (i=0;i<knum;i++)
    {
        for (j=0;j<Number_of_columns_in_transfer_file;j++)
      {
        fprintf(f,"%.10e\t",OutputTable[i][j]);
      }
      fprintf(f,"\n");
    }
    fclose(f);
    if (verb>1) printf("Written file %s\n",file);
    deallocate_matrix(OutputTable,knum,Number_of_columns_in_transfer_file);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Creating the file %s\n"
           "raised a fatal exception\n",file);
    frame(error);
    exit(-1);
  }
}

void print_ngenic_transfer(char file[],int knum,double *k,double *Db,
                        double *Dc,double *Dn,double *Dm,double *P0)
{
  // In this case there are no velocities and transfers are defined
  // in another way

  // 1. Create power spectra for each species
  double Pb[knum];
  double Pc[knum];
  double Pn[knum];
  double Pm[knum];
  int i;
  for (i=0;i<knum;i++)
  {
    Pb[i] = Db[i]*Db[i]*P0[i];
    Pc[i] = Dc[i]*Dc[i]*P0[i];
    Pn[i] = Dn[i]*Dn[i]*P0[i];
    Pm[i] = Dm[i]*Dm[i]*P0[i];
  }
  // 2. Compute camb normalization and use it to compute transfers
  double Tb[knum];
  double Tc[knum];
  double Tn[knum];
  double Tm[knum];
  double Tcb;
  for (i=0;i<knum;i++)
  {
    Tb[i] = sqrt(Pb[i]/Pm[i]);
    Tc[i] = sqrt(Pc[i]/Pm[i]);
    Tn[i] = sqrt(Pn[i]/Pm[i]);
    Tm[i] = sqrt(Pm[i]/Pm[i]);
  }

  // IMPORTANT !!!!
  // The following is relative to NGenIC
  // Hopefully, these things don't change too much among different versions
  // of NGenIC, but SOMETIMES they do.
  // Please check what your code really wants as an input and, in case,
  // modify this part accordingly (then 'make clean' and 'make' again).
  // I hope this is general enough, so it should be simple to adjust to
  // everyone's needs.

  // Specify the total number of columns you want your file to have
  // Irrelevant columns will be set to 0
  int Number_of_columns_in_transfer_file = 7;

  // Now specify in which column you want each output
  // (Columns are numbered starting from 0)
  int k_col = 0;
  int neutrino_transfer_col = 5;
  int matter_transfer_col = 6;
  int cdm_baryon_transfer_col_1 = 1;
  int cdm_baryon_transfer_col_2 = 2;

  double **OutputTable = allocate_matrix(knum,Number_of_columns_in_transfer_file);
  int j;
  for (i=0;i<Number_of_columns_in_transfer_file;i++)
  {
    if (i==k_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = k[j];
      }
    }
    else if(i==neutrino_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tn[j];
      }
    }
    else if(i==matter_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tm[j];
      }
    }
    else if(i==cdm_baryon_transfer_col_1 || i==cdm_baryon_transfer_col_2)
    {
      for (j=0;j<knum;j++)
      {
        Tcb = (OB0/(OB0+OC0))*Tb[j] + (OC0/(OB0+OC0))*Tc[j];
        OutputTable[j][i] = Tcb;
      }
    }
    else
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = 0.0;
      }
    }
  }

  FILE * f = fopen(file,"w");
  if(f!=NULL)
  {
    for (i=0;i<knum;i++)
    {
        for (j=0;j<Number_of_columns_in_transfer_file;j++)
      {
        fprintf(f,"%.10e\t",OutputTable[i][j]);
      }
      fprintf(f,"\n");
    }
    fclose(f);
    if (verb>1) printf("Written file %s\n",file);
    deallocate_matrix(OutputTable,knum,Number_of_columns_in_transfer_file);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Creating the file %s\n"
           "raised a fatal exception\n",file);
    frame(error);
    exit(-1);
  }
}

void print_ngenic_old_transfer(char file[],int knum,double *k,double *Db,
                               double *Dc,double *Dn,double *Dm,double *P0)
{
  // In this case there are no velocities and transfers are defined
  // in another way

  // 1. Create power spectra for each species
  double Pb[knum];
  double Pc[knum];
  double Pn[knum];
  double Pm[knum];
  int i;
  for (i=0;i<knum;i++)
  {
    Pb[i] = Db[i]*Db[i]*P0[i];
    Pc[i] = Dc[i]*Dc[i]*P0[i];
    Pn[i] = Dn[i]*Dn[i]*P0[i];
    Pm[i] = Dm[i]*Dm[i]*P0[i];
  }
  // 2. Compute camb normalization and use it to compute transfers
  double Tb[knum];
  double Tc[knum];
  double Tn[knum];
  double Tm[knum];
  for (i=0;i<knum;i++)
  {
    Tb[i] = sqrt(Pb[i]/Pm[i]);
    Tc[i] = sqrt(Pc[i]/Pm[i]);
    Tn[i] = sqrt(Pn[i]/Pm[i]);
    Tm[i] = sqrt(Pm[i]/Pm[i]);
  }

  // IMPORTANT !!!!
  // The following is relative to NGenIC
  // Hopefully, these things don't change too much among different versions
  // of NGenIC, but SOMETIMES they do.
  // Please check what your code really wants as an input and, in case,
  // modify this part accordingly (then 'make clean' and 'make' again).
  // I hope this is general enough, so it should be simple to adjust to
  // everyone's needs.

  // Specify the total number of columns you want your file to have
  // Irrelevant columns will be set to 0
  int Number_of_columns_in_transfer_file = 7;

  // Now specify in which column you want each output
  // (Columns are numbered starting from 0)
  int k_col = 0;
  int neutrino_transfer_col = 5;
  int matter_transfer_col = 6;
  int cdm_transfer_col = 2;
  int baryon_transfer_col = 1;

  double **OutputTable = allocate_matrix(knum,Number_of_columns_in_transfer_file);
  int j;
  for (i=0;i<Number_of_columns_in_transfer_file;i++)
  {
    if (i==k_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = k[j];
      }
    }
    else if(i==neutrino_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tn[j];
      }
    }
    else if(i==matter_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tm[j];
      }
    }
    else if(i==cdm_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tc[j];
      }
    }
    else if(i==baryon_transfer_col)
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = Tb[j];
      }
    }
    else
    {
      for (j=0;j<knum;j++)
      {
        OutputTable[j][i] = 0.0;
      }
    }
  }

  FILE * f = fopen(file,"w");
  if(f!=NULL)
  {
    for (i=0;i<knum;i++)
    {
        for (j=0;j<Number_of_columns_in_transfer_file;j++)
      {
        fprintf(f,"%.10e\t",OutputTable[i][j]);
      }
      fprintf(f,"\n");
    }
    fclose(f);
    if (verb>1) printf("Written file %s\n",file);
    deallocate_matrix(OutputTable,knum,Number_of_columns_in_transfer_file);
  }
  else
  {
    char error[2000];
    sprintf(error,"Error! Creating the file %s\n"
           "raised a fatal exception\n",file);
    frame(error);
    exit(-1);
  }
}

void retrieve_Pm_z0(int knum,double *ktrue, double *k, double *P)
{
  if (compute_Pk_0=='T')
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
  else
  {
    int lines_tmp = count_lines(file_Pk_0_in);
    int headerlines_tmp = count_header_lines(file_Pk_0_in);
    int knum_tmp = lines_tmp-headerlines_tmp;
    double k_tmp[knum_tmp];
    double P_tmp[knum_tmp];
    read_power_spectrum(file_Pk_0_in,knum_tmp,k_tmp,P_tmp);
    int i;
    for(i=0;i<knum;i++)
    {
      k[i] = ktrue[i];
      P[i] = lin_interp_at(ktrue[i],k_tmp,P_tmp,knum_tmp);
    }
  }
}

void write_output(int knum, double *k,
             double **Delta_b,double **Delta_c,double **Delta_n,double **Delta_m,
             double **growth_b,double **growth_c,double **growth_n,double **growth_m)
{
  double ktmp[knum];
  double Pmz0[knum];
  retrieve_Pm_z0(knum,k,ktmp,Pmz0);

  char currentfile[1000];

  // int i;
  // for (i=0;i<knum;i++) printf("%e   %e\n",k[i],Delta_c[output_number-1][i]);

  int index_out;

  if(strcmp(output_format,"camb")==0)
  {
    for (index_out=0; index_out<output_number; index_out++)
    {
      sprintf(currentfile,"%s_rescaled_transfer_z%.4lf.txt",outputfile,z_output[index_out]);
      print_camblike_transfer(currentfile,knum,k,Delta_b[index_out],Delta_c[index_out],
                              Delta_n[index_out],Delta_m[index_out],growth_b[index_out],
                              growth_c[index_out],growth_n[index_out],
                              growth_m[index_out],Pmz0,z_output[index_out]);
    }
  }
  else if(strcmp(output_format,"ngenic_old")==0)
  {
    for (index_out=0; index_out<output_number; index_out++)
    {
      sprintf(currentfile,"%s_rescaled_transfer_z%.4lf.txt",outputfile,z_output[index_out]);
      print_ngenic_old_transfer(currentfile,knum,k,Delta_b[index_out],Delta_c[index_out],
                              Delta_n[index_out],Delta_m[index_out],Pmz0);
      sprintf(currentfile,"%s_Pm_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_m[index_out],Pmz0);

      sprintf(currentfile,"%s_fb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_b[index_out]);

      sprintf(currentfile,"%s_fc_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_c[index_out]);

      double growth_cb[knum];
      double Delta_cb,Delta_b_Delta_cb,Delta_c_Delta_cb;
      int i;
      for(i=0;i<knum; i++)
      {
        Delta_cb = (OB0/(OB0+OC0))*Delta_b[index_out][i] + (OC0/(OB0+OC0))*Delta_c[index_out][i];
        Delta_b_Delta_cb = Delta_b[index_out][i] / Delta_cb;
        Delta_c_Delta_cb = Delta_c[index_out][i] / Delta_cb;
        growth_cb[i] = Delta_b_Delta_cb*(OB0/(OB0+OC0))*growth_b[index_out][i] +
                       Delta_c_Delta_cb*(OC0/(OB0+OC0))*growth_c[index_out][i];
      }

      sprintf(currentfile,"%s_fcb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_cb);

      sprintf(currentfile,"%s_fn_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_n[index_out]);

      sprintf(currentfile,"%s_fm_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_m[index_out]);
    }
  }
  else if(strcmp(output_format,"ngenic")==0)
  {
    for (index_out=0; index_out<output_number; index_out++)
    {
      sprintf(currentfile,"%s_rescaled_transfer_z%.4lf.txt",outputfile,z_output[index_out]);
      print_ngenic_transfer(currentfile,knum,k,Delta_b[index_out],Delta_c[index_out],
                              Delta_n[index_out],Delta_m[index_out],Pmz0);

      sprintf(currentfile,"%s_Pm_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_m[index_out],Pmz0);

      sprintf(currentfile,"%s_Pb_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_b[index_out],Pmz0);

      sprintf(currentfile,"%s_Pc_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_c[index_out],Pmz0);

      sprintf(currentfile,"%s_Pn_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_n[index_out],Pmz0);

      double Delta_cb[knum];
      int i;
      for(i=0;i<knum; i++)
      Delta_cb[i] = (OB0/(OB0+OC0))*Delta_b[index_out][i] +
                    (OC0/(OB0+OC0))*Delta_c[index_out][i] ;
      sprintf(currentfile,"%s_Pcb_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_cb,Pmz0);

      sprintf(currentfile,"%s_fb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_b[index_out]);

      sprintf(currentfile,"%s_fc_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_c[index_out]);

      double growth_cb[knum];
      double Delta_b_Delta_cb,Delta_c_Delta_cb;
      for(i=0;i<knum; i++)
      {
        Delta_b_Delta_cb = Delta_b[index_out][i] / Delta_cb[i];
        Delta_c_Delta_cb = Delta_c[index_out][i] / Delta_cb[i];
        growth_cb[i] = Delta_b_Delta_cb*(OB0/(OB0+OC0))*growth_b[index_out][i] +
                       Delta_c_Delta_cb*(OC0/(OB0+OC0))*growth_c[index_out][i];
      }

      sprintf(currentfile,"%s_fcb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_cb);

      sprintf(currentfile,"%s_fn_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_n[index_out]);

      sprintf(currentfile,"%s_fm_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_m[index_out]);
    }
  }
  else
  {
    for (index_out=0; index_out<output_number; index_out++)
    {
      int i;
      sprintf(currentfile,"%s_Pb_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_b[index_out],Pmz0);

      sprintf(currentfile,"%s_Pc_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_c[index_out],Pmz0);

      sprintf(currentfile,"%s_Pn_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_n[index_out],Pmz0);

      sprintf(currentfile,"%s_Pm_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_m[index_out],Pmz0);

      double Delta_cb[knum];
      for(i=0;i<knum; i++)
      Delta_cb[i] = (OB0/(OB0+OC0))*Delta_b[index_out][i] +
                     (OC0/(OB0+OC0))*Delta_c[index_out][i];
      sprintf(currentfile,"%s_Pcb_rescaled_z%.4lf.txt",outputfile,z_output[index_out]);
      print_power_spectrum(currentfile,knum,k,Delta_cb,Pmz0);

      sprintf(currentfile,"%s_fb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_b[index_out]);

      sprintf(currentfile,"%s_fc_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_c[index_out]);

      sprintf(currentfile,"%s_fn_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_n[index_out]);

      sprintf(currentfile,"%s_fm_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_m[index_out]);

      double growth_cb[knum];
      double Delta_b_Delta_cb,Delta_c_Delta_cb;
      for(i=0;i<knum; i++)
      {
        Delta_b_Delta_cb = Delta_b[index_out][i] / Delta_cb[i];
        Delta_c_Delta_cb = Delta_c[index_out][i] / Delta_cb[i];
        growth_cb[i] = Delta_b_Delta_cb*(OB0/(OB0+OC0))*growth_b[index_out][i] +
                       Delta_c_Delta_cb*(OC0/(OB0+OC0))*growth_c[index_out][i];
      }

      sprintf(currentfile,"%s_fcb_z%.4lf.txt",outputfile,z_output[index_out]);
      print_growth_rate(currentfile,knum,k,growth_cb);
    }
  }
}
