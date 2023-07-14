#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_global.h"
#include "background.h"

void create_boltzmann_ini_file (char dir_chain[])
{
  double bc_zz[50];
  if (mode==1)
  {
    int i_bc;
    double zmin = z_initial-2.;
    double zmax = z_initial+2.;

    double xmin = log(1.+zmin);
    double xmax = log(1.+zmax);
    double bc_step = (xmax-xmin)/(50.-1.);

    for (i_bc=0; i_bc<50; i_bc++)
      bc_zz[i_bc] = exp(xmin + i_bc*bc_step) - 1.;
  }

  if (strcmp(boltzmann_code,"camb")==0)
  {
    char fileopen[200];

    if (mode==0)
    {
      printf("Creating %s input_file\n   %s/PK_TABS/power.ini\n",boltzmann_code,dir_chain);
      sprintf(fileopen,"%s/PK_TABS/power.ini",dir_chain);
    }
    if (mode==1 || mode==2)
    {
      printf("Creating %s input_file\n   %s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini\n",boltzmann_code,dir_chain);
      sprintf(fileopen,"%s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini",dir_chain);
    }
    if (mode==3)
    {
      sprintf(fileopen,"%s/PK_TABS/power_norm.ini",dir_chain);
    }

    double Nur,N_eigenstates,N_degeneracies;
    if (M_nu==0.)
    {
      Nur = 3.046;
      N_eigenstates = 0.0;
      N_degeneracies = 0.0;
    }
    else
    {
      Nur = 3.0 - (double)N_nu + 0.046;
      N_eigenstates = 1.0;
      N_degeneracies = (double)N_nu;
    }

    int i;
    char command_mkdir[200];

    FILE *out = fopen(fileopen,"w");
    if (out == NULL)
    {
      if (mode==0 || mode==3)
      {
        printf("Does the folder PK_TABS exist? Creating...\n");
        sprintf(command_mkdir,"mkdir %s/PK_TABS",dir_chain);
        system(command_mkdir);
        out = fopen(fileopen,"w");
        if (out==NULL)
        {
          printf("Error creating params file...\n");
          exit(1);
        }
      }
      else
      {
        printf("Does the folder BOUNDARY_CONDITIONS_MODULE/tabs exist? Creating...\n");
        sprintf(command_mkdir,"mkdir %s/BOUNDARY_CONDITIONS_MODULE",dir_chain);
        system(command_mkdir);
        sprintf(command_mkdir,"mkdir %s/BOUNDARY_CONDITIONS_MODULE/tabs",dir_chain);
        system(command_mkdir);
        out = fopen(fileopen,"w");
        if (out==NULL)
        {
          printf("Error creating params file...\n");
          exit(1);
        }
      }
    }

    int znum=0;

    if (mode==0) {znum=output_number; fprintf(out,"output_root = %s/PK_TABS/power\n",dir_chain);}
    if (mode==1) {znum=50; fprintf(out,"output_root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power\n",dir_chain);}
    if (mode==2) {znum=1; fprintf(out,"output_root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power_zin\n",dir_chain);}
    if (mode==3) {znum=2; fprintf(out,"output_root = %s/PK_TABS/power_norm\n",dir_chain);}

    fprintf(out,
    "get_scalar_cls = F \n"
    "get_vector_cls = F\n"
    "get_tensor_cls = F\n"
    "get_transfer = T\n"
    "do_lensing= F\n"
    "do_nonlinear = 0\n"
    "l_max_scalar	   = 2000\n"
    "k_eta_max_scalar  = 4000\n"
    "l_max_tensor	   = 1500\n"
    "k_eta_max_tensor  = 3000\n"
    "												     \n"
    "#Main cosmological parameters, neutrino masses are assumed degenerate\n"
    "# If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k\n"
    "use_physical	= T\n"
    "												     \n"
    "omk= 0.0\n"
    "hubble= %lf\n"
    "w= %lf\n"
    "wa = %lf\n"
    "cs2_lam= 1\n"
    "ombh2= %.10lf\n"
    "omch2= %.10lf\n"
    "omega_lambda= %.10lf\n"
    "omnuh2 = %.10lf\n"
    "												     \n"
    "temp_cmb= 2.7255\n"
    "helium_fraction= 0.24\n"
    "massless_neutrinos = %lf\n"
    "massive_neutrinos  = %i\n"
    "												     \n"
    "												     \n"
    "share_delta_neff = F \n"
    "nu_mass_eigenstates = %i\n"
    "nu_mass_degeneracies = %i\n"
    "nu_mass_fractions = 1\n"
    "												     \n"
    "initial_power_num = 1\n"
    "pivot_scalar = 0.05\n"
    "pivot_tensor = 0.05\n"
    "scalar_amp(1) = %e\n"
    "scalar_spectral_index(1) = %.8lf\n"
    "scalar_nrun(1) = 0\n"
    "tensor_spectral_index(1)  = 0\n"
    "initial_ratio(1)= 1\n"
    "												     \n"
    "reionization= T\n"
    "												     \n"
    "re_use_optical_depth = T\n"
    "re_optical_depth= 0.0925\n"
    "re_redshift= 11\n"
    "re_delta_redshift= 0.5\n"
    "re_ionization_frac= -1\n"
    "												     \n"
    "RECFAST_fudge = 1.14\n"
    "RECFAST_fudge_He = 0.86\n"
    "RECFAST_Heswitch = 6\n"
    "RECFAST_Hswitch  = T\n"
    "												     \n"
    "initial_condition= 1\n"
    "initial_vector = -1 0 0 0 0\n"
    "												     \n"
    "vector_mode = 0\n"
    "												     \n"
    "COBE_normalize = F\n"
    "CMB_outputscale = 7.42835025e12 \n"
    "												     \n"
    "transfer_power_var = 7\n"
    "transfer_high_precision = T\n"
    "transfer_kmax = 50\n"
    "transfer_k_per_logint= %i\n"
    "												     \n"
    "transfer_num_redshifts = %i\n"
    "												     \n"
    "transfer_interp_matterpower = F\n"
    "												     \n"
    "scalar_output_file = scalCls.dat\n"
    "vector_output_file = vecCls.dat\n"
    "tensor_output_file = tensCls.dat\n"
    "total_output_file  = totCls.dat\n"
    "lensed_output_file = lensedCls.dat\n"
    "lensed_total_output_file  =lensedtotCls.dat\n"
    "lens_potential_output_file = lenspotentialCls.dat\n"
    "FITS_filename      = scalCls.fits\n"
    "												     \n"
    "feedback_level = 1\n"
    "												     \n"
    "lensing_method = 1\n"
    "accurate_BB = F\n"
    "												     \n"
    "massive_nu_approx = 1\n"
    "												     \n"
    "accurate_polarization   = T\n"
    "accurate_reionization   = T\n"
    "do_tensor_neutrinos	 = F\n"
    "												     \n"
    "#Whether to turn off small-scale late time radiation hierarchies (save time,v. accurate)\n"
    "do_late_rad_truncation = F\n"
    "high_accuracy_default = F\n"
    "use_spline_template = F\n"
    "												     \n"
    "number_of_threads= 0\n"
    "accuracy_boost= 2\n"
    "l_accuracy_boost= 2\n"
    "l_sample_boost= 2\n",
    h*100.0,
    w0,
    wa,
    OB0*h*h,
    OC0*h*h,
    1.0-OB0-OC0-ONE2(1.0),
    ONE2(1.)*h*h,
    Nur,
    (int)round(N_nu),
    (int)round(N_eigenstates),
    (int)round(N_degeneracies),
    As,
    ns,
    k_per_logint_camb,
    znum);

    if (mode==0) {
    for (i=0; i< output_number; i++)
    {
    	fprintf(out,
    	"transfer_redshift(%i) = %lf \n"
    	"transfer_filename(%i) = z%i_tk.dat \n"
    	"transfer_matterpower(%i) = z%i_pk.dat \n",
    	i+1,
    	z_output[output_number-1-i],
    	i+1,
    	output_number-i,
    	i+1,
    	output_number-i);
    } }

    if (mode==1) {
    for (i=0; i<50; i++)
    {
    	fprintf(out,
    	"transfer_redshift(%i) = %lf \n"
    	"transfer_filename(%i) = z%i_tk.dat \n"
    	"transfer_matterpower(%i) = z%i_pk.dat \n",
    	i+1,
    	bc_zz[49-i],
    	i+1,
    	50-i,
    	i+1,
    	50-i);
    } }

    if (mode==2) {
    fprintf(out,
    "transfer_redshift(%i) = %lf \n"
    "transfer_filename(%i) = tk.dat \n"
    "transfer_matterpower(%i) = pk.dat \n",
    1,z_initial,1,1);
    }

    if (mode==3) {
    fprintf(out,
    "transfer_redshift(%i) = %lf \n"
    "transfer_filename(%i) = z%i_tk.dat \n"
    "transfer_matterpower(%i) = z%i_pk.dat \n"
    "transfer_redshift(%i) = %lf \n"
    "transfer_filename(%i) = z%i_tk.dat \n"
    "transfer_matterpower(%i) = z%i_pk.dat \n",
    	1,z_initial,1,2,1,2,
      2,z_final,2,1,2,1);
    }
    fclose(out);
  }
  else
  {
    char fileopen[200];

    if (mode==0)
    {
      printf("Creating %s input_file\n   %s/PK_TABS/power.ini\n",boltzmann_code,dir_chain);
      sprintf(fileopen,"%s/PK_TABS/power.ini",dir_chain);
    }
    if (mode==1 || mode==2)
    {
      printf("Creating %s input_file\n   %s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini\n",boltzmann_code,dir_chain);
      sprintf(fileopen,"%s/BOUNDARY_CONDITIONS_MODULE/tabs/power.ini",dir_chain);
    }
    if (mode==3)
    {
      sprintf(fileopen,"%s/PK_TABS/power_norm.ini",dir_chain);
    }

    double Nur;
    if (M_nu==0.) Nur = 3.046;
    else Nur = Neff;

    int i;
    char command_mkdir[200];

    FILE *out = fopen(fileopen,"w");
    if (out == NULL)
    {
      if (mode==0 || mode==3)
      {
        printf("Does the folder PK_TABS exist? Creating...\n");
        sprintf(command_mkdir,"mkdir %s/PK_TABS",dir_chain);
        system(command_mkdir);
        out = fopen(fileopen,"w");
        if (out==NULL)
        {
          printf("Error creating params file...\n");
          exit(1);
        }
      }
      else
      {
        printf("Does the folder BOUNDARY_CONDITIONS_MODULE/tabs exist? Creating...\n");
        sprintf(command_mkdir,"mkdir %s/BOUNDARY_CONDITIONS_MODULE",dir_chain);
        system(command_mkdir);
        sprintf(command_mkdir,"mkdir %s/BOUNDARY_CONDITIONS_MODULE/tabs",dir_chain);
        system(command_mkdir);
        out = fopen(fileopen,"w");
        if (out==NULL)
        {
          printf("Error creating params file...\n");
          exit(1);
        }
      }
    }

    if (strcmp(class_base_par_file,"none")==0)
    {
      fprintf(out,
              "H0 = %.8lf                      \n"
              "T_cmb = 2.7255                      \n"
              "Omega_b = %.8lf                     \n"
              "N_ur = %.8lf                      \n"
              "Omega_cdm = %.8lf                     \n"
              "Omega_dcdmdr = 0.0                    \n"
              "Gamma_dcdm = 0.0                    \n"
              "N_ncdm = %lf                      \n",
              h * 100.,
              OB0,
              Nur,
              OC0,
              N_nu);

      fprintf(out,
      "m_ncdm = %lf",M_nu/N_nu);
      for(i=1;i<N_nu;i++) fprintf(out,", %lf",M_nu/N_nu);
      fprintf(out,"\n");

      fprintf(out,
      "Omega_ncdm =                        \n"
      "T_ncdm = ");
      if(M_nu!=0)
      {
        fprintf(out,"0.71611");
        for(i=1;i<N_nu;i++) fprintf(out,", 0.71611");
      }
      fprintf(out,"\n");
      if(M_nu!=0)
      {
        printf("deg_ncdm = 1");
        for(i=1;i<N_nu;i++) fprintf(out,", 1");
      }
      fprintf(out,"\n");

      fprintf(out,
      "Omega_k = 0.                      \n"
      "#Omega_Lambda = 0.7                     \n"
      "                        \n"
      "attractor_ic_scf = yes                    \n"
      "#scf_parameters = [scf_lambda, scf_alpha, scf_A, scf_B, phi, phi_prime]       \n"
      "scf_parameters = 10.0, 0.0, 0.0, 0.0, 100.0, 0.0            \n"
      "scf_tuning_index = 0                    \n"
      "                        \n"
      "YHe = BBN                       \n"
      "recombination = RECFAST                   \n"
      "reio_parametrization = reio_camb                \n"
      "tau_reio = %.8lf                    \n"
      "reionization_exponent = 1.5                   \n"
      "reionization_width = 0.5                  \n"
      "helium_fullreio_redshift = 3.5                  \n"
      "helium_fullreio_width = 0.5                   \n"
      "binned_reio_num = 3                     \n"
      "binned_reio_z = 8,12,16                   \n"
      "binned_reio_xe = 0.8,0.2,0.1                  \n"
      "binned_reio_step_sharpness = 0.3                \n"
      "annihilation = 0.                     \n"
      "annihilation_variation = 0.                   \n"
      "annihilation_z = 1000                     \n"
      "annihilation_zmax = 2500                  \n"
      "annihilation_zmin = 30                    \n"
      "annihilation_f_halo= 20                   \n"
      "annihilation_z_halo= 8                    \n"
      "on the spot = yes                     \n"
      "decay = 0.                      \n"
      "                        \n"
      "output = mPk,mTk                    \n"
      "                        \n"
      "non linear =                        \n"
      "                        \n"
      "modes = s                       \n"
      "                        \n"
      "lensing = no                      \n"
      "                        \n"
      "ic = ad                       \n"
      "gauge = synchronous                     \n"
      "                        \n"
      "P_k_ini type = analytic_Pk                  \n"
      "k_pivot = 0.05                      \n"
      "A_s = %.8e                          \n"
      "n_s = %.8lf                       \n"
      "alpha_s = 0.                      \n"
      "                        \n"
      "f_bi = 1.                       \n"
      "n_bi = 1.5                      \n"
      "f_cdi=1.                      \n"
      "f_nid=1.                      \n"
      "n_nid=2.                      \n"
      "alpha_nid= 0.01                     \n"
      "                        \n"
      "c_ad_bi = 0.5                       \n"
      "c_ad_cdi = -1.                      \n"
      "c_bi_nid = 1.                       \n"
      "                        \n"
      "r = 1.                        \n"
      "n_t = scc                       \n"
      "alpha_t = scc                       \n"
      "                        \n"
      "potential = polynomial                    \n"
      "V_0=1.e-13                      \n"
      "V_1=-1.e-14                       \n"
      "V_2=7.e-14                      \n"
      "V_3=                        \n"
      "V_4=                        \n"
      "                        \n"
      "H_0=1.e-13                      \n"
      "H_1=-1.e-14                       \n"
      "H_2=7.e-14                      \n"
      "H_3=                        \n"
      "H_4=                        \n"
      "                        \n"
      "phi_end =                       \n"
      "Vparam0 =                       \n"
      "Vparam1 =                       \n"
      "Vparam2 =                       \n"
      "Vparam3 =                       \n"
      "Vparam4 =                       \n"
      "ln_aH_ratio = 50                    \n"
      "                        \n"
      "k1=0.002                      \n"
      "k2=0.1                        \n"
      "                        \n"
      "P_{RR}^1 = 2.3e-9                     \n"
      "P_{RR}^2 = 2.3e-9                     \n"
      "P_{II}^1 = 1.e-11                     \n"
      "P_{II}^2 = 1.e-11                     \n"
      "P_{RI}^1 = -1.e-13                    \n"
      "|P_{RI}^2| = 1.e-13                     \n"
      "                        \n"
      "special_iso =                       \n"
      "                        \n"
      "command = cat external_Pk/Pk_example.dat              \n"
      "                        \n"
      "custom1 = 0.05     # In the example command: k_pivot            \n"
      "custom2 = 2.215e-9 # In the example command: A_s            \n"
      "custom3 = 0.9624   # In the example command: n_s            \n"
      "custom4 = 2e-10    # In the example (with tensors) command: A_t         \n"
      "custom5 = -0.1     # In the example (with tensors) command: n_t         \n"
      "#custom6 = 0                      \n"
      "#custom7 = 0                      \n"
      "#custom8 = 0                      \n"
      "#custom9 = 0                      \n"
      "#custom10 = 0                       \n"
      "                        \n"
      "P_k_max_h/Mpc = %.1lf                     \n"
      "z_pk = ",
      //H0,OB,Nur,OC,MNU/3.,MNU/3.,MNU/3.,tau,AS,NS,Z
      tau_reio,
      As,
      ns,
      kmax);

      if (mode==0)
      {
        for(i=0; i<(output_number-1); i++) fprintf(out,"%.4lf,",z_output[i]);
        i = output_number-1; fprintf(out,"%.4lf\n",z_output[i]);
      }
      if (mode==1)
      {
        for(i=0;i<49;i++) fprintf(out,"%.4lf,",bc_zz[i]);
        fprintf(out,"%.4lf\n",bc_zz[49]);
      }
      if (mode==2)
      {
        fprintf(out,"%.4lf\n",z_initial);
      }
      if (mode==3)
      {
        fprintf(out,"%.4lf, %.4lf\n",0.,z_initial);
      }

      fprintf(out,
      "                        \n"
      "selection=gaussian                    \n"
      "selection_mean = 0.98,0.99,1.0,1.1,1.2                \n"
      "selection_width = 0.1                     \n"
      "non_diagonal=4                      \n"
      "                        \n"
      "dNdz_selection =                    \n"
      "dNdz_evolution =                    \n"
      "bias = 1.                       \n"
      "                        \n");

      if (mode==0) fprintf(out,"root = %s/PK_TABS/power_\n",dir_chain);
      if (mode==1) fprintf(out,"root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power_\n",dir_chain);
      if (mode==2) fprintf(out,"root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power_zin_\n",dir_chain);
      if (mode==3) fprintf(out,"root = %s/PK_TABS/power_norm_\n",dir_chain);

      fprintf(out,
      "headers = no                      \n"
      "format = class                      \n"
      "                        \n"
      "write background = no                     \n"
      "write thermodynamics = no                   \n"
      "write primordial = no                     \n"
      "write parameters = no                           \n"
      "write warnings =                    \n"
      "                        \n"
      "input_verbose = 1                     \n"
      "background_verbose = 1                    \n"
      "thermodynamics_verbose = 1                  \n"
      "perturbations_verbose = 1                   \n"
      "transfer_verbose = 1                    \n"
      "primordial_verbose = 1                    \n"
      "spectra_verbose = 1                     \n"
      "nonlinear_verbose = 1                     \n"
      "lensing_verbose = 1                     \n"
      "output_verbose = 1                    \n");
      fclose(out);
    }
    else
    {
      fprintf(out, "headers = no                      \n"
                   "format = class                      \n"
                   "output = mPk,mTk                    \n"
                   "P_k_max_h/Mpc = %.1lf                     \n"
                   "z_pk = ",
              kmax);
      if (mode==0)
      {
        for(i=0; i<(output_number-1); i++) fprintf(out,"%.4lf,",z_output[i]);
        i = output_number-1; fprintf(out,"%.4lf\n",z_output[i]);
      }
      if (mode==1)
      {
        for(i=0;i<49;i++) fprintf(out,"%.4lf,",bc_zz[i]);
        fprintf(out,"%.4lf\n",bc_zz[49]);
      }
      if (mode==2)
      {
        fprintf(out,"%.4lf\n",z_initial);
      }
      if (mode==3)
      {
        fprintf(out,"%.4lf, %.4lf\n",0.,z_initial);
      }
      if (mode==0) fprintf(out,"root = %s/PK_TABS/power_\n",dir_chain);
      if (mode==1) fprintf(out,"root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power_\n",dir_chain);
      if (mode==2) fprintf(out,"root = %s/BOUNDARY_CONDITIONS_MODULE/tabs/power_zin_\n",dir_chain);
      if (mode==3) fprintf(out,"root = %s/PK_TABS/power_norm_\n",dir_chain);
      fclose(out);
      char buf[1024];
      sprintf(buf, "cat %s >> %s", class_base_par_file, fileopen);
      system(buf);
    }
  }
}
