/******************************************************************************/
/*    INITIAL SETTINGS - GLOBAL VARIABLES                                     */
/******************************************************************************/
char input_file[200];
char outputfile[200];
char boltzmann_code[200];
char boltzmann_folder[500];

int wrong_nu = 0;
char do_rescaled_ps;
char print_hubble;

double z_initial,z_final;

double *z_output;

int output_number;
int bin_beta_in,bin_fc_in,bin_fn_in;
int bin_out;
int mode = 0;

double **matrix_M;
double **matrix_M2;

double H0 = 100.;
double OM0,OB0,OC0,OX0,OG0,OR0,h,M_nu,tau_reio,ns,As,kmax,N_nu,Neff;
double w0,wa;

double Tcmb_0 = 2.7255;      //K
double Kb = 8.617342e-5;     //eV K^-1
double Gamma_nu = 0.71611;   //Tnu0/Tcmb0

double *ytab;
double *FFtab;
double *GGtab;
double ytab_max;
double ytab_min;
double F_inf,F_0;
double G_inf,G_0;
double ytabstep;
double ps_norm_at_z;
int ntab;

//Hard coded boundary condition: the amplitude of the total matter
//growth factor at redshift z0 (which is usually z0=0) is fixed to
//be 1.
double alpha = 1.;
