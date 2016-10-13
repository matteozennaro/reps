/******************************************************************************/
/*    DECLARATION OF FUNCTIONS                                                */
/******************************************************************************/

/**   GENERAL PURPOSE *********************************************************/
extern void fscanf_error(int n);
extern double det(double **a_in,int n);
extern double **allocate_matrix(int rows, int columns);
extern double *allocate_double_vec(int n_elems);
extern void deallocate_matrix(double **m, int rows, int columns);
extern int find_z_bin(double Z, double *vec, int n_elems);
extern void reallocate_matrix(double **m,int rows,int newcols);
extern int count_lines(char file[]);
extern double lin_interp_between(double x, double x0, double x1, double y0, double y1);

/**   READ INI FILE ***********************************************************/
extern void read_parameter_file(char parfile[]);

/**   RUNGE-KUTTA SOLVER ******************************************************/
extern void RK (int k_num, double *k, double *BETAB, double *BETANU, double *FB, double *FC, double *FN);

/**   RESCALE PS **************************************************************/
extern void read_D(char filename[],int n, double *dc, double *dn, double *dm);
extern void rescale_camb_ps(int knum, double *k);
extern void rescale_class_ps(int knum, double *k);

/**   BOLTZMANN SOLVER ********************************************************/
extern void create_boltzmann_ini_file (char dir_chain[]);

/**  BACKGROUND FUNCTIONS *****************************************************/
extern double func_1_3w(double A);
extern double A_func (double A, double OR_rk, double OCB_rk, double OX_rk,
                      double rk_1_plus_3w, double ON_rk, double E2_rk);
extern double E2(double A, double ON_CURRENT);
extern double OCB(double A, double E2);
extern double set_ON0();
extern double ONE2(double A);
extern double ON(double ON_e2, double E2);
extern double OR(double A, double E2);
extern double OX(double A, double E2);
extern void print_hubble_table();

/**   NEUTRINO DISTRIBUTION FUNCTION ******************************************/
extern void read_GG_FF_tabs();
extern double FF(double Y);
extern double GG(double Y);

/******************************************************************************/
