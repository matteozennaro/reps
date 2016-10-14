/**   GENERAL PURPOSE *********************************************************/
extern void frame(char str[]);
extern void fscanf_error(int n);
extern double det(double **a_in,int n);
extern double **allocate_matrix(int rows, int columns);
extern double *allocate_double_vec(int n_elems);
extern void deallocate_matrix(double **m, int rows, int columns);
extern int find_z_bin(double Z, double *vec, int n_elems);
extern void reallocate_matrix(double **m,int rows,int newcols);
extern void sort_double_vec(int n_elems, double * vec);
extern int count_lines(char file[]);
extern int count_header_lines(char file[]);
extern int count_number_of_columns(char file[], int number_of_header_lines);
extern double lin_interp_between(double x, double x0, double x1, double y0, double y1);
/******************************************************************************/
