/**   RUNGE-KUTTA SOLVER ******************************************************/
extern void RK (int k_num, double *k,
      double *BETAB, double *BETANU,
      double *FB, double *FC, double *FN,
      double **Delta_b, double **Delta_c, double **Delta_n, double **Delta_m,
      double **growth_b, double **growth_c, double **growth_n, double **growth_m);
/******************************************************************************/
