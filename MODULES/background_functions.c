#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>

#include "include_extern.h"
extern double FF(double y);
extern double GG(double y);

/******************************************************************************/
/*    BACKGROUND                                                              */
/******************************************************************************/
double A_func (double A, double OR_rk, double OCB_rk, double OX_rk,
               double rk_1_plus_3w, double ON_rk, double E2_rk)
{
  if (wrong_nu==1)
  {
    // In case05, being the background exact, we need to use
    // the correct value of w_nu and Omgega_nu in the computation
    // of A(z). I will then overwrite the approximated values with
    // the correct ones just in this function.
    double Gnu4,pi4,a4,y;
    if (M_nu==0.) ON_rk= 0.0;
    else
    {
      Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
      pi4 = M_PI*M_PI*M_PI*M_PI;
      a4 = A*A*A*A;
      y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
      ON_rk = ((15.*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
      ON_rk /= E2_rk;
    }
    if (M_nu==0.) rk_1_plus_3w = 1.;
    else rk_1_plus_3w = 2.0 - y*(GG(y)/FF(y));
  }
  return 0.5*(OCB_rk+2.*OR_rk+(1.+3.*(w0+wa*(1.-A)))*OX_rk+rk_1_plus_3w*ON_rk-2.);
}

double E2(double A, double ONE2_CURRENT)
{
  if (wrong_nu==1)
  {
    // In case05 we need the exact H(z), so I'm substituting the
    // approximated neutrino fraction with the correct one
    // only here in E2(z)
    if (M_nu==0.) ONE2_CURRENT= 0.0;
    else
    {
      double Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
      double pi4 = M_PI*M_PI*M_PI*M_PI;
      double a4 = A*A*A*A;
      double y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
      ONE2_CURRENT = ((15.*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
    }
  }
  return (OR0/(A*A*A*A) + (OC0+OB0)/(A*A*A) +
          OX0*(pow(A,-3.0*(1.0+w0+wa))*exp(3.0*wa*(A-1.0))) + ONE2_CURRENT);
}

double func_1_3w(double A)
{
  if (wrong_nu != 0)
    return 1.;
  else
  {
    if (M_nu==0.) return 1.;
    double y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
    return y*(GG(y)/FF(y));
  }
}

double OCB(double A, double E2)
{
  return ((OC0+OB0)/(A*A*A))/E2;
}

double set_ON0()
{
  if (M_nu==0.) return 0.0;
  else
  {
    if (wrong_nu==2)
      return (M_nu/(93.14*h*h));
    else
    {
      double A,Gnu4,pi4,a4,y;
      A = 1.0;
      Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
      pi4 = M_PI*M_PI*M_PI*M_PI;
      a4 = A*A*A*A;
      y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
      return ((15.*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
    }
  }
}

double ONE2(double A)
{
  if (wrong_nu==0)
  {
    if (M_nu==0.) return 0.;
    double Gnu4 = Gamma_nu*Gamma_nu*Gamma_nu*Gamma_nu;
    double pi4 = M_PI*M_PI*M_PI*M_PI;
    double a4 = A*A*A*A;
    double y = (M_nu*A)/(Gamma_nu*Kb*N_nu*Tcmb_0);
    return ((15.*Gnu4*N_nu*(2.469e-05/(h*h)))/(a4*pi4))*FF(y);
  }
  else
  return (M_nu/(93.14*h*h))/(A*A*A);
}

double ON(double ON_e2, double E2)
{
  return ON_e2/E2;
}

double OR(double A, double E2)
{
    return (OR0/(A*A*A*A))/E2;
}

double OX(double A, double E2)
{
  return (OX0*(pow(A,-3.0*(1.0+w0+wa))*exp(3.0*wa*(A-1.0))))/E2;
}

void print_hubble_table()
{
  char openf[200];
  sprintf(openf,"%s_hubble.txt",outputfile);

  printf("\nHubble table written to %s\n",openf);

  FILE *hub = fopen(openf,"w");

  int z_nstep=1000;
  double z_step=(double)((log(1.0+z_initial)-log(1.0+z_final))/(z_nstep-1));

  double zz,e2;

  int i;
  for (i=0; i<z_nstep; i++)
  {
    zz = exp(log(1.+z_initial)+i*z_step);
    e2 = E2(1.0/zz,ONE2(1.0/zz));

    fprintf(hub,"%.10e\t%.10e\n",zz-1.,h*H0*sqrt(e2));
  }
}
