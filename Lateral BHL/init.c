/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/* * x1 = r, x2 = phi, x3 = z (ignored in 2D)
 *********************************************************************** */
{
  double v_inf = g_inputParam[V_INF];
  double mach  = g_inputParam[MACH_INF];
  double gamma = g_inputParam[GAMMA];
  
  /* Calculate sound speed from Mach and V_inf */
  double cs  = v_inf / mach;

  /* BACKGROUND AMBIENT MEDIUM 
   * Set density to 0.0108 as requested.
   * Calculate pressure based on this density to maintain Mach consistency.
   */
  v[RHO] = 0.0108; 
  v[VX1] = 0.0; 
  v[VX2] = 0.0; 
  v[VX3] = 0.0;
  v[PRS] = (0.0108 * cs * cs) / gamma; 
  v[TRC] = 0.0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid) {}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid) {}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  int  i, j, k;
  double *phi = grid->x[JDIR];

  double v_inf = g_inputParam[V_INF];
  double mach  = g_inputParam[MACH_INF];
  double gamma = g_inputParam[GAMMA];
  
  /* Derived pressure for the incoming gas (rho = 0.5) */
  double cs    = v_inf / mach;
  double prs_in = (0.5 * cs * cs) / gamma;

  /* Target the outer boundary */
  if (side == X1_END){  
    
    /* Check variable centering for cell-centered variables */
    if (box->vpos == CENTER) {
      
      BOX_LOOP(box,k,j,i){
        
        double phi_min = 265.0 * CONST_PI / 180.0;
        double phi_max = 270.0 * CONST_PI / 180.0;

        if (phi[j] >= phi_min && phi[j] <= phi_max) {
          /* INFLOW SLIT: Incoming gas density = 0.5 */
          d->Vc[RHO][k][j][i] = 0.5;
          d->Vc[VX1][k][j][i] = v_inf * cos(phi[j]);
          d->Vc[VX2][k][j][i] = -v_inf * sin(phi[j]);
          d->Vc[PRS][k][j][i] = prs_in;
        } else {
          /* OUTFLOW: Copy from last interior cell */
          d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][box->iend];
          d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][box->iend];
          d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][box->iend];
          d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][box->iend];
        }

        /* --- PRESSURE FLOOR --- 
         * Numerical safety to prevent NaN crashes
         */
        if (d->Vc[PRS][k][j][i] < 1.e-12) d->Vc[PRS][k][j][i] = 1.e-12;
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  g[IDIR] = 0.0; g[JDIR] = 0.0; g[KDIR] = 0.0;
}

/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
{
  double r = (x1 < 1.01) ? 1.01 : x1;
  return -0.5 / (r - 1.0); 
}
#endif
