/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/* ********************************************************************* */
{
  double r = x1;
  double z = x2;
  double dist = sqrt(r*r + z*z);

  /* Apply absorbing sphere logic at t=0 to prevent initial NaNs */
  if (dist < g_inputParam[R_SINK]) {
    v[RHO] = 1.e-6;
    v[PRS] = 1.e-6;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
  } else {
    /* These names refer to the indices 0, 1, 2, 3 from definitions.h */
    double v_inf = g_inputParam[V_INF];
    double mach  = g_inputParam[MACH_INF];
    double gamma = g_inputParam[GAMMA];
    
    /* Calculate p_inf based on Mach number: cs = v/M -> P = rho*v^2 / (gamma * M^2) */
    double p_inf = (v_inf * v_inf) / (gamma * mach * mach);

    v[RHO] = 1.0;
    v[VX1] = 0.0;
    v[VX2] = v_inf;
    v[VX3] = 0.0;
    v[PRS] = p_inf;
  }
}

/* MANDATORY EMPTY FUNCTIONS */
void InitDomain (Data *d, Grid *grid) { }
void Analysis (const Data *d, Grid *grid) { }

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* ********************************************************************* */
{
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];

  /* ================================================================
     side == 0: Internal boundary (Spherical Absorbing Sink at origin)
     ================================================================ */
  if (side == 0) {    
    BOX_LOOP(box, k, j, i) {
      double r = x1[i];
      double z = x2[j];
      double dist = sqrt(r*r + z*z);
      
      if (dist < g_inputParam[R_SINK]) {
        d->Vc[RHO][k][j][i] = 1.e-6; 
        d->Vc[PRS][k][j][i] = 1.e-6; 
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
      }
    }
  }

  /* ================================================================
     X1_BEG: Custom axis boundary at r = 0
     - Near the BH (|z| < R_SINK): outflow-style zero-gradient copy,
       so it behaves like the absorbing sphere's extension.
     - Far from BH: proper axisymmetric mirror (symmetry axis).
     ================================================================ */
  if (side == X1_BEG) {
    if (box->vpos == CENTER) {
      BOX_LOOP(box, k, j, i) {
        double z = x2[j];

        if (fabs(z) < g_inputParam[R_SINK]) {
          /* Zero-gradient outflow copy from first interior cell */
          d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
          d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
          d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IBEG];
          d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IBEG];
          d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][IBEG];
        } else {
          /* Axisymmetric mirror: scalars copied, radial velocity flipped */
          int i_mirror = 2*IBEG - i - 1;
          d->Vc[RHO][k][j][i] =  d->Vc[RHO][k][j][i_mirror];
          d->Vc[PRS][k][j][i] =  d->Vc[PRS][k][j][i_mirror];
          d->Vc[VX1][k][j][i] = -d->Vc[VX1][k][j][i_mirror]; /* flip r-velocity */
          d->Vc[VX2][k][j][i] =  d->Vc[VX2][k][j][i_mirror]; /* z-velocity unchanged */
          d->Vc[VX3][k][j][i] =  d->Vc[VX3][k][j][i_mirror];
        }
      }
    }
  }

  /* ================================================================
     X2_BEG: Upstream wind injector
     ================================================================ */
  if (side == X2_BEG) {  
    if (box->vpos == CENTER) {
      double v_inf = g_inputParam[V_INF];
      double mach  = g_inputParam[MACH_INF];
      double gamma = g_inputParam[GAMMA];
      double p_inf = (v_inf * v_inf) / (gamma * mach * mach);
      BOX_LOOP(box, k, j, i) {
        d->Vc[RHO][k][j][i] = 1.0;
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = v_inf;
        d->Vc[PRS][k][j][i] = p_inf;
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/* ********************************************************************* */
{
  double dist = sqrt(x1*x1 + x2*x2);
  /* Paczyński-Wiita Potential with safety cap at 1.01rg */
  if (dist < 1.01) dist = 1.01; 
  return -0.5 / (dist - 1.0); 
}

void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  g[IDIR] = g[JDIR] = g[KDIR] = 0.0;
}
#endif
