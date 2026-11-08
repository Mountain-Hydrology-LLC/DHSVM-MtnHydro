/*
 * DESCRIPTION:  Soil temperature profile with freezing and thawing in one
 *               pass per time step, as in CLM 5 (Lawrence et al. 2018,
 *               "Soil and Snow Temperatures", "Phase Change"; Oleson et al.
 *               2013 ch. 6): a Crank-Nicolson finite volume heat equation
 *               on exponentially spaced nodes (VIC 5, Cherkauer and
 *               Lettenmaier 1999) with the start-of-step conductivity and
 *               heat capacity and no phase change, the surface node taking
 *               a flux linear in its temperature and the bottom node at a
 *               fixed temperature; then each node freezes or thaws with the
 *               heat of its departure from 0 C, limited by its liquid water
 *               above the unfrozen water content or by its ice.
 */
#include <math.h>
#include <stdio.h>
#include "settings.h"
#include "constants.h"
#include "soiltemp.h"

#define PROFILE_ALPHA 0.5       /* Weight of the end of the step in the conduction (Crank-Nicolson) */

double SoilNodeDepth[SOIL_NODES];
static double SoilNodeDz[SOIL_NODES];    /* Thickness of each node's volume (m) */

/*****************************************************************************
  InitSoilNodes()

  Node depths z_k = exp(B k) - 1, from the surface to SOIL_DAMPING_DEPTH,
  B = ln(SOIL_DAMPING_DEPTH + 1) / (SOIL_NODES - 1), and the volumes of the
  nodes above the bottom one: from the surface or the midpoint above to the
  midpoint below
*****************************************************************************/
void InitSoilNodes(void)
{
  double B;
  int k;

  B = log(SOIL_DAMPING_DEPTH + 1.) / (double) (SOIL_NODES - 1);
  for (k = 0; k < SOIL_NODES; k++)
    SoilNodeDepth[k] = exp(B * k) - 1.;
  SoilNodeDepth[SOIL_NODES - 1] = SOIL_DAMPING_DEPTH;
  SoilNodeDz[0] = SoilNodeDepth[1] / 2.;
  for (k = 1; k < SOIL_NODES - 1; k++)
    SoilNodeDz[k] = (SoilNodeDepth[k + 1] - SoilNodeDepth[k - 1]) / 2.;
  SoilNodeDz[SOIL_NODES - 1] = 0.;
}

/*****************************************************************************
  Tridiagonal()

  Solves a[i] x[i-1] + b[i] x[i] + c[i] x[i+1] = r[i], i = 0..n-1, into r
  (Thomas algorithm); c is overwritten
*****************************************************************************/
static void Tridiagonal(double *a, double *b, double *c, double *r, int n)
{
  int i;
  double Factor;

  c[0] /= b[0];
  r[0] /= b[0];
  for (i = 1; i < n; i++) {
    Factor = b[i] - a[i] * c[i - 1];
    c[i] /= Factor;
    r[i] = (r[i] - a[i] * r[i - 1]) / Factor;
  }
  for (i = n - 2; i >= 0; i--)
    r[i] -= c[i] * r[i + 1];
}

/*****************************************************************************
  SolveProfile()

  Advances the node temperatures TNode (C) and ice IceNode (m3/m3) by Dt
  (s). The surface node takes the flux into the ground H0 + Beta (T0 - T0
  at the start) (W/m2, Beta <= 0); the bottom node stays at TBot. Nodes
  holds the node water and the start-of-step properties (NodeMoisture()).
  Conductance between two nodes: the harmonic mean of their conductivities
  over their distance. Returns the flux into the ground (W/m2, at the
  surface node's temperature of the conduction solve); *Storage receives
  the change of sensible and latent heat of the column and *Rest the
  energy balance error, flux in - flux out at the bottom - Storage (W/m2).
*****************************************************************************/
double SolveProfile(float *TNode, float *IceNode, float TBot, SOILNODES *Nodes,
                    double Dt, double H0, double Beta, double *Storage,
                    double *Rest)
{
  const int n = SOIL_NODES - 1;   /* Nodes solved: 0..n-1 */
  const double Alpha = PROFILE_ALPHA;
  double a[SOIL_NODES], b[SOIL_NODES], c[SOIL_NODES], r[SOIL_NODES];
  double Cap[SOIL_NODES], Cond[SOIL_NODES], Flux[SOIL_NODES];
  double T, Ice, Liquid, Excess, Change, Latent = 0., Sensible = 0.;
  double Into, Bottom;
  double *z = SoilNodeDepth;
  int k;

  /* Heat capacity of each node's volume per step (W/(m2 C)), conductance
     between node k and k + 1 (W/(m2 C)) and the flux down between them at
     the start of the step (W/m2) */
  for (k = 0; k < n; k++) {
    Cap[k] = Nodes->Cs[k] * SoilNodeDz[k] / Dt;
    Cond[k] = 2. * Nodes->Kappa[k] * Nodes->Kappa[k + 1] /
      (Nodes->Kappa[k] + Nodes->Kappa[k + 1]) / (z[k + 1] - z[k]);
    Flux[k] = Cond[k] * (TNode[k] - TNode[k + 1]);
  }

  /* Conduction: Crank-Nicolson; the surface flux implicit */
  for (k = 0; k < n; k++) {
    a[k] = (k > 0) ? -Alpha * Cond[k - 1] : 0.;
    c[k] = -Alpha * Cond[k];
    b[k] = Cap[k] + Alpha * Cond[k] + ((k > 0) ? Alpha * Cond[k - 1] : -Beta);
    r[k] = Cap[k] * TNode[k] - (1. - Alpha) * Flux[k] +
      ((k > 0) ? (1. - Alpha) * Flux[k - 1] : H0 - Beta * TNode[0]);
  }
  r[n - 1] += Alpha * Cond[n - 1] * TBot;
  Tridiagonal(a, b, c, r, n);
  Into = H0 + Beta * (r[0] - TNode[0]);
  Bottom = (1. - Alpha) * Flux[n - 1] + Alpha * Cond[n - 1] * (r[n - 1] - TBot);

  /* Phase change: a node below 0 C freezes its liquid water above the
     unfrozen water content, a node above 0 C thaws its ice, with at most
     the heat of the departure from 0 C; the rest stays in its temperature */
  for (k = 0; k < n; k++) {
    T = r[k];
    Ice = IceNode[k];
    Liquid = Nodes->Moist[k] - Ice;
    Change = 0.;
    if (T < 0.) {
      Excess = Liquid - UnfrozenWater(T, Nodes->Porosity[k], Nodes->Press[k],
                                      Nodes->PoreDist[k]);
      if (Excess > 0.) {
        Change = -Nodes->Cs[k] * T / (WATER_DENSITY * LF);
        if (Change > Excess) {
          Change = Excess;
          T += WATER_DENSITY * LF * Change / Nodes->Cs[k];
        }
        else
          T = 0.;
      }
    }
    else if (T > 0. && Ice > 0.) {
      Change = -Nodes->Cs[k] * T / (WATER_DENSITY * LF);
      if (Change < -Ice) {
        Change = -Ice;
        T += WATER_DENSITY * LF * Change / Nodes->Cs[k];
      }
      else
        T = 0.;
    }
    Sensible += Cap[k] * (T - TNode[k]);
    Latent -= WATER_DENSITY * LF * Change * SoilNodeDz[k] / Dt;
    TNode[k] = (float) T;
    IceNode[k] = (float) (Ice + Change);
  }
  TNode[n] = TBot;
  IceNode[n] = (float) NodeIce(TBot, Nodes->Moist[n], Nodes->Porosity[n],
                               Nodes->Press[n], Nodes->PoreDist[n]);

  *Storage = Sensible + Latent;
  *Rest = Into - Bottom - *Storage;
  return Into;
}

/*****************************************************************************
  FrostZones()

  The frozen zones (runs of frozen nodes: below 0 C, or at 0 C holding ice)
  whose tops lie above SoilDepth, shallowest first: their tops TopZ and
  bottoms BottomZ (m, bottoms at most SoilDepth); returns their number (at
  most SOIL_FROST_ZONES). Zone edges are the 0 C points between nodes; a
  zone starting in a node at 0 C starts at the top of that node's volume; a
  zone ending in a node at 0 C (partly frozen) ends within that node's
  volume, at the fraction of its water (Water, m3/m3) that is ice, counted
  from the top of the volume; a zone reaching the bottom node ends there.
*****************************************************************************/
int FrostZones(float *TNode, float *IceNode, double *Water, float SoilDepth,
               float *TopZ, float *BottomZ)
{
  double Top, Bottom, VolTop, VolBottom, Frac;
  int k, Start, n = 0;
  double *z = SoilNodeDepth;

#define FROZEN(k) (TNode[k] < 0. || (TNode[k] <= 0. && IceNode[k] > 0.))
  k = 0;
  while (k < SOIL_NODES) {
    if (!FROZEN(k)) {
      k++;
      continue;
    }
    Start = k;
    while (k < SOIL_NODES && FROZEN(k))
      k++;
    if (Start == 0)
      Top = 0.;
    else if (TNode[Start] < 0.)
      Top = z[Start - 1] + (0. - TNode[Start - 1]) *
        (z[Start] - z[Start - 1]) / (TNode[Start] - TNode[Start - 1]);
    else
      Top = (z[Start - 1] + z[Start]) / 2.;
    if (k == SOIL_NODES)
      Bottom = z[SOIL_NODES - 1];
    else if (TNode[k - 1] < 0.)
      Bottom = z[k - 1] + (0. - TNode[k - 1]) * (z[k] - z[k - 1]) /
        (TNode[k] - TNode[k - 1]);
    else {
      VolTop = (k - 1 == 0) ? 0. : (z[k - 2] + z[k - 1]) / 2.;
      VolBottom = (z[k - 1] + z[k]) / 2.;
      Frac = (Water[k - 1] > IceNode[k - 1]) ? IceNode[k - 1] / Water[k - 1] : 1.;
      Bottom = VolTop + Frac * (VolBottom - VolTop);
    }
    if (Top < SoilDepth) {
      TopZ[n] = (float) Top;
      BottomZ[n] = (float) ((Bottom < SoilDepth) ? Bottom : SoilDepth);
      n++;
    }
  }
#undef FROZEN
  return n;
}

/*****************************************************************************
  NodeMean()

  Mean over Top..Bottom (m) of node values f, linear between nodes and
  constant below the bottom node
*****************************************************************************/
static double NodeMean(double *f, double Top, double Bottom)
{
  double *z = SoilNodeDepth;
  double Sum = 0., a, b, fa, fb;
  int k;

  for (k = 0; k < SOIL_NODES - 1; k++) {
    a = (Top > z[k]) ? Top : z[k];
    b = (Bottom < z[k + 1]) ? Bottom : z[k + 1];
    if (b <= a)
      continue;
    fa = f[k] + (f[k + 1] - f[k]) * (a - z[k]) / (z[k + 1] - z[k]);
    fb = f[k] + (f[k + 1] - f[k]) * (b - z[k]) / (z[k + 1] - z[k]);
    Sum += (b - a) * (fa + fb) / 2.;
  }
  a = (Top > z[SOIL_NODES - 1]) ? Top : z[SOIL_NODES - 1];
  if (Bottom > a)
    Sum += (Bottom - a) * f[SOIL_NODES - 1];
  return Sum / (Bottom - Top);
}

/*****************************************************************************
  LayerAverages()

  Mean temperature TLayer (C) and ice Ice (m3/m3, at most Moist) of the
  root zone layers and the deep layer (NLayers + 1 values) from the node
  temperatures and ice
*****************************************************************************/
void LayerAverages(float *TNode, float *IceNode, int NLayers,
                   float *LayerThickness, float SoilDepth, float *Moist,
                   float *TLayer, float *Ice)
{
  double T[SOIL_NODES], I[SOIL_NODES];
  double Top = 0., Bottom, Mean;
  int k, l;

  for (k = 0; k < SOIL_NODES; k++) {
    T[k] = TNode[k];
    I[k] = IceNode[k];
  }
  for (l = 0; l <= NLayers; l++) {
    Bottom = (l < NLayers) ? Top + LayerThickness[l] : SoilDepth;
    TLayer[l] = NodeMean(T, Top, Bottom);
    Mean = NodeMean(I, Top, Bottom);
    Ice[l] = (Mean < Moist[l]) ? Mean : Moist[l];
    Top = Bottom;
  }
}
