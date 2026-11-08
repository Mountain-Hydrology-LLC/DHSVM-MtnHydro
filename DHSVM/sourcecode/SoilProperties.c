/*
 * DESCRIPTION:  Soil thermal properties for frozen soil, as in VIC 5
 *               (soil_conduction.c): unfrozen water from freezing point
 *               depression, Johansen thermal conductivity with an organic
 *               fraction (Farouki 1981) and volumetric heat capacity, and the
 *               moisture and properties of the soil temperature nodes.
 */
#include <math.h>
#include "settings.h"
#include "constants.h"
#include "soiltemp.h"

/*****************************************************************************
  UnfrozenWater()

  Largest liquid water content (m3/m3) at temperature T (C) of a soil with
  Porosity, Brooks-Corey bubbling pressure Press (m) and pore size
  distribution index PoreDist
*****************************************************************************/
double UnfrozenWater(double T, double Porosity, double Press, double PoreDist)
{
  double Unfrozen;

  if (T >= 0.)
    return Porosity;

  Unfrozen = Porosity * pow((-LF * T) / SOIL_TKTRIP / (G * Press), -PoreDist);
  if (Unfrozen > Porosity)
    Unfrozen = Porosity;
  if (Unfrozen < 0.)
    Unfrozen = 0.;
  return Unfrozen;
}

/*****************************************************************************
  NodeIce()

  Ice content (m3/m3, as water) in equilibrium at T (C) with total water
  Moist (m3/m3)
*****************************************************************************/
double NodeIce(double T, double Moist, double Porosity, double Press,
               double PoreDist)
{
  double Ice;

  if (T >= 0.)
    return 0.;
  Ice = Moist - UnfrozenWater(T, Porosity, Press, PoreDist);
  return (Ice > 0.) ? Ice : 0.;
}

/*****************************************************************************
  SoilConductivity()

  Thermal conductivity (W/(m C)), Johansen's method (Farouki 1986). Moist:
  total water, Liquid: liquid water (m3/m3), Quartz: quartz fraction of the
  mineral solids, Organic: organic fraction of the solids by volume. The
  mineral bulk density follows the porosity.
*****************************************************************************/
double SoilConductivity(double Moist, double Liquid, double Porosity,
                        double Quartz, double Organic)
{
  double BulkDensity = (1. - Porosity) * SOIL_PARTICLE_DENSITY;
  double Kdry, Ks, Ksat, Ke, Sr, K;

  Kdry = (0.135 * BulkDensity + 64.7) / (SOIL_PARTICLE_DENSITY - 0.947 * BulkDensity);
  Kdry = (1. - Organic) * Kdry + Organic * SOIL_KDRY_ORG;

  if (Moist <= 0.)
    return Kdry;

  /* The geometric means as exponentials of weighted logarithms (the
     logarithms of constants fold at compile time) */
  Sr = Moist / Porosity;
  Ks = exp(Quartz * log(7.7) + (1. - Quartz) * ((Quartz < 0.2) ? log(3.0) : log(2.2)));
  Ks = (1. - Organic) * Ks + Organic * SOIL_KS_ORG;

  if (Liquid == Moist) {
    Ksat = exp((1. - Porosity) * log(Ks) + Porosity * log(SOIL_K_WATER));
    Ke = 0.7 * log10(Sr) + 1.;
  }
  else {
    Ksat = exp((1. - Porosity) * log(Ks) + (Porosity - Liquid) * log(SOIL_K_ICE) +
               Liquid * log(SOIL_K_WATER));
    Ke = Sr;
  }

  K = (Ksat - Kdry) * Ke + Kdry;
  return (K > Kdry) ? K : Kdry;
}

/*****************************************************************************
  SoilHeatCapacity()

  Volumetric heat capacity (J/(m3 C)) of a soil with Porosity, Liquid water
  and Ice (m3/m3) and an Organic fraction of the solids
*****************************************************************************/
double SoilHeatCapacity(double Porosity, double Liquid, double Ice,
                        double Organic)
{
  double Solids = 1. - Porosity;

  return SOIL_CS_MIN * Solids * (1. - Organic) + SOIL_CS_ORG * Solids * Organic +
    SOIL_CS_WATER * Liquid + SOIL_CS_ICE * Ice +
    SOIL_CS_AIR * (1. - (Solids + Liquid + Ice));
}

/*****************************************************************************
  NodeLayer()

  The layer holding node k (0..NLayers-1 the root zone layers, NLayers the
  deep layer down to SoilDepth, NLayers + 1 below the soil); *n = 2 if the
  node lies on the boundary of that layer and the one below, else 1
*****************************************************************************/
static int NodeLayer(int k, int NLayers, float *LayerThickness, float SoilDepth,
                     int *n)
{
  double Bottom = 0.;
  int l;

  for (l = 0; l <= NLayers; l++) {
    Bottom += (l < NLayers) ? LayerThickness[l] : 0.;
    if (l == NLayers)
      Bottom = SoilDepth;
    if (SoilNodeDepth[k] <= Bottom)
      break;
  }
  *n = (k > 0 && l <= NLayers && SoilNodeDepth[k] == Bottom) ? 2 : 1;
  return l;
}

/*****************************************************************************
  NodeWater()

  Total water (liquid and ice, m3/m3) of each temperature node, as in
  NodeMoisture()
*****************************************************************************/
void NodeWater(int NLayers, float *LayerThickness, float SoilDepth,
               float *Moist, float *Porosity, double *Water)
{
  double Pore;
  int k, l, m, n, Side;

  for (k = 0; k < SOIL_NODES; k++) {
    l = NodeLayer(k, NLayers, LayerThickness, SoilDepth, &n);
    Water[k] = Pore = 0.;
    for (Side = 0; Side < n; Side++) {
      m = l + Side;
      Water[k] += ((m > NLayers) ? Porosity[NLayers] : Moist[m]) / n;
      Pore += Porosity[(m > NLayers) ? NLayers : m] / n;
    }
    if (Water[k] > Pore)
      Water[k] = Pore;
  }
}

/*****************************************************************************
  NodeMoisture()

  Fills Nodes with the water content and properties of each temperature
  node. Layers 0..NLayers-1 are the root zone layers (LayerThickness), layer
  NLayers the deep layer down to SoilDepth; below it the soil is saturated
  with the deep layer's properties. A node on a layer boundary takes the mean
  of the two layers. Moist, Porosity, Press, PoreDist, Quartz and Organic
  have NLayers + 1 values. IceNode (m3/m3) is capped at the node's water;
  conductivity and heat capacity use that ice.
*****************************************************************************/
void NodeMoisture(int NLayers, float *LayerThickness, float SoilDepth,
                  float *Moist, float *Porosity, float *Press, float *PoreDist,
                  float *Quartz, float *Organic, float *IceNode,
                  SOILNODES *Nodes)
{
  int k, l, m, n, Side;
  double Liquid;

  for (k = 0; k < SOIL_NODES; k++) {
    l = NodeLayer(k, NLayers, LayerThickness, SoilDepth, &n);

    /* On a boundary, average layer l and the one below */
    Nodes->Moist[k] = Nodes->Porosity[k] = Nodes->Press[k] = 0.;
    Nodes->PoreDist[k] = Nodes->Quartz[k] = Nodes->Organic[k] = 0.;
    for (Side = 0; Side < n; Side++) {
      m = l + Side;
      if (m > NLayers) {
        Nodes->Moist[k] += Porosity[NLayers] / n;
        m = NLayers;
      }
      else
        Nodes->Moist[k] += Moist[m] / n;
      Nodes->Porosity[k] += Porosity[m] / n;
      Nodes->Press[k] += Press[m] / n;
      Nodes->PoreDist[k] += PoreDist[m] / n;
      Nodes->Quartz[k] += Quartz[m] / n;
      Nodes->Organic[k] += Organic[m] / n;
    }
    if (Nodes->Moist[k] > Nodes->Porosity[k])
      Nodes->Moist[k] = Nodes->Porosity[k];

    if (IceNode[k] > Nodes->Moist[k])
      IceNode[k] = Nodes->Moist[k];
    Liquid = Nodes->Moist[k] - IceNode[k];
    if (k > 0 && IceNode[k] == IceNode[k - 1] && Nodes->Moist[k] == Nodes->Moist[k - 1] &&
        Nodes->Porosity[k] == Nodes->Porosity[k - 1] && Nodes->Quartz[k] == Nodes->Quartz[k - 1] &&
        Nodes->Organic[k] == Nodes->Organic[k - 1]) {
      /* As the node above (the saturated nodes below the soil) */
      Nodes->Kappa[k] = Nodes->Kappa[k - 1];
      Nodes->Cs[k] = Nodes->Cs[k - 1];
      continue;
    }
    Nodes->Kappa[k] = SoilConductivity(Nodes->Moist[k], Liquid,
      Nodes->Porosity[k], Nodes->Quartz[k], Nodes->Organic[k]);
    Nodes->Cs[k] = SoilHeatCapacity(Nodes->Porosity[k], Liquid, IceNode[k],
      Nodes->Organic[k]);
  }
}
