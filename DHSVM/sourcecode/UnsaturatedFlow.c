
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"
#include "slopeaspect.h"

/*****************************************************************************
Function name: UnsaturatedFlow()

Purpose      : Calculate the unsaturated flow in the soil column, and
adjust the moisture in each soil layer

Sources:
Bras, R. A., Hydrology, an introduction to hydrologic science, Addisson
Wesley, Inc., Reading, etc., 1990.
Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed
hydrology-vegetation model for complex terrain, Water Resour. Res.,
30(6), 1665-1679, 1994.

This function is based on Wigmosta et al [1994], and assumes a unit
hydraulic gradient in the unsaturated zone.  This implies a
steady-state situation and uniform moisture distribution.  This is not
unreasonable for situations where the groundwater is fairly deep
[Bras, 1990], but may be unrealistic for the boreal forest system, or
other similar situations where the groundwater table is relatively close
to the surface. The current formulation does not allow for upward
movement of water in the unsaturated zone, by capillary and matrix
forces.  No unsaturated flux is assumed to occur if the water content
drops below the field capacity.

The unsaturated hydraulic conductivity is calculated using the
Brooks-Corey equation see for example Wigmosta et al. [1994], eq.41

The calculated amount of drainage is averaged with the amount calculated
for the previous timestep, see eq. 42 [Wigmosta et al., 1994].

The residual moisture content is assumed to be zero (easy to adapt if
needed).

If the amount of soil moisture in a layer after drainage still
exceeds the porosity for that layer, this additional amount of water is
added to the drainage to the next layer.

CHANGES:

Changes have been made to account for the potental loss of soil storage
in a grid cell due to a channel.  Correction coefficents are
calculated in AdjustStorage() and CutBankGeometry()

Eli Boardman added unsaturated lateral gravimetric interflow (2026)

With soil freezing the ice of each layer (SoilMap Ice) stays in place: the
conductivity follows the liquid water, Ks ((Moist - Ice)/Porosity)^Exponent,
and drainage stops when the liquid water reaches the field capacity. The
water above the porosity still spills to the layer below. The table depth
counts the ice as water filling the pores.
*****************************************************************************/
void UnsaturatedFlow(OPTIONSTRUCT *Options, int Dt, float DX, float DY, float Infiltration,
  int NSoilLayers,
  float TotalDepth, float Area, float *RootDepth, float *Ks, float KsAnisotropy,
  float *PoreDist, float *Porosity, float *FCap,
  float *Perc, float *PercArea, float *Adjust,
  int CutBankZone, float BankHeight, float *TableDepth,
  float *IExcess, float *Moist, int InfiltOption,
  VEGPIX **VegMap, NETSTRUCT **Network, SOILPIX **SoilMap,
  float cosTheta, float sinTheta, TOPOPIX **TopoMap,
  MAPSIZE *Map, int y, int x)
{
  float LayerBottomDepth, NextLayerDepth;
  float DeepLayerDepth;		/* depth of the layer below the deepest root layer */
  float Drainage;		    /* amount of water drained from each soil
                               layer during the current timestep */
  float InterFlow, InterFlowK; /* amount of water drained laterally */
  float Exponent;		    /* Brooks-Corey exponent */
  float Multiplier; /* Brooks-Corey multiplicative factor on Ks */
  float FieldCapacity;		/* amount of water in soil at field capacity (m) */
  float MaxSoilWater;		/* maximum allowable amount of soil moiture in each layer (m) */
  float SoilWater;		    /* amount of water in each soil layer (m) */
  float LateralFrac; /* Fraction of total percolation that is routed laterally */
  float DownhillCapacity; /* Amount of water that neighbor cell's layer can accept */
  int i, k, nx, ny;			        /* counter */
  float *Ice;               /* Ice of each layer (m3/m3; NULL without soil freezing) */
  float Liquid;             /* Liquid water content of a layer (m3/m3) */

  Ice = Options->SoilFreezing ? SoilMap[y][x].Ice : NULL;

  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NSoilLayers; i++)
    DeepLayerDepth -= RootDepth[i];
  
  if (*TableDepth <= 0) { /* watertable above surface */
    *IExcess += Infiltration;
    if (InfiltOption == DYNAMIC) 
      Infiltration = 0.;
  }
  else {
    Moist[0] += Infiltration / (RootDepth[0] * Adjust[0]);
  }

  /* From top to bottom soil layer */
  LayerBottomDepth = 0.0;
  for (i = 0; i < NSoilLayers; i++) {
    
    LayerBottomDepth += RootDepth[i];
    Liquid = (Ice != NULL) ? Moist[i] - Ice[i] : Moist[i];

    /* No movement if soil moisture is below field capacity (with ice, the
       liquid water, unless the layer holds more than its porosity) */
    if (Liquid > FCap[i] || (Ice != NULL && Moist[i] > Porosity[i])) {
      Exponent = 2.0 / PoreDist[i] + 3.0;

      if (Liquid > Porosity[i]) {
        /* This can happen because the moisture content can exceed the
        porosity the way the algorithm is implemented */
        Drainage = Ks[i] * cosTheta;
        InterFlow = Ks[i] * KsAnisotropy * sinTheta;
      }
      else if (Liquid > FCap[i]) {
        Multiplier = pow((double)(Liquid/Porosity[i]), (double)Exponent);
        Drainage = Ks[i] * Multiplier * cosTheta;
        InterFlow = Ks[i] * Multiplier * KsAnisotropy * sinTheta;
      }
      else {
        Drainage = 0.0;
        InterFlow = 0.0;
      }
      
      /* Only permit interflow in soil layers that do NOT contain the water table
         to prevent conflicting routing directions and redundancy
         between unsaturated and saturated flow */
      if (!(Options->UseInterflow) ||
          (*TableDepth < LayerBottomDepth)) {
        InterFlow = 0.0;
      }
      
      Drainage += InterFlow;
      if (Drainage > 1e-10) {
        LateralFrac = InterFlow / Drainage;
      } else {
        LateralFrac = 0.0;
      }
      
      /* Average across timesteps and convert to m */
      Perc[i] = 0.5 * (Perc[i] + Drainage * Dt) * PercArea[i];

      MaxSoilWater = RootDepth[i] * Porosity[i] * Adjust[i];
      SoilWater = RootDepth[i] * Moist[i] * Adjust[i];
      FieldCapacity = RootDepth[i] * FCap[i] * Adjust[i];
      if (Ice != NULL)
        FieldCapacity += RootDepth[i] * Ice[i] * Adjust[i];

      /* No unsaturated flow if the moisture content drops below field
         capacity (with ice, the liquid water) */
      if ((SoilWater - Perc[i]) < FieldCapacity)
        Perc[i] = SoilWater - FieldCapacity;
      if (Ice != NULL && Perc[i] < 0.0)
        Perc[i] = 0.0;

      /* If moisture > porosity add the additional to percolation */
      SoilWater -= Perc[i];
      if (SoilWater > MaxSoilWater)
        Perc[i] += SoilWater - MaxSoilWater;

      /* Loop through adjacent cells and update the interflow (current layer) in down-slope cells */
      InterFlow = 0.0;
      for (k = 0; k < NDIRS; k++) {
        if (TopoMap[y][x].Dir[k] > 0) {
          nx = xdirection[k] + x;
          ny = ydirection[k] + y;
          if (valid_cell(Map, nx, ny) && INBASIN(TopoMap[ny][nx].Mask)) {
            
            /* Follows the topographic flow conventions of RouteSurface() */
            InterFlowK = Perc[i] * LateralFrac * ((float) TopoMap[y][x].Dir[k] / (float) TopoMap[y][x].TotalDir);
            
            /* Only allow enough lateral flow to raise neighboring layer to saturation */
            DownhillCapacity = (SoilMap[ny][nx].Porosity[i] - SoilMap[ny][nx].Moist[i] - SoilMap[ny][nx].InterFlow[i]) *
                                SoilMap[ny][nx].LayerThickness[i] * Network[ny][nx].Adjust[i];
            if (DownhillCapacity < 0.0)
              DownhillCapacity = 0.0;
            if (InterFlowK > DownhillCapacity) {
              InterFlowK = DownhillCapacity;
              Perc[i] -= (InterFlowK - DownhillCapacity);
            }
            InterFlow += InterFlowK;
            SoilMap[ny][nx].InterFlow[i] += InterFlowK / (SoilMap[ny][nx].LayerThickness[i] * Network[ny][nx].Adjust[i]);
          }
        }
      }
      
      /* Update the moisture content in the current layer and the layer immediately below it */
      
      if (i < (NSoilLayers - 1))
        NextLayerDepth = RootDepth[i + 1];
      else
        NextLayerDepth = DeepLayerDepth;
      
      Moist[i] -= Perc[i] / (RootDepth[i] * Adjust[i]);
      Moist[i+1] += (Perc[i] - InterFlow) / (NextLayerDepth * Adjust[i+1]);
    }
    else
      Perc[i] = 0.0;

    /* Convert back to straight 1-d flux */
    Perc[i] /= PercArea[i];
  }
  
  /* Calculate the depth of the water table based on the soil moisture
  profile and adjust the soil moisture profile, to assure that the soil
  moisture is never more than the maximum allowed soil moisture amount,
  i.e. the porosity.  A negative water table depth means that the water is
  ponding on the surface.  This amount of water becomes surface IExcess */

  *TableDepth = WaterTableDepth(NSoilLayers, TotalDepth, RootDepth, Porosity, FCap, Adjust, Moist);
  
  if (*TableDepth < 0.0) {
    *IExcess += -(*TableDepth);
    if (InfiltOption == DYNAMIC) {
      if (Infiltration > -(*TableDepth))
        Infiltration += *TableDepth;
      else
        Infiltration = 0.;
    }

    *TableDepth = 0.0;
  }
}
