
/*
 When calculating lateral outflow, available water was calculated using all 3 root zone layers 
 and the deep layer underneath but outflow was extracted from the bottom deep soil layer. This 
 might lead to negative soil moisture in deep soil while the layer above it remains saturated. This tends
 to happen more often in dry climate. In order to avoid negative deep soil moisture, here I add a loop 
 to redistribution of water extraction. The outflow starts from the (top) water table layer, extract the 
 excess water larger than field capacity. While there's not enough water from the current layer, extract 
 water from one layer below it until it reaches the bottom layer.
 Added 06/09/2016 by Zhuoran Duan(zhuoran.duan@pnnl.gov)
 
 Further modified by Eli Boardman to subtract water from deep infiltration
 (essentially a leaky bottom to the watershed)

 Ice (the layer ice, m3/m3; NULL without soil freezing) takes no part in
 the outflow: a layer holding ice gives only its liquid water above the
 field capacity, and none if that is negative. The inflow fills the pore
 space the ice leaves.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "settings.h"
#include "functions.h"
#include "soilmoisture.h"

/*****************************************************************************/
void DistributeSatflow(int Dt, float DX, float DY, float SatFlow,
  int NSoilLayers, float TotalDepth, float *RootDepth,
  float *Porosity, float *FCap, float *Adjust,
  float *TableDepth, float *IExcess, float *Moist, float *Ice)

{
  float DeepLayerDepth;		/* depth of the layer below the deepest root layer */
  int i;			        /* counter */
  float DeepFCap;		/* field capacity of the layer below the deepest root layer */
  float DeepPorosity;		/* porosity of the layer below the deepest root layer */
  float DeepAvaWater;
  float AvaWater;
  float DeepWaterGap;
  float WaterGap;
  float ExtracWater;
  float DeepExtracWater;
  float Depth;

  DeepPorosity = Porosity[NSoilLayers];
  DeepFCap = FCap[NSoilLayers];
  
  DeepLayerDepth = TotalDepth;
  for (i = 0; i < NSoilLayers; i++)
    DeepLayerDepth -= RootDepth[i];
  
  /* New algorithm for SatFlow, remove water from top layer to bottom layer */
  if (SatFlow < 0.0) {
    
    Depth = 0.0;
    for (i = 0; i < NSoilLayers && Depth < TotalDepth; i++) {
      AvaWater = 0.0;
      ExtracWater = 0.0;
      if (RootDepth[i] < (TotalDepth - Depth))
        Depth += RootDepth[i];
      else
        Depth = TotalDepth;

      /* If water table is in the ith root zone layer */
      if (Depth > *TableDepth) {
        /* Fully saturated in the current layer */
        if ((Depth - *TableDepth) > RootDepth[i])
          AvaWater = (Porosity[i]-FCap[i]) * RootDepth[i] * Adjust[i];
        else
          AvaWater = (Moist[i]-FCap[i]) * RootDepth[i] * Adjust[i];
        if (Ice != NULL && Ice[i] > 0.0) {
          AvaWater -= Ice[i] * RootDepth[i] * Adjust[i];
          if (AvaWater < 0.0)
            AvaWater = 0.0;
        }
      }
      ExtracWater = (-SatFlow > AvaWater) ? -AvaWater : SatFlow;

      Moist[i] += ExtracWater / (RootDepth[i] * Adjust[i]);

      SatFlow -= ExtracWater;
      if (fequal(SatFlow, 0.0))
        break;
    }

    if (SatFlow < 0.0) {
      DeepAvaWater = 0.0;
      DeepExtracWater = 0.0;
      if (Depth < TotalDepth) {
        Depth = TotalDepth;
        if ((Depth - *TableDepth) > DeepLayerDepth)
          DeepAvaWater = (DeepPorosity - DeepFCap) * DeepLayerDepth * Adjust[NSoilLayers];
        else
          DeepAvaWater = (Moist[NSoilLayers] - DeepFCap) * DeepLayerDepth * Adjust[NSoilLayers];;
        if (Ice != NULL && Ice[NSoilLayers] > 0.0) {
          DeepAvaWater -= Ice[NSoilLayers] * DeepLayerDepth * Adjust[NSoilLayers];
          if (DeepAvaWater < 0.0)
            DeepAvaWater = 0.0;
        }
      }

      DeepExtracWater = (-SatFlow > DeepAvaWater) ? -DeepAvaWater : SatFlow;

      Moist[NSoilLayers] += DeepExtracWater / (DeepLayerDepth * Adjust[NSoilLayers]);

      SatFlow -= DeepExtracWater;
    }
  } else if (SatFlow > 0.0) {

    Depth = DeepLayerDepth;
    DeepWaterGap = 0.0;
    DeepExtracWater = 0.0;
    if (Depth > (TotalDepth - *TableDepth)) {
      DeepWaterGap = (DeepPorosity - Moist[NSoilLayers]) * DeepLayerDepth * Adjust[NSoilLayers];
      DeepExtracWater = (SatFlow > DeepWaterGap) ? DeepWaterGap : SatFlow;
      SatFlow -= DeepExtracWater;
      Moist[NSoilLayers] += DeepExtracWater / (DeepLayerDepth * Adjust[NSoilLayers]);
    }

    if (SatFlow > 0.0) {
      for (i = NSoilLayers - 1; i >= 0; i--) {
        WaterGap = 0.0;
        ExtracWater = 0.0;
        Depth += RootDepth[i];
        if (Depth > (TotalDepth - *TableDepth)) {
          WaterGap = (Porosity[i] - Moist[i]) * RootDepth[i] * Adjust[i];
          ExtracWater = (SatFlow > WaterGap) ? WaterGap : SatFlow;
          SatFlow -= ExtracWater;
          Moist[i] += ExtracWater / (RootDepth[i] * Adjust[i]);
        }
        if (fequal(SatFlow, 0.0))
          break;
      }
    }
  }

  if (SatFlow > 0.0)
    *IExcess += SatFlow;
}
