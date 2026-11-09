
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "functions.h"

/*****************************************************************************
EvapoTranspiration()

Evaporation of the interception and transpiration of the vegetation layer
at position Layer of the cell, with the parameters of its story's class
(Story: the overstory or the understory class) and the stomatal scheme
(Canopy, which also receives the scheme's canopy resistance, net
assimilation per cell area and soil moisture factor)
*****************************************************************************/
void EvapoTranspiration(int Layer, int Dt, PIXMET *Met,
  float NetRad, float Rp, STORYTABLE *Story, SOILTABLE *SType,
  float MoistureFlux, float *Moist, float *SoilTemp, float *Int,
  float *EPot, float *EInt, float **ESoil, float *EAct, float *ETot,
  float *Adjust, float *LayerThickness, float *FCap, float Ra, VEGPIX *LocalVeg,
  CANOPYSCHEME *Canopy, float Overflow, float *Throughfall)
{
  float *Rc;			/* canopy resistance associated with
                        conditions in each soil layer (s/m) */
  float *Share;			/* each soil layer's share of the transpiration */
  float DryEvapTime;	/* amount of time remaining during a timestep
                        after the interception storage is depleted (sec) */
  float F;			    /* Fractional coverage by vegetation layer */
  float SoilMoisture;	/* Amount of water in each soil layer (m) */
  float Available;		/* Water of the layer above the wilting point (m) */
  float WetArea;		/* relative leaf area wetted by interception storage */
  float WetEvapRate;	/* evaporation rate from wetted fraction per unit ground area (m/s) */
  float WetEvapTime;	/* amount of time needed to evaporate the amount of water
                        in interception storage (sec) */
  float MaxInt;			/* Maximum interception storage as if the entire pixel
                        is covered (m) */
  float Wet;            /* Water the wet leaves can evaporate during the step:
                        the storage plus the step's rain its full storage
                        passed on (m) */
  float MaxRefill;      /* Most rain the storage can take back in the step (m) */
  float Refill;         /* Rain taken back from the throughfall to refill the
                        storage as it evaporates (m) */
  int i;			    /* counter */


  F = LocalVeg->Fract[Layer];
  NetRad /= F;

  /* Convert the water amounts related to partial canopy cover to a pixel depth
  as if the entire pixel is covered. These depths will be converted back later on. */
  *Int /= F;
  MoistureFlux /= F;
  MaxInt = LocalVeg->MaxInt[Layer] / F;

  /* allocate memory for the canopy resistance and share arrays */
  if (!(Rc = (float *)calloc(2 * SType->NLayers, sizeof(float))))
    ReportError("EvapoTranspiration()", 1);
  Share = Rc + SType->NLayers;

  /* Calculate the evaporation rate in m/s */
  EPot[Layer] = (Met->Slope * NetRad + Met->AirDens * CP * Met->Vpd / Ra) /
    (WATER_DENSITY * Met->Lv * (Met->Slope + Met->Gamma));

  /* The potential evaporation rate accounts for the amount of moisture that
  the atmosphere can absorb. If we do not account for the amount of
  evaporation from overlying evaporation, we can end up with the situation
  that all vegetation layers and the soil layer transpire/evaporate at the
  potential rate, resulting in an overprediction of the actual evaporation
  rate. Thus we subtract the amount of evaporation that has already
  been calculated for overlying layers from the potential evaporation.
  Another mechanism that could be used to account for this would be to
  decrease the vapor pressure deficit while going down through the canopy
  (not implemented here) */
  EPot[Layer] -= MoistureFlux / Dt;

  if (EPot[Layer] < 0)
    EPot[Layer] = 0;

  /* WetArea = pow(*Int/MaxInt, (double) 2.0/3.0); */

  WetArea = cbrt(*Int / MaxInt);
  WetArea = MIN(WetArea, 1);
  WetArea = WetArea * WetArea;
  
  /* calculate the amount of water that can evaporate from the interception
  storage.  Given this evaporation rate, calculate the amount of time it
  will take to evaporate the entire amount of intercepted water.  If this
  time period is shorter than the length of the time interval, the
  previously wetted leaves can transpire during the remaining part of the
  time interval */

  /* The step's rain falls throughout the step: what the full storage passed
  on (Overflow, at most the throughfall left) refills the storage as it
  evaporates, at most one storage per hour of the step including the storage
  itself, so the wet leaves can evaporate the storage plus that rain, and the
  storage ends at most full. The rain taken back leaves the throughfall.
  Without overflow this is the storage alone. */
  if (Overflow > *Throughfall)
    Overflow = *Throughfall;
  if (Overflow < 0.0)
    Overflow = 0.0;
  Overflow /= F;
  MaxRefill = MaxInt * ((float) Dt / SECPHOUR - 1.0);
  if (MaxRefill < 0.0)
    MaxRefill = 0.0;
  if (Overflow > MaxRefill)
    Overflow = MaxRefill;
  Wet = *Int + Overflow;

  WetEvapRate = WetArea * EPot[Layer];
  if (WetEvapRate > 0) {
    WetEvapTime = Wet / WetEvapRate;
    if (WetEvapTime > Dt) {
      WetEvapTime = Dt;
    }
  }
  else if (*Int > 0)
    WetEvapTime = Dt;
  else
    WetEvapTime = 0;

  if (WetEvapRate > 0) {
    if (WetEvapTime < Dt) {
      EInt[Layer] = Wet;
      Refill = Overflow;
      *Int = 0.0;
      DryEvapTime = Dt - WetEvapTime;
    }
    else {
      EInt[Layer] = Dt * WetEvapRate;
      Refill = MIN(EInt[Layer], Overflow);
      *Int += Refill - EInt[Layer];
      DryEvapTime = 0.0;
    }
    *Throughfall -= Refill * F;
  }
  else {
    EInt[Layer] = 0.0;
    if (*Int > 0)
      DryEvapTime = 0.0;
    else
      DryEvapTime = Dt;
  }

  /* Correct the evaporation from interception and the interception storage for
  the fractional overstory coverage */
  EInt[Layer] *= F;
  *ETot += EInt[Layer];
  *Int *= F;

  /* calculate the canopy resistances associated with the conditions in
  each of the soil layers, and the layers' shares of the transpiration */
  CanopyResistanceLayers(Layer, Met, Rp, Story, SType, FCap, Moist, SoilTemp, LocalVeg,
    Canopy, Rc, Share);
  Canopy->An *= F;

  /* calculate the transpiration rate for the current vegetation layer,
  and adjust the soil moisture content in each of the soil layers */
  for (i = 0; i < SType->NLayers; i++) {
    ESoil[Layer][i] = (Met->Slope + Met->Gamma) /
      (Met->Slope + Met->Gamma * (1 + Rc[i] / Ra)) * Share[i] *
      EPot[Layer] * Adjust[i];

    /* calculate the amounts of water transpirated during each timestep based
    on the evaporation and transpiration rates.  While there is still water
    in interception storage only the area that is not covered by intercep-
    ted water will transpire.  When all of the interception storage has
    disappeared all leaves will contribute to the transpiration */
    ESoil[Layer][i] *= WetEvapTime * (1 - WetArea) + DryEvapTime;

    /* correct the transpiration for the fractional coverage, at most the
    layer's water above the wilting point, and update the soil moisture */
    ESoil[Layer][i] *= F;
    SoilMoisture = Moist[i] * LayerThickness[i] * Adjust[i];
    Available = (Moist[i] - SType->WP[i]) * LayerThickness[i] * Adjust[i];
    if (Available < 0)
      Available = 0;
    if (ESoil[Layer][i] > Available)
      ESoil[Layer][i] = Available;
    SoilMoisture -= ESoil[Layer][i];

    Moist[i] = SoilMoisture / (LayerThickness[i] * Adjust[i]);

  }

  for (i = 0, EAct[Layer] = 0; i < SType->NLayers; i++) {
    *ETot += ESoil[Layer][i];
    EAct[Layer] += ESoil[Layer][i];
  }
  
  free(Rc);
}

