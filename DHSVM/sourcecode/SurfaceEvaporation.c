
/*
 SUMMARY:      SurfaceEvaporation.c - Calculate evaporation from various surfaces
 
 Note this is adapted (by Eli) from SoilEvaporation.c originally by Bart Nijssen
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"
#include "settings.h"
#include "DHSVMChannel.h"

/*****************************************************************************
 PondEvaporation()
*****************************************************************************/
float PondEvaporation(int Dt, float DXDY, float ChannelArea,
                      float Temp, float Slope, float Gamma,
                      float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                      float Evapotranspiration, float *IExcess)
{
  float EPot;			/* Potential evaporation from open water during timestep (m) */
  float PondEvap;	/* Amount of evaporation from surface ponding (m) */
  
  EPot = (Slope * NetRad + AirDens * CP * Vpd / LowerRa) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPot -= Evapotranspiration;
  
  /* EPot must only include non-channel area, since channel evaporation is separate */
  if (ChannelArea > 0.0)
    EPot *= 1.0 - ChannelArea / DXDY;
  
  if (EPot < 0.0)
    EPot = 0.0;
  
  /* Water evaporates from surface ponding up to IExcess or EPot */
  if (EPot > *IExcess)
    PondEvap = *IExcess;
  else
    PondEvap = EPot;
  
  if (PondEvap < 0.0)
    PondEvap = 0.0;
  
  *IExcess -= PondEvap;
  
  return PondEvap;
}

/*****************************************************************************
 ChannelEvaporation()
 *****************************************************************************/
void ChannelEvaporation(int Dt, float DXDY,
                        float Temp, float Slope, float Gamma,
                        float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                        float Evapotranspiration, int x, int y, CHANNEL *ChannelData)
{
  float EPot;			/* Potential evaporation rate during timestep (m) */
  float EPotCell;  /* Unsatisfied atmospheric demand during timestep (m) */
  
  /* Maximum atmospheric demand if the whole grid cell is covered by water */
  EPot = (Slope * NetRad + AirDens * CP * Vpd / LowerRa) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPotCell = EPot;
  EPotCell -= Evapotranspiration;
  
  /* Water evaporates from each channel segment at the potential rate
   up to the limit imposed by the total grid cell EPot */
  
  if (EPotCell > 0.0)
    channel_grid_calc_evaporation(ChannelData->stream_map, x, y, EPot, (EPotCell * DXDY));
  else
    channel_grid_calc_evaporation(ChannelData->stream_map, x, y, 0.0, 0.0);
}

/*****************************************************************************
 ChannelSoilEvaporation()
 *****************************************************************************/
float ChannelSoilEvaporation(int Dt, float DXDY,
                      float Temp, float Slope, float Gamma, float Lv,
                      float AirDens, float Vpd, float NetRad, float RaSoil,
                      float Evapotranspiration, float *Porosity, float *FCap, float *Ks,
                      float *Press, float *m, float LayerThickness,
                      float *MoistContent, float *Adjust,
                      int x, int y, CHANNEL *ChannelData, int CutBankZone)
{
  float EPot;			/* Potential evaporation rate during timestep (m) */
  float EPotCell;  /* Unsatisfied atmospheric demand during timestep (m) */
  float SoilEvap;		/* Amount of evaporation directly from the soil (m) */
  
  EPot = (Slope * NetRad + AirDens * CP * Vpd / RaSoil) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;
  
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;
  
  EPotCell = EPot;
  EPotCell -= Evapotranspiration;
  
  if (EPotCell > 0.0)
    SoilEvap = channel_grid_dry_evaporation(ChannelData->stream_map, x, y,
                                            EPot, EPotCell, DXDY,
                                            Dt, Porosity, FCap, Ks,
                                            Press, m, LayerThickness,
                                            MoistContent, Adjust, CutBankZone);
  else
    SoilEvap = 0.0;
  
  return SoilEvap;
}

/*****************************************************************************
 SoilEvaporation()

 Evaporation from the soil of the cell's bare fraction BareFract (outside
 the litter, under and outside the understory); *PotEvap returns the
 potential evaporation after the overlying demand (m per time step). Two
 stages, per unit bare ground area, which take BareFract times their depth
 from the layer's mean moisture:
 1. Above the field capacity the soil evaporates at the potential rate, down
    to the field capacity, and *DesorptionAcc is 0.
 2. At or below the field capacity, the rest of the step's potential
    evaporation is limited by the desorption capacity of the step, following
    the cumulative desorption S sqrt(t) of a soil that has delivered
    *DesorptionAcc since stage 1 ended: sqrt(E^2 + D^2) - E, with D the
    capacity of a freshly wetted soil (Desorption()) at the moisture where
    stage 2 starts, down to half the wilting point. The stage 2 evaporation
    is added to *DesorptionAcc.
 *****************************************************************************/
float SoilEvaporation(int Dt,
                      float Temp, float Slope, float Gamma, float Lv,
                      float AirDens, float Vpd, float NetRad, float RaSoil,
                      float Evapotranspiration, float Porosity, float FCap, float Ks,
                      float Press, float m, float WP, float RootDepth,
                      float *MoistContent, float Adjust, float BareFract, float *PotEvap,
                      float *DesorptionAcc)
{
  float DesorptionVolume;	/* Amount of water the soil can deliver to the atmosphere during a timestep (m) */
  float EPot;			/* Potential evaporation from soil during timestep (m) */
  float GroundEvap;		/* Evaporation per unit ground area (m) */
  float StageEvap;		/* Evaporation per unit ground area of one stage (m) */
  float Floor;			/* Lowest moisture content of stage 2 (0-1) */
  double Acc;

  /* Eq.4 Wigmosta et al [1994] */

  EPot = (Slope * NetRad + AirDens * CP * Vpd / RaSoil) /
    (WATER_DENSITY * Lv * (Slope + Gamma)) * Dt;

  /* The potential evaporation is reduced by the evaporation already
     calculated for the overlying layers (prior snow vapor deposition does
     not increase it) */
  if (Evapotranspiration < 0.0)
    Evapotranspiration = 0.0;

  EPot -= Evapotranspiration;
  if (EPot < 0.0)
    EPot = 0.0;
  *PotEvap = EPot;

  GroundEvap = 0.0;

  /* Stage 1 */
  if (*MoistContent > FCap) {
    *DesorptionAcc = 0.0;
    StageEvap = (*MoistContent - FCap) * RootDepth;
    if (EPot < StageEvap) {
      StageEvap = EPot;
      *MoistContent -= BareFract * StageEvap / RootDepth;
    }
    else
      *MoistContent = FCap + (1. - BareFract) * (*MoistContent - FCap);
    GroundEvap = StageEvap;
    EPot -= StageEvap;
  }

  /* Stage 2 */
  Floor = 0.5 * WP;
  if (EPot > 0.0 && *MoistContent > Floor) {

    DesorptionVolume = Desorption(Dt, *MoistContent, Porosity, Ks, Press, m);

    if (*DesorptionAcc > 0.0) {
      Acc = *DesorptionAcc;
      DesorptionVolume = sqrt(Acc * Acc + (double) DesorptionVolume * DesorptionVolume) - Acc;
    }

    /* Eq.8 Wigmosta et al [1994] */
    StageEvap = MIN(EPot, DesorptionVolume);

    if (StageEvap >= (*MoistContent - Floor) * RootDepth) {
      StageEvap = (*MoistContent - Floor) * RootDepth;
      *MoistContent = Floor + (1. - BareFract) * (*MoistContent - Floor);
    }
    else
      *MoistContent -= BareFract * StageEvap / RootDepth;

    *DesorptionAcc += StageEvap;
    GroundEvap += StageEvap;
  }

  return GroundEvap * BareFract * Adjust;
}
