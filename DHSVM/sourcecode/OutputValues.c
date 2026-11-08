/*
 * DESCRIPTION:  Output values: the value code of each output variable of
 *               the registry, the value of a variable layer in one cell
 *               (point and map outputs) and of a channel segment; the snow
 *               cold content derived from the snow state.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "output.h"
#include "registry.h"
#include "photosynthesis.h"
#include "snow.h"

static const struct {
  const char *Name;
  int Code;
} CodeNames[] = {
  {"Met.Precipitation", V_MET_PRECIPITATION},
  {"Met.Snowfall", V_MET_SNOWFALL},
  {"Met.Rainfall", V_MET_RAINFALL},
  {"Met.AirTemperature", V_MET_AIRTEMPERATURE},
  {"Met.WindSpeed", V_MET_WINDSPEED},
  {"Met.RelativeHumidity", V_MET_RELATIVEHUMIDITY},
  {"Met.PrecipitationWaterYear", V_MET_PRECIPITATIONWATERYEAR},
  {"Rad.ShortwaveInObserved", V_RAD_SHORTWAVEINOBSERVED},
  {"Rad.ShortwaveIn", V_RAD_SHORTWAVEIN},
  {"Rad.NetShortwave", V_RAD_NETSHORTWAVE},
  {"Rad.LongwaveIn", V_RAD_LONGWAVEIN},
  {"Rad.CanopyNetShortwave", V_RAD_CANOPYNETSHORTWAVE},
  {"Rad.CanopyLongwaveIn", V_RAD_CANOPYLONGWAVEIN},
  {"Rad.NetRadiation", V_RAD_NETRADIATION},
  {"Snow.Covered", V_SNOW_COVERED},
  {"Snow.DaysSinceSnowfall", V_SNOW_DAYSSINCESNOWFALL},
  {"Snow.SWE", V_SNOW_SWE},
  {"Snow.Outflow", V_SNOW_OUTFLOW},
  {"Snow.PackLiquid", V_SNOW_PACKLIQUID},
  {"Snow.PackTemperature", V_SNOW_PACKTEMPERATURE},
  {"Snow.SurfaceLiquid", V_SNOW_SURFACELIQUID},
  {"Snow.SurfaceTemperature", V_SNOW_SURFACETEMPERATURE},
  {"Snow.ColdContent", V_SNOW_COLDCONTENT},
  {"Snow.Albedo", V_SNOW_ALBEDO},
  {"Snow.Depth", V_SNOW_DEPTH},
  {"Snow.Vapor", V_SNOW_VAPOR},
  {"Snow.NetShortwave", V_SNOW_NETSHORTWAVE},
  {"Snow.NetLongwave", V_SNOW_NETLONGWAVE},
  {"Snow.SensibleHeat", V_SNOW_SENSIBLEHEAT},
  {"Snow.LatentHeat", V_SNOW_LATENTHEAT},
  {"Snow.RainHeat", V_SNOW_RAINHEAT},
  {"Snow.MeltEnergy", V_SNOW_MELTENERGY},
  {"Snow.PeakSWE", V_SNOW_PEAKSWE},
  {"Snow.PeakSWEDay", V_SNOW_PEAKSWEDAY},
  {"Snow.MeltOutDay", V_SNOW_MELTOUTDAY},
  {"Snow.AccumulationWaterYear", V_SNOW_ACCUMULATIONWATERYEAR},
  {"Snow.OutflowWaterYear", V_SNOW_OUTFLOWWATERYEAR},
  {"Canopy.InterceptedRain", V_CANOPY_INTERCEPTEDRAIN},
  {"Canopy.InterceptedSnow", V_CANOPY_INTERCEPTEDSNOW},
  {"Canopy.TemporaryStorage", V_CANOPY_TEMPORARYSTORAGE},
  {"Canopy.Temperature", V_CANOPY_TEMPERATURE},
  {"Canopy.Storage", V_CANOPY_STORAGE},
  {"Canopy.SnowVapor", V_CANOPY_SNOWVAPOR},
  {"Canopy.Resistance", V_CANOPY_RESISTANCE},
  {"Canopy.Assimilation", V_CANOPY_ASSIMILATION},
  {"Canopy.MoistureStress", V_CANOPY_MOISTURESTRESS},
  {"Canopy.Acclimation", V_CANOPY_ACCLIMATION},
  {"Canopy.AcclimationTemperature", V_CANOPY_ACCLIMATIONTEMPERATURE},
  {"Canopy.GrowthTemperature", V_CANOPY_GROWTHTEMPERATURE},
  {"Evap.Total", V_EVAP_TOTAL},
  {"Evap.Potential", V_EVAP_POTENTIAL},
  {"Evap.Actual", V_EVAP_ACTUAL},
  {"Evap.Interception", V_EVAP_INTERCEPTION},
  {"Evap.Transpiration", V_EVAP_TRANSPIRATION},
  {"Evap.Soil", V_EVAP_SOIL},
  {"Evap.Channel", V_EVAP_CHANNEL},
  {"Soil.Moisture", V_SOIL_MOISTURE},
  {"Soil.Percolation", V_SOIL_PERCOLATION},
  {"Soil.Interflow", V_SOIL_INTERFLOW},
  {"Soil.TableDepth", V_SOIL_TABLEDEPTH},
  {"Soil.SatFlow", V_SOIL_SATFLOW},
  {"Soil.WaterLevel", V_SOIL_WATERLEVEL},
  {"Soil.SurfaceWater", V_SOIL_SURFACEWATER},
  {"Soil.InfiltrationAccumulated", V_SOIL_INFILTRATIONACCUMULATED},
  {"Soil.DesorptionAccumulated", V_SOIL_DESORPTIONACCUMULATED},
  {"Soil.Water", V_SOIL_WATER},
  {"Soil.DeepFlux", V_SOIL_DEEPFLUX},
  {"Soil.LatentHeat", V_SOIL_LATENTHEAT},
  {"Soil.Ice", V_SOIL_ICE},
  {"Soil.Temperature", V_SOIL_TEMPERATURE},
  {"Soil.NodeTemperature", V_SOIL_NODETEMPERATURE},
  {"Soil.NodeIce", V_SOIL_NODEICE},
  {"Soil.SurfaceTemperature", V_SOIL_SURFACETEMPERATURE},
  {"Soil.FrostDepth", V_SOIL_FROSTDEPTH},
  {"Soil.FrostTop", V_SOIL_FROSTTOP},
  {"Soil.GroundHeatFlux", V_SOIL_GROUNDHEATFLUX},
  {"Soil.SensibleHeatFlux", V_SOIL_SENSIBLEHEATFLUX},
  {"Surface.Detention", V_SURFACE_DETENTION},
  {"Surface.Runoff", V_SURFACE_RUNOFF},
  {"Stream.Interception", V_STREAM_INTERCEPTION},
  {"Stream.Infiltration", V_STREAM_INFILTRATION},
  {"Lake.Storage", V_LAKE_STORAGE},
  {"Lake.Level", V_LAKE_LEVEL},
  {"Balance.WaterInput", V_BALANCE_WATERINPUT},
  {"Balance.Storage", V_BALANCE_STORAGE},
  {"Balance.Error", V_BALANCE_ERROR},
  {"Balance.CumulativeError", V_BALANCE_CUMULATIVEERROR},
  {"Balance.RainOrMelt", V_BALANCE_RAINORMELT},
  {"Soil.SaturatedFraction", V_SOIL_SATURATEDFRACTION},
  {"Channel.NetworkLateralInflow", V_CHANNEL_NETWORKLATERALINFLOW},
  {"Channel.NetworkOutflow", V_CHANNEL_NETWORKOUTFLOW},
  {"Channel.NetworkStorageChange", V_CHANNEL_NETWORKSTORAGECHANGE},
  {"Channel.NetworkInfiltration", V_CHANNEL_NETWORKINFILTRATION},
  {"Channel.NetworkEvaporation", V_CHANNEL_NETWORKEVAPORATION},
  {"Channel.NetworkError", V_CHANNEL_NETWORKERROR},
  {"Channel.NetworkStorage", V_CHANNEL_NETWORKSTORAGE},
  {"Channel.Inflow", V_CHANNEL_INFLOW},
  {"Channel.LateralInflow", V_CHANNEL_LATERALINFLOW},
  {"Channel.Outflow", V_CHANNEL_OUTFLOW},
  {"Channel.Storage", V_CHANNEL_STORAGE},
  {"Channel.StorageChange", V_CHANNEL_STORAGECHANGE},
  {"Channel.Infiltration", V_CHANNEL_INFILTRATION},
  {"Channel.Evaporation", V_CHANNEL_EVAPORATION}
};

static int VarOfCode[NVALUECODES];	/* Registry row of each value code */
static int *CodeOfVar = NULL;		/* Value code of each registry row, -1 if none */

/*****************************************************************************
  InitValueCodes()

  Binds the value codes to the registry rows: every code names a registry
  variable, and every variable that an output kind accepts has a code
*****************************************************************************/
void InitValueCodes(void)
{
  int NCodes = (int) (sizeof(CodeNames) / sizeof(CodeNames[0]));
  int i, Var;

  if (NCodes != NVALUECODES)
    ReportError("InitValueCodes", 26);
  if (!(CodeOfVar = (int *) malloc(RegistryCount() * sizeof(int))))
    ReportError("InitValueCodes", 1);
  for (Var = 0; Var < RegistryCount(); Var++)
    CodeOfVar[Var] = -1;
  for (i = 0; i < NCodes; i++) {
    Var = RegistryFind(CodeNames[i].Name);
    if (Var < 0)
      ReportError((char *) CodeNames[i].Name, 26);
    VarOfCode[CodeNames[i].Code] = Var;
    CodeOfVar[Var] = CodeNames[i].Code;
  }
  for (Var = 0; Var < RegistryCount(); Var++) {
    if (RegistryKinds(Var) != 0 && CodeOfVar[Var] < 0)
      ReportError((char *) RegistryRow(Var)->Name, 26);
  }
}

/* Value code of a registry row (-1 if it has none), registry row of a code */
int ValueCode(int Var)
{
  return CodeOfVar[Var];
}

int CodeVar(int Code)
{
  return VarOfCode[Code];
}

/*****************************************************************************
  CellValue()

  Value of layer Layer of the variable with value code Code in mask cell
  (y, x). *Has is FALSE (and the value NA) where the cell does not have the
  layer: a canopy layer the cell lacks, the lake level off the
  lakes, the canopy resistance and moisture factor of a buried understory,
  a canopy temperature memory not yet set (the start of a cold start).
*****************************************************************************/
float CellValue(int Code, int Layer, MODELVIEW *Model, int y, int x, int *Has)
{
  EVAPPIX *Evap = &(Model->Evap[y][x]);
  PRECIPPIX *Precip = &(Model->Precip[y][x]);
  PIXRAD *Rad = &(Model->Rad[y][x]);
  SNOWPIX *Snow = &(Model->Snow[y][x]);
  SOILPIX *Soil = &(Model->Soil[y][x]);
  VEGPIX *Veg = &(Model->Veg[y][x]);
  NETSTRUCT *Network = &(Model->Network[y][x]);
  int L = Model->Map->SoilLayers;
  int Index;
  int i;
  float Value;
  float DeepDepth;

  *Has = TRUE;

  /* Position of a canopy layer among the cell's vegetation layers */
  switch (Code) {
  case V_RAD_CANOPYNETSHORTWAVE:
  case V_RAD_CANOPYLONGWAVEIN:
  case V_CANOPY_INTERCEPTEDRAIN:
  case V_CANOPY_INTERCEPTEDSNOW:
  case V_EVAP_INTERCEPTION:
  case V_CANOPY_RESISTANCE:
  case V_CANOPY_ASSIMILATION:
  case V_CANOPY_MOISTURESTRESS:
    Index = CanopyIndex(Veg, Layer);
    break;
  case V_EVAP_POTENTIAL:
  case V_EVAP_ACTUAL:
    Index = (Layer == 2) ? Veg->NVegLayers : CanopyIndex(Veg, Layer);
    break;
  case V_EVAP_TRANSPIRATION:
    Index = CanopyIndex(Veg, Layer / L);
    break;
  default:
    Index = 0;
    break;
  }
  if (Index < 0) {
    *Has = FALSE;
    return NA;
  }

  switch (Code) {
  case V_MET_PRECIPITATION: return Precip->Precip;
  case V_MET_SNOWFALL: return Precip->SnowFall;
  case V_MET_RAINFALL: return Precip->RainFall;
  case V_MET_AIRTEMPERATURE: return Rad->Tair;
  case V_MET_WINDSPEED: return Rad->Wind;
  case V_MET_RELATIVEHUMIDITY: return Rad->Rh;
  case V_MET_PRECIPITATIONWATERYEAR: return Precip->SumPrecip;
  case V_RAD_SHORTWAVEINOBSERVED: return Rad->ObsShortIn;
  case V_RAD_SHORTWAVEIN: return Rad->BeamIn + Rad->DiffuseIn;
  case V_RAD_NETSHORTWAVE: return Rad->PixelNetShort;
  case V_RAD_LONGWAVEIN: return Rad->PixelLongIn;
  case V_RAD_CANOPYNETSHORTWAVE: return Rad->NetShort[Index];
  case V_RAD_CANOPYLONGWAVEIN: return Rad->LongIn[Index];
  case V_RAD_NETRADIATION: return Rad->NetRadiation[0] + Rad->NetRadiation[1];
  case V_SNOW_COVERED: return (float) Snow->HasSnow;
  case V_SNOW_DAYSSINCESNOWFALL: return Snow->LastSnow;
  case V_SNOW_SWE: return Snow->Swq;
  case V_SNOW_OUTFLOW: return Snow->Outflow;
  case V_SNOW_PACKLIQUID: return Snow->PackWater;
  case V_SNOW_PACKTEMPERATURE: return Snow->TPack;
  case V_SNOW_SURFACELIQUID: return Snow->SurfWater;
  case V_SNOW_SURFACETEMPERATURE: return Snow->TSurf;
  case V_SNOW_COLDCONTENT: return SnowColdContent(Snow);
  case V_SNOW_ALBEDO: return Snow->Albedo;
  case V_SNOW_DEPTH: return SnowDepth(Snow->Swq);
  case V_SNOW_VAPOR: return Snow->VaporMassFlux;
  case V_SNOW_NETSHORTWAVE: return Snow->Qsw;
  case V_SNOW_NETLONGWAVE: return Snow->Qlw;
  case V_SNOW_SENSIBLEHEAT: return Snow->Qs;
  case V_SNOW_LATENTHEAT: return Snow->Qe;
  case V_SNOW_RAINHEAT: return Snow->Qp;
  case V_SNOW_MELTENERGY: return Snow->MeltEnergy;
  case V_SNOW_PEAKSWE: return Snow->MaxSwe;
  case V_SNOW_PEAKSWEDAY: return (float) Snow->MaxSweDate;
  case V_SNOW_MELTOUTDAY: return (float) Snow->MeltOutDate;
  case V_SNOW_ACCUMULATIONWATERYEAR: return Precip->SnowAccum;
  case V_SNOW_OUTFLOWWATERYEAR: return Precip->SnowMelt;
  case V_CANOPY_INTERCEPTEDRAIN: return Precip->IntRain[Index];
  case V_CANOPY_INTERCEPTEDSNOW: return Precip->IntSnow[Index];
  case V_CANOPY_TEMPORARYSTORAGE: return Precip->TempIntStorage;
  case V_CANOPY_TEMPERATURE: return Veg->Tcanopy;
  case V_CANOPY_STORAGE:
    Value = 0.0;
    for (i = 0; i < Veg->NVegLayers; i++)
      Value += Precip->IntRain[i] + Precip->IntSnow[i];
    return Value;
  case V_CANOPY_SNOWVAPOR: return Snow->CanopyVaporMassFlux;
  case V_CANOPY_RESISTANCE:
    *Has = (Evap->Rc[Index] != NA);
    return Evap->Rc[Index];
  case V_CANOPY_ASSIMILATION: return Evap->An[Index];
  case V_CANOPY_MOISTURESTRESS:
    *Has = (Evap->Beta[Index] != NA);
    return Evap->Beta[Index];
  case V_CANOPY_ACCLIMATION: return PhotoAcclimation(Veg->TAcclim);
  case V_CANOPY_ACCLIMATIONTEMPERATURE:
    *Has = (Veg->TAcclim > (float) PHOTO_UNSET);
    return *Has ? Veg->TAcclim : NA;
  case V_CANOPY_GROWTHTEMPERATURE:
    *Has = (Veg->TGrowth > (float) PHOTO_UNSET);
    return *Has ? Veg->TGrowth : NA;
  case V_EVAP_TOTAL: return Evap->ETot;
  case V_EVAP_POTENTIAL: return Evap->EPot[Index] * Model->Dt;
  case V_EVAP_ACTUAL: return Evap->EAct[Index];
  case V_EVAP_INTERCEPTION: return Evap->EInt[Index];
  case V_EVAP_TRANSPIRATION: return Evap->ESoil[Index][Layer % L];
  case V_EVAP_SOIL: return Evap->EvapSoil;
  case V_EVAP_CHANNEL: return Evap->EvapChannel;
  case V_SOIL_MOISTURE: return Soil->Moist[Layer];
  case V_SOIL_PERCOLATION: return Soil->Perc[Layer];
  case V_SOIL_INTERFLOW: return Soil->InterFlow[Layer];
  case V_SOIL_TABLEDEPTH: return Soil->TableDepth;
  case V_SOIL_SATFLOW: return Soil->SatFlow;
  case V_SOIL_WATERLEVEL: return Soil->WaterLevel;
  case V_SOIL_SURFACEWATER: return Soil->IExcess;
  case V_SOIL_INFILTRATIONACCUMULATED: return Soil->InfiltAcc;
  case V_SOIL_DESORPTIONACCUMULATED: return Soil->DesorptionAcc;
  case V_SOIL_WATER:
    Value = 0.0;
    DeepDepth = 0.0;
    for (i = 0; i < L; i++) {
      Value += Soil->Moist[i] * Soil->LayerThickness[i] * Network->Adjust[i];
      DeepDepth += Soil->LayerThickness[i];
    }
    Value += Soil->Moist[L] * (Soil->Depth - DeepDepth) * Network->Adjust[L];
    return Value;
  case V_SOIL_DEEPFLUX: return Soil->DeepFlux;
  case V_SOIL_LATENTHEAT: return Soil->Qe;
  case V_SOIL_ICE: return Soil->Ice[Layer];
  case V_SOIL_TEMPERATURE: return Soil->TLayer[Layer];
  case V_SOIL_NODETEMPERATURE: return Soil->TNode[Layer];
  case V_SOIL_NODEICE: return Soil->IceNode[Layer];
  case V_SOIL_SURFACETEMPERATURE: return Soil->TSurf;
  case V_SOIL_FROSTDEPTH: return Soil->FrostDepth;
  case V_SOIL_FROSTTOP: return (Soil->NFrost > 0) ? Soil->FrostTop[0] : 0.0;
  case V_SOIL_GROUNDHEATFLUX: return Soil->Qg;
  case V_SOIL_SENSIBLEHEATFLUX: return Soil->Qs;
  case V_SURFACE_DETENTION: return Soil->DetentionStorage;
  case V_SURFACE_RUNOFF: return Soil->Runoff;
  case V_STREAM_INTERCEPTION: return Soil->ChannelInt;
  case V_STREAM_INFILTRATION: return Soil->ChannelInfiltration;
  case V_LAKE_STORAGE: return (Model->Topo[y][x].LakeID > 0) ? Soil->IExcess : 0.0;
  case V_LAKE_LEVEL:
    if (Model->Topo[y][x].LakeID > 0)
      return Model->LakeLevel[Model->Topo[y][x].LakeID - 1];
    *Has = FALSE;
    return NA;
  case V_BALANCE_WATERINPUT:
    Value = (Precip->Precip + (Snow->OldSwq - Snow->Swq)) + Snow->VaporMassFlux;
    if (fabs(Value) <= 1.e-12)
      Value = 0.;
    return Value;
  default:
    *Has = FALSE;
    return NA;
  }
}

/*****************************************************************************
  SnowColdContent()

  Cold content of the snow pack (J/m2) from the cell's snow state, formed as
  SnowMelt forms it at the start of a step
*****************************************************************************/
float SnowColdContent(SNOWPIX *Snow)
{
  double Ice = (double) Snow->Swq - Snow->PackWater - Snow->SurfWater;
  double SurfaceSwq = (Ice > MAX_SURFACE_SWE) ? MAX_SURFACE_SWE : Ice;
  double PackSwq = Ice - SurfaceSwq;
  float SurfaceCC = CH_ICE * SurfaceSwq * Snow->TSurf;
  float PackCC = CH_ICE * PackSwq * Snow->TPack;

  return SurfaceCC + PackCC;
}

/*****************************************************************************
  SegmentValue()

  Value of a segment variable: volumes per step as rates (m3/s), the storage
  in m3
*****************************************************************************/
double SegmentValue(int Code, Channel *Segment, int Dt)
{
  float Change;

  switch (Code) {
  case V_CHANNEL_INFLOW: return (double) Segment->inflow / Dt;
  case V_CHANNEL_LATERALINFLOW: return (double) Segment->lateral_inflow / Dt;
  case V_CHANNEL_OUTFLOW: return (double) Segment->outflow / Dt;
  case V_CHANNEL_STORAGE: return (double) Segment->storage;
  case V_CHANNEL_STORAGECHANGE:
    Change = Segment->storage - Segment->last_storage;
    return (double) Change / Dt;
  case V_CHANNEL_INFILTRATION: return (double) Segment->infiltration / Dt;
  case V_CHANNEL_EVAPORATION: return (double) Segment->evaporation / Dt;
  default:
    return NA;
  }
}
