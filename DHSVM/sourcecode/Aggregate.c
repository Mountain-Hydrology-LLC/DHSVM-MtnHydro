/*
 * DESCRIPTION:  Basin and class means. One pass over the mask cells (north
 *               to south, west to east) adds each cell's values to running
 *               totals in float: the basin's and, if requested, those of
 *               the cell's overstory and understory class combination. Most variables add one value per cell;
 *               canopy storage, soil water and incoming shortwave are formed
 *               term by term. Means are the totals divided by the number of
 *               cells (cells with snow for the snow pack variables). The pass
 *               also sets the root zone moisture <= 0 to 0 before it is added
 *               to the moisture and soil water totals, and forms the level of
 *               each lake (the mean surface water of its cells).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "output.h"
#include "registry.h"
#include "photosynthesis.h"
#include "snow.h"

static int Base[NVALUECODES];	/* First slot of a value code with means, -1 if none */
static int IsSnowMean[NVALUECODES];	/* TRUE: mean over the cells with snow */
static int NLayers[NVALUECODES];	/* Slots of a value code */
static int SlotDiffuse;		/* Diffuse shortwave (incoming shortwave is beam + diffuse) */
static int NSlots;

static void InitAggregate(AGGREGATE *Aggregate);
static int CellClass(TOTALS *Totals, VEGPIX *Cell);
static void AddCell(AGGREGATE **Target, int NTarget, MODELVIEW *Model, int y, int x);
static void FinishAggregate(AGGREGATE *Aggregate, int Dt);

/*****************************************************************************
  InitTotals()

  Slots of the variables with basin and class means, the basin aggregate and
  one aggregate per class combination of [OUTPUT] Classes
*****************************************************************************/
void InitTotals(MAPSIZE *Map, TOPOPIX **TopoMap, VEGPIX **VegMap,
                DUMPSTRUCT *Dump, TOTALS *Totals)
{
  const char *Aggregation;
  int Code, Var, i;

  NSlots = 0;
  for (Code = 0; Code < NVALUECODES; Code++) {
    Var = CodeVar(Code);
    Aggregation = RegistryRow(Var)->Aggregation;
    Base[Code] = -1;
    NLayers[Code] = 0;
    IsSnowMean[Code] = (strcmp(Aggregation, "SnowMean") == 0);
    if (strcmp(Aggregation, "Mean") == 0 || IsSnowMean[Code]) {
      Base[Code] = NSlots;
      NLayers[Code] = RegistryLayerCount(Var, Map->SoilLayers);
      NSlots += NLayers[Code];
    }
  }
  SlotDiffuse = NSlots++;

  memset(Totals, 0, sizeof(TOTALS));
  InitAggregate(&(Totals->Basin));
  Totals->NClasses = Dump->NClasses;
  if (Totals->NClasses > 0) {
    if (!(Totals->Class = (AGGREGATE *) calloc(Totals->NClasses, sizeof(AGGREGATE))) ||
        !(Totals->ClassOver = (int *) calloc(Totals->NClasses, sizeof(int))) ||
        !(Totals->ClassUnder = (int *) calloc(Totals->NClasses, sizeof(int))))
      ReportError("InitTotals", 1);
    for (i = 0; i < Totals->NClasses; i++) {
      InitAggregate(&(Totals->Class[i]));
      Totals->ClassOver[i] = Dump->Classes[i].Over;
      Totals->ClassUnder[i] = Dump->Classes[i].Under;
    }
  }
}

static void InitAggregate(AGGREGATE *Aggregate)
{
  if (!(Aggregate->Sum = (float *) calloc(NSlots, sizeof(float))) ||
      !(Aggregate->Mean = (float *) calloc(NSlots, sizeof(float))))
    ReportError("InitAggregate", 1);
  Aggregate->NCells = 0;
  Aggregate->NSnow = 0;
  Aggregate->OldSwq = 0.0;
}

/*****************************************************************************
  Aggregate()

  Basin and class means of the current state and of the time step's fluxes,
  the fraction of mask cells whose water table is above 90% of the soil
  depth, and the lake levels
*****************************************************************************/
void Aggregate(MODELVIEW *Model, TOTALS *Totals)
{
  MAPSIZE *Map = Model->Map;
  AGGREGATE *Target[2];
  SOILPIX *Soil;
  float mgrid;
  int NTarget;
  int Class;
  int Lake;
  int x, y, i;

  for (i = -1; i < Totals->NClasses; i++) {
    AGGREGATE *A = (i < 0) ? &(Totals->Basin) : &(Totals->Class[i]);
    memset(A->Sum, 0, NSlots * sizeof(float));
    A->NCells = 0;
    A->NSnow = 0;
  }
  Totals->NSaturated = 0;
  for (i = 0; i < Model->NumLakes; i++) {
    Model->LakeLevel[i] = 0.0;
    Model->LakeCells[i] = 0;
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(Model->Topo[y][x].Mask)) {
        NTarget = 0;
        Target[NTarget++] = &(Totals->Basin);
        Class = CellClass(Totals, &(Model->Veg[y][x]));
        if (Class >= 0)
          Target[NTarget++] = &(Totals->Class[Class]);
        AddCell(Target, NTarget, Model, y, x);

        /* Water table above MTHRESH of the soil depth */
        Soil = &(Model->Soil[y][x]);
        mgrid = (Soil->Depth - Soil->TableDepth) / Soil->Depth;
        if (mgrid > MTHRESH)
          Totals->NSaturated++;

        Lake = Model->Topo[y][x].LakeID;
        if (Lake > 0) {
          Model->LakeLevel[Lake - 1] += Soil->IExcess;
          Model->LakeCells[Lake - 1]++;
        }
      }
    }
  }
  for (i = 0; i < Model->NumLakes; i++) {
    if (Model->LakeCells[i] > 0)
      Model->LakeLevel[i] /= Model->LakeCells[i];
  }

  FinishAggregate(&(Totals->Basin), Model->Dt);
  for (i = 0; i < Totals->NClasses; i++)
    FinishAggregate(&(Totals->Class[i]), Model->Dt);
  Totals->SaturatedFraction = (float) Totals->NSaturated / (float) Totals->Basin.NCells;
}

/* The class combination of a cell's overstory and understory classes, -1
   if it has no means */
static int CellClass(TOTALS *Totals, VEGPIX *Cell)
{
  int i;

  for (i = 0; i < Totals->NClasses; i++) {
    if (Totals->ClassOver[i] == Cell->OverClass && Totals->ClassUnder[i] == Cell->UnderClass)
      return i;
  }
  return -1;
}

/* Adds Value to slot Layer of value code Code of every target */
#define ADD(Code, Layer, Value) \
  for (t = 0; t < NTarget; t++) Target[t]->Sum[Base[Code] + (Layer)] += (Value)

/*****************************************************************************
  AddCell()

  Adds the values of mask cell (y, x) to the totals of the targets; the root
  zone moisture is set to 0 where it is <= 0 before it is added
*****************************************************************************/
static void AddCell(AGGREGATE **Target, int NTarget, MODELVIEW *Model, int y, int x)
{
  EVAPPIX *Evap = &(Model->Evap[y][x]);
  PRECIPPIX *Precip = &(Model->Precip[y][x]);
  PIXRAD *Rad = &(Model->Rad[y][x]);
  SNOWPIX *Snow = &(Model->Snow[y][x]);
  SOILPIX *Soil = &(Model->Soil[y][x]);
  VEGPIX *Veg = &(Model->Veg[y][x]);
  NETSTRUCT *Network = &(Model->Network[y][x]);
  int L = Model->Map->SoilLayers;
  int NVeg = Veg->NVegLayers;
  int Index[2];
  int c, i, t;
  float Term;
  float DeepDepth;

  for (t = 0; t < NTarget; t++)
    Target[t]->NCells++;

  Index[0] = CanopyIndex(Veg, 0);
  Index[1] = CanopyIndex(Veg, 1);

  /* Meteorology and radiation */
  ADD(V_MET_PRECIPITATION, 0, Precip->Precip);
  ADD(V_MET_SNOWFALL, 0, Precip->SnowFall);
  ADD(V_MET_RAINFALL, 0, Precip->RainFall);
  ADD(V_MET_AIRTEMPERATURE, 0, Rad->Tair);
  ADD(V_MET_WINDSPEED, 0, Rad->Wind);
  ADD(V_MET_RELATIVEHUMIDITY, 0, Rad->Rh);
  ADD(V_MET_PRECIPITATIONWATERYEAR, 0, Precip->SumPrecip);
  ADD(V_RAD_SHORTWAVEINOBSERVED, 0, Rad->ObsShortIn);
  ADD(V_RAD_SHORTWAVEIN, 0, Rad->BeamIn);
  for (t = 0; t < NTarget; t++)
    Target[t]->Sum[SlotDiffuse] += Rad->DiffuseIn;
  ADD(V_RAD_NETSHORTWAVE, 0, Rad->PixelNetShort);
  ADD(V_RAD_LONGWAVEIN, 0, Rad->PixelLongIn);
  Term = Rad->NetRadiation[0] + Rad->NetRadiation[1];
  ADD(V_RAD_NETRADIATION, 0, Term);

  /* Vegetation layers, by canopy */
  for (c = 0; c < 2; c++) {
    if (Index[c] < 0)
      continue;
    ADD(V_RAD_CANOPYNETSHORTWAVE, c, Rad->NetShort[Index[c]]);
    ADD(V_RAD_CANOPYLONGWAVEIN, c, Rad->LongIn[Index[c]]);
    ADD(V_CANOPY_INTERCEPTEDRAIN, c, Precip->IntRain[Index[c]]);
    ADD(V_CANOPY_INTERCEPTEDSNOW, c, Precip->IntSnow[Index[c]]);
    ADD(V_EVAP_POTENTIAL, c, Evap->EPot[Index[c]]);
    ADD(V_EVAP_ACTUAL, c, Evap->EAct[Index[c]]);
    ADD(V_EVAP_INTERCEPTION, c, Evap->EInt[Index[c]]);
    ADD(V_CANOPY_ASSIMILATION, c, Evap->An[Index[c]]);
    for (i = 0; i < L; i++)
      ADD(V_EVAP_TRANSPIRATION, c * L + i, Evap->ESoil[Index[c]][i]);
  }
  ADD(V_EVAP_POTENTIAL, 2, Evap->EPot[NVeg]);
  ADD(V_EVAP_ACTUAL, 2, Evap->EAct[NVeg]);

  /* Canopy storage, layer by layer in position order */
  for (i = 0; i < NVeg; i++) {
    Term = Precip->IntRain[i] + Precip->IntSnow[i];
    ADD(V_CANOPY_STORAGE, 0, Term);
  }
  ADD(V_CANOPY_TEMPORARYSTORAGE, 0, Precip->TempIntStorage);
  ADD(V_CANOPY_SNOWVAPOR, 0, Snow->CanopyVaporMassFlux);
  Term = PhotoAcclimation(Veg->TAcclim);
  ADD(V_CANOPY_ACCLIMATION, 0, Term);

  /* Evapotranspiration */
  ADD(V_EVAP_TOTAL, 0, Evap->ETot);
  ADD(V_EVAP_SOIL, 0, Evap->EvapSoil);
  ADD(V_EVAP_CHANNEL, 0, Evap->EvapChannel);

  /* Snow: the pack variables over the cells with snow */
  if (Snow->HasSnow) {
    for (t = 0; t < NTarget; t++)
      Target[t]->NSnow++;
    ADD(V_SNOW_ALBEDO, 0, Snow->Albedo);
    ADD(V_SNOW_DEPTH, 0, SnowDepth(Snow->Swq));
    ADD(V_SNOW_DAYSSINCESNOWFALL, 0, Snow->LastSnow);
    ADD(V_SNOW_PACKLIQUID, 0, Snow->PackWater);
    ADD(V_SNOW_PACKTEMPERATURE, 0, Snow->TPack);
    ADD(V_SNOW_SURFACELIQUID, 0, Snow->SurfWater);
    ADD(V_SNOW_SURFACETEMPERATURE, 0, Snow->TSurf);
    Term = SnowColdContent(Snow);
    ADD(V_SNOW_COLDCONTENT, 0, Term);
    ADD(V_SNOW_LATENTHEAT, 0, Snow->Qe);
    ADD(V_SNOW_SENSIBLEHEAT, 0, Snow->Qs);
    ADD(V_SNOW_NETSHORTWAVE, 0, Snow->Qsw);
    ADD(V_SNOW_NETLONGWAVE, 0, Snow->Qlw);
    ADD(V_SNOW_RAINHEAT, 0, Snow->Qp);
    ADD(V_SNOW_MELTENERGY, 0, Snow->MeltEnergy);
  }
  Term = (float) Snow->HasSnow;
  ADD(V_SNOW_COVERED, 0, Term);
  ADD(V_SNOW_SWE, 0, Snow->Swq);
  ADD(V_SNOW_OUTFLOW, 0, Snow->Outflow);
  ADD(V_SNOW_VAPOR, 0, Snow->VaporMassFlux);
  ADD(V_SNOW_ACCUMULATIONWATERYEAR, 0, Precip->SnowAccum);
  ADD(V_SNOW_OUTFLOWWATERYEAR, 0, Precip->SnowMelt);

  /* Soil: root zone layers (moisture set >= 0, then added), then the deep
     layer below the root zone */
  DeepDepth = 0.0;
  for (i = 0; i < L; i++) {
    if (Soil->Moist[i] <= 0.0)
      Soil->Moist[i] = 0.0;
    ADD(V_SOIL_MOISTURE, i, Soil->Moist[i]);
    ADD(V_SOIL_INTERFLOW, i, Soil->InterFlow[i]);
    ADD(V_SOIL_PERCOLATION, i, Soil->Perc[i]);
    Term = Soil->Moist[i] * Soil->LayerThickness[i] * Network->Adjust[i];
    ADD(V_SOIL_WATER, 0, Term);
    DeepDepth += Soil->LayerThickness[i];
  }
  ADD(V_SOIL_MOISTURE, L, Soil->Moist[L]);
  for (i = 0; i <= L; i++) {
    ADD(V_SOIL_ICE, i, Soil->Ice[i]);
    ADD(V_SOIL_TEMPERATURE, i, Soil->TLayer[i]);
  }
  ADD(V_SOIL_SURFACETEMPERATURE, 0, Soil->TSurf);
  ADD(V_SOIL_FROSTDEPTH, 0, Soil->FrostDepth);
  ADD(V_SOIL_FROSTTOP, 0, (Soil->NFrost > 0) ? Soil->FrostTop[0] : 0.0);
  ADD(V_SOIL_GROUNDHEATFLUX, 0, Soil->Qg);
  ADD(V_SOIL_SENSIBLEHEATFLUX, 0, Soil->Qs);
  Term = Soil->Moist[L] * (Soil->Depth - DeepDepth) * Network->Adjust[L];
  ADD(V_SOIL_WATER, 0, Term);
  ADD(V_SOIL_TABLEDEPTH, 0, Soil->TableDepth);
  ADD(V_SOIL_SATFLOW, 0, Soil->SatFlow);
  ADD(V_SOIL_DEEPFLUX, 0, Soil->DeepFlux);
  ADD(V_SOIL_LATENTHEAT, 0, Soil->Qe);
  ADD(V_SOIL_SURFACEWATER, 0, Soil->IExcess);
  Term = (Model->Topo[y][x].LakeID > 0) ? Soil->IExcess : 0.0;
  ADD(V_LAKE_STORAGE, 0, Term);
  ADD(V_SURFACE_DETENTION, 0, Soil->DetentionStorage);
  ADD(V_SOIL_INFILTRATIONACCUMULATED, 0, Soil->InfiltAcc);
  ADD(V_STREAM_INTERCEPTION, 0, Soil->ChannelInt);
  ADD(V_STREAM_INFILTRATION, 0, Soil->ChannelInfiltration);
}

/*****************************************************************************
  FinishAggregate()

  Means from the totals: divided by the number of cells, or by the number of
  cells with snow (at least 1) for the snow pack variables; potential
  evapotranspiration from m/s to m per step after the division; incoming
  shortwave as the mean beam plus the mean diffuse
*****************************************************************************/
static void FinishAggregate(AGGREGATE *Aggregate, int Dt)
{
  float Beam, Diffuse;
  int NSnow = (Aggregate->NSnow < 1) ? 1 : Aggregate->NSnow;
  int Code, k, s;

  /* A class without cells (vegetation updates): no means (NA) */
  if (Aggregate->NCells == 0)
    return;

  for (Code = 0; Code < NVALUECODES; Code++) {
    for (k = 0; k < NLayers[Code]; k++) {
      s = Base[Code] + k;
      if (IsSnowMean[Code])
        Aggregate->Mean[s] = Aggregate->Sum[s] / NSnow;
      else
        Aggregate->Mean[s] = Aggregate->Sum[s] / Aggregate->NCells;
    }
  }
  for (k = 0; k < NLayers[V_EVAP_POTENTIAL]; k++)
    Aggregate->Mean[Base[V_EVAP_POTENTIAL] + k] *= Dt;
  Beam = Aggregate->Sum[Base[V_RAD_SHORTWAVEIN]] / Aggregate->NCells;
  Diffuse = Aggregate->Sum[SlotDiffuse] / Aggregate->NCells;
  Aggregate->Mean[Base[V_RAD_SHORTWAVEIN]] = Beam + Diffuse;
}

/*****************************************************************************
  AggregateValue(), SetAggregateValue()

  Mean of layer Layer of a value code with means
*****************************************************************************/
float AggregateValue(AGGREGATE *Aggregate, int Code, int Layer)
{
  return Aggregate->Mean[Base[Code] + Layer];
}

void SetAggregateValue(AGGREGATE *Aggregate, int Code, int Layer, float Value)
{
  Aggregate->Mean[Base[Code] + Layer] = Value;
}
