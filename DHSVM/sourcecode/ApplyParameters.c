/*
 * DESCRIPTION:  The [PARAMETERS] values applied to the class and lake tables
 *               as read (ApplyTableParameters), then to the input maps, once,
 *               at initialization, each by its one formula, in the order
 *               numbered here:
 *                 T1 DeepFlux = (soil table DeepFlux + offset) / (365 x steps
 *                   per day), every soil class
 *                 T2 RsMin = MinResistance x scale, every row of the
 *                   overstory and the understory table (JARVIS)
 *                 T3 SnowIntEff = min(overstory table SnowInterceptionEff
 *                   x scale, 1), every overstory class
 *                 T4 PowLawScale = lake table OutflowScale x scale, every
 *                   lake (Lake Dynamics)
 *                 T5 PowLawExponent = lake table OutflowExponent x scale,
 *                   every lake (Lake Dynamics)
 *                 T6 MoistThres = MoistureThreshold x scale (relative: 0 at
 *                   the wilting point, 1 at field capacity), every row of
 *                   the overstory and the understory table
 *                 T7 G1 = StomatalSlope x scale, every row of the overstory
 *                   and the understory table (MEDLYN)
 *                 1 Depth = max(Soil.Depth + offset, sum of
 *                   Soil.LayerThickness + 0.01)
 *                 2 KsLat = Soil.LateralConductivity x scale
 *                 3 KsLatExp = Soil.ExponentialDecrease x scale
 *                   then the vertical conductivity, the maximum infiltration
 *                   rate and the deep flux of each cell
 *                 4 FCap[l] = max(Soil.FieldCapacity[l] x scale, WP[l] + 0.01)
 *                 5 Porosity[l] = max(Soil.Porosity[l] x scale, FCap[l] + 0.01)
 *                 6 overstory cover = min(Veg.OverstoryCover x scale, 1),
 *                   understory cover = min(Veg.UnderstoryCover x scale, 1)
 *                 7 overstory LAI of a month = Veg.OverstoryLAI x the
 *                   overstory class fraction of the month x scale,
 *                   understory LAI of a month = Veg.UnderstoryLAI x the
 *                   understory class fraction of the month x scale
 *                 8-15 the snow parameters of every cell, with the forest
 *                   offset of the albedo decay rates
 *               and the melt scale of every cell from Snow.MeltScale (1
 *               without it).
 *               Only cells inside the basin mask are set. The input maps
 *               were checked before (Validate.c), so the bounds of rows 1,
 *               4, 5 and 6 act only through the parameters.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "Calendar.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

#define FCAP_SUCTION 3.3   /* Suction at field capacity (m) */
#define WP_SUCTION 150.0   /* Suction at the wilting point (m) */

static void ApplySoilParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, SOILPIX ***SoilMap, int Dt);
static void ApplyVegParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, VEGTABLES *Veg, VEGPIX ***VegMap);
static void ApplySnowParameters(MAPSIZE *Map, INPUTMAPS *Raw, TOPOPIX **TopoMap,
  VEGPIX **VegMap, SNOWPIX **SnowMap);

/*****************************************************************************
  ApplyTableParameters()

  Rows T1-T7, on the class and lake tables as read, before the maps: the
  deep flux offset of every soil class (m/yr, then m per time step), the
  minimum stomatal resistance of every overstory and understory class
  (JARVIS), the snow interception efficiency of every overstory class (at
  most 1), the outflow scale and exponent of every lake, the soil moisture
  threshold and the stomatal slope (MEDLYN) of every overstory and
  understory class
*****************************************************************************/
void ApplyTableParameters(OPTIONSTRUCT *Options, int NDaySteps, SOILTABLE *SType, int NSoils,
  VEGTABLES *Veg, LAKETABLE *LType, int NLakes)
{
  int i;

  /* T1 Deep Flux Offset, then m/yr to m per time step */
  for (i = 0; i < NSoils; i++) {
    SType[i].DeepFlux += SOIL_DEEPFLUX_ADJ;
    SType[i].DeepFlux /= (DAYPYEAR * NDaySteps);
  }

  for (i = 0; i < Veg->NOver; i++) {
    /* T2 Min Stomatal Resistance Scale (the maximum unchanged) */
    if (Options->StomatalScheme == JARVIS)
      Veg->Over[i].Story.RsMin *= VEG_RSMIN_ADJ;

    /* T3 Snow Interception Efficiency Scale: at most 1 */
    Veg->Over[i].SnowIntEff *= VEG_SNOWINT_ADJ;
    if (Veg->Over[i].SnowIntEff > 1.0)
      Veg->Over[i].SnowIntEff = 1.0;
  }
  for (i = 0; i < Veg->NUnder; i++) {
    /* T2 Min Stomatal Resistance Scale (the maximum unchanged) */
    if (Options->StomatalScheme == JARVIS)
      Veg->Under[i].Story.RsMin *= VEG_RSMIN_ADJ;
  }

  /* T4 Lake Outflow Rate Scale, T5 Lake Outflow Exponent Scale */
  if (Options->LakeDynamics) {
    for (i = 0; i < NLakes; i++) {
      LType[i].PowLawScale *= LAKE_OUTFLOW_ADJ;
      LType[i].PowLawExponent *= LAKE_EXPONENT_ADJ;
    }
  }

  /* T6 Moisture Threshold Scale, T7 Stomatal Slope Scale */
  for (i = 0; i < Veg->NOver; i++) {
    Veg->Over[i].Story.MoistThres *= VEG_MOIST_ADJ;
    if (Options->StomatalScheme == MEDLYN)
      Veg->Over[i].Story.G1 *= VEG_G1_ADJ;
  }
  for (i = 0; i < Veg->NUnder; i++) {
    Veg->Under[i].Story.MoistThres *= VEG_MOIST_ADJ;
    if (Options->StomatalScheme == MEDLYN)
      Veg->Under[i].Story.G1 *= VEG_G1_ADJ;
  }
}

/*****************************************************************************
  ApplyParameters()

  Sets up the soil and vegetation maps of every cell from the input maps and
  the [PARAMETERS] values, and the snow parameters of every cell. Frees the
  input maps it uses.
*****************************************************************************/
void ApplyParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, VEGTABLES *Veg,
  SOILPIX ***SoilMap, VEGPIX ***VegMap, SNOWPIX **SnowMap, int Dt)
{
  Message("\nApplying parameters\n");

  ApplySoilParameters(Options, Map, Raw, TopoMap, Soil, SType, SoilMap, Dt);
  ApplyVegParameters(Options, Map, Raw, TopoMap, Veg, VegMap);
  ApplySnowParameters(Map, Raw, TopoMap, *VegMap, SnowMap);
}

/*****************************************************************************
  ApplySoilParameters()

  Rows 1-5 and the values derived from them
*****************************************************************************/
static void ApplySoilParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, SOILPIX ***SoilMap, int Dt)
{
  const char *Routine = "ApplySoilParameters";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NCells = Map->NY * Map->NX;
  int NLayers = Map->SoilLayers;	/* root zone layers; the deep layer is index NLayers */
  int NSet, VSet;
  int sidx;
  SOILPIX *Cell;
  float LayerDepth, Transmissivity, KsVertCalc;
  double MaxSeepage;		/* Largest seepage loss per time step (negative) */

  if (!(*SoilMap = (SOILPIX **)calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *)calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      Cell = &((*SoilMap)[y][x]);
      Cell->Soil = Raw->SoilClass[i];
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      sidx = Cell->Soil - 1;

      /* allocate memory for the root zone layers and the deep layer */
      if (!(Cell->KsVert = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->FCap = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Porosity = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->PoreDist = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Press = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Moist = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Perc = (float *)calloc(NLayers, sizeof(float))) ||
          !(Cell->InterFlow = (float *)calloc(NLayers + 1, sizeof(float))) ||
          !(Cell->Temp = (float *)calloc(NLayers, sizeof(float))) ||
          !(Cell->LayerThickness = (float *)calloc(NLayers, sizeof(float))))
        ReportError((char *)Routine, 1);
      if (Options->Infiltration == DYNAMIC)
        Cell->InfiltAcc = 0.;
      Cell->MoistInit = 0.;

      /* Root zone layer thicknesses and their sum */
      Cell->TotalThickness = 0.0;
      for (NSet = 0; NSet < NLayers; NSet++) {
        Cell->LayerThickness[NSet] = Raw->LayerThickness[NSet * NCells + i];
        Cell->TotalThickness += Cell->LayerThickness[NSet];
      }

      /* 1 Soil Depth Offset: at least 0.01 m deeper than the root zone */
      Cell->Depth = Raw->Depth[i];
      Cell->Depth += SOIL_DEPTH_ADJ;
      if (Cell->Depth < (Cell->TotalThickness + 0.01))
        Cell->Depth = (Cell->TotalThickness + 0.01);

      /* 2 Lateral Conductivity Scale */
      Cell->KsLat = Raw->LateralConductivity[i];
      Cell->KsLat *= SOIL_KSAT_ADJ;

      /* 3 Exponential Decrease Scale */
      Cell->KsLatExp = Raw->ExponentialDecrease[i];
      Cell->KsLatExp *= SOIL_EXPDEC_ADJ;

      /* Vertical conductivity of each root zone layer: the transmissivity of
         the layer over its thickness and the anisotropy, or the soil table */
      for (NSet = 0; NSet < NLayers; NSet++) {
        if (Options->UseKsatAnisotropy) {
          /* Get total depth of current soil layer bottom */
          LayerDepth = 0.0;
          for (VSet = 0; VSet <= NSet; VSet++)
            LayerDepth += Cell->LayerThickness[VSet];

          /* Calculate effective lateral conductivity of current soil layer */
          Transmissivity =
            CalcTransmissivity(LayerDepth, LayerDepth - Cell->LayerThickness[NSet],
                               Cell->KsLat, Cell->KsLatExp, SType[sidx].DepthThresh);
          KsVertCalc = Transmissivity / Cell->LayerThickness[NSet];

          /* Account for user-supplied vertical anisotropy */
          KsVertCalc /= SType[sidx].KsAnisotropy;

          Cell->KsVert[NSet] = KsVertCalc;
        }
        else
          Cell->KsVert[NSet] = SType[sidx].Ks[NSet];
      }

      /* Vertical conductivity of the deep layer */
      NSet = NLayers;
      if (Options->UseKsatAnisotropy) {
        /* Transmissivity of the deep layer, from the bottom of the root zone
           to the soil depth */
        Transmissivity =
          CalcTransmissivity(Cell->Depth, Cell->TotalThickness,
                             Cell->KsLat, Cell->KsLatExp, SType[sidx].DepthThresh);
        KsVertCalc = Transmissivity / (Cell->Depth - Cell->TotalThickness);

        /* Account for user-supplied vertical anisotropy */
        KsVertCalc /= SType[sidx].KsAnisotropy;

        Cell->KsVert[NSet] = KsVertCalc;
      }
      else
        Cell->KsVert[NSet] = SType[sidx].Ks[NSet];

      /* Maximum infiltration rate: surface layer KsLat over the anisotropy,
         or the soil table */
      if (Options->UseKsatAnisotropy) {
        Cell->MaxInfiltrationRate = Cell->KsLat;
        Cell->MaxInfiltrationRate /= SType[sidx].KsAnisotropy;
      }
      else
        Cell->MaxInfiltrationRate = SType[sidx].MaxInfiltrationRate;

      /* Deep flux of the soil class; a seepage loss (negative deep flux) is
         limited to the deep layer's vertical conductivity over one time step */
      Cell->MaxDeepFlux = SType[sidx].DeepFlux;
      if (Cell->MaxDeepFlux < 0.0) {
        MaxSeepage = Cell->KsVert[NLayers] * -1.0;
        MaxSeepage *= Dt;
        if (Cell->MaxDeepFlux < MaxSeepage)
          Cell->MaxDeepFlux = MaxSeepage;
      }

      /* 4 Field Capacity Scale: at least 0.01 above the wilting point */
      for (NSet = 0; NSet <= NLayers; NSet++) {
        Cell->FCap[NSet] = Raw->FieldCapacity[NSet * NCells + i] * SOIL_FIELDCAP_ADJ;
        if (Cell->FCap[NSet] < (SType[sidx].WP[NSet] + 0.01))
          Cell->FCap[NSet] = SType[sidx].WP[NSet] + 0.01;
      }

      /* 5 Porosity Scale: at least 0.01 above the field capacity */
      for (NSet = 0; NSet <= NLayers; NSet++) {
        Cell->Porosity[NSet] = Raw->Porosity[NSet * NCells + i] * SOIL_POROSITY_ADJ;
        if (Cell->Porosity[NSet] < (Cell->FCap[NSet] + 0.01))
          Cell->Porosity[NSet] = Cell->FCap[NSet] + 0.01;
      }

      /* Brooks-Corey curve theta / Porosity = (Press / suction)^PoreDist
         through the field capacity at FCAP_SUCTION and the wilting point at
         WP_SUCTION */
      for (NSet = 0; NSet <= NLayers; NSet++) {
        Cell->PoreDist[NSet] = log(Cell->FCap[NSet] / SType[sidx].WP[NSet]) / log(WP_SUCTION / FCAP_SUCTION);
        Cell->Press[NSet] = FCAP_SUCTION * pow(Cell->FCap[NSet] / Cell->Porosity[NSet], 1.0 / Cell->PoreDist[NSet]);
      }
    }
  }

  free(Raw->SoilClass);
  free(Raw->Depth);
  free(Raw->LateralConductivity);
  free(Raw->ExponentialDecrease);
  free(Raw->FieldCapacity);
  free(Raw->Porosity);
  free(Raw->LayerThickness);
  Raw->SoilClass = NULL;
  Raw->Depth = Raw->LateralConductivity = Raw->ExponentialDecrease = NULL;
  Raw->FieldCapacity = Raw->Porosity = Raw->LayerThickness = NULL;
}

/*****************************************************************************
  ApplyVegParameters()

  Rows 6 and 7 on every mask cell from the eight vegetation maps: the
  classes, the overstory and understory cover and LAI, the heights and the
  wind profile and aerodynamic resistances
*****************************************************************************/
static void ApplyVegParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, VEGTABLES *Veg, VEGPIX ***VegMap)
{
  const char *Routine = "ApplyVegParameters";
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  VEGPIX *Cell;
  VEGMAPS *V = &(Raw->Veg);

  if (!(*VegMap = (VEGPIX **)calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *)calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      Cell = &((*VegMap)[y][x]);
      Cell->OverClass = V->OverstoryClass[i];
      Cell->UnderClass = V->UnderstoryClass[i];
      Cell->Tcanopy = 0.0;
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;

      /* Room for both canopies, whatever the classes (vegetation updates
         change them) */
      if (!(Cell->Fract = (float *)calloc(MAXVEGLAYERS, sizeof(float))) ||
          !(Cell->LAIFull = (float *)calloc(MAXVEGLAYERS, sizeof(float))) ||
          !(Cell->LAI = (float *)calloc(MAXVEGLAYERS, sizeof(float))) ||
          !(Cell->MaxInt = (float *)calloc(MAXVEGLAYERS, sizeof(float))) ||
          !(Cell->Height = (float *)calloc(MAXVEGLAYERS, sizeof(float))))
        ReportError((char *)Routine, 1);

      SetCellVegetation(Options, Veg, Cell, V, i, y, x);
    }
  }

  FreeVegMaps(V);
}

/*****************************************************************************
  SetCellVegetation()

  Rows 6 and 7 on one mask cell from the values of cell i of the vegetation
  maps V: the overstory and understory classes and the layers present, the
  cover of each layer, the full-leaf LAI of each layer, the
  heights, the wind profile and aerodynamic resistances (layers in position
  order, overstory first)
*****************************************************************************/
void SetCellVegetation(OPTIONSTRUCT *Options, VEGTABLES *Veg, VEGPIX *Cell, VEGMAPS *V,
                       size_t i, int y, int x)
{
  OVERTABLE *Over;
  int Lower;			/* Position of the understory layer */

  /* Classes and the layers present: an understory class with cover 0 has
     no understory layer (the cover before its scale) */
  Cell->OverClass = V->OverstoryClass[i];
  Cell->UnderClass = V->UnderstoryClass[i];
  Cell->OverStory = (Cell->OverClass > 0) ? TRUE : FALSE;
  Cell->UnderStory = (Cell->UnderClass > 0 && V->UnderstoryCover[i] > 0.0) ? TRUE : FALSE;
  Cell->NVegLayers = Cell->OverStory + Cell->UnderStory;
  Over = OverstoryRow(Veg, Cell);
  Lower = CanopyIndex(Cell, 1);

  /* 6 Overstory and Understory Cover Scale: each cover times its scale, at
     most 1 */
  if (Cell->OverStory == TRUE) {
    Cell->Fract[0] = V->OverstoryCover[i];
    Cell->Fract[0] *= VEG_FC_ADJ;
    /* Limit to maximum of 100% cover */
    if (Cell->Fract[0] > 1.0)
      Cell->Fract[0] = 1.0;
  }
  Cell->UnderstoryCover = 0.0;
  if (Cell->UnderStory == TRUE) {
    Cell->UnderstoryCover = V->UnderstoryCover[i];
    Cell->UnderstoryCover *= VEG_UNDER_FC_ADJ;
    if (Cell->UnderstoryCover > 1.0)
      Cell->UnderstoryCover = 1.0;
    Cell->Fract[Lower] = Cell->UnderstoryCover;
  }

  /* 7 Overstory and Understory LAI Scale: each LAI map times its scale, the
     full-leaf LAI (SetCellLAI forms the LAI of the day) */
  if (Cell->OverStory == TRUE) {
    Cell->LAIFull[0] = V->OverstoryLAI[i];
    Cell->LAIFull[0] *= VEG_LAI_ADJ;
  }
  if (Cell->UnderStory == TRUE) {
    Cell->LAIFull[Lower] = V->UnderstoryLAI[i];
    Cell->LAIFull[Lower] *= VEG_UNDER_LAI_ADJ;
  }

  /* Vegetation heights of the cell from the overstory and understory height
     maps */
  if (Cell->OverStory == TRUE)
    Cell->Height[0] = V->OverstoryHeight[i];
  if (Cell->UnderStory == TRUE)
    Cell->Height[Lower] = V->UnderstoryHeight[i];

  /* Wind profile and aerodynamic resistances of the cell's layers,
     normalized for a reference height wind speed of 1 m/s (scaled by the
     actual wind each time step); the overstory's wind attenuation and trunk
     space, the understory cover */
  if (Cell->OverStory == TRUE)
    CalcAerodynamic(y, x, Cell->NVegLayers, Cell->OverStory, Cell->UnderstoryCover,
                    Over->Cn, Cell->Height, Over->Trunk, Cell->U, &(Cell->USnow),
                    Cell->Ra, &(Cell->RaSnow));
  else
    CalcAerodynamic(y, x, Cell->NVegLayers, Cell->OverStory, Cell->UnderstoryCover,
                    NOT_APPLICABLE, Cell->Height, NOT_APPLICABLE, Cell->U, &(Cell->USnow),
                    Cell->Ra, &(Cell->RaSnow));
}

/*****************************************************************************
  GroundAlbedo()

  The albedo of the ground under the snowpack of a cell in a month: the
  understory albedo over the understory cover, the soil albedo over the
  rest (the soil albedo without an understory layer)
*****************************************************************************/
float GroundAlbedo(VEGTABLES *Veg, VEGPIX *Cell, float SoilAlbedo, int Month)
{
  float Fu = Cell->UnderstoryCover;

  if (Cell->UnderStory != TRUE)
    return SoilAlbedo;
  return Fu * UnderstoryRow(Veg, Cell)->Story.AlbedoMonthly[Month - 1] +
    (1 - Fu) * SoilAlbedo;
}

/*****************************************************************************
  OverstoryRow(), UnderstoryRow()

  The overstory and understory table rows of a cell's classes, NULL for
  class 0. A cell with an understory class and cover 0 has an understory row
  (its impervious fraction) but no understory layer.
*****************************************************************************/
OVERTABLE *OverstoryRow(VEGTABLES *Veg, VEGPIX *Cell)
{
  return (Cell->OverClass > 0) ? &(Veg->Over[Cell->OverClass - 1]) : NULL;
}

UNDERTABLE *UnderstoryRow(VEGTABLES *Veg, VEGPIX *Cell)
{
  return (Cell->UnderClass > 0) ? &(Veg->Under[Cell->UnderClass - 1]) : NULL;
}

/*****************************************************************************
  ImperviousFraction()

  The impervious fraction of a cell: its understory class's, 0 without an
  understory class
*****************************************************************************/
float ImperviousFraction(VEGTABLES *Veg, VEGPIX *Cell)
{
  return (Cell->UnderClass > 0) ? Veg->Under[Cell->UnderClass - 1].ImpervFrac : 0.0;
}

/*****************************************************************************
  ApplySnowParameters()

  Rows 8-15: the rain and snow thresholds, the fresh snow albedo, the albedo
  decay rates (with the forest offset on overstory cells, after row 6) and
  minima of every cell; the melt scale of every mask cell from its map, or 1
*****************************************************************************/
static void ApplySnowParameters(MAPSIZE *Map, INPUTMAPS *Raw, TOPOPIX **TopoMap,
  VEGPIX **VegMap, SNOWPIX **SnowMap)
{
  int x, y, i;

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN(TopoMap[y][x].Mask))
        SnowMap[y][x].MeltScale = (Raw->MeltScale != NULL) ? Raw->MeltScale[i] : 1.0;
      SnowMap[y][x].Tr = MIN_RAIN_TEMP;
      SnowMap[y][x].Ts = MAX_SNOW_TEMP;
      SnowMap[y][x].amax = ALB_MAX;
      SnowMap[y][x].LamdaAcc = ALB_ACC_LAMBDA;
      SnowMap[y][x].LamdaMelt = ALB_MELT_LAMBDA;
      SnowMap[y][x].AccMin = ALB_ACC_MIN;
      SnowMap[y][x].MeltMin = ALB_MELT_MIN;
    }
  }

  /* Adjust albedo melt/acc lambda in forested regions */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask) && VegMap[y][x].OverStory == TRUE) {
        SnowMap[y][x].LamdaAcc += (LAMBDA_FOREST_OFFSET * VegMap[y][x].Fract[0]);
        SnowMap[y][x].LamdaMelt += (LAMBDA_FOREST_OFFSET * VegMap[y][x].Fract[0]);
      }
    }
  }

  free(Raw->MeltScale);
  Raw->MeltScale = NULL;
}
