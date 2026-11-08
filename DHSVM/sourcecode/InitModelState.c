
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "soilmoisture.h"
#include "channel_grid.h"
#include "Calendar.h"
#include "mapio.h"
#include "photosynthesis.h"
#include "snow.h"

 /*****************************************************************************
   Function name: InitModelState()

   Purpose      : Initialize the state of the model: a cold start
                  (Initial State = NONE) or the snapshot valid at Start
                  (StoreModelState(), ReadModelState())

   Comments :
     Cold start: no interception, the canopy temperature memories unset (the
     first step sets them to its air temperature), snow only where the
     initial snow water
     equivalent map has some (its liquid water at capacity), every soil
     layer at its porosity with the
     water table at the surface and no soil evaporation accumulated since a
     wetting, every channel segment filled to the initial
     channel depth. The groundwater spin-up, if any, follows the cold start.
     With Soil Freezing, every soil temperature node and the soil surface
     start at the cell's bottom temperature, with the node ice in
     equilibrium with it.
     Snapshot: the values as read; with Recompute Derived State, the states
     that follow from other states are then rebuilt from them
     (RecomputeDerivedState()).
     With Soil Freezing, after either: the layer temperatures, layer ice
     and frost depth from the node state (UpdateSoilThermalDerived()); at a
     cold start or with Recompute Derived State the table depth becomes the
     deeper of the water table and the frost depth.
 *****************************************************************************/
static void ColdStartCells(MAPSIZE *Map, OPTIONSTRUCT *Options,
                           PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                           SOILPIX **SoilMap, VEGPIX **VegMap,
                           TOPOPIX **TopoMap);
static void RecomputeDerivedState(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **SnowMap,
                                  SOILPIX **SoilMap, SOILTABLE *SType,
                                  NETSTRUCT **Network, CHANNEL *ChannelData);

void InitModelState(TIMESTRUCT *Time, MAPSIZE *Map, OPTIONSTRUCT *Options,
  PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, SOILPIX **SoilMap, LAYER Soil,
  SOILTABLE *SType, VEGPIX **VegMap, VEGTABLES *Veg, char *Path,
  TOPOPIX **TopoMap, NETSTRUCT **Network, CHANNEL *ChannelData, LAKETABLE *LType)
{
  const char *Routine = "InitModelState";
  int i, j;		         /* counter */
  int x;				 /* counter */
  int y;				 /* counter */

  /* Channel segments and cells */
  if (Options->ColdStart) {
    Message("Initializing the model state (cold start)\n");
#ifndef SNOW_ONLY
    ColdStartChannelState(Options->InitChannelDepth, ChannelData->streams);
#endif
    ColdStartCells(Map, Options, PrecipMap, SnowMap, SoilMap, VegMap, TopoMap);
  }
  else {
    ReadModelState(Path, &(Time->Start), Time->Dt, Map, Options->SoilFreezing, TopoMap,
                   PrecipMap, SnowMap, VegMap, SoilMap, ChannelData, LType);
    if (Options->RecomputeDerived)
      RecomputeDerivedState(Map, TopoMap, SnowMap, SoilMap, SType, Network, ChannelData);
  }

  /* Albedo of the ground or understory under the snow for the month of
     Start; at a cold start also the snow albedo */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        SnowMap[y][x].AlbedoGround = GroundAlbedo(Veg, &(VegMap[y][x]),
          SType[SoilMap[y][x].Soil - 1].Albedo, Time->Start.Month);

        if (Options->ColdStart) {
          if (SnowMap[y][x].HasSnow)
            SnowMap[y][x].Albedo = CalcSnowAlbedo(&(SnowMap[y][x]), 0.0);
          else
            SnowMap[y][x].Albedo = 0;
        }
      }
    }
  }

  if (Options->ColdStart) {
    /* Water table at the surface; the channel cells take the table depth of
       their cell */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        SoilMap[y][x].SatFlow = 0.0;
        SoilMap[y][x].TableDepth = 0.0;
        if (INBASIN(TopoMap[y][x].Mask)) {
          SoilMap[y][x].WaterLevel = TopoMap[y][x].Dem - SoilMap[y][x].TableDepth;
          SoilMap[y][x].WaterLevelLast = SoilMap[y][x].WaterLevel;
        }
        channel_grid_init_table(ChannelData->stream_map, x, y, SoilMap[y][x].TableDepth);
      }
    }

    /* Spinup groundwater state (optional, not included in mass balance) */
    if (Options->GW_SPINUP == TRUE) {

      Options->GW_SPINUP_RECHARGE /= DAYPYEAR; /* Convert m/yr to m/day */

      float **SubFlowGrad;	        /* Magnitude of subsurface flow gradient slope * width */
      unsigned char ***SubDir;      /* Fraction of flux moving in each direction*/
      unsigned int **SubTotalDir;	/* Sum of Dir array */

      if (!(SubFlowGrad = (float **)calloc((unsigned int) Map->NY, sizeof(float *))))
        ReportError((char *) Routine, 1);
      for(i=0; i<Map->NY; i++) {
        if (!(SubFlowGrad[i] = (float *)calloc(Map->NX, sizeof(float))))
          ReportError((char *) Routine, 1);
      }

      if (!((SubDir) = (unsigned char ***) calloc((unsigned int) Map->NY, sizeof(unsigned char **))))
        ReportError((char *) Routine, 1);
      for (i=0; i<Map->NY; i++) {
        if (!((SubDir)[i] = (unsigned char **) calloc(Map->NX, sizeof(unsigned char*))))
          ReportError((char *) Routine, 1);
        for (j=0; j<Map->NX; j++) {
          if (!(SubDir[i][j] = (unsigned char *)calloc(NDIRS, sizeof(unsigned char ))))
            ReportError((char *) Routine, 1);
        }
      }

      if (!(SubTotalDir = (unsigned int **)calloc((unsigned int) Map->NY, sizeof(unsigned int *))))
        ReportError((char *) Routine, 1);
      for (i=0; i<Map->NY; i++) {
        if (!(SubTotalDir[i] = (unsigned int *)calloc(Map->NX, sizeof(unsigned int))))
          ReportError((char *) Routine, 1);
      }

      /* Each spin-up day adds one day of deep flux */
      for (i = 0; i < Options->GW_SPINUP_YRS; i++) {
        Message("Groundwater spinup: year %d\n",i+1);
        for (j = 0; j < DAYPYEAR; j++) {
          RouteSubSurfaceSpinup(Time->Dt * Time->NDaySteps, Time->NDaySteps, Map, TopoMap,
                                VegMap, Network, SType, SoilMap, Options,
                                SubFlowGrad, SubDir, SubTotalDir);
        }
      }
      Message("Groundwater spinup complete\n\n");

      for (i = 0; i < Map->NY; i++) {
        for (j = 0; j < Map->NX; j++)
          free(SubDir[i][j]);
        free(SubDir[i]);
        free(SubFlowGrad[i]);
        free(SubTotalDir[i]);
      }
      free(SubDir);
      free(SubFlowGrad);
      free(SubTotalDir);

      /* The water level, and its time-averaging memory, from the table
         depth reached by the spin-up */
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          SoilMap[y][x].SatFlow = 0.0;
          if (INBASIN(TopoMap[y][x].Mask)) {
            SoilMap[y][x].WaterLevel = TopoMap[y][x].Dem - SoilMap[y][x].TableDepth;
            SoilMap[y][x].WaterLevelLast = SoilMap[y][x].WaterLevel;
            channel_grid_init_table(ChannelData->stream_map, x, y, SoilMap[y][x].TableDepth);
          }
        }
      }
    }
  }

  /* Soil temperature nodes: cold start, then the layer values and frost
     depth that follow from them */
  if (Options->SoilFreezing) {
    if (Options->ColdStart) {
      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN(TopoMap[y][x].Mask))
            ColdStartSoilThermal(&(SoilMap[y][x]), &(SType[SoilMap[y][x].Soil - 1]));
        }
      }
    }
    UpdateSoilThermalDerived(Map, TopoMap, SoilMap, SType);
  }

  /* Flood detention storage of impervious cells: empty at a cold start */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (Options->ColdStart)
        SoilMap[y][x].DetentionStorage = 0.0;
      SoilMap[y][x].DetentionIn = 0.0;
      SoilMap[y][x].DetentionOut = 0.0;
    }
  }
}

/*****************************************************************************
  ColdStartCells()

  Cold start of every mask cell (the soil at its porosity in every layer and
  freshly wetted, snow from the optional initial snow water equivalent map,
  the liquid water of both snow layers at capacity)
*****************************************************************************/
static void ColdStartCells(MAPSIZE *Map, OPTIONSTRUCT *Options,
                           PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                           SOILPIX **SoilMap, VEGPIX **VegMap,
                           TOPOPIX **TopoMap)
{
  MAPFILE File;
  float *SWE = NULL;
  int NVeg;
  int i, x, y;
  double Ice;                  /* Ice of the initial snow (m) */
  double SurfaceIce;           /* Ice of its surface layer (m) */

  if (!IsEmptyStr(Options->InitialSWEPath)) {
    if (!(SWE = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
      ReportError("ColdStartCells", 1);
    OpenInputMap(&File, Options->InitialSWEPath, "Snow.InitialSWE", NC_FLOAT, "m",
                 MAP_ONE, Map);
    ReadMapLayer(&File, 0, SWE, Map);
    CloseMap(&File);
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        /* Interception */
        NVeg = VegMap[y][x].NVegLayers;
        for (i = 0; i < NVeg; i++) {
          PrecipMap[y][x].IntRain[i] = 0.0;
          PrecipMap[y][x].IntSnow[i] = 0.0;
        }
        PrecipMap[y][x].TempIntStorage = 0.0;
        VegMap[y][x].Tcanopy = 0.0;
        VegMap[y][x].TAcclim = (float) PHOTO_UNSET;
        VegMap[y][x].TGrowth = (float) PHOTO_UNSET;

        /* Snow */
        SnowMap[y][x].Swq = (SWE != NULL) ? SWE[y * Map->NX + x] : 0.0;
        SnowMap[y][x].HasSnow = (SnowMap[y][x].Swq > 0.0) ? TRUE : FALSE;
        SnowMap[y][x].LastSnow = Options->InitDaysSinceSnowfall;
        /* Each layer holds the capacity times its ice, so the ice is
           Swq / (1 + capacity), the surface layer up to its maximum */
        Ice = SnowMap[y][x].Swq / (1.0 + LIQUID_WATER_CAPACITY);
        SurfaceIce = (Ice > MAX_SURFACE_SWE) ? MAX_SURFACE_SWE : Ice;
        SnowMap[y][x].PackWater = LIQUID_WATER_CAPACITY * (Ice - SurfaceIce);
        SnowMap[y][x].TPack = 0.0;
        SnowMap[y][x].SurfWater = LIQUID_WATER_CAPACITY * SurfaceIce;
        SnowMap[y][x].TSurf = 0.0;
        SnowMap[y][x].AccumSeason = FALSE;

        /* Soil: every layer, the deep one included, at its porosity */
        for (i = 0; i <= Map->SoilLayers; i++) {
          SoilMap[y][x].Moist[i] = SoilMap[y][x].Porosity[i];
          SoilMap[y][x].InterFlow[i] = 0.0;
        }
        SoilMap[y][x].IExcess = 0.0;
        SoilMap[y][x].DesorptionAcc = 0.0;
      }
    }
  }

  if (SWE != NULL)
    free(SWE);
}

/*****************************************************************************
  RecomputeDerivedState()

  After a snapshot is read, the states that are functions of other states
  are rebuilt from them, for every mask cell in row-major order:
  - snow cover where the snow water equivalent is > 0; where it is not,
    the pack and surface liquid water and temperatures are 0
    (as SnowMelt() leaves them after melt-out);
  - the water table depth from the soil moisture as UnsaturatedFlow() ends
    a time step (moisture above the porosity moves up; ponding goes to the
    surface water and the table depth to 0);
  - the water level (elevation minus table depth), its averaging memory
    (the water level) and the table depth of the cell's channel records.
  Then the top and bottom water depths of every channel segment from its
  storage, as the cold start computes them.
*****************************************************************************/
static void RecomputeDerivedState(MAPSIZE *Map, TOPOPIX **TopoMap, SNOWPIX **SnowMap,
                                  SOILPIX **SoilMap, SOILTABLE *SType,
                                  NETSTRUCT **Network, CHANNEL *ChannelData)
{
  SNOWPIX *Snow;
  SOILPIX *Soil;
  int x, y;
#ifndef SNOW_ONLY
  Channel *Segment;
#endif

  Message("Rebuilding the derived state from the snapshot\n");

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;

      Snow = &(SnowMap[y][x]);
      if (Snow->Swq > 0.0) {
        Snow->HasSnow = TRUE;
      }
      else {
        Snow->HasSnow = FALSE;
        Snow->PackWater = 0.0;
        Snow->SurfWater = 0.0;
        Snow->TPack = 0.0;
        Snow->TSurf = 0.0;
      }

      Soil = &(SoilMap[y][x]);
      Soil->TableDepth = WaterTableDepth(SType[Soil->Soil - 1].NLayers, Soil->Depth,
                                         Soil->LayerThickness, Soil->Porosity,
                                         Soil->FCap, Network[y][x].Adjust, Soil->Moist);
      if (Soil->TableDepth < 0.0) {
        Soil->IExcess += -(Soil->TableDepth);
        Soil->TableDepth = 0.0;
      }
      Soil->WaterLevel = TopoMap[y][x].Dem - Soil->TableDepth;
      Soil->WaterLevelLast = Soil->WaterLevel;
      channel_grid_init_table(ChannelData->stream_map, x, y, Soil->TableDepth);
    }
  }

#ifndef SNOW_ONLY
  for (Segment = ChannelData->streams; Segment != NULL; Segment = Segment->next) {
    Segment->top_water_depth = Segment->storage / (Segment->class2->width * Segment->length);
    Segment->bottom_water_depth = Segment->top_water_depth;
  }
#endif
}
