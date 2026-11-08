/* The Distributed Hydrology Soil Vegetation Model (DHSVM)
 * Custom version maintained by Eli Boardman | Mountain Hydrology LLC
 * https://mountainhydrology.com
 * Contact: eli.boardman@mountainhydrology.com
 * License: CC BY-NC-SA 4.0
*/

/******************************************************************************/
/* INCLUDES */
/******************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "fileio.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "channel.h"
#include "output.h"

/******************************************************************************/
/* GLOBAL VARIABLES */
/******************************************************************************/

char *version = "Version " DHSVM_VERSION;    /* store version string */
char commandline[BUFSIZE + 1] = "";	/* store command line */
char errorstr[BUFSIZ + 1] = "";			/* error message */

/******************************************************************************/
/* MAIN */
/******************************************************************************/

int main(int argc, char **argv) {
  float **PrismMap = NULL;
  float **SnowPatternMap = NULL;
  float **SnowPatternMapBase = NULL;
  int MaxStreamID;
  time_t WallStart;
  int i, x, y;
  int NStats;
  uchar ***MetWeights = NULL;

  TOTALS Totals;			/* Basin and class means, water balance */
  MODELVIEW Model;			/* The model state the outputs read */
  CHANNEL ChannelData = {NULL, NULL, NULL};
  CHANNELTABLES ChannelTables;
  DUMPSTRUCT Dump;
  EVAPPIX **EvapMap = NULL;
  INPUTFILES InFiles;
  INPUTMAPS Raw;
  LAYER Soil;
  VEGTABLES Veg;			/* Overstory and understory tables */
  LISTPTR Input = NULL;
  MAPSIZE Map;
  METLOCATION *Stat = NULL;
  OPTIONSTRUCT Options;
  PIXMET LocalMet;
  PRECIPPIX **PrecipMap = NULL;
  PIXRAD **RadiationMap = NULL;
  NETSTRUCT **Network	= NULL;
  SNOWPIX **SnowMap = NULL;
  SOILPIX **SoilMap = NULL;
  SOILTABLE *SType = NULL;
  SOLARGEOMETRY SolarGeo;
  SHADE Shade;
  double ShadeFactorCell = 0.0;
  float SkyViewCell = 1.0;
  TIMESTRUCT Time;
  TOPOPIX **TopoMap = NULL;
  LAKETABLE *LType = NULL;
  VEGPIX **VegMap = NULL;
    
/******************************************************************************/
/* INITIALIZATION */
/******************************************************************************/

  if (argc != 2)
    ReportErrorCode(101, "usage: DHSVM_%s <config file>", DHSVM_VERSION);

  snprintf(commandline, sizeof(commandline), "%s %s", argv[0], argv[1]);
  if (strlen(argv[1]) > BUFSIZE)
    ReportErrorCode(101, "%s: config file path longer than %d characters", argv[1], BUFSIZE);
  strcpy(InFiles.Const, argv[1]);

  WallStart = time(NULL);

  /* Config (every key), time steps, output requests and output directory */
  InitValueCodes();
  ReadInitFile(InFiles.Const, &Input);
  InitConstants(Input, &Options, &Map, &SolarGeo, &Time, &Dump);

  Message("%s \n", commandline);
  Message("\nRunning DHSVM %s\n", version);
  Message("\nMountain Hydrology Version Copyright Eli Boardman\n");
  Message("\nLICENSE: CC BY-NC-SA 4.0 (Non-Commercial Use Only!)\n");
#ifdef SNOW_ONLY
  Message("----------------------------------\n");
  Message("WARNING: USING SNOW ONLY MODULES (prescribed in makefile)!\n");
  Message("----------------------------------\n");
#endif
  Message("\nSTARTING INITIALIZATION PROCEDURES\n\n");
  PrintRunSettings(&Options);

  /* Input files: maps, tables (with their [PARAMETERS] rows T1-T5), the maps
     formed from the class tables, station data */
  CheckInputFiles(Input, &Options, OutputListNamesAny(&(Dump.List[OUT_POINT])));
  ReadInputMaps(Input, &Options, &Map, &Raw);
  InitTables(Time.NDaySteps, Input, &Options, &Map, &SType, &Soil, &Veg, &LType, &Time);
  FormClassMaps(&Options, &Map, &Soil, SType, &Veg, &Raw);
  ApplyTableParameters(&Options, Time.NDaySteps, SType, Soil.NTypes, &Veg, LType, Map.NumLakes);
  ReadChannelTables(Input, &ChannelTables);
  ReadPointTable(Input, &Dump);
  InitMetSources(Input, &Options, &Map, &Time, &NStats, &Stat);

  /* Vegetation updates checked; the latest before Start gives the
     vegetation maps */
  InitVegUpdates(&Options, &Map, &Time, &Raw, &Veg, &Dump);

  /* The input maps as given, then the parameters */
  ValidateInputs(Input, &Options, &Map, &Raw, &Soil, SType, &Veg, NStats, Stat);
  InitTopoMap(&Options, &Map, &Raw, &TopoMap, LType);
  InitSnowMap(&Map, &SnowMap, &Time);
  ApplyParameters(&Options, &Map, &Raw, TopoMap, &Soil, SType, &Veg, &SoilMap, &VegMap,
                  SnowMap, Time.Dt);
  ValidateParameters(&Options, &Map, TopoMap, SoilMap);
  PrintClassSummary(&Options, &Veg, Soil, SType, &Map, TopoMap, VegMap, SoilMap);

  /* Channels bound to the cells, storage adjustment */
  ValidateChannels(&Options, &Map, TopoMap, SoilMap, &ChannelTables, LType);
  InitChannel(&Map, Time.Dt, &ChannelData, &ChannelTables, SType, SoilMap,
              LType, TopoMap, &MaxStreamID, &Options);
  InitNetwork(Map.NY, Map.NX, Map.DX, Map.DY, TopoMap, SoilMap,
	      Map.SoilLayers, &Network, &ChannelData);

  /* Station weights, meteorological maps, outputs */
  InitMetMaps(&Map, &Options, &Raw,
	      &PrismMap, &SnowPatternMap, &SnowPatternMapBase,
	      &EvapMap, &PrecipMap, &RadiationMap, SoilMap, &Soil, VegMap, TopoMap);
  InitInterpolationWeights(&Map, &Options, TopoMap, &MetWeights, Stat, NStats);
  ValidateWeights(&Map, TopoMap, MetWeights, NStats);
  if (Options.SoilFreezing)
    InitSoilThermal(&Map, TopoMap, SoilMap, Stat, NStats, MetWeights);

  /* Impervious surface routing, horizon angles of the shading */
  if (Options.Impervious)
    InitSurfaceRouting(&Map, TopoMap, &ChannelData);
  if (Options.Shading)
    InitShading(&Options, &Map, TopoMap, &Raw, Time.Dt, &Shade);

  /* Output requests checked, output files created */
  ValidateRequests(Input, &Options, &Map, TopoMap, VegMap, &Veg, ChannelData.streams, &Dump);
  CreateOutputs(&Map, &Time, ChannelData.streams, &Dump);

  /* Done with initialization, warn about keys nothing read and delete the
     list with input strings */
  WarnUnreadKeys(Input);
  DeleteList(Input);

  InitModelState(&Time, &Map, &Options, PrecipMap, SnowMap, SoilMap,
		 Soil, SType, VegMap, &Veg, Dump.InitStatePath,
		 TopoMap, Network, &ChannelData, LType);
  InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, SnowPatternMap, SnowPatternMapBase,
	       &InFiles, &Veg, NStats, Stat, Dump.InitStatePath, &VegMap, SnowMap,
	       SType, SoilMap);
  InitNewDay(Time.Current.JDay, &SolarGeo);
  InitNewDayLAI(&Map, TopoMap, &Veg, VegMap, Time.Current.JDay);

  /* What the outputs read */
  Model.Map = &Map;
  Model.Topo = TopoMap;
  Model.Soil = SoilMap;
  Model.SType = SType;
  Model.Veg = VegMap;
  Model.Snow = SnowMap;
  Model.Precip = PrecipMap;
  Model.Evap = EvapMap;
  Model.Rad = RadiationMap;
  Model.Network = Network;
  Model.Streams = ChannelData.streams;
  Model.Dt = Time.Dt;
  Model.SoilFreezing = Options.SoilFreezing;
  Model.NumLakes = Map.NumLakes;
  Model.LakeLevel = NULL;
  Model.LakeCells = NULL;
  if (Model.NumLakes > 0 &&
      (!(Model.LakeLevel = (float *) calloc(Model.NumLakes, sizeof(float))) ||
       !(Model.LakeCells = (int *) calloc(Model.NumLakes, sizeof(int)))))
    ReportError("MainDHSVM", 1);

  /* Basin storage and means at Start, the series values at Start, the Start
     snapshots and map layers */
  InitTotals(&Map, TopoMap, VegMap, &Dump, &Totals);
  Aggregate(&Model, &Totals);
  StartBalance(&Totals);
#ifndef SNOW_ONLY
  StartNetworkStorage(&Model, Totals.Basin.NCells, &Totals);
#endif
  HoldStartValues(&Model, &Totals, &Dump);
  StoreSnapshots(0, &Map, &Dump, Options.SoilFreezing, TopoMap, PrecipMap, SnowMap,
                 VegMap, SoilMap, &ChannelData, LType);
  WriteMapLayers(0, &Model, &Dump);

/******************************************************************************/
/* CALCULATIONS */
/******************************************************************************/

  while (Time.Step < Time.NSteps) {

    ApplyVegUpdate(&Options, &Time, &Map, TopoMap, &Veg, VegMap, SType, SoilMap,
                   SnowMap, PrecipMap);
    if (IsNewWaterYear(&(Time.Current)))
      InitNewWaterYear(&Time, &Options, &Map, TopoMap, SnowMap, PrecipMap);
    if (IsNewMonth(&(Time.Current), Time.Dt))
      InitNewMonth(&Time, &Options, &Map, TopoMap, PrismMap, SnowPatternMap, SnowPatternMapBase,
		   &InFiles, &Veg, NStats, Stat, Dump.InitStatePath, &VegMap, SnowMap,
		   SType, SoilMap);
    if (IsNewDay(Time.DayStep)) {
      InitNewDay(Time.Current.JDay, &SolarGeo);
      InitNewDayLAI(&Map, TopoMap, &Veg, VegMap, Time.Current.JDay);
      Message("%04d-%02d-%02d\n", Time.Current.Year, Time.Current.Month, Time.Current.Day);
    }
    InitNewStep(&InFiles, &Map, &Time, &Options, NStats, Stat,
                &SolarGeo, TopoMap, SoilMap);
    if (Options.StomatalScheme == MEDLYN)
      ATMOS_CO2 = YearCO2(Time.Current.Year);
    channel_step_initialize_network(ChannelData.streams);
    if (Options.Shading)
      ShadeStep(&Shade, &SolarGeo, Time.DayStep, Time.Dt);

    for (y = 0; y < Map.NY; y++) {
      for (x = 0; x < Map.NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask)) {
          for (i = 0; i <= SType[SoilMap[y][x].Soil - 1].NLayers; i++)
            SoilMap[y][x].InterFlow[i] = 0.0;
          SoilMap[y][x].ChannelInt = 0.0;
          SoilMap[y][x].ChannelInfiltration = 0.0;
          SnowMap[y][x].OldSwq = SnowMap[y][x].Swq;
        }
      }
    }
    
    for (y = 0; y < Map.NY; y++) {
      for (x = 0; x < Map.NX; x++) {
  	    if (INBASIN(TopoMap[y][x].Mask)) {
  	      /* Beam factor of the terrain shading and the sky view of the cell */
  	      if (Options.Shading) {
  	        ShadeFactorCell = ShadeFactor(&Shade, &(TopoMap[y][x]), y, x);
  	        SkyViewCell = Shade.SkyView[y * Map.NX + x];
  	      }
  	      LocalMet =
  	        MakeLocalMetData(y, x, &Map, Time.DayStep, Time.NDaySteps, &Options, NStats,
                            Stat, MetWeights[y][x], TopoMap[y][x].Dem,
                            &(RadiationMap[y][x]), &(PrecipMap[y][x]),
                            PrismMap, SnowPatternMap, &(SnowMap[y][x]),
                            PRECIP_MULTIPLIER, Time.Current.Month, ShadeFactorCell,
                            SkyViewCell, TopoMap[y][x].WindScale);
  	      
  		    /* Get surface temperature of each soil layer */
    		  for (i = 0; i < Soil.MaxLayers; i++) {
    	        SoilMap[y][x].Temp[i] = LocalMet.Tair;
    		  }
    		  
    		  MassEnergyBalance(&Options, y, x, SolarGeo.SineSolarAltitude,
                          Map.DX, Map.DY, Time.Dt,
                          Options.Infiltration, Soil.MaxLayers,
                          MAXVEGLAYERS, &LocalMet, &(Network[y][x]),
                          &(PrecipMap[y][x]),
                          OverstoryRow(&Veg, &(VegMap[y][x])),
                          UnderstoryRow(&Veg, &(VegMap[y][x])), &(VegMap[y][x]),
                          &(SType[SoilMap[y][x].Soil - 1]), &(SoilMap[y][x]),
                          &(SnowMap[y][x]), &(RadiationMap[y][x]),
                          &(EvapMap[y][x]), &ChannelData,
                          VegMap, Network, SoilMap,
                          &(TopoMap[y][x]), TopoMap, &Map);
          
  		    PrecipMap[y][x].SumPrecip += PrecipMap[y][x].Precip;
  		    PrecipMap[y][x].SnowAccum += PrecipMap[y][x].SnowFall;
  		    PrecipMap[y][x].SnowMelt += SnowMap[y][x].Outflow;
		    }
	    }
    }
    
#ifndef SNOW_ONLY
    
    RouteSubSurface(Time.Dt, &Map, TopoMap, VegMap, Network,
		    SType, SoilMap, &ChannelData, &Time, &Options);

    RouteChannel(&ChannelData, &Time, &Map, TopoMap, SoilMap,
		 &Options, Network, SType, VegMap, EvapMap, LType);

    RouteSurface(&Map, &Time, TopoMap, SoilMap, &Options,
      &Dump, VegMap, &Veg, LType, SType, &ChannelData, LocalMet.Tair, LocalMet.Rh);

#endif

    /* Layer ice and temperatures and frost depth from the soil temperature
       nodes */
    if (Options.SoilFreezing)
      UpdateSoilThermalDerived(&Map, TopoMap, SoilMap, SType);

    /* Basin and class means (and the root zone moisture >= 0), snow
       statistics, water balance, outputs */
    Aggregate(&Model, &Totals);

    SnowStats(&(Time.Current), &Map, &Options, TopoMap, SnowMap, Time.Dt);

    MassBalance(&Totals);
#ifndef SNOW_ONLY
    NetworkSums(&Model, Totals.Basin.NCells, &Totals);
#endif

    ExecDump(&Map, &(Time.Current), Time.Step + 1, &Model, &Totals, &Dump, PrecipMap,
             SnowMap, &ChannelData, LType);

    IncreaseTime(&Time);
  } /* End of calculation loop over time steps */

  CloseOutputs(&Dump);
  Message("\nEND OF MODEL RUN\n");
#ifndef SNOW_ONLY
  FinalMassBalance(&Totals);
#endif

  /* Wall clock time of the run in whole seconds, also with Verbose = FALSE */
  fflush(stderr);
  printf("Wall clock time (s)              %.0f\n", difftime(time(NULL), WallStart));
  fflush(stdout);

  return EXIT_SUCCESS;
}
