
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "data.h"
#include "channel.h"
#include "DHSVMChannel.h"

void CalcAerodynamic(int y, int x, int NVegLayers, unsigned char OverStory,
		     float UnderCover, float n, float *Height, float Trunk, float *U,
		     float *U2mSnow, float *Ra, float *RaSnow);

double CalcDistance(COORD *LocA, COORD *LocB);

float CalcSnowAlbedo(SNOWPIX *LocalSnow, float NewSnow);

float CalcTransmissivity(float SoilDepth, float WaterTable, float LateralKs,
			 float KsExponent, float DepthThresh);

void CalcWeights(METLOCATION *Station, int NStats, int NX, int NY,
		 uchar **BasinMask, uchar ****WeightArray,
		 OPTIONSTRUCT *Options);

void ApplyTableParameters(OPTIONSTRUCT *Options, int NDaySteps, SOILTABLE *SType, int NSoils,
  VEGTABLES *Veg, LAKETABLE *LType, int NLakes);

void ApplyParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX **TopoMap, LAYER *Soil, SOILTABLE *SType, VEGTABLES *Veg,
  SOILPIX ***SoilMap, VEGPIX ***VegMap, SNOWPIX **SnowMap, int Dt);

void SetCellVegetation(OPTIONSTRUCT *Options, VEGTABLES *Veg, VEGPIX *Cell, VEGMAPS *V,
  size_t i, int y, int x);

OVERTABLE *OverstoryRow(VEGTABLES *Veg, VEGPIX *Cell);

UNDERTABLE *UnderstoryRow(VEGTABLES *Veg, VEGPIX *Cell);

float ImperviousFraction(VEGTABLES *Veg, VEGPIX *Cell);

float GroundAlbedo(VEGTABLES *Veg, VEGPIX *Cell, float SoilAlbedo, int Month);

void SetCellLAI(VEGTABLES *Veg, VEGPIX *Cell, int Day);
float StoryLAIFraction(STORYTABLE *Story, int Day);
void InitNewDayLAI(MAPSIZE *Map, TOPOPIX **TopoMap, VEGTABLES *Veg, VEGPIX **VegMap, int Day);

void InitVegUpdateTimes(LISTPTR Input, OPTIONSTRUCT *Options, TIMESTRUCT *Time);

void InitVegUpdates(OPTIONSTRUCT *Options, MAPSIZE *Map, TIMESTRUCT *Time, INPUTMAPS *Raw,
  VEGTABLES *Veg, DUMPSTRUCT *Dump);

void ApplyVegUpdate(OPTIONSTRUCT *Options, TIMESTRUCT *Time, MAPSIZE *Map, TOPOPIX **TopoMap,
  VEGTABLES *Veg, VEGPIX **VegMap, SOILTABLE *SType, SOILPIX **SoilMap,
  SNOWPIX **SnowMap, PRECIPPIX **PrecipMap);

void FreeVegMaps(VEGMAPS *V);

float **VegMapValues(VEGMAPS *V, int m);

void ValidateVegMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, unsigned char *Mask, VEGMAPS *V,
  VEGTABLES *Veg, const char *Suffix);

void CheckInputFiles(LISTPTR Input, OPTIONSTRUCT *Options, int PointTable);

unsigned char dequal(double a, double b);

void ExecDump(MAPSIZE *Map, DATE *Current, int StepsDone, MODELVIEW *Model, TOTALS *Totals,
              DUMPSTRUCT *Dump, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, CHANNEL *ChannelData,
              LAKETABLE *LType);

unsigned char fequal(float a, float b);

float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
             TOPOPIX **TopoMap, SOILTABLE *SType);

void GetMetData(OPTIONSTRUCT *Options, TIMESTRUCT *Time,
		int NStats, float SunMax, METLOCATION *Stat);

void InitCharArray(char *Array, int Size);

void InitConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
	SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time, DUMPSTRUCT *Dump);

void PrintRunSettings(OPTIONSTRUCT *Options);

void InitEvapMap(MAPSIZE *Map, EVAPPIX ***EvapMap, SOILPIX **SoilMap,
		 LAYER *Soil, VEGPIX **VegMap, TOPOPIX **TopoMap);

void InitInterpolationWeights(MAPSIZE *Map, OPTIONSTRUCT *Options,
			      TOPOPIX **TopoMap, uchar ****MetWeights,
			      METLOCATION *Stats, int NStats);

void InitMetMaps(MAPSIZE *Map, OPTIONSTRUCT *Options, INPUTMAPS *Raw,
                 float ***PrismMap, float ***SnowPatternMap, float ***SnowPatternMapBase,
                 EVAPPIX ***EvapMap, PRECIPPIX ***PrecipMap, PIXRAD ***RadMap,
                 SOILPIX **SoilMap, LAYER *Soil, VEGPIX **VegMap,
                 TOPOPIX **TopoMap);

void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
            TIMESTRUCT *Time, int *NStats, METLOCATION **Stat);

void InitModelState(TIMESTRUCT *Time, MAPSIZE *Map, OPTIONSTRUCT *Options,
		    PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
		    SOILPIX **SoilMap, LAYER Soil, SOILTABLE *SType,
		    VEGPIX **VegMap, VEGTABLES *Veg, char *Path,
		    TOPOPIX **TopoMap,
		    NETSTRUCT **Network, CHANNEL *ChannelData, LAKETABLE *LType);

void InitNetwork(int NY, int NX, float DX, float DY, TOPOPIX **TopoMap,
		 SOILPIX **SoilMap, int SoilLayers,
		 NETSTRUCT ***Network, CHANNEL *ChannelData);

void InitNewDay(int DayOfYear, SOLARGEOMETRY *SolarGeo);

void InitNewMonth(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
		  TOPOPIX **TopoMap, float **PrismMap, float **SnowPatternMap, float **SnowPatternMapBase,
		  INPUTFILES *InFiles, VEGTABLES *Veg, int NStats,
		  METLOCATION *Stat, char *Path, VEGPIX ***VegMap, SNOWPIX **SnowMap,
		  SOILTABLE *SType, SOILPIX **SoilMap);

void InitNewStep(INPUTFILES *InFiles, MAPSIZE *Map, TIMESTRUCT *Time,
		 OPTIONSTRUCT *Options, int NStats,
		 METLOCATION *Stat, SOLARGEOMETRY *SolarGeo, 
		 TOPOPIX **TopoMap, SOILPIX **SoilMap);

void InitNewWaterYear(TIMESTRUCT *Time, OPTIONSTRUCT *Options, MAPSIZE *Map,
                TOPOPIX **TopoMap, SNOWPIX **SnowMap, PRECIPPIX **PrecipMap);

void InitPrismMap(int NY, int NX, float ***PrismMap);

void InitSnowPatternMap(INPUTMAPS *Raw, float ***SnowPatternMap, float ***SnowPatternMapBase,
                        MAPSIZE *Map);

void InitShading(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap, INPUTMAPS *Raw,
                 int Dt, SHADE *Shade);

void InitSurfaceRouting(MAPSIZE *Map, TOPOPIX **TopoMap, CHANNEL *ChannelData);

void ShadeStep(SHADE *Shade, SOLARGEOMETRY *SolarGeo, int DayStep, int Dt);

double ShadeFactor(SHADE *Shade, TOPOPIX *Topo, int y, int x);

void InitPrecipMap(MAPSIZE *Map, PRECIPPIX ***PrecipMap, VEGPIX **VegMap,
		   TOPOPIX **TopoMap);

void InitRadMap(MAPSIZE *Map, PIXRAD ***RadMap);

void InitSatVaporTable(void);

void InitSnowMap(MAPSIZE *Map, SNOWPIX ***SnowMap, TIMESTRUCT *Time);

int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE **SType, const char *Path,
                  LAYER *Soil, int SoilLayers, TIMESTRUCT *Time);

void InitStateDump(LISTPTR Input, TIMESTRUCT *Time, DUMPSTRUCT *Dump);

float YearCO2(int Year);
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLES *Veg,
  LAKETABLE **LType, TIMESTRUCT *Time);

void InitTopoMap(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX ***TopoMap, LAKETABLE *LType);

int InitLakeTable(LAKETABLE **LType, const char *Path, OPTIONSTRUCT *Options);

float evalexpint(int n, float x);

uchar IsStationLocation(COORD *Loc, int NStats, METLOCATION *Station,
			int *WhichStation);

float LapseT(float Temp, float FromElev, float ToElev, float LapseRate);
float AirPressure(float Elev, float LapseRate);

PIXMET MakeLocalMetData(int y, int x, MAPSIZE *Map, int DayStep, int NDaySteps,
			OPTIONSTRUCT *Options, int NStats, METLOCATION *Stat,
      uchar *MetWeights, float LocalElev, PIXRAD *RadMap,
			PRECIPPIX *PrecipMap,
			float **PrismMap, float **SnowPatternMap, SNOWPIX *LocalSnow,
			float precipMultiplier, int Month, double ShadeFactor, float SkyView,
			float WindScale);

void MassEnergyBalance(OPTIONSTRUCT *Options, int y, int x,
                       float SineSolarAltitude, float DX, float DY,
                       int Dt, int InfiltOption, int MaxSoilLayers, int MaxVegLayers, PIXMET *LocalMet,
                       NETSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip,
                       OVERTABLE *Over, UNDERTABLE *Under, VEGPIX *LocalVeg, SOILTABLE *SType,
                       SOILPIX *LocalSoil, SNOWPIX *LocalSnow, PIXRAD *LocalRad,
                       EVAPPIX *LocalEvap, CHANNEL *ChannelData,
                       VEGPIX **VegMap, NETSTRUCT **Network, SOILPIX **SoilMap,
                       TOPOPIX *LocalTopo, TOPOPIX **TopoMap, MAPSIZE *Map);

double pow (double a, double b);

void quick(ITEM *OrderedCells, int count);

void qs(ITEM *OrderedCells, int left, int right);

void ColdStartChannelState(float Depth, Channel *Head);

void ReadChannelState(char *FileName, int deltat, Channel *Head);

void PrintClassSummary(OPTIONSTRUCT *Options, VEGTABLES *Veg, LAYER Soil,
  SOILTABLE *SType, MAPSIZE *Map, TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap);

void ReadInputMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw);

void FormClassMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, LAYER *Soil, SOILTABLE *SType,
  VEGTABLES *Veg, INPUTMAPS *Raw);

void FormVegMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, VEGTABLES *Veg, VEGMAPS *V);

void ReadMetRecord(DATE *Current, METLOCATION *Station);

int Round(double x);

void RouteSubSurface(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
		     VEGPIX **VegMap,
		     NETSTRUCT **Network, SOILTABLE *SType,
		     SOILPIX **SoilMap, CHANNEL *ChannelData, 
		     TIMESTRUCT *Time, OPTIONSTRUCT *Options);

void RouteSubSurfaceSpinup(int Dt, int DeepFluxSteps, MAPSIZE *Map, TOPOPIX **TopoMap,
                           VEGPIX **VegMap,
                           NETSTRUCT **Network, SOILTABLE *SType,
                           SOILPIX **SoilMap, OPTIONSTRUCT *Options,
                           float **SubFlowGrad, unsigned char ***SubDir, unsigned int **SubTotalDir);

void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLES *Veg, LAKETABLE *LType, SOILTABLE *SType, CHANNEL *ChannelData,
  float Tair, float Rh);

float SatVaporPressure(float Temperature);

void StoreChannelState(char *FileName, Channel *Head);

void StoreModelState(char *Path, DATE *Valid, MAPSIZE *Map, TOPOPIX **TopoMap,
                     PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, VEGPIX **VegMap,
                     SOILPIX **SoilMap, CHANNEL *ChannelData,
                     LAKETABLE *LType);

void StoreSnapshots(int StepsDone, MAPSIZE *Map, DUMPSTRUCT *Dump,
                    TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                    VEGPIX **VegMap, SOILPIX **SoilMap,
                    CHANNEL *ChannelData, LAKETABLE *LType);

void ReadModelState(char *Path, DATE *Start, int Dt, MAPSIZE *Map,
                    TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                    VEGPIX **VegMap, SOILPIX **SoilMap,
                    CHANNEL *ChannelData, LAKETABLE *LType);

void SnapshotFile(char *Buffer, size_t Size, char *Dir, DATE *Valid, const char *Suffix);

int CanopyIndex(VEGPIX *Cell, int Canopy);

void SnowStats(DATE *Now, MAPSIZE *Map, OPTIONSTRUCT *Options, 
        TOPOPIX **TopoMap, SNOWPIX **Snow, int Dt);

void ValidateInputs(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  LAYER *Soil, SOILTABLE *SType, VEGTABLES *Veg, int NStats,
  METLOCATION *Stat);

void ValidateParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap,
  SOILPIX **SoilMap);

void ValidateWeights(MAPSIZE *Map, TOPOPIX **TopoMap, uchar ***MetWeights, int NStats);

float viscosity(float Tair, float Rh);

#endif
