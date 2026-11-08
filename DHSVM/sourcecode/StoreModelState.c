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
#include "channel_grid.h"
#include "mapio.h"

 /*****************************************************************************
   Snapshots

   A snapshot holds every value the model carries from one time step to the
   next, valid at one instant T, in three files named by the stamp
   YYYY-MM-DD_HHMM of T:

     State.<stamp>.map               the cell values, one FLOAT layer each
                                     (SnapshotLayers() gives the order), -9999
                                     outside the mask and where a cell does
                                     not have the layer
     State.<stamp>.Segments.csv      the channel segments (ChannelState.c)
     State.<stamp>.ChannelCells.csv  the channel cells, one row per row of the
                                     stream map, in its order

   Values are stored exactly (float); flags as 0.0 and 1.0.
 *****************************************************************************/

#define CELL_COLUMNS "Row,Column,Segment,Channel.CellTableDepth,Channel.CellAvailableWater"
#define CELL_UNITS "-,-,-,m,m"

/* What a snapshot layer holds */
enum SNAPSHOTSLOT {
  SNAP_INTRAIN, SNAP_INTSNOW, SNAP_TEMPINT, SNAP_TCANOPY,
  SNAP_HASSNOW, SNAP_LASTSNOW, SNAP_SWQ, SNAP_PACKWATER, SNAP_TPACK,
  SNAP_SURFWATER, SNAP_TSURF, SNAP_ACCUMSEASON, SNAP_ALBEDO,
  SNAP_MOIST, SNAP_IEXCESS, SNAP_PERC, SNAP_SATFLOW, SNAP_TABLEDEPTH,
  SNAP_WATERLEVEL, SNAP_WATERLEVELLAST, SNAP_INFILTACC, SNAP_MOISTINIT,
  SNAP_PRECIPSTART, SNAP_RUNOFF, SNAP_STARTRUNOFF, SNAP_DETENTION,
  SNAP_LAKEOUTFLOW, SNAP_SUMPRECIP, SNAP_SNOWACCUM, SNAP_SNOWMELT,
  SNAP_MAXSWE, SNAP_MAXSWEDATE, SNAP_MELTOUTDATE, SNAP_DESORPTIONACC,
  SNAP_TACCLIM, SNAP_TGROWTH, SNAP_TNODE, SNAP_ICENODE,
  SNAP_SOILTSURF
};

typedef struct {
  int Slot;			/* SNAPSHOTSLOT */
  int Sub;			/* canopy (0 overstory, 1 understory), soil layer or node */
} SNAPLAYER;

static int SnapshotLayers(int L, int SoilFreezing, SNAPLAYER **Layers,
                          char (**Names)[MAP_NAMESIZE]);
static int CellHasLayer(SNAPLAYER *Layer, int y, int x, TOPOPIX **TopoMap,
                        VEGPIX **VegMap, LAKETABLE *LType);
static float GetLayerValue(SNAPLAYER *Layer, int y, int x, TOPOPIX **TopoMap,
                           PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, VEGPIX **VegMap,
                           SOILPIX **SoilMap, LAKETABLE *LType);
static void SetLayerValue(SNAPLAYER *Layer, int y, int x, float Value,
                          TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                          VEGPIX **VegMap, SOILPIX **SoilMap,
                          LAKETABLE *LType);

/*****************************************************************************
  SnapshotFile()

  Path of a snapshot file: <Dir>State.<stamp><Suffix>
*****************************************************************************/
void SnapshotFile(char *Buffer, size_t Size, char *Dir, DATE *Valid, const char *Suffix)
{
  char Stamp[NAMESIZE + 1];

  SPrintIsoStamp(Valid, Stamp);
  MakePath(Buffer, Size, "%sState.%s%s", Dir, Stamp, Suffix);
}

/*****************************************************************************
  CanopyIndex()

  Position of the overstory (Canopy 0) or the understory (Canopy 1) among
  the vegetation layers of a cell (0 = the top layer present), -1 if the
  cell lacks it
*****************************************************************************/
int CanopyIndex(VEGPIX *Cell, int Canopy)
{
  if (Canopy == 0)
    return (Cell->OverStory == TRUE) ? 0 : -1;
  if (Cell->UnderStory != TRUE)
    return -1;
  return (Cell->OverStory == TRUE) ? 1 : 0;
}

/*****************************************************************************
  SnapshotLayers()

  The layers of a snapshot map in file order, with their names; L is the
  number of root zone layers. Returns the number of layers: 37 + 2 L, and
  with SoilFreezing 2 SOIL_NODES + 1 more (the node temperatures, the node
  ice, the soil surface temperature).
*****************************************************************************/
static void AddLayer(SNAPLAYER *Layers, char (*Names)[MAP_NAMESIZE], int *n,
                     int Slot, int Sub, const char *Name)
{
  Layers[*n].Slot = Slot;
  Layers[*n].Sub = Sub;
  snprintf(Names[*n], MAP_NAMESIZE, "%s", Name);
  (*n)++;
}

static int SnapshotLayers(int L, int SoilFreezing, SNAPLAYER **Layers,
                          char (**Names)[MAP_NAMESIZE])
{
  char Name[MAP_NAMESIZE];
  int NLayers = 37 + 2 * L + (SoilFreezing ? 2 * SOIL_NODES + 1 : 0);
  int n = 0;
  int l;

  if (!(*Layers = (SNAPLAYER *) calloc(NLayers, sizeof(SNAPLAYER))) ||
      !(*Names = calloc(NLayers, MAP_NAMESIZE)))
    ReportError("SnapshotLayers", 1);

  AddLayer(*Layers, *Names, &n, SNAP_INTRAIN, 0, "Canopy.InterceptedRain.Overstory");
  AddLayer(*Layers, *Names, &n, SNAP_INTRAIN, 1, "Canopy.InterceptedRain.Understory");
  AddLayer(*Layers, *Names, &n, SNAP_INTSNOW, 0, "Canopy.InterceptedSnow.Overstory");
  AddLayer(*Layers, *Names, &n, SNAP_INTSNOW, 1, "Canopy.InterceptedSnow.Understory");
  AddLayer(*Layers, *Names, &n, SNAP_TEMPINT, 0, "Canopy.TemporaryStorage");
  AddLayer(*Layers, *Names, &n, SNAP_TCANOPY, 0, "Canopy.Temperature");
  AddLayer(*Layers, *Names, &n, SNAP_HASSNOW, 0, "Snow.Covered");
  AddLayer(*Layers, *Names, &n, SNAP_LASTSNOW, 0, "Snow.DaysSinceSnowfall");
  AddLayer(*Layers, *Names, &n, SNAP_SWQ, 0, "Snow.SWE");
  AddLayer(*Layers, *Names, &n, SNAP_PACKWATER, 0, "Snow.PackLiquid");
  AddLayer(*Layers, *Names, &n, SNAP_TPACK, 0, "Snow.PackTemperature");
  AddLayer(*Layers, *Names, &n, SNAP_SURFWATER, 0, "Snow.SurfaceLiquid");
  AddLayer(*Layers, *Names, &n, SNAP_TSURF, 0, "Snow.SurfaceTemperature");
  AddLayer(*Layers, *Names, &n, SNAP_ACCUMSEASON, 0, "Snow.AccumulationSeason");
  AddLayer(*Layers, *Names, &n, SNAP_ALBEDO, 0, "Snow.Albedo");
  for (l = 0; l < L; l++) {
    snprintf(Name, sizeof(Name), "Soil.Moisture.%d", l + 1);
    AddLayer(*Layers, *Names, &n, SNAP_MOIST, l, Name);
  }
  AddLayer(*Layers, *Names, &n, SNAP_MOIST, L, "Soil.Moisture.Deep");
  AddLayer(*Layers, *Names, &n, SNAP_IEXCESS, 0, "Soil.SurfaceWater");
  for (l = 0; l < L; l++) {
    snprintf(Name, sizeof(Name), "Soil.Percolation.%d", l + 1);
    AddLayer(*Layers, *Names, &n, SNAP_PERC, l, Name);
  }
  AddLayer(*Layers, *Names, &n, SNAP_SATFLOW, 0, "Soil.SatFlow");
  AddLayer(*Layers, *Names, &n, SNAP_TABLEDEPTH, 0, "Soil.TableDepth");
  AddLayer(*Layers, *Names, &n, SNAP_WATERLEVEL, 0, "Soil.WaterLevel");
  AddLayer(*Layers, *Names, &n, SNAP_WATERLEVELLAST, 0, "Soil.WaterLevelLast");
  AddLayer(*Layers, *Names, &n, SNAP_INFILTACC, 0, "Soil.InfiltrationAccumulated");
  AddLayer(*Layers, *Names, &n, SNAP_MOISTINIT, 0, "Soil.InfiltrationStartMoisture");
  AddLayer(*Layers, *Names, &n, SNAP_PRECIPSTART, 0, "Met.PrecipitationStart");
  AddLayer(*Layers, *Names, &n, SNAP_RUNOFF, 0, "Surface.Runoff");
  AddLayer(*Layers, *Names, &n, SNAP_STARTRUNOFF, 0, "Surface.RunoffStart");
  AddLayer(*Layers, *Names, &n, SNAP_DETENTION, 0, "Surface.Detention");
  AddLayer(*Layers, *Names, &n, SNAP_LAKEOUTFLOW, 0, "Lake.Outflow");
  AddLayer(*Layers, *Names, &n, SNAP_SUMPRECIP, 0, "Met.PrecipitationWaterYear");
  AddLayer(*Layers, *Names, &n, SNAP_SNOWACCUM, 0, "Snow.AccumulationWaterYear");
  AddLayer(*Layers, *Names, &n, SNAP_SNOWMELT, 0, "Snow.OutflowWaterYear");
  AddLayer(*Layers, *Names, &n, SNAP_MAXSWE, 0, "Snow.PeakSWE");
  AddLayer(*Layers, *Names, &n, SNAP_MAXSWEDATE, 0, "Snow.PeakSWEDay");
  AddLayer(*Layers, *Names, &n, SNAP_MELTOUTDATE, 0, "Snow.MeltOutDay");
  AddLayer(*Layers, *Names, &n, SNAP_DESORPTIONACC, 0, "Soil.DesorptionAccumulated");
  AddLayer(*Layers, *Names, &n, SNAP_TACCLIM, 0, "Canopy.AcclimationTemperature");
  AddLayer(*Layers, *Names, &n, SNAP_TGROWTH, 0, "Canopy.GrowthTemperature");
  if (SoilFreezing) {
    for (l = 0; l < SOIL_NODES; l++) {
      snprintf(Name, sizeof(Name), "Soil.NodeTemperature.%d", l + 1);
      AddLayer(*Layers, *Names, &n, SNAP_TNODE, l, Name);
    }
    for (l = 0; l < SOIL_NODES; l++) {
      snprintf(Name, sizeof(Name), "Soil.NodeIce.%d", l + 1);
      AddLayer(*Layers, *Names, &n, SNAP_ICENODE, l, Name);
    }
    AddLayer(*Layers, *Names, &n, SNAP_SOILTSURF, 0, "Soil.SurfaceTemperature");
  }

  return n;
}

/*****************************************************************************
  CellHasLayer()

  TRUE if a mask cell has the layer: a canopy layer only if the cell has
  that canopy, the lake outflow only on lake cells
*****************************************************************************/
static int CellHasLayer(SNAPLAYER *Layer, int y, int x, TOPOPIX **TopoMap,
                        VEGPIX **VegMap, LAKETABLE *LType)
{
  switch (Layer->Slot) {
  case SNAP_INTRAIN:
  case SNAP_INTSNOW:
    return CanopyIndex(&(VegMap[y][x]), Layer->Sub) >= 0;
  case SNAP_LAKEOUTFLOW:
    return (TopoMap[y][x].LakeID > 0 && LType != NULL);
  default:
    return TRUE;
  }
}

/*****************************************************************************
  GetLayerValue(), SetLayerValue()

  The value of a snapshot layer at a mask cell that has the layer
*****************************************************************************/
static float GetLayerValue(SNAPLAYER *Layer, int y, int x, TOPOPIX **TopoMap,
                           PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, VEGPIX **VegMap,
                           SOILPIX **SoilMap, LAKETABLE *LType)
{
  PRECIPPIX *Precip = &(PrecipMap[y][x]);
  SNOWPIX *Snow = &(SnowMap[y][x]);
  SOILPIX *Soil = &(SoilMap[y][x]);
  int i;

  switch (Layer->Slot) {
  case SNAP_INTRAIN:
    i = CanopyIndex(&(VegMap[y][x]), Layer->Sub);
    return Precip->IntRain[i];
  case SNAP_INTSNOW:
    i = CanopyIndex(&(VegMap[y][x]), Layer->Sub);
    return Precip->IntSnow[i];
  case SNAP_TEMPINT: return Precip->TempIntStorage;
  case SNAP_TCANOPY: return VegMap[y][x].Tcanopy;
  case SNAP_HASSNOW: return (float) Snow->HasSnow;
  case SNAP_LASTSNOW: return Snow->LastSnow;
  case SNAP_SWQ: return Snow->Swq;
  case SNAP_PACKWATER: return Snow->PackWater;
  case SNAP_TPACK: return Snow->TPack;
  case SNAP_SURFWATER: return Snow->SurfWater;
  case SNAP_TSURF: return Snow->TSurf;
  case SNAP_ACCUMSEASON: return (float) Snow->AccumSeason;
  case SNAP_ALBEDO: return Snow->Albedo;
  case SNAP_MOIST: return Soil->Moist[Layer->Sub];
  case SNAP_IEXCESS: return Soil->IExcess;
  case SNAP_PERC: return Soil->Perc[Layer->Sub];
  case SNAP_SATFLOW: return Soil->SatFlow;
  case SNAP_TABLEDEPTH: return Soil->TableDepth;
  case SNAP_WATERLEVEL: return Soil->WaterLevel;
  case SNAP_WATERLEVELLAST: return Soil->WaterLevelLast;
  case SNAP_INFILTACC: return Soil->InfiltAcc;
  case SNAP_MOISTINIT: return Soil->MoistInit;
  case SNAP_PRECIPSTART: return (float) Precip->PrecipStart;
  case SNAP_RUNOFF: return Soil->Runoff;
  case SNAP_STARTRUNOFF: return Soil->startRunoff;
  case SNAP_DETENTION: return Soil->DetentionStorage;
  case SNAP_LAKEOUTFLOW: return LType[TopoMap[y][x].LakeID - 1].Outflow;
  case SNAP_SUMPRECIP: return Precip->SumPrecip;
  case SNAP_SNOWACCUM: return Precip->SnowAccum;
  case SNAP_SNOWMELT: return Precip->SnowMelt;
  case SNAP_MAXSWE: return Snow->MaxSwe;
  case SNAP_MAXSWEDATE: return (float) Snow->MaxSweDate;
  case SNAP_MELTOUTDATE: return (float) Snow->MeltOutDate;
  case SNAP_DESORPTIONACC: return Soil->DesorptionAcc;
  case SNAP_TACCLIM: return VegMap[y][x].TAcclim;
  case SNAP_TGROWTH: return VegMap[y][x].TGrowth;
  case SNAP_TNODE: return Soil->TNode[Layer->Sub];
  case SNAP_ICENODE: return Soil->IceNode[Layer->Sub];
  case SNAP_SOILTSURF: return Soil->TSurf;
  default:
    return NA;
  }
}

static void SetLayerValue(SNAPLAYER *Layer, int y, int x, float Value,
                          TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                          VEGPIX **VegMap, SOILPIX **SoilMap,
                          LAKETABLE *LType)
{
  PRECIPPIX *Precip = &(PrecipMap[y][x]);
  SNOWPIX *Snow = &(SnowMap[y][x]);
  SOILPIX *Soil = &(SoilMap[y][x]);
  int i;

  switch (Layer->Slot) {
  case SNAP_INTRAIN:
    i = CanopyIndex(&(VegMap[y][x]), Layer->Sub);
    Precip->IntRain[i] = Value;
    break;
  case SNAP_INTSNOW:
    i = CanopyIndex(&(VegMap[y][x]), Layer->Sub);
    Precip->IntSnow[i] = Value;
    break;
  case SNAP_TEMPINT: Precip->TempIntStorage = Value; break;
  case SNAP_TCANOPY: VegMap[y][x].Tcanopy = Value; break;
  case SNAP_HASSNOW: Snow->HasSnow = (unsigned char) Value; break;
  case SNAP_LASTSNOW: Snow->LastSnow = Value; break;
  case SNAP_SWQ: Snow->Swq = Value; break;
  case SNAP_PACKWATER: Snow->PackWater = Value; break;
  case SNAP_TPACK: Snow->TPack = Value; break;
  case SNAP_SURFWATER: Snow->SurfWater = Value; break;
  case SNAP_TSURF: Snow->TSurf = Value; break;
  case SNAP_ACCUMSEASON: Snow->AccumSeason = (unsigned char) Value; break;
  case SNAP_ALBEDO: Snow->Albedo = Value; break;
  case SNAP_MOIST:
    Soil->Moist[Layer->Sub] = Value;
    Soil->InterFlow[Layer->Sub] = 0.0;
    break;
  case SNAP_IEXCESS: Soil->IExcess = Value; break;
  case SNAP_PERC: Soil->Perc[Layer->Sub] = Value; break;
  case SNAP_SATFLOW: Soil->SatFlow = Value; break;
  case SNAP_TABLEDEPTH: Soil->TableDepth = Value; break;
  case SNAP_WATERLEVEL: Soil->WaterLevel = Value; break;
  case SNAP_WATERLEVELLAST: Soil->WaterLevelLast = Value; break;
  case SNAP_INFILTACC: Soil->InfiltAcc = Value; break;
  case SNAP_MOISTINIT: Soil->MoistInit = Value; break;
  case SNAP_PRECIPSTART: Precip->PrecipStart = (int) Value; break;
  case SNAP_RUNOFF: Soil->Runoff = Value; break;
  case SNAP_STARTRUNOFF: Soil->startRunoff = Value; break;
  case SNAP_DETENTION: Soil->DetentionStorage = Value; break;
  case SNAP_LAKEOUTFLOW: LType[TopoMap[y][x].LakeID - 1].Outflow = Value; break;
  case SNAP_SUMPRECIP: Precip->SumPrecip = Value; break;
  case SNAP_SNOWACCUM: Precip->SnowAccum = Value; break;
  case SNAP_SNOWMELT: Precip->SnowMelt = Value; break;
  case SNAP_MAXSWE: Snow->MaxSwe = Value; break;
  case SNAP_MAXSWEDATE: Snow->MaxSweDate = (unint) Value; break;
  case SNAP_MELTOUTDATE: Snow->MeltOutDate = (unint) Value; break;
  case SNAP_DESORPTIONACC: Soil->DesorptionAcc = Value; break;
  case SNAP_TACCLIM: VegMap[y][x].TAcclim = Value; break;
  case SNAP_TGROWTH: VegMap[y][x].TGrowth = Value; break;
  case SNAP_TNODE: Soil->TNode[Layer->Sub] = Value; break;
  case SNAP_ICENODE: Soil->IceNode[Layer->Sub] = Value; break;
  case SNAP_SOILTSURF: Soil->TSurf = Value; break;
  default:
    break;
  }
}

/*****************************************************************************
  StoreModelState()

  Writes the snapshot valid at instant Valid into directory Path
*****************************************************************************/
void StoreModelState(char *Path, DATE *Valid, MAPSIZE *Map, int SoilFreezing,
                     TOPOPIX **TopoMap,
                     PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, VEGPIX **VegMap,
                     SOILPIX **SoilMap, CHANNEL *ChannelData,
                     LAKETABLE *LType)
{
  char FileName[BUFSIZE * 2 + 1];
  char Time[NAMESIZE + 1];
  char (*Names)[MAP_NAMESIZE];
  SNAPLAYER *Layers;
  ChannelMapPtr Cell;
  FILE *OutFile;
  float *Array;
  int NLayers;
  int k, x, y;

  if (DEBUG)
    Message("Storing model state\n");

  /* Cell values */
  NLayers = SnapshotLayers(Map->SoilLayers, SoilFreezing, &Layers, &Names);
  if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError("StoreModelState", 1);
  SPrintIsoDate(Valid, Time);
  SnapshotFile(FileName, sizeof(FileName), Path, Valid, ".map");
  OutFile = CreateMap(FileName, "State", "snapshot", NLayers, Names, Time, Map);
  for (k = 0; k < NLayers; k++) {
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask) &&
            CellHasLayer(&(Layers[k]), y, x, TopoMap, VegMap, LType))
          Array[y * Map->NX + x] =
            GetLayerValue(&(Layers[k]), y, x, TopoMap, PrecipMap, SnowMap, VegMap,
                          SoilMap, LType);
        else
          Array[y * Map->NX + x] = NA;
      }
    }
    WriteMapLayer(OutFile, Array, Map, FileName);
  }
  CloseMapFile(OutFile, FileName);
  free(Array);
  free(Layers);
  free(Names);

  /* Channel segments */
  SnapshotFile(FileName, sizeof(FileName), Path, Valid, ".Segments.csv");
  StoreChannelState(FileName, ChannelData->streams);

  /* Channel cells, in the order of the stream map */
  SnapshotFile(FileName, sizeof(FileName), Path, Valid, ".ChannelCells.csv");
  if (!(OutFile = fopen(FileName, "wb")))
    ReportErrorCode(113, "%s: cannot write", FileName);
  fprintf(OutFile, "%s\n%s\n", CELL_COLUMNS, CELL_UNITS);
  for (Cell = channel_grid_first_row(); Cell != NULL; Cell = Cell->next_row)
    fprintf(OutFile, "%d,%d,%d,%.9g,%.9g\n", Cell->row + 1, Cell->col + 1,
            (int) Cell->channel->id, Cell->table_depth, Cell->avail_water);
  if (ferror(OutFile) || fclose(OutFile) != 0)
    ReportErrorCode(113, "%s: cannot write", FileName);
}

/*****************************************************************************
  ReadModelState()

  Reads the snapshot valid at Start from directory Path: the cell values,
  the channel segments and the channel cells. Every file must match the
  model's grid, mask, classes and channel network.
*****************************************************************************/
void ReadModelState(char *Path, DATE *Start, int Dt, MAPSIZE *Map, int SoilFreezing,
                    TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap,
                    VEGPIX **VegMap, SOILPIX **SoilMap,
                    CHANNEL *ChannelData, LAKETABLE *LType)
{
  const char *Suffix[3] = { ".map", ".Segments.csv", ".ChannelCells.csv" };
  char FileName[3][BUFSIZE * 2 + 1];
  char Time[NAMESIZE + 1];
  char Line[BUFSIZE * 2 + 1];
  char *Fields[6];
  char (*Names)[MAP_NAMESIZE];
  SNAPLAYER *Layers;
  ChannelMapPtr Cell;
  MAPFILE File;
  FILE *InFile;
  float *Array;
  float Value[2];
  long Number[3];
  int NLayers;
  int NoData, FirstY = 0, FirstX = 0;
  int LineNo;
  int i, k, x, y;

  Message("Reading the model state\n");

  /* The three files exist */
  SPrintIsoDate(Start, Time);
  for (i = 0; i < 3; i++) {
    SnapshotFile(FileName[i], sizeof(FileName[i]), Path, Start, Suffix[i]);
    if (!(InFile = fopen(FileName[i], "rb")))
      ReportErrorCode(112, "%s: not found (Initial State for Start %s)", FileName[i], Time);
    fclose(InFile);
  }

  /* Cell values */
  NLayers = SnapshotLayers(Map->SoilLayers, SoilFreezing, &Layers, &Names);
  CheckMap(&File, FileName[0], "State", NC_FLOAT, "snapshot", NLayers, Names, Time, Map, 112);
  if (!(Array = (float *) calloc(Map->NY * Map->NX, sizeof(float))))
    ReportError("ReadModelState", 1);
  for (k = 0; k < NLayers; k++) {
    ReadMapLayer(&File, k, Array, Map);
    NoData = 0;
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask) &&
            CellHasLayer(&(Layers[k]), y, x, TopoMap, VegMap, LType) &&
            Array[y * Map->NX + x] == NA) {
          if (NoData == 0) {
            FirstY = y;
            FirstX = x;
          }
          NoData++;
        }
      }
    }
    if (NoData > 0)
      ReportErrorCode(112, "%s: layer %s is no-data in %d mask cells; first at row %d, column %d",
                      FileName[0], Names[k], NoData, FirstY + 1, FirstX + 1);
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask) &&
            CellHasLayer(&(Layers[k]), y, x, TopoMap, VegMap, LType))
          SetLayerValue(&(Layers[k]), y, x, Array[y * Map->NX + x], TopoMap,
                        PrecipMap, SnowMap, VegMap, SoilMap, LType);
      }
    }
  }
  CloseMap(&File);
  free(Array);
  free(Layers);
  free(Names);

  /* Channel segments (no channel routing without the soil modules) */
#ifndef SNOW_ONLY
  ReadChannelState(FileName[1], Dt, ChannelData->streams);
#endif

  /* Channel cells, in the order of the stream map */
  if (!(InFile = fopen(FileName[2], "rb")))
    ReportErrorCode(112, "%s: not found", FileName[2]);
  if (!ReadTextLine(InFile, Line, sizeof(Line), FileName[2], 1, 112) ||
      strcmp(Line, CELL_COLUMNS) != 0)
    ReportErrorCode(112, "%s: line 1: the columns must be %s", FileName[2], CELL_COLUMNS);
  if (!ReadTextLine(InFile, Line, sizeof(Line), FileName[2], 2, 112) ||
      strcmp(Line, CELL_UNITS) != 0)
    ReportErrorCode(112, "%s: line 2: the units must be %s", FileName[2], CELL_UNITS);
  LineNo = 2;
  for (Cell = channel_grid_first_row(); Cell != NULL; Cell = Cell->next_row) {
    LineNo++;
    if (!ReadTextLine(InFile, Line, sizeof(Line), FileName[2], LineNo, 112) ||
        IsBlankLine(Line))
      ReportErrorCode(112, "%s: line %d: no row, expected (%d, %d, %d)", FileName[2],
                      LineNo, Cell->row + 1, Cell->col + 1, (int) Cell->channel->id);
    if (SplitFields(Line, Fields, 6) != 5)
      ReportErrorCode(112, "%s: line %d: expected 5 fields", FileName[2], LineNo);
    for (i = 0; i < 3; i++) {
      if (!TextToInteger(Fields[i], &(Number[i])))
        Number[i] = -1;
    }
    if (Number[0] != Cell->row + 1 || Number[1] != Cell->col + 1 ||
        Number[2] != (long) Cell->channel->id)
      ReportErrorCode(112, "%s: line %d: (%s, %s, %s), expected (%d, %d, %d)", FileName[2],
                      LineNo, Fields[0], Fields[1], Fields[2], Cell->row + 1,
                      Cell->col + 1, (int) Cell->channel->id);
    for (i = 0; i < 2; i++) {
      if (!TextToFloat(Fields[i + 3], &(Value[i])))
        ReportErrorCode(112, "%s: line %d: column %d: \"%s\" is not a number", FileName[2],
                        LineNo, i + 4, Fields[i + 3]);
    }
    Cell->table_depth = Value[0];
    Cell->avail_water = Value[1];
  }
  while (ReadTextLine(InFile, Line, sizeof(Line), FileName[2], ++LineNo, 112)) {
    if (!IsBlankLine(Line))
      ReportErrorCode(112, "%s: line %d: a row after the last channel cell", FileName[2],
                      LineNo);
  }
  fclose(InFile);
}
