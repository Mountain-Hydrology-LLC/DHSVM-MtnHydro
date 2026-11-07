/*
 * DESCRIPTION:  Vegetation updates: at each [INPUTS] Vegetation Update <n>
 *               instant the eight vegetation maps of every mask cell are
 *               replaced by <Vegetation Update Directory>/<Name>.<stamp>.map
 *               (Veg.OverstoryClass, Veg.UnderstoryClass, Veg.OverstoryCover,
 *               Veg.OverstoryLAI, Veg.OverstoryHeight, Veg.UnderstoryCover,
 *               Veg.UnderstoryLAI, Veg.UnderstoryHeight), before the time
 *               step starting then; a map without an [INPUTS] key is formed
 *               from its class table through the update's class map, and
 *               has no file. The latest update before Start gives
 *               the vegetation at Start. The water of a canopy the cell
 *               loses, and the detention storage of a cell without
 *               impervious fraction, becomes surface water of the cell.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "mapio.h"
#include "DHSVMChannel.h"

/*****************************************************************************
  InitVegUpdateTimes()

  The update instants [INPUTS] Vegetation Update 1, 2, ... (the list ends at
  the first missing number), on the step grid and increasing, with
  Vegetation Update Directory; any instant before, during or after the run
*****************************************************************************/
void InitVegUpdateTimes(LISTPTR Input, OPTIONSTRUCT *Options, TIMESTRUCT *Time)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  char KeyName[BUFSIZE + 1];
  char Value[BUFSIZE + 1];
  long long Seconds;
  int n, k;

  for (n = 0; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "Vegetation Update %d", n + 1);
    if (!GetConfigValue(Input, "INPUTS", KeyName, Value))
      break;
  }
  U->N = n;
  U->Next = 0;
  U->Time = NULL;
  U->Step = NULL;
  if (n > 0 && IsEmptyStr(U->Directory))
    ReportErrorCode(103, "[INPUTS] Vegetation Update 1 needs Vegetation Update Directory");
  if (n == 0 && !IsEmptyStr(U->Directory))
    ReportErrorCode(103, "[INPUTS] Vegetation Update Directory needs Vegetation Update 1");
  if (n == 0)
    return;

  if (!(U->Time = (DATE *) calloc(n, sizeof(DATE))) ||
      !(U->Step = (long long *) calloc(n, sizeof(long long))))
    ReportError("InitVegUpdateTimes", 1);

  for (k = 0; k < n; k++) {
    snprintf(KeyName, sizeof(KeyName), "Vegetation Update %d", k + 1);
    RequireConfigValue(Input, "INPUTS", KeyName, Value);
    if (!SScanIsoDate(Value, &(U->Time[k])))
      ReportErrorCode(103, "[INPUTS] %s = %s: expected YYYY-MM-DD HH:MM", KeyName, Value);
    Seconds = (MinutesOf(&(U->Time[k])) - MinutesOf(&(Time->Start))) * SECPMIN;
    if (Seconds % Time->Dt != 0)
      ReportErrorCode(110, "[INPUTS] %s = %s: not on the time step grid", KeyName, Value);
    U->Step[k] = Seconds / Time->Dt;
    if (k > 0 && U->Step[k] <= U->Step[k - 1])
      ReportErrorCode(110, "[INPUTS] %s = %s: not after Vegetation Update %d", KeyName, Value, k);
  }
}

/* Reads one map of update k into Values */
static void ReadUpdateMap(VEGUPDATES *U, int k, const char *Name, int NumberType,
                          const char *Units, MAPSIZE *Map, void *Values)
{
  char Path[BUFSIZE * 2 + 1];
  char Stamp[NAMESIZE + 1];
  MAPFILE File;

  SPrintIsoStamp(&(U->Time[k]), Stamp);
  if (snprintf(Path, sizeof(Path), "%s%s.%s.map", U->Directory, Name, Stamp) >= (int) sizeof(Path))
    ReportErrorCode(104, "%s%s.%s.map: path too long", U->Directory, Name, Stamp);
  OpenInputMap(&File, Path, Name, NumberType, Units, MAP_ONE, Map);
  ReadMapLayer(&File, 0, Values, Map);
  CloseMap(&File);
}

/*****************************************************************************
  ReadVegUpdate()

  Reads and checks the maps of update k: the two class maps and the maps
  whose [INPUTS] keys are given; the others are formed from the class tables
  through the update's class maps (the checks of the input vegetation maps)
*****************************************************************************/
static void ReadVegUpdate(OPTIONSTRUCT *Options, VEGUPDATES *U, int k, MAPSIZE *Map,
                          unsigned char *Mask, VEGTABLES *Veg, VEGMAPS *V)
{
  size_t NCells = (size_t) Map->NY * Map->NX;
  char Suffix[NAMESIZE + 2];
  float **Values;
  int m;

  memset(V, 0, sizeof(VEGMAPS));
  if (!(V->OverstoryClass = (unsigned char *) calloc(NCells, sizeof(unsigned char))) ||
      !(V->UnderstoryClass = (unsigned char *) calloc(NCells, sizeof(unsigned char))))
    ReportError("ReadVegUpdate", 1);
  ReadUpdateMap(U, k, "Veg.OverstoryClass", NC_BYTE, "-", Map, V->OverstoryClass);
  ReadUpdateMap(U, k, "Veg.UnderstoryClass", NC_BYTE, "-", Map, V->UnderstoryClass);
  for (m = CM_OVERCOVER; m <= CM_UNDERHEIGHT; m++) {
    if (Options->FromTable[m])
      continue;
    Values = VegMapValues(V, m);
    if (!(*Values = (float *) calloc(NCells, sizeof(float))))
      ReportError("ReadVegUpdate", 1);
    ReadUpdateMap(U, k, ClassMaps[m].Name, NC_FLOAT, ClassMaps[m].Units, Map, *Values);
  }
  FormVegMaps(Options, Map, Veg, V);

  Suffix[0] = '.';
  SPrintIsoStamp(&(U->Time[k]), Suffix + 1);
  ValidateVegMaps(Options, Map, Mask, V, Veg, Suffix);
}

/*****************************************************************************
  VegMapValues()

  The array of map m of ClassMaps[] (CM_OVERCOVER..CM_UNDERHEIGHT) in V
*****************************************************************************/
float **VegMapValues(VEGMAPS *V, int m)
{
  switch (m) {
  case CM_OVERCOVER:
    return &(V->OverstoryCover);
  case CM_OVERLAI:
    return &(V->OverstoryLAI);
  case CM_OVERHEIGHT:
    return &(V->OverstoryHeight);
  case CM_UNDERCOVER:
    return &(V->UnderstoryCover);
  case CM_UNDERLAI:
    return &(V->UnderstoryLAI);
  default:
    return &(V->UnderstoryHeight);
  }
}

/*****************************************************************************
  FreeVegMaps()

  Frees the eight vegetation maps of V
*****************************************************************************/
void FreeVegMaps(VEGMAPS *V)
{
  free(V->OverstoryClass);
  free(V->UnderstoryClass);
  free(V->OverstoryCover);
  free(V->OverstoryLAI);
  free(V->OverstoryHeight);
  free(V->UnderstoryCover);
  free(V->UnderstoryLAI);
  free(V->UnderstoryHeight);
  memset(V, 0, sizeof(VEGMAPS));
}

/* Marks the class combinations of [OUTPUT] Classes that the update maps V
   give to a mask cell */
static void MarkClasses(VEGMAPS *V, MAPSIZE *Map, unsigned char *Mask, DUMPSTRUCT *Dump)
{
  size_t NCells = (size_t) Map->NY * Map->NX;
  size_t i;
  int j;

  for (j = 0; j < Dump->NClasses; j++) {
    for (i = 0; i < NCells && !Dump->Classes[j].Present; i++) {
      if (Mask[i] == 1 && V->OverstoryClass[i] == Dump->Classes[j].Over &&
          V->UnderstoryClass[i] == Dump->Classes[j].Under)
        Dump->Classes[j].Present = TRUE;
    }
  }
}

/*****************************************************************************
  InitVegUpdates()

  Reads and checks every update before End (so that a bad map stops the run
  before it starts) and marks the class combinations of [OUTPUT] Classes
  they give to mask cells; the latest one before Start replaces the
  vegetation maps as read, which then give the vegetation at Start
*****************************************************************************/
void InitVegUpdates(OPTIONSTRUCT *Options, MAPSIZE *Map, TIMESTRUCT *Time, INPUTMAPS *Raw,
                    VEGTABLES *Veg, DUMPSTRUCT *Dump)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  char Stamp[NAMESIZE + 1];
  VEGMAPS V;
  int Latest = -1;
  int k;

  if (U->N == 0)
    return;

  Message("\nChecking the vegetation updates\n");

  /* The vegetation maps as given are checked even when an update replaces
     them */
  ValidateVegMaps(Options, Map, Raw->Mask, &(Raw->Veg), Veg, "");

  U->Next = U->N;
  for (k = 0; k < U->N && U->Step[k] < Time->NSteps; k++) {
    if (U->Step[k] < 0)
      Latest = k;
    else if (U->Next == U->N)
      U->Next = k;
    ReadVegUpdate(Options, U, k, Map, Raw->Mask, Veg, &V);
    MarkClasses(&V, Map, Raw->Mask, Dump);
    FreeVegMaps(&V);
  }

  if (Latest >= 0) {
    ReadVegUpdate(Options, U, Latest, Map, Raw->Mask, Veg, &V);
    FreeVegMaps(&(Raw->Veg));
    Raw->Veg = V;
    SPrintIsoStamp(&(U->Time[Latest]), Stamp);
    Message("Vegetation at Start from the update of %s\n", Stamp);
  }
}

/*****************************************************************************
  ApplyVegUpdate()

  At the start of a time step with an update: the vegetation of every mask
  cell from the update's maps (rows 6 and 7, the albedo decay rates of the snow, the LAI of the
  month and the ground albedo under the snow); intercepted water moves with
  its canopy, a canopy the cell loses leaves its water on the ground
*****************************************************************************/
void ApplyVegUpdate(OPTIONSTRUCT *Options, TIMESTRUCT *Time, MAPSIZE *Map, TOPOPIX **TopoMap,
                    VEGTABLES *Veg, VEGPIX **VegMap, SOILTABLE *SType,
                    SOILPIX **SoilMap, SNOWPIX **SnowMap, PRECIPPIX **PrecipMap)
{
  VEGUPDATES *U = &(Options->VegUpdates);
  unsigned char *Mask;
  char Stamp[NAMESIZE + 1];
  VEGMAPS V;
  PRECIPPIX *Precip;
  VEGPIX *Cell;
  float Rain[2], Snow[2];
  float Surface;
  int x, y, i, c, j;

  if (U->Next >= U->N || U->Step[U->Next] != (long long) Time->Step)
    return;

  SPrintIsoStamp(&(U->Time[U->Next]), Stamp);
  Message("Updating the vegetation (%s)\n", Stamp);

  if (!(Mask = (unsigned char *) calloc((size_t) Map->NY * Map->NX, sizeof(unsigned char))))
    ReportError("ApplyVegUpdate", 1);
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      Mask[i] = INBASIN(TopoMap[y][x].Mask) ? 1 : 0;
  ReadVegUpdate(Options, U, U->Next, Map, Mask, Veg, &V);

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (!Mask[i])
        continue;
      Cell = &(VegMap[y][x]);
      Precip = &(PrecipMap[y][x]);

      /* Intercepted water by canopy (overstory, understory) of the old
         vegetation */
      for (c = 0; c < 2; c++) {
        j = CanopyIndex(Cell, c);
        Rain[c] = (j >= 0) ? Precip->IntRain[j] : 0.0;
        Snow[c] = (j >= 0) ? Precip->IntSnow[j] : 0.0;
      }

      /* The cell's new vegetation */
      SetCellVegetation(Options, Veg, Cell, &V, i, y, x);
      SetCellLAI(Veg, Cell, Time->Current.JDay);

      /* The water of a canopy the cell lost becomes surface water; the
         water of the others moves to their new positions */
      Surface = 0.0;
      if (Cell->OverStory != TRUE) {
        Surface += Rain[0] + Snow[0] + Precip->TempIntStorage;
        Precip->TempIntStorage = 0.0;
      }
      if (Cell->UnderStory != TRUE)
        Surface += Rain[1] + Snow[1];
      for (j = 0; j < MAXVEGLAYERS; j++) {
        Precip->IntRain[j] = 0.0;
        Precip->IntSnow[j] = 0.0;
      }
      for (c = 0; c < 2; c++) {
        j = CanopyIndex(Cell, c);
        if (j >= 0) {
          Precip->IntRain[j] = Rain[c];
          Precip->IntSnow[j] = Snow[c];
        }
      }

      /* Detention storage only drains from an impervious fraction */
      if (!(ImperviousFraction(Veg, Cell) > 0.0)) {
        Surface += SoilMap[y][x].DetentionStorage;
        SoilMap[y][x].DetentionStorage = 0.0;
      }
      SoilMap[y][x].IExcess += Surface;

      /* Albedo decay rates of the snow (the forest offset on overstory
         cells) and the ground albedo under the snow, as at Start */
      SnowMap[y][x].LamdaAcc = ALB_ACC_LAMBDA;
      SnowMap[y][x].LamdaMelt = ALB_MELT_LAMBDA;
      if (Cell->OverStory == TRUE) {
        SnowMap[y][x].LamdaAcc += (LAMBDA_FOREST_OFFSET * Cell->Fract[0]);
        SnowMap[y][x].LamdaMelt += (LAMBDA_FOREST_OFFSET * Cell->Fract[0]);
      }
      SnowMap[y][x].AlbedoGround = GroundAlbedo(Veg, Cell,
        SType[SoilMap[y][x].Soil - 1].Albedo, Time->Current.Month);
    }
  }

  FreeVegMaps(&V);
  free(Mask);
  U->Next++;
}
