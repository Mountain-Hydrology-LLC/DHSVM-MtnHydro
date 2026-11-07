
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
#include "slopeaspect.h"
#include "mapio.h"

static void *ReadMap(LISTPTR Input, const char *Name, int NumberType, const char *Units,
                     int LayerSet, MAPSIZE *Map);
static float *ReadGivenMap(LISTPTR Input, OPTIONSTRUCT *Options, int m, MAPSIZE *Map);
static float ClassValue(int m, int Layer, int Row, SOILTABLE *SType, VEGTABLES *Veg);
static float *FormMap(int m, const unsigned char *Class, int NClasses, MAPSIZE *Map,
                      SOILTABLE *SType, VEGTABLES *Veg);

/* The soil and vegetation maps a class table can give, CM_* order: map
   name ([INPUTS] key), class table ([INPUTS] key), class column, units,
   layer set */
const CLASSMAP ClassMaps[NCLASSMAPS] = {
  {"Soil.Depth", "Soil Table", "Depth", "m", MAP_ONE},
  {"Soil.LateralConductivity", "Soil Table", "LateralConductivity", "m/s", MAP_ONE},
  {"Soil.ExponentialDecrease", "Soil Table", "ExponentialDecrease", "1/m", MAP_ONE},
  {"Soil.Porosity", "Soil Table", "Porosity", "m3/m3", MAP_SOIL},
  {"Soil.FieldCapacity", "Soil Table", "FieldCapacity", "m3/m3", MAP_SOIL},
  {"Veg.OverstoryCover", "Overstory Table", "Cover", "fraction", MAP_ONE},
  {"Veg.OverstoryLAI", "Overstory Table", "LAI", "m2/m2", MAP_ONE},
  {"Veg.OverstoryHeight", "Overstory Table", "Height", "m", MAP_ONE},
  {"Veg.UnderstoryCover", "Understory Table", "Cover", "fraction", MAP_ONE},
  {"Veg.UnderstoryLAI", "Understory Table", "LAI", "m2/m2", MAP_ONE},
  {"Veg.UnderstoryHeight", "Understory Table", "Height", "m", MAP_ONE}
};

/*****************************************************************************
  ReadMap()

  Allocates and reads every layer of the input map Name ([INPUTS] key);
  the header is checked first
*****************************************************************************/
static void *ReadMap(LISTPTR Input, const char *Name, int NumberType, const char *Units,
                     int LayerSet, MAPSIZE *Map)
{
  void *Values;
  size_t NValues = (size_t) Map->NY * Map->NX * MapLayerCount(LayerSet, Map->SoilLayers);

  if (!(Values = calloc(NValues, SizeOfNumberType(NumberType))))
    ReportError("ReadMap", 1);
  ReadInputMap(Input, Name, NumberType, Units, LayerSet, Map, Values);
  return Values;
}

/* Map m of ClassMaps[] as read, NULL when FormClassMaps() forms it */
static float *ReadGivenMap(LISTPTR Input, OPTIONSTRUCT *Options, int m, MAPSIZE *Map)
{
  if (Options->FromTable[m])
    return NULL;
  return ReadMap(Input, ClassMaps[m].Name, NC_FLOAT, ClassMaps[m].Units, ClassMaps[m].LayerSet,
                 Map);
}

/*****************************************************************************
  ReadInputMaps()

  Reads every input map the run uses, as given (the values are checked by
  ValidateInputs(), then ApplyParameters() applies the parameters); the soil
  and vegetation maps without an [INPUTS] key are formed by FormClassMaps().
  The monthly precipitation pattern stays open: its layers are read by
  ValidateInputs() and, a month at a time, by InitNewMonth().
*****************************************************************************/
void ReadInputMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw)
{
  char Path[BUFSIZE + 1];

  Message("\nReading input maps\n");

  memset(Raw, 0, sizeof(INPUTMAPS));
  Raw->Elevation = ReadMap(Input, "Terrain.Elevation", NC_FLOAT, "m", MAP_ONE, Map);
  Raw->Mask = ReadMap(Input, "Terrain.Mask", NC_BYTE, "-", MAP_ONE, Map);

  /* Horizon angles: the map is open until InitShading() has read it */
  Options->HorizonAngle.Ptr = NULL;
  if (Options->Shading) {
    RequireConfigValue(Input, "INPUTS", "Terrain.HorizonAngle", Path);
    OpenInputMap(&(Options->HorizonAngle), Path, "Terrain.HorizonAngle", NC_FLOAT,
                 "degrees", MAP_AZIMUTH, Map);
    Raw->SkyView = ReadMap(Input, "Terrain.SkyView", NC_FLOAT, "fraction", MAP_ONE, Map);
  }

  Raw->SoilClass = ReadMap(Input, "Soil.Class", NC_BYTE, "-", MAP_ONE, Map);
  Raw->Depth = ReadGivenMap(Input, Options, CM_DEPTH, Map);
  Raw->LateralConductivity = ReadGivenMap(Input, Options, CM_KSLAT, Map);
  Raw->ExponentialDecrease = ReadGivenMap(Input, Options, CM_EXPDEC, Map);
  Raw->Porosity = ReadGivenMap(Input, Options, CM_POROSITY, Map);
  Raw->FieldCapacity = ReadGivenMap(Input, Options, CM_FCAP, Map);
  Raw->LayerThickness = ReadMap(Input, "Soil.LayerThickness", NC_FLOAT, "m", MAP_ROOT, Map);
  Raw->Veg.OverstoryClass = ReadMap(Input, "Veg.OverstoryClass", NC_BYTE, "-", MAP_ONE, Map);
  Raw->Veg.UnderstoryClass = ReadMap(Input, "Veg.UnderstoryClass", NC_BYTE, "-", MAP_ONE, Map);
  Raw->Veg.OverstoryCover = ReadGivenMap(Input, Options, CM_OVERCOVER, Map);
  Raw->Veg.OverstoryLAI = ReadGivenMap(Input, Options, CM_OVERLAI, Map);
  Raw->Veg.OverstoryHeight = ReadGivenMap(Input, Options, CM_OVERHEIGHT, Map);
  Raw->Veg.UnderstoryCover = ReadGivenMap(Input, Options, CM_UNDERCOVER, Map);
  Raw->Veg.UnderstoryLAI = ReadGivenMap(Input, Options, CM_UNDERLAI, Map);
  Raw->Veg.UnderstoryHeight = ReadGivenMap(Input, Options, CM_UNDERHEIGHT, Map);

  Options->PrecipPattern.Ptr = NULL;
  if (Options->Prism) {
    RequireConfigValue(Input, "INPUTS", "Met.PrecipPattern", Path);
    OpenInputMap(&(Options->PrecipPattern), Path, "Met.PrecipPattern", NC_FLOAT,
                 "-", MAP_MONTH, Map);
  }
  if (Options->SnowPattern)
    Raw->SnowPattern = ReadMap(Input, "Met.SnowPattern", NC_FLOAT, "-", MAP_ONE, Map);
  if (Options->MeltScale)
    Raw->MeltScale = ReadMap(Input, "Snow.MeltScale", NC_FLOAT, "-", MAP_ONE, Map);
  if (Options->WindScale)
    Raw->WindScale = ReadMap(Input, "Met.WindScale", NC_FLOAT, "-", MAP_ONE, Map);
  if (Options->LakeDynamics)
    Raw->LakeID = ReadMap(Input, "Lake.ID", NC_BYTE, "-", MAP_ONE, Map);
  if (!IsEmptyStr(Options->InitialSWEPath))
    Raw->InitialSWE = ReadMap(Input, "Snow.InitialSWE", NC_FLOAT, "m", MAP_ONE, Map);
}

/* Value of map m of ClassMaps[] (soil layer Layer) of class row Row */
static float ClassValue(int m, int Layer, int Row, SOILTABLE *SType, VEGTABLES *Veg)
{
  switch (m) {
  case CM_DEPTH:
    return SType[Row].Depth;
  case CM_KSLAT:
    return SType[Row].KsLat;
  case CM_EXPDEC:
    return SType[Row].KsLatExp;
  case CM_POROSITY:
    return SType[Row].Porosity[Layer];
  case CM_FCAP:
    return SType[Row].FCap[Layer];
  case CM_OVERCOVER:
    return Veg->Over[Row].Story.Cover;
  case CM_OVERLAI:
    return Veg->Over[Row].Story.LAI;
  case CM_OVERHEIGHT:
    return Veg->Over[Row].Story.Height;
  case CM_UNDERCOVER:
    return Veg->Under[Row].Story.Cover;
  case CM_UNDERLAI:
    return Veg->Under[Row].Story.LAI;
  default:
    return Veg->Under[Row].Story.Height;
  }
}

/* Map m of ClassMaps[] formed from its class column: on every layer the
   value of the cell's class where the class map holds a row 1..NClasses,
   0 elsewhere (a vegetation class 0 has no story and leaves it unused; a
   soil class 0, a class outside the table or no-data is reported by the
   class checks) */
static float *FormMap(int m, const unsigned char *Class, int NClasses, MAPSIZE *Map,
                      SOILTABLE *SType, VEGTABLES *Veg)
{
  size_t NCells = (size_t) Map->NY * Map->NX;
  int NLayers = MapLayerCount(ClassMaps[m].LayerSet, Map->SoilLayers);
  float *Values;
  size_t i;
  int k;

  if (!(Values = (float *) calloc(NCells * NLayers, sizeof(float))))
    ReportError("FormMap", 1);
  for (k = 0; k < NLayers; k++) {
    for (i = 0; i < NCells; i++) {
      if (Class[i] >= 1 && Class[i] <= NClasses)
        Values[k * NCells + i] = ClassValue(m, k, Class[i] - 1, SType, Veg);
    }
  }
  return Values;
}

/*****************************************************************************
  FormVegMaps()

  The vegetation maps of V without an [INPUTS] key, formed from the
  overstory and understory tables through the class maps of V
*****************************************************************************/
void FormVegMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, VEGTABLES *Veg, VEGMAPS *V)
{
  int m;

  for (m = CM_OVERCOVER; m <= CM_OVERHEIGHT; m++) {
    if (Options->FromTable[m])
      *VegMapValues(V, m) = FormMap(m, V->OverstoryClass, Veg->NOver, Map, NULL, Veg);
  }
  for (m = CM_UNDERCOVER; m <= CM_UNDERHEIGHT; m++) {
    if (Options->FromTable[m])
      *VegMapValues(V, m) = FormMap(m, V->UnderstoryClass, Veg->NUnder, Map, NULL, Veg);
  }
}

/*****************************************************************************
  FormClassMaps()

  The soil and vegetation input maps without an [INPUTS] key, formed from
  the class tables through the class maps (after ReadInputMaps() and
  InitTables(); checked with the maps as read by ValidateInputs())
*****************************************************************************/
void FormClassMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, LAYER *Soil, SOILTABLE *SType,
                   VEGTABLES *Veg, INPUTMAPS *Raw)
{
  if (Options->FromTable[CM_DEPTH])
    Raw->Depth = FormMap(CM_DEPTH, Raw->SoilClass, Soil->NTypes, Map, SType, NULL);
  if (Options->FromTable[CM_KSLAT])
    Raw->LateralConductivity = FormMap(CM_KSLAT, Raw->SoilClass, Soil->NTypes, Map, SType, NULL);
  if (Options->FromTable[CM_EXPDEC])
    Raw->ExponentialDecrease = FormMap(CM_EXPDEC, Raw->SoilClass, Soil->NTypes, Map, SType, NULL);
  if (Options->FromTable[CM_POROSITY])
    Raw->Porosity = FormMap(CM_POROSITY, Raw->SoilClass, Soil->NTypes, Map, SType, NULL);
  if (Options->FromTable[CM_FCAP])
    Raw->FieldCapacity = FormMap(CM_FCAP, Raw->SoilClass, Soil->NTypes, Map, SType, NULL);
  FormVegMaps(Options, Map, Veg, &(Raw->Veg));
}

/*****************************************************************************
  InitTopoMap()

  Elevation and basin mask of every cell, lake ID of the cells inside the
  mask; lake areas; slope, aspect and flow directions
*****************************************************************************/
void InitTopoMap(OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  TOPOPIX ***TopoMap, LAKETABLE *LType)
{
  const char *Routine = "InitTopoMap";
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  float Area; /* Area of a lake */

  Message("\nInitializing terrain maps\n");

  if (!(*TopoMap = (TOPOPIX **)calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *)calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *)Routine, 1);
  }

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      (*TopoMap)[y][x].Dem = Raw->Elevation[i];
      (*TopoMap)[y][x].Mask = Raw->Mask[i];
      (*TopoMap)[y][x].LakeID = 0;
      if (Options->LakeDynamics && INBASIN((*TopoMap)[y][x].Mask))
        (*TopoMap)[y][x].LakeID = Raw->LakeID[i];
      (*TopoMap)[y][x].WindScale = (Raw->WindScale != NULL) ? Raw->WindScale[i] : 1.0;
    }
  }
  free(Raw->Elevation);
  free(Raw->Mask);
  free(Raw->WindScale);
  Raw->Elevation = NULL;
  Raw->Mask = NULL;
  Raw->WindScale = NULL;
  if (Raw->LakeID != NULL) {
    free(Raw->LakeID);
    Raw->LakeID = NULL;
  }

  /* Calculate the area of each lake */
  if (Options->LakeDynamics) {
    Message("\nSetting up lakes\n");
    for (i = 0; i < Map->NumLakes; i++) {
      Area = 0.0;

      for (y = 0; y < Map->NY; y++) {
        for (x = 0; x < Map->NX; x++) {
          if (INBASIN((*TopoMap)[y][x].Mask) && (*TopoMap)[y][x].LakeID == i+1)
            Area += 1.0;
        }
      }

      LType[i].Area = Area * (Map->DX * Map->DY);
      LType[i].Storage = 0.0;
      LType[i].Inflow = 0.0;
      LType[i].Outflow = 0.0;

      Message("%20s: %4.2f km^2\t\tQ ~ %3.2f*L^%3.2f\n",
             LType[i].Name, LType[i].Area/(1000*1000),
             LType[i].PowLawScale, LType[i].PowLawExponent);
    }
    Message("Total: %d lakes\n\n",Map->NumLakes);
  }

  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and
     fraction of flow flowing in each direction based on the land surface
     slope. */
  ElevationSlopeAspect(Map, *TopoMap, Options->MultiFlowDir);
}
