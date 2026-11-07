/*
 * DESCRIPTION:  Checks of the inputs before the run starts. The run stops
 *               at the first failure with one line on stderr and the exit
 *               code of the check:
 *               104 an input file or directory cannot be opened
 *               108 cell values: the input maps and the vegetation
 *                   updates as given, before any parameter applies; the
 *                   porosity after the parameters; the station weights
 *               109 the channel tables and the lakes
 *               Cell messages name the map, the condition, the number of
 *               cells and the first one (rows and columns from 1, row after
 *               row from north to south, each row from west to east).
 *               The summary of the classes in the basin is printed here.
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
#include "mapio.h"
#include "DHSVMChannel.h"

static const char *MonthNames[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

/* Cells failing one condition: how many, and the first one */
typedef struct {
  long N;
  int Row;
  int Column;
  double Value[2];
} FAILURES;

static void Fail(FAILURES *F, int y, int x, double Value1, double Value2)
{
  if (F->N == 0) {
    F->Row = y + 1;
    F->Column = x + 1;
    F->Value[0] = Value1;
    F->Value[1] = Value2;
  }
  F->N++;
}

/* Stops the run if a condition failed: "<What>: <Condition> in <n> cells;
   first at row <r>, column <c> (<Name1> = <v>[, <Name2> = <v>])" */
static void ReportCells(FAILURES *F, const char *What, const char *Condition,
                        const char *Name1, const char *Name2)
{
  if (F->N == 0)
    return;
  if (Name2 == NULL)
    ReportErrorCode(108, "%s: %s in %ld cells; first at row %d, column %d (%s = %g)",
                    What, Condition, F->N, F->Row, F->Column, Name1, F->Value[0]);
  else
    ReportErrorCode(108, "%s: %s in %ld cells; first at row %d, column %d (%s = %g, %s = %g)",
                    What, Condition, F->N, F->Row, F->Column, Name1, F->Value[0], Name2,
                    F->Value[1]);
}

/* Name of a soil layer of a map: <Name>.1 .. <Name>.L, <Name>.Deep */
static void SoilLayerName(char *Buffer, size_t Size, const char *Name, int Layer,
                          int SoilLayers)
{
  if (Layer < SoilLayers)
    snprintf(Buffer, Size, "%s.%d", Name, Layer + 1);
  else
    snprintf(Buffer, Size, "%s.Deep", Name);
}

/* Name of map m of ClassMaps[] (soil layer Layer) in messages:
   <Name>[.<layer>]<Suffix>, followed by " (<table key> <column>[_<layer>])"
   when the map is formed from its class table */
static void ClassMapName(char *Buffer, size_t Size, OPTIONSTRUCT *Options, int m, int Layer,
                         int SoilLayers, const char *Suffix)
{
  const CLASSMAP *C = &(ClassMaps[m]);
  char Token[NAMESIZE + 2] = "";

  if (C->LayerSet == MAP_SOIL && Layer < SoilLayers)
    snprintf(Token, sizeof(Token), "%d", Layer + 1);
  else if (C->LayerSet == MAP_SOIL)
    strcpy(Token, "Deep");
  if (!Options->FromTable[m])
    snprintf(Buffer, Size, "%s%s%s%s", C->Name, (Token[0] != '\0') ? "." : "", Token, Suffix);
  else
    snprintf(Buffer, Size, "%s%s%s%s (%s %s%s%s)", C->Name, (Token[0] != '\0') ? "." : "",
             Token, Suffix, C->Table, C->Column, (Token[0] != '\0') ? "_" : "", Token);
}

/*****************************************************************************
  CheckInputFiles()

  Every input file the run uses can be opened, and the station data
  directory exists (the station data files are opened with the station
  table). PointTable: TRUE when [OUTPUT] Point names a variable (a missing
  Point Table key is an output request error)
*****************************************************************************/
void CheckInputFiles(LISTPTR Input, OPTIONSTRUCT *Options, int PointTable)
{
  const char *Always[] = {
    "Terrain.Elevation", "Terrain.Mask", "Soil.Class", "Soil.LayerThickness",
    "Veg.OverstoryClass", "Veg.UnderstoryClass", "Soil Table", "Overstory Table",
    "Understory Table", "Channel Class Table", "Channel Segment Table", "Channel Cell Table",
    "Station Table"
  };
  const char *Used[NCLASSMAPS + 12];
  char Path[BUFSIZE + 1];
  FILE *File;
  int n = 0;
  int i;

  for (i = 0; i < NCLASSMAPS; i++) {
    if (!Options->FromTable[i])
      Used[n++] = ClassMaps[i].Name;
  }
  if (Options->Shading) {
    Used[n++] = "Terrain.HorizonAngle";
    Used[n++] = "Terrain.SkyView";
  }
  if (Options->Prism)
    Used[n++] = "Met.PrecipPattern";
  if (Options->SnowPattern)
    Used[n++] = "Met.SnowPattern";
  if (Options->MeltScale)
    Used[n++] = "Snow.MeltScale";
  if (Options->WindScale)
    Used[n++] = "Met.WindScale";
  if (Options->LakeDynamics) {
    Used[n++] = "Lake.ID";
    Used[n++] = "Lake Table";
  }
  if (Options->StomatalScheme == MEDLYN)
    Used[n++] = "CO2 File";
  if (PointTable && GetConfigValue(Input, "INPUTS", "Point Table", Path) && !IsEmptyStr(Path))
    Used[n++] = "Point Table";

  for (i = 0; i < (int) (sizeof(Always) / sizeof(Always[0])) + n; i++) {
    RequireConfigValue(Input, "INPUTS",
                       (i < (int) (sizeof(Always) / sizeof(Always[0]))) ? Always[i] :
                       Used[i - (int) (sizeof(Always) / sizeof(Always[0]))], Path);
    if (!(File = fopen(Path, "rb")))
      ReportErrorCode(104, "%s: cannot open", Path);
    fclose(File);
  }
  if (!IsEmptyStr(Options->InitialSWEPath)) {
    if (!(File = fopen(Options->InitialSWEPath, "rb")))
      ReportErrorCode(104, "%s: cannot open", Options->InitialSWEPath);
    fclose(File);
  }

  RequireConfigValue(Input, "INPUTS", "Station Data Directory", Path);
  if (!IsDirectory(Path))
    ReportErrorCode(104, "%s: cannot open", Path);
}

/*****************************************************************************
  ValidateVegMaps()

  The eight vegetation maps as given or formed from the class tables (the
  input maps or one vegetation update): the classes of every mask cell rows
  of their tables or 0; on overstory cells (overstory class > 0) no no-data
  in the overstory maps, the overstory cover in (0, 1], the overstory LAI
  and height positive; on cells with an understory class no no-data in the
  understory cover, which is in [0, 1]; where that cover is positive no
  no-data in the understory LAI and height, both positive. Suffix follows
  the map names in messages
*****************************************************************************/
void ValidateVegMaps(OPTIONSTRUCT *Options, MAPSIZE *Map, unsigned char *Mask, VEGMAPS *V,
  VEGTABLES *Veg, const char *Suffix)
{
  int x, y, i;
  char Name[BUFSIZE + 1];
  char What[BUFSIZE + 1];
  FAILURES F;

#define VMASK(i) (Mask[i] == 1)
#define VOVER(i) (V->OverstoryClass[i] > 0)
#define VUNDER(i) (V->UnderstoryClass[i] > 0)
#define VCOVERED(i) (V->UnderstoryClass[i] > 0 && V->UnderstoryCover[i] > 0.0)
#define VFAIL(Array, Where, Condition)                              \
  memset(&F, 0, sizeof(F));                                         \
  for (y = 0, i = 0; y < Map->NY; y++)                              \
    for (x = 0; x < Map->NX; x++, i++)                              \
      if (VMASK(i) && (Where) && (Condition))                       \
        Fail(&F, y, x, (Array)[i], 0);
#define CLASSCHECK(Array, MapName, Condition, Text)                 \
  VFAIL(Array, TRUE, Condition)                                     \
  snprintf(What, sizeof(What), "%s%s", MapName, Suffix);            \
  ReportCells(&F, What, Text, MapName, NULL);
#define VCHECK(Array, Where, m, Condition, Text, ValueName)          \
  VFAIL(Array, Where, Condition)                                    \
  ClassMapName(What, sizeof(What), Options, m, 0, 0, Suffix);       \
  ReportCells(&F, What, Text, ValueName, NULL);

  /* Classes: no no-data, rows of their tables or 0 */
  CLASSCHECK(V->OverstoryClass, "Veg.OverstoryClass", V->OverstoryClass[i] == MAP_NODATA_BYTE,
             "no-data");
  snprintf(Name, sizeof(Name), "class not in 0..%d (the overstory table rows)", Veg->NOver);
  CLASSCHECK(V->OverstoryClass, "Veg.OverstoryClass", V->OverstoryClass[i] > Veg->NOver, Name);
  CLASSCHECK(V->UnderstoryClass, "Veg.UnderstoryClass",
             V->UnderstoryClass[i] == MAP_NODATA_BYTE, "no-data");
  snprintf(Name, sizeof(Name), "class not in 0..%d (the understory table rows)", Veg->NUnder);
  CLASSCHECK(V->UnderstoryClass, "Veg.UnderstoryClass", V->UnderstoryClass[i] > Veg->NUnder,
             Name);

  /* No no-data where a value is used */
  VCHECK(V->OverstoryCover, VOVER(i), CM_OVERCOVER,
         V->OverstoryCover[i] == MAP_NODATA_FLOAT, "no-data", "Veg.OverstoryCover");
  VCHECK(V->OverstoryLAI, VOVER(i), CM_OVERLAI,
         V->OverstoryLAI[i] == MAP_NODATA_FLOAT, "no-data", "Veg.OverstoryLAI");
  VCHECK(V->OverstoryHeight, VOVER(i), CM_OVERHEIGHT,
         V->OverstoryHeight[i] == MAP_NODATA_FLOAT, "no-data", "Veg.OverstoryHeight");
  VCHECK(V->UnderstoryCover, VUNDER(i), CM_UNDERCOVER,
         V->UnderstoryCover[i] == MAP_NODATA_FLOAT, "no-data", "Veg.UnderstoryCover");
  VCHECK(V->UnderstoryLAI, VCOVERED(i), CM_UNDERLAI,
         V->UnderstoryLAI[i] == MAP_NODATA_FLOAT, "no-data", "Veg.UnderstoryLAI");
  VCHECK(V->UnderstoryHeight, VCOVERED(i), CM_UNDERHEIGHT,
         V->UnderstoryHeight[i] == MAP_NODATA_FLOAT, "no-data", "Veg.UnderstoryHeight");

  /* Overstory cover in (0, 1], overstory LAI and height positive; the
     understory cover in [0, 1], the understory LAI and height positive
     where it covers */
  VCHECK(V->OverstoryCover, VOVER(i), CM_OVERCOVER,
         !(V->OverstoryCover[i] > 0.0 && V->OverstoryCover[i] <= 1.0),
         "cover not in (0, 1] on overstory cells", "cover");
  VCHECK(V->OverstoryLAI, VOVER(i), CM_OVERLAI, !(V->OverstoryLAI[i] > 0.0),
         "LAI <= 0 on overstory cells", "LAI");
  VCHECK(V->OverstoryHeight, VOVER(i), CM_OVERHEIGHT, !(V->OverstoryHeight[i] > 0.0),
         "height <= 0 on overstory cells", "height");
  VCHECK(V->UnderstoryCover, VUNDER(i), CM_UNDERCOVER,
         !(V->UnderstoryCover[i] >= 0.0 && V->UnderstoryCover[i] <= 1.0),
         "cover not in [0, 1] on understory cells", "cover");
  VCHECK(V->UnderstoryLAI, VCOVERED(i), CM_UNDERLAI, !(V->UnderstoryLAI[i] > 0.0),
         "LAI <= 0 on understory cells with cover > 0", "LAI");
  VCHECK(V->UnderstoryHeight, VCOVERED(i), CM_UNDERHEIGHT,
         !(V->UnderstoryHeight[i] > 0.0), "height <= 0 on understory cells with cover > 0",
         "height");

#undef VMASK
#undef VOVER
#undef VUNDER
#undef VCOVERED
#undef VFAIL
#undef CLASSCHECK
#undef VCHECK
}

/*****************************************************************************
  ValidateInputs()

  The input maps as given (before the parameters) and the station normals
*****************************************************************************/
void ValidateInputs(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, INPUTMAPS *Raw,
  LAYER *Soil, SOILTABLE *SType, VEGTABLES *Veg, int NStats,
  METLOCATION *Stat)
{
  int NCells = Map->NY * Map->NX;
  int L = Map->SoilLayers;
  int x, y, i, k, Ones;

  float TotalRootDepth, Fc, Wp, Por;
  float *Layer;
  char Name[BUFSIZE + 1];
  char Path[BUFSIZE + 1];
  FAILURES F;

  Message("\nChecking the inputs\n");

#define INMASK(i) (Raw->Mask[i] == 1)
#define NODATA(v) ((v) == MAP_NODATA_FLOAT)
#define SOILOK(i) (Raw->SoilClass[i] >= 1 && Raw->SoilClass[i] <= Soil->NTypes)

  /* The mask is 0 or 1 everywhere, with at least one 1 */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0, Ones = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Raw->Mask[i] != 0 && Raw->Mask[i] != 1)
        Fail(&F, y, x, Raw->Mask[i], 0);
      if (Raw->Mask[i] == 1)
        Ones++;
    }
  }
  ReportCells(&F, "Terrain.Mask", "value other than 0 or 1", "Terrain.Mask", NULL);
  if (Ones == 0)
    ReportErrorCode(108, "Terrain.Mask: no cell has the value 1");

  /* The elevation has a value everywhere */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (NODATA(Raw->Elevation[i]))
        Fail(&F, y, x, Raw->Elevation[i], 0);
  ReportCells(&F, "Terrain.Elevation", "no-data", "Terrain.Elevation", NULL);

  /* No no-data where a value is used: every mask cell, the overstory
     maps on overstory cells */
#define CHECK_NODATA(Array, Where, MapName)                         \
  memset(&F, 0, sizeof(F));                                         \
  for (y = 0, i = 0; y < Map->NY; y++)                              \
    for (x = 0; x < Map->NX; x++, i++)                              \
      if ((Where) && NODATA((Array)[i]))                            \
        Fail(&F, y, x, (Array)[i], 0);                              \
  ReportCells(&F, MapName, "no-data", MapName, NULL);
  if (Options->HorizonAngle.Ptr != NULL) {
    if (!(Layer = calloc(NCells, sizeof(float))))
      ReportError("ValidateInputs", 1);
    for (k = 0; k < MAP_NAZIMUTHS; k++) {
      ReadMapLayer(&(Options->HorizonAngle), k, Layer, Map);
      snprintf(Name, sizeof(Name), "Terrain.HorizonAngle.Az%03d", k * 360 / MAP_NAZIMUTHS);
      CHECK_NODATA(Layer, INMASK(i), Name);
    }
    free(Layer);
  }
  if (Raw->SkyView != NULL) {
    CHECK_NODATA(Raw->SkyView, INMASK(i), "Terrain.SkyView");
  }
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && Raw->SoilClass[i] == MAP_NODATA_BYTE)
        Fail(&F, y, x, Raw->SoilClass[i], 0);
  ReportCells(&F, "Soil.Class", "no-data", "Soil.Class", NULL);
  ClassMapName(Name, sizeof(Name), Options, CM_DEPTH, 0, L, "");
  CHECK_NODATA(Raw->Depth, INMASK(i), Name);
  ClassMapName(Name, sizeof(Name), Options, CM_KSLAT, 0, L, "");
  CHECK_NODATA(Raw->LateralConductivity, INMASK(i), Name);
  ClassMapName(Name, sizeof(Name), Options, CM_EXPDEC, 0, L, "");
  CHECK_NODATA(Raw->ExponentialDecrease, INMASK(i), Name);
  for (k = 0; k <= L; k++) {
    Layer = Raw->Porosity + k * NCells;
    ClassMapName(Name, sizeof(Name), Options, CM_POROSITY, k, L, "");
    CHECK_NODATA(Layer, INMASK(i), Name);
  }
  for (k = 0; k <= L; k++) {
    Layer = Raw->FieldCapacity + k * NCells;
    ClassMapName(Name, sizeof(Name), Options, CM_FCAP, k, L, "");
    CHECK_NODATA(Layer, INMASK(i), Name);
  }
  for (k = 0; k < L; k++) {
    Layer = Raw->LayerThickness + k * NCells;
    SoilLayerName(Name, sizeof(Name), "Soil.LayerThickness", k, L);
    CHECK_NODATA(Layer, INMASK(i), Name);
  }
  if (Options->Prism) {
    if (!(Layer = calloc(NCells, sizeof(float))))
      ReportError("ValidateInputs", 1);
    for (k = 0; k < 12; k++) {
      ReadMapLayer(&(Options->PrecipPattern), k, Layer, Map);
      snprintf(Name, sizeof(Name), "Met.PrecipPattern.%s", MonthNames[k]);
      CHECK_NODATA(Layer, INMASK(i), Name);
    }
    free(Layer);
  }
  if (Options->SnowPattern) {
    CHECK_NODATA(Raw->SnowPattern, INMASK(i), "Met.SnowPattern");
  }
  if (Options->MeltScale) {
    CHECK_NODATA(Raw->MeltScale, INMASK(i), "Snow.MeltScale");
  }
  if (Options->WindScale) {
    CHECK_NODATA(Raw->WindScale, INMASK(i), "Met.WindScale");
  }
  if (Options->LakeDynamics) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        if (INMASK(i) && Raw->LakeID[i] == MAP_NODATA_BYTE)
          Fail(&F, y, x, Raw->LakeID[i], 0);
    ReportCells(&F, "Lake.ID", "no-data", "Lake.ID", NULL);
  }
  if (Raw->InitialSWE != NULL) {
    CHECK_NODATA(Raw->InitialSWE, INMASK(i), "Snow.InitialSWE");
    free(Raw->InitialSWE);
    Raw->InitialSWE = NULL;
  }

  /* Classes in the tables */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !SOILOK(i))
        Fail(&F, y, x, Raw->SoilClass[i], 0);
  snprintf(Name, sizeof(Name), "class not in 1..%d (the soil table rows)", Soil->NTypes);
  ReportCells(&F, "Soil.Class", Name, "Soil.Class", NULL);

  /* The vegetation maps */
  ValidateVegMaps(Options, Map, Raw->Mask, &(Raw->Veg), Veg, "");

  /* Root zone layers thicker than 0.001 m */
  for (k = 0; k < L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        if (INMASK(i) && !(Raw->LayerThickness[k * NCells + i] > 0.001))
          Fail(&F, y, x, Raw->LayerThickness[k * NCells + i], 0);
    SoilLayerName(Name, sizeof(Name), "Soil.LayerThickness", k, L);
    ReportCells(&F, Name, "thickness <= 0.001", "thickness", NULL);
  }

  /* The depth map at least 0.01 m below the root zone */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INMASK(i)) {
        TotalRootDepth = 0.0;
        for (k = 0; k < L; k++)
          TotalRootDepth += Raw->LayerThickness[k * NCells + i];
        if (Raw->Depth[i] < (TotalRootDepth + 0.01))
          Fail(&F, y, x, Raw->Depth[i], TotalRootDepth);
      }
    }
  }
  ClassMapName(Name, sizeof(Name), Options, CM_DEPTH, 0, L, "");
  ReportCells(&F, Name, "depth < root zone depth + 0.01", "depth", "root zone depth");

  /* Conductivity and exponential decrease, and their scales */
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !(Raw->LateralConductivity[i] > 0.0))
        Fail(&F, y, x, Raw->LateralConductivity[i], 0);
  ClassMapName(Name, sizeof(Name), Options, CM_KSLAT, 0, L, "");
  ReportCells(&F, Name, "value <= 0", "Soil.LateralConductivity", NULL);
  if (!(SOIL_KSAT_ADJ > 0.0))
    ReportErrorCode(108, "[PARAMETERS] Lateral Conductivity Scale = %g: must be > 0",
                    SOIL_KSAT_ADJ);
  memset(&F, 0, sizeof(F));
  for (y = 0, i = 0; y < Map->NY; y++)
    for (x = 0; x < Map->NX; x++, i++)
      if (INMASK(i) && !(Raw->ExponentialDecrease[i] >= 0.0))
        Fail(&F, y, x, Raw->ExponentialDecrease[i], 0);
  ClassMapName(Name, sizeof(Name), Options, CM_EXPDEC, 0, L, "");
  ReportCells(&F, Name, "value < 0", "Soil.ExponentialDecrease", NULL);
  if (!(SOIL_EXPDEC_ADJ >= 0.0))
    ReportErrorCode(108, "[PARAMETERS] Exponential Decrease Scale = %g: must be >= 0",
                    SOIL_EXPDEC_ADJ);

  /* Field capacity map at least 0.01 above the wilting point */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Fc = Raw->FieldCapacity[k * NCells + i];
          Wp = SType[Raw->SoilClass[i] - 1].WP[k];
          if (Fc < (Wp + 0.01))
            Fail(&F, y, x, Fc, Wp);
        }
      }
    }
    ClassMapName(Name, sizeof(Name), Options, CM_FCAP, k, L, "");
    ReportCells(&F, Name, "field capacity < wilting point + 0.01", "field capacity",
                "wilting point");
  }

  /* Porosity map at least 0.01 above the field capacity map */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Por = Raw->Porosity[k * NCells + i];
          Fc = Raw->FieldCapacity[k * NCells + i];
          if (Por < (Fc + 0.01))
            Fail(&F, y, x, Por, Fc);
        }
      }
    }
    ClassMapName(Name, sizeof(Name), Options, CM_POROSITY, k, L, "");
    ReportCells(&F, Name, "porosity < field capacity + 0.01", "porosity", "field capacity");
  }

  /* Porosity map at most 1 */
  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (INMASK(i)) {
          Por = Raw->Porosity[k * NCells + i];
          if (Por > 1.0)
            Fail(&F, y, x, Por, 0);
        }
      }
    }
    ClassMapName(Name, sizeof(Name), Options, CM_POROSITY, k, L, "");
    ReportCells(&F, Name, "porosity > 1", "porosity", NULL);
  }

  /* The scales of the overstory and the understory cover positive, the
     leafless LAI fraction of every deciduous overstory and understory class
     positive */
  if (!(VEG_FC_ADJ > 0.0))
    ReportErrorCode(108, "[PARAMETERS] Overstory Cover Scale = %g: must be > 0", VEG_FC_ADJ);
  if (!(VEG_UNDER_FC_ADJ > 0.0))
    ReportErrorCode(108, "[PARAMETERS] Understory Cover Scale = %g: must be > 0", VEG_UNDER_FC_ADJ);
  RequireConfigValue(Input, "INPUTS", "Overstory Table", Path);
  for (i = 0; i < Veg->NOver; i++) {
    if (!(Veg->Over[i].Story.MinLAIFraction > 0.0))
      ReportErrorCode(108, "%s: line %d: MinLAIFraction = %g must be > 0", Path, i + 3,
                      Veg->Over[i].Story.MinLAIFraction);
  }
  RequireConfigValue(Input, "INPUTS", "Understory Table", Path);
  for (i = 0; i < Veg->NUnder; i++) {
    if (!(Veg->Under[i].Story.MinLAIFraction > 0.0))
      ReportErrorCode(108, "%s: line %d: MinLAIFraction = %g must be > 0", Path, i + 3,
                      Veg->Under[i].Story.MinLAIFraction);
  }

  /* Lake IDs in the lake table (Lake Dynamics) */
  if (Options->LakeDynamics) {
    memset(&F, 0, sizeof(F));
    for (y = 0, i = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++, i++)
        if (INMASK(i) && Raw->LakeID[i] > Map->NumLakes)
          Fail(&F, y, x, Raw->LakeID[i], 0);
    snprintf(Name, sizeof(Name), "lake not in 0..%d (the lake table rows)", Map->NumLakes);
    ReportCells(&F, "Lake.ID", Name, "Lake.ID", NULL);
  }

  /* Station normals positive */
  RequireConfigValue(Input, "INPUTS", "Station Table", Path);
  for (i = 0; i < NStats; i++) {
    for (k = 0; Options->Prism && k < 12; k++) {
      if (!(Stat[i].PrismPrecip[k] > 0.0))
        ReportErrorCode(108, "%s: line %d: PrecipNormal_%s = %g must be > 0", Path, i + 3,
                        MonthNames[k], Stat[i].PrismPrecip[k]);
    }
    if (Options->SnowPattern && !(Stat[i].SnowPatternBase > 0.0))
      ReportErrorCode(108, "%s: line %d: SnowNormal = %g must be > 0", Path, i + 3,
                      Stat[i].SnowPatternBase);
  }

#undef CHECK_NODATA
#undef SOILOK
#undef NODATA
#undef INMASK
}

/*****************************************************************************
  ValidateParameters()

  The porosity after the field capacity and porosity parameters at
  most 1
*****************************************************************************/
void ValidateParameters(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap,
  SOILPIX **SoilMap)
{
  int x, y, k;
  int L = Map->SoilLayers;
  char Name[BUFSIZE + 1];
  FAILURES F;

  for (k = 0; k <= L; k++) {
    memset(&F, 0, sizeof(F));
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        if (INBASIN(TopoMap[y][x].Mask) && (double) SoilMap[y][x].Porosity[k] > 1.0)
          Fail(&F, y, x, SoilMap[y][x].Porosity[k], SoilMap[y][x].FCap[k]);
      }
    }
    ClassMapName(Name, sizeof(Name), Options, CM_POROSITY, k, L, "");
    strcat(Name, " x Porosity Scale");
    ReportCells(&F, Name, "porosity > 1", "porosity", "field capacity");
  }
}

/*****************************************************************************
  ValidateWeights()

  The station weights of every mask cell (bytes) sum to more than 0
*****************************************************************************/
void ValidateWeights(MAPSIZE *Map, TOPOPIX **TopoMap, uchar ***MetWeights, int NStats)
{
  int x, y, i;
  long Sum;
  FAILURES F;

  memset(&F, 0, sizeof(F));
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        for (i = 0, Sum = 0; i < NStats; i++)
          Sum += MetWeights[y][x][i];
        if (Sum <= 0)
          Fail(&F, y, x, Sum, 0);
      }
    }
  }
  ReportCells(&F, "station weights", "weights sum to 0", "sum", NULL);
}

/*****************************************************************************
  PrintClassSummary()

  The overstory, understory and soil classes in the basin and their
  fractions
*****************************************************************************/
void PrintClassSummary(OPTIONSTRUCT *Options, VEGTABLES *Veg, LAYER Soil,
  SOILTABLE *SType, MAPSIZE *Map, TOPOPIX **TopoMap, VEGPIX **VegMap, SOILPIX **SoilMap)
{
  int y, x, i;
  int *ocount = NULL, *ucount = NULL, *scount = NULL;
  double *usnow = NULL;
  int npixels;

  if (!(ocount = calloc(Veg->NOver + 1, sizeof(int))) ||
      !(ucount = calloc(Veg->NUnder + 1, sizeof(int))) ||
      !(scount = calloc(Soil.NTypes, sizeof(int))) ||
      !(usnow = calloc(Veg->NOver + 1, sizeof(double))))
    ReportError("PrintClassSummary", 1);

  npixels = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
        ocount[VegMap[y][x].OverClass]++;
        ucount[VegMap[y][x].UnderClass]++;
        usnow[VegMap[y][x].OverClass] += VegMap[y][x].USnow;
        scount[SoilMap[y][x].Soil - 1]++;
        npixels++;
      }
    }
  }
  Message("\nBasin has %d active pixels \n", npixels);

  Message("\nThe following OVERSTORY types are in the current basin \n");
  for (i = 0; i <= Veg->NOver; i++) {
    if (ocount[i] > 0)
      Message("Class # %d of Type: %s has fraction basin area: %5.3f\n",
        i, (i > 0) ? Veg->Over[i - 1].Desc : "(no overstory)", (float)ocount[i] / (float)npixels);
  }

  Message("\nThe following UNDERSTORY types are in the current basin \n");
  for (i = 0; i <= Veg->NUnder; i++) {
    if (ucount[i] > 0)
      Message("Class # %d of Type: %s has fraction basin area: %5.3f\n",
        i, (i > 0) ? Veg->Under[i - 1].Desc : "(no understory)",
        (float)ucount[i] / (float)npixels);
  }

  Message("\nThe following SOIL types are in the current basin \n");
  for (i = 0; i < Soil.NTypes; i++)
    if (scount[i] > 0)
      Message("Class # %d of Type: %s has fraction basin area: %5.3f\n",
        i + 1, SType[i].Desc, (float)scount[i] / (float)npixels);

  Message("\nSome estimates for current vegetation specification\n");
  for (i = 0; i <= Veg->NOver; i++) {
    if (ocount[i] > 0) {
      Message("\nOverstory Type: %s\n", (i > 0) ? Veg->Over[i - 1].Desc : "(no overstory)");
      Message("2meter    wind speed fraction of ref level %1.3f (mean of the cells)\n",
              usnow[i] / ocount[i]);
    }
  }

  if (Veg->NOver > 0) {
    Message("\nShortwave extinction coeff. per unit LAI for the overstory types: \n\n");
    for (i = 0; i < Veg->NOver; i++)
      Message("%s: %.3f\n", Veg->Over[i].Desc, Veg->Over[i].ExtnCoeff);
    Message("\n");
  }

  free(ocount);
  free(ucount);
  free(scount);
  free(usnow);
}

/*****************************************************************************
  ValidateChannels()

  The channel class, segment and cell tables and the lakes, after the
  parameters (the channel cuts are compared with the soil depth)
*****************************************************************************/
void ValidateChannels(OPTIONSTRUCT *Options, MAPSIZE *Map, TOPOPIX **TopoMap,
                      SOILPIX **SoilMap, CHANNELTABLES *Tables, LAKETABLE *LType)
{
  int i, j, k, MaxOrder, Found, InMask;
  int x, y;
  const char *c;

  /* Segment and channel class IDs positive and unique */
  for (i = 0; i < Tables->NSegments; i++) {
    if (Tables->SegmentID[i] < 1 || Tables->SegmentID[i] > 65535)
      ReportErrorCode(109, "%s: line %d: Segment %d is not in 1..65535", Tables->SegmentPath,
                      i + 3, Tables->SegmentID[i]);
    for (j = 0; j < i; j++) {
      if (Tables->SegmentID[j] == Tables->SegmentID[i])
        ReportErrorCode(109, "%s: line %d: Segment %d repeated", Tables->SegmentPath, i + 3,
                        Tables->SegmentID[i]);
    }
  }
  for (i = 0; i < Tables->NClasses; i++) {
    if (Tables->ClassID[i] < 1 || Tables->ClassID[i] > 65535)
      ReportErrorCode(109, "%s: line %d: Class %d is not in 1..65535", Tables->ClassPath,
                      i + 3, Tables->ClassID[i]);
    for (j = 0; j < i; j++) {
      if (Tables->ClassID[j] == Tables->ClassID[i])
        ReportErrorCode(109, "%s: line %d: Class %d repeated", Tables->ClassPath, i + 3,
                        Tables->ClassID[i]);
    }
  }

  /* The class of every segment is in the class table */
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0, Found = FALSE; j < Tables->NClasses && !Found; j++)
      Found = (Tables->ClassID[j] == Tables->Class[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Class %d not in %s", Tables->SegmentPath, i + 3,
                      Tables->Class[i], Tables->ClassPath);
  }

  /* Every outlet 0 or a segment */
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0, Found = (Tables->Outlet[i] == 0); j < Tables->NSegments && !Found; j++)
      Found = (Tables->SegmentID[j] == Tables->Outlet[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Outlet %d is not a segment", Tables->SegmentPath,
                      i + 3, Tables->Outlet[i]);
  }

  /* Orders 1..max all present, every segment drains to a higher order */
  for (i = 0, MaxOrder = 0; i < Tables->NSegments; i++) {
    if (Tables->Order[i] < 1)
      ReportErrorCode(109, "%s: line %d: Order %d must be >= 1", Tables->SegmentPath, i + 3,
                      Tables->Order[i]);
    if (Tables->Order[i] > MaxOrder)
      MaxOrder = Tables->Order[i];
  }
  for (k = 1; k <= MaxOrder; k++) {
    for (i = 0, Found = FALSE; i < Tables->NSegments && !Found; i++)
      Found = (Tables->Order[i] == k);
    if (!Found)
      ReportErrorCode(109, "%s: no segment of order %d", Tables->SegmentPath, k);
  }
  for (i = 0; i < Tables->NSegments; i++) {
    for (j = 0; j < Tables->NSegments && Tables->Outlet[i] != 0; j++) {
      if (Tables->SegmentID[j] == Tables->Outlet[i] && Tables->Order[j] <= Tables->Order[i])
        ReportErrorCode(109, "%s: line %d: segment %d (order %d) drains to segment %d (order %d)",
                        Tables->SegmentPath, i + 3, Tables->SegmentID[i], Tables->Order[i],
                        Tables->SegmentID[j], Tables->Order[j]);
    }
  }

  /* Positive slopes, lengths, widths, bank heights, roughness */
  for (i = 0; i < Tables->NSegments; i++) {
    if (!(Tables->Slope[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Slope = %g must be > 0", Tables->SegmentPath, i + 3,
                      Tables->Slope[i]);
    if (!(Tables->Length[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Length = %g must be > 0", Tables->SegmentPath, i + 3,
                      Tables->Length[i]);
  }
  for (i = 0; i < Tables->NClasses; i++) {
    if (!(Tables->Width[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Width = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->Width[i]);
    if (!(Tables->BankHeight[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: BankHeight = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->BankHeight[i]);
    if (!(Tables->Friction[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: ManningsN = %g must be > 0", Tables->ClassPath, i + 3,
                      Tables->Friction[i]);
  }

  /* Channel cells inside the grid and the mask */
  for (i = 0; i < Tables->NCells; i++) {
    y = Tables->Row[i] - 1;
    x = Tables->Column[i] - 1;
    InMask = (y >= 0 && y < Map->NY && x >= 0 && x < Map->NX && INBASIN(TopoMap[y][x].Mask));
    if (!InMask)
      ReportErrorCode(109, "%s: line %d: row %d, column %d is outside the mask", Tables->CellPath,
                      i + 3, Tables->Row[i], Tables->Column[i]);
  }

  /* The segments of the channel cells exist, lengths positive (the lateral
     inflow of a cell is split by length), cuts not negative */
  for (i = 0; i < Tables->NCells; i++) {
    for (j = 0, Found = FALSE; j < Tables->NSegments && !Found; j++)
      Found = (Tables->SegmentID[j] == Tables->Segment[i]);
    if (!Found)
      ReportErrorCode(109, "%s: line %d: Segment %d is not in %s", Tables->CellPath, i + 3,
                      Tables->Segment[i], Tables->SegmentPath);
    if (!(Tables->CellLength[i] > 0.0))
      ReportErrorCode(109, "%s: line %d: Length = %g must be > 0", Tables->CellPath, i + 3,
                      Tables->CellLength[i]);
    if (!(Tables->CutWidth[i] >= 0.0))
      ReportErrorCode(109, "%s: line %d: CutWidth = %g must be >= 0", Tables->CellPath, i + 3,
                      Tables->CutWidth[i]);
    if (!(Tables->CutHeight[i] >= 0.0))
      ReportErrorCode(109, "%s: line %d: CutHeight = %g must be >= 0", Tables->CellPath, i + 3,
                      Tables->CutHeight[i]);
  }

  /* Channel cuts not deeper than the soil (after the depth offset) */
  for (i = 0; i < Tables->NCells; i++) {
    y = Tables->Row[i] - 1;
    x = Tables->Column[i] - 1;
    if (Tables->CutHeight[i] > SoilMap[y][x].Depth)
      ReportErrorCode(109, "%s: line %d: CutHeight %g is deeper than the soil (%g m) at row %d, column %d",
                      Tables->CellPath, i + 3, Tables->CutHeight[i], SoilMap[y][x].Depth,
                      Tables->Row[i], Tables->Column[i]);
  }

  /* Save names of letters, digits, _ . - and unique */
  for (i = 0; i < Tables->NSegments; i++) {
    if (Tables->SaveName[i][0] == '\0')
      continue;
    for (c = Tables->SaveName[i]; *c != '\0'; c++) {
      if (!((*c >= 'A' && *c <= 'Z') || (*c >= 'a' && *c <= 'z') || (*c >= '0' && *c <= '9') ||
            *c == '_' || *c == '.' || *c == '-'))
        ReportErrorCode(109, "%s: line %d: SaveName \"%s\" has a character outside A-Z a-z 0-9 _ . -",
                        Tables->SegmentPath, i + 3, Tables->SaveName[i]);
    }
    for (j = 0; j < i; j++) {
      if (strcmp(Tables->SaveName[j], Tables->SaveName[i]) == 0)
        ReportErrorCode(109, "%s: line %d: SaveName \"%s\" repeated", Tables->SegmentPath, i + 3,
                        Tables->SaveName[i]);
    }
  }

  /* Lake outlets are segments, every lake has a cell in the mask */
  if (Options->LakeDynamics) {
    for (i = 0; i < Map->NumLakes; i++) {
      for (j = 0, Found = FALSE; j < Tables->NSegments && !Found; j++)
        Found = (Tables->SegmentID[j] == LType[i].Outlet);
      if (!Found)
        ReportErrorCode(109, "%s: line %d: OutletSegment %d is not a segment",
                        Options->LakeTablePath, i + 3, LType[i].Outlet);
    }
    for (i = 0; i < Map->NumLakes; i++) {
      for (y = 0, Found = FALSE; y < Map->NY && !Found; y++)
        for (x = 0; x < Map->NX && !Found; x++)
          Found = (INBASIN(TopoMap[y][x].Mask) && TopoMap[y][x].LakeID == i + 1);
      if (!Found)
        ReportErrorCode(109, "Lake.ID: lake %d has no cell in the mask", i + 1);
    }
  }
}
