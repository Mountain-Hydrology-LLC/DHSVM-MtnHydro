/*
 * DESCRIPTION:  The output requests of [OUTPUT]: variable lists expanded to
 *               columns, and their checks (exit 111): names, the output
 *               kinds that accept them, options, classes, points, saved
 *               segments.
 *
 *               A list holds names separated by blanks: a variable (all its
 *               layers), a variable layer (<Name>.<token>) or a group
 *               wildcard <Group>.* (every variable of the group that the
 *               output kind accepts and the options make available). The
 *               list of a kind is its keys (<Key>, <Key> 2, ...) followed by
 *               the request table's variables flagged for it. Columns keep
 *               the order of the names; a variable's layers follow its layer
 *               set. A column already in the list is not added again, and a
 *               wildcard that matches nothing adds nothing.
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "output.h"
#include "registry.h"

static const char *Unavailable(int Var, OPTIONSTRUCT *Options);
static void AddColumns(OUTLIST *List, int *Capacity, int Var, int First, int Count,
                       int SoilLayers);

/*****************************************************************************
  Unavailable()

  What a variable needs that the run does not have, NULL if it is available
*****************************************************************************/
static const char *Unavailable(int Var, OPTIONSTRUCT *Options)
{
  const char *Name = RegistryRow(Var)->Name;

#ifdef SNOW_ONLY
  if (strncmp(Name, "Stream.", 7) == 0 || strncmp(Name, "Balance.", 8) == 0 ||
      strncmp(Name, "Channel.", 8) == 0 || strncmp(Name, "Lake.", 5) == 0)
    return "the channel routing, which the snow-only build does not run";
#endif
  if (strcmp(Name, "Soil.InfiltrationAccumulated") == 0 && Options->Infiltration != DYNAMIC)
    return "Infiltration = DYNAMIC";
  if (strncmp(Name, "Lake.", 5) == 0 && !Options->LakeDynamics)
    return "Lake Dynamics = TRUE";
  if ((strcmp(Name, "Soil.Ice") == 0 || strcmp(Name, "Soil.Temperature") == 0 ||
       strcmp(Name, "Soil.NodeTemperature") == 0 || strcmp(Name, "Soil.NodeIce") == 0 ||
       strcmp(Name, "Soil.SurfaceTemperature") == 0 || strcmp(Name, "Soil.FrostDepth") == 0 ||
       strcmp(Name, "Soil.FrostTop") == 0 ||
       strcmp(Name, "Soil.GroundHeatFlux") == 0 || strcmp(Name, "Soil.SensibleHeatFlux") == 0) &&
      !Options->SoilFreezing)
    return "Soil Freezing = TRUE";
  if ((strcmp(Name, "Canopy.Assimilation") == 0 || strcmp(Name, "Canopy.Acclimation") == 0) &&
      Options->StomatalScheme != MEDLYN)
    return "Stomatal Scheme = MEDLYN";
  return NULL;
}

/* Adds layers First .. First + Count - 1 of a variable to the columns, except
   the columns already in the list */
static void AddColumns(OUTLIST *List, int *Capacity, int Var, int First, int Count,
                       int SoilLayers)
{
  char Name[OUT_NAMESIZE];
  int i, k;

  for (k = First; k < First + Count; k++) {
    RegistryLayerName(Var, k, SoilLayers, Name, sizeof(Name));
    for (i = 0; i < List->N; i++) {
      if (strcmp(List->Vars[i].Name, Name) == 0)
        break;
    }
    if (i < List->N)
      continue;
    if (List->N == *Capacity) {
      *Capacity = 2 * *Capacity + 16;
      if (!(List->Vars = (OUTVAR *) realloc(List->Vars, *Capacity * sizeof(OUTVAR))))
        ReportError("AddColumns", 1);
    }
    List->Vars[List->N].Var = Var;
    List->Vars[List->N].Layer = k;
    List->Vars[List->N].Code = ValueCode(Var);
    List->Vars[List->N].Instant = !RegistryIsInterval(Var);
    strcpy(List->Vars[List->N].Name, Name);
    List->N++;
  }
}

/*****************************************************************************
  ExpandList()

  The columns of a variable list for an output kind (Kind, a KIND_* bit;
  KindName for messages). Check selects the one check whose failures stop
  the run: 0 none (names that fail are left out), 1 the name exists, 2 the
  kind accepts it, 3 the options make it available. Returns the number of
  columns.
*****************************************************************************/
int ExpandList(OUTLIST *List, int Kind, const char *KindName, OPTIONSTRUCT *Options,
               int SoilLayers, int Check)
{
  char Value[BUFSIZE + 1];
  char *Token;
  const char *Need;
  size_t Length;
  int Capacity = 0;
  int Var, Layer, First, Count;
  int i;

  free(List->Vars);
  List->Vars = NULL;
  List->N = 0;
  for (i = 0; i < List->NKeys; i++) {
    strcpy(Value, List->Value[i]);
    for (Token = strtok(Value, " \t"); Token != NULL; Token = strtok(NULL, " \t")) {
      Length = strlen(Token);

      /* A group wildcard */
      if (Length > 2 && strcmp(Token + Length - 2, ".*") == 0) {
        for (Var = 0; Var < RegistryCount(); Var++) {
          if (RegistryGroupMatches(Var, Token, Length - 2) && (RegistryKinds(Var) & Kind) &&
              Unavailable(Var, Options) == NULL)
            AddColumns(List, &Capacity, Var, 0, RegistryLayerCount(Var, SoilLayers),
                       SoilLayers);
        }
        continue;
      }

      /* A variable, or one layer of a variable */
      Var = RegistryFind(Token);
      if (Var >= 0) {
        First = 0;
        Count = RegistryLayerCount(Var, SoilLayers);
      }
      else if (RegistryFindLayer(Token, SoilLayers, &Var, &Layer)) {
        First = Layer;
        Count = 1;
      }
      else {
        if (Check == 1)
          ReportErrorCode(111, "%s: %s is not a DHSVM variable", List->Key[i], Token);
        continue;
      }
      if (!(RegistryKinds(Var) & Kind)) {
        if (Check == 2)
          ReportErrorCode(111, "%s: %s cannot be written to %s", List->Key[i], Token,
                          KindName);
        continue;
      }
      Need = Unavailable(Var, Options);
      if (Need != NULL) {
        if (Check == 3)
          ReportErrorCode(111, "%s: %s needs %s", List->Key[i], Token, Need);
        continue;
      }
      AddColumns(List, &Capacity, Var, First, Count, SoilLayers);
    }
  }
  return List->N;
}

/*****************************************************************************
  ValidateRequests()

  The checks of the output requests in order: checks 1-3 of ExpandList()
  over every variable list, then the classes, the points and the saved
  segments. Leaves the expanded lists and the points' cells in Dump.
*****************************************************************************/
void ValidateRequests(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
                      TOPOPIX **TopoMap, VEGPIX **VegMap, VEGTABLES *Veg, Channel *Streams,
                      DUMPSTRUCT *Dump)
{
  const int Kinds[NOUTLISTS] = { KIND_BASIN, KIND_POINT, KIND_CLASS, KIND_CHANNELS, KIND_MAP };
  const char *KindNames[NOUTLISTS] = { "Basin", "Point", "Class", "Channels", "a map" };
  char Path[BUFSIZE + 1];
  const char *c;
  Channel *Segment;
  CLASSOUT *Class;
  int Found;
  int Check, i, j, k, x, y;
  double Row, Col;

  for (Check = 1; Check <= 3; Check++) {
    for (k = 0; k < NOUTLISTS; k++)
      ExpandList(&(Dump->List[k]), Kinds[k], KindNames[k], Options, Map->SoilLayers, Check);
  }

  /* Class combinations, paired with a Class list: each class a row of its
     table or 0, each combination on a mask cell at Start or given to one by
     a vegetation update before End */
  if (Dump->List[OUT_CLASS].N > 0 && IsEmptyStr(Dump->ClassesText))
    ReportErrorCode(111, "[OUTPUT] Class needs Classes");
  if (Dump->List[OUT_CLASS].N == 0 && !IsEmptyStr(Dump->ClassesText))
    ReportErrorCode(111, "[OUTPUT] Classes needs Class");
  for (j = 0; j < Dump->NClasses; j++) {
    Class = &(Dump->Classes[j]);
    if (Class->Over < 0 || Class->Over > Veg->NOver) {
      GetConfigValue(Input, "INPUTS", "Overstory Table", Path);
      ReportErrorCode(111, "[OUTPUT] Classes: %d_%d: overstory class %d is not in 0..%d "
                      "(the rows of %s)", Class->Over, Class->Under, Class->Over, Veg->NOver,
                      Path);
    }
    if (Class->Under < 0 || Class->Under > Veg->NUnder) {
      GetConfigValue(Input, "INPUTS", "Understory Table", Path);
      ReportErrorCode(111, "[OUTPUT] Classes: %d_%d: understory class %d is not in 0..%d "
                      "(the rows of %s)", Class->Over, Class->Under, Class->Under, Veg->NUnder,
                      Path);
    }
    Found = Class->Present;
    for (y = 0; y < Map->NY && !Found; y++)
      for (x = 0; x < Map->NX && !Found; x++)
        if (INBASIN(TopoMap[y][x].Mask) && VegMap[y][x].OverClass == Class->Over &&
            VegMap[y][x].UnderClass == Class->Under)
          Found = TRUE;
    if (!Found)
      ReportErrorCode(111, "[OUTPUT] Classes: %d_%d has no cell in the mask", Class->Over,
                      Class->Under);
  }

  /* The points: a point table, names, cells inside the mask. Names are
     unique ignoring case, since they name files. The row and column are
     tested before rounding: a cell index rounds into the grid exactly when
     -0.5 < index < size - 0.5 */
  if (Dump->List[OUT_POINT].N > 0) {
    if (!GetConfigValue(Input, "INPUTS", "Point Table", Path) || IsEmptyStr(Path))
      ReportErrorCode(111, "[OUTPUT] Point needs [INPUTS] Point Table");
    for (i = 0; i < Dump->NPoints; i++) {
      POINTOUT *Point = &(Dump->Points[i]);
      for (c = Point->Name; *c != '\0'; c++) {
        if (!(isalnum((unsigned char) *c) || *c == '_' || *c == '-'))
          ReportErrorCode(111, "%s: line %d: point name \"%s\" has a character outside "
                          "A-Z a-z 0-9 _ -", Path, Point->Line, Point->Name);
      }
      for (j = 0; j < i; j++)
        if (IsChoice(Dump->Points[j].Name, Point->Name))
          ReportErrorCode(111, "%s: line %d: point %s repeated", Path, Point->Line, Point->Name);
      Row = ((Map->Yorig - 0.5 * Map->DY) - Point->North) / Map->DY;
      Col = (Point->East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX;
      if (!(Row > -0.5 && Row < Map->NY - 0.5 && Col > -0.5 && Col < Map->NX - 0.5))
        ReportErrorCode(111, "%s: line %d: point %s is outside the mask", Path, Point->Line,
                        Point->Name);
      Point->y = Round(Row);
      Point->x = Round(Col);
      if (!INBASIN(TopoMap[Point->y][Point->x].Mask))
        ReportErrorCode(111, "%s: line %d: point %s is outside the mask", Path, Point->Line,
                        Point->Name);
    }
  }

  /* Saved segments for Streamflow.csv; channel outputs need the
     channel routing */
  if (Dump->Streamflow) {
#ifdef SNOW_ONLY
    ReportErrorCode(111, "[OUTPUT] Streamflow = TRUE needs the channel routing, which the "
                    "snow-only build does not run");
#endif
    for (Segment = Streams; Segment != NULL; Segment = Segment->next)
      if (Segment->record)
        break;
    if (Segment == NULL)
      ReportErrorCode(111, "[OUTPUT] Streamflow = TRUE but no segment has a SaveName");
  }
#ifdef SNOW_ONLY
  if (Dump->NetworkSteps > 0)
    ReportErrorCode(111, "[OUTPUT] Network Dynamics Interval needs the channel routing, which "
                    "the snow-only build does not run");
#endif
}
