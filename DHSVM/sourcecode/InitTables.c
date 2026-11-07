
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "Calendar.h"
#include "constants.h"
#include "fileio.h"
#include "getinit.h"
#include "mapio.h"
#include "tables.h"

static float *TableFloats(TABLE *T, int Row, int *Columns, int N, float *Values);
static int InitOverstoryTable(OVERTABLE **OType, const char *Path, OPTIONSTRUCT *Options,
                              int SoilLayers);
static int InitUnderstoryTable(UNDERTABLE **UType, const char *Path, OPTIONSTRUCT *Options,
                               int SoilLayers);
static void InitCO2Table(const char *Path, TIMESTRUCT *Time);

/* Atmospheric CO2 of the years of the run (MEDLYN) */
static int CO2FirstYear = 0;
static float *CO2OfYear = NULL;

/*******************************************************************************/
/*				  InitTables()                                 */
/*******************************************************************************/
void InitTables(int StepsPerDay, LISTPTR Input, OPTIONSTRUCT *Options,
  MAPSIZE *Map, SOILTABLE **SType, LAYER *Soil, VEGTABLES *Veg,
  LAKETABLE **LType, TIMESTRUCT *Time)
{
  char Path[BUFSIZE + 1];

  Message("Initializing tables\n");

  RequireConfigValue(Input, "INPUTS", "Soil Table", Path);
  Soil->NTypes = InitSoilTable(Options, SType, Path, Soil, Map->SoilLayers, Time);

  RequireConfigValue(Input, "INPUTS", "Overstory Table", Path);
  Veg->NOver = InitOverstoryTable(&(Veg->Over), Path, Options, Map->SoilLayers);

  RequireConfigValue(Input, "INPUTS", "Understory Table", Path);
  Veg->NUnder = InitUnderstoryTable(&(Veg->Under), Path, Options, Map->SoilLayers);

  if (Options->LakeDynamics) {
    RequireConfigValue(Input, "INPUTS", "Lake Table", Path);
    Map->NumLakes = InitLakeTable(LType, Path, Options);
  }
  else {
    GetConfigValue(Input, "INPUTS", "Lake Table", Path);
    Map->NumLakes = 0;
  }

  if (Options->StomatalScheme == MEDLYN) {
    RequireConfigValue(Input, "INPUTS", "CO2 File", Path);
    InitCO2Table(Path, Time);
  }

  InitSatVaporTable();
}

/* Values of the columns Columns[0..N-1] of a row */
static float *TableFloats(TABLE *T, int Row, int *Columns, int N, float *Values)
{
  int k;

  for (k = 0; k < N; k++)
    Values[k] = TableFloat(T, Row, Columns[k]);
  return Values;
}

/* A value of row i, NOT_APPLICABLE where the column is not read (-1) */
static float OptionalValue(TABLE *T, int i, int Column)
{
  return (Column >= 0) ? TableFloat(T, i, Column) : NOT_APPLICABLE;
}

/* The columns of the class column of map m of ClassMaps[] (one column, or
   the soil layers 1..L and Deep): read when the map is formed from the
   table, otherwise known and not read (-1) */
static void ClassColumns(TABLE *T, OPTIONSTRUCT *Options, int m, int SoilLayers, int *Columns)
{
  const CLASSMAP *C = &(ClassMaps[m]);
  int k;

  if (C->LayerSet == MAP_SOIL) {
    if (Options->FromTable[m])
      TableLayerColumns(T, C->Column, MAP_SOIL, SoilLayers, C->Units, TRUE, Columns);
    else {
      TableKnownLayers(T, C->Column, MAP_SOIL, SoilLayers);
      for (k = 0; k <= SoilLayers; k++)
        Columns[k] = -1;
    }
  }
  else if (Options->FromTable[m])
    Columns[0] = TableColumn(T, C->Column, C->Units, TRUE);
  else {
    TableKnown(T, C->Column);
    Columns[0] = -1;
  }
}

/********************************************************************************
Function Name: InitSoilTable()

Purpose      : Initialize the soil lookup table from the soil class table
               (one row per class, Class = row number). Every class has the
               L root zone layers and the deep layer. The class columns of
               the soil maps without an [INPUTS] key are read (the others
               are known and not read).

Returns      : Number of soil classes
********************************************************************************/
int InitSoilTable(OPTIONSTRUCT *Options, SOILTABLE **SType, const char *Path,
  LAYER *Soil, int SoilLayers, TIMESTRUCT *Time)
{
  const char *Routine = "InitSoilTable";
  TABLE T;
  int i;
  int NSoils;
  int NLayers = SoilLayers + 1;	/* root zone layers and the deep layer */
  int CName, CDepthThresh, CAnisotropy, CMaxInfilt, CCapillary, CDeepFlux;
  int CAlbedo, CManning, CDepth, CKsLat, CKsLatExp;
  int *CWP, *CKs, *CPorosity, *CFCap;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", 1, TABLE_MAXROWS);
  NSoils = T.NRows;

  CWP = calloc(NLayers, sizeof(int));
  CKs = calloc(NLayers, sizeof(int));
  CPorosity = calloc(NLayers, sizeof(int));
  CFCap = calloc(NLayers, sizeof(int));
  if (!CWP || !CKs || !CPorosity || !CFCap)
    ReportError((char *)Routine, 1);

  CName = TableColumn(&T, "Name", "-", TRUE);
  CDepthThresh = TableColumn(&T, "DepthThreshold", "m", TRUE);
  CAnisotropy = TableColumn(&T, "VerticalAnisotropy", "-", TRUE);
  CDeepFlux = TableColumn(&T, "DeepFlux", "m/yr", TRUE);
  CAlbedo = TableColumn(&T, "Albedo", "fraction", TRUE);
  TableLayerColumns(&T, "WiltingPoint", MAP_SOIL, SoilLayers, "m3/m3", TRUE, CWP);
  CMaxInfilt = CCapillary = CManning = -1;
  if (!Options->UseKsatAnisotropy) {
    CMaxInfilt = TableColumn(&T, "MaxInfiltration", "m/s", TRUE);
    TableLayerColumns(&T, "VerticalConductivity", MAP_SOIL, SoilLayers, "m/s", TRUE, CKs);
  }
  else {
    TableKnown(&T, "MaxInfiltration");
    TableKnownLayers(&T, "VerticalConductivity", MAP_SOIL, SoilLayers);
  }
  if (Options->Infiltration == DYNAMIC)
    CCapillary = TableColumn(&T, "CapillaryDrive", "m", TRUE);
  else
    TableKnown(&T, "CapillaryDrive");
  if (Options->Routing)
    CManning = TableColumn(&T, "ManningsN", "-", TRUE);
  else
    TableKnown(&T, "ManningsN");
  ClassColumns(&T, Options, CM_DEPTH, SoilLayers, &CDepth);
  ClassColumns(&T, Options, CM_KSLAT, SoilLayers, &CKsLat);
  ClassColumns(&T, Options, CM_EXPDEC, SoilLayers, &CKsLatExp);
  ClassColumns(&T, Options, CM_POROSITY, SoilLayers, CPorosity);
  ClassColumns(&T, Options, CM_FCAP, SoilLayers, CFCap);
  TableWarnUnknown(&T);

  if (!(Soil->NLayers = (int *)calloc(NSoils, sizeof(int))))
    ReportError((char *)Routine, 1);

  if (!(*SType = (SOILTABLE *)calloc(NSoils, sizeof(SOILTABLE))))
    ReportError((char *)Routine, 1);

  Soil->MaxLayers = SoilLayers;

  for (i = 0; i < NSoils; i++) {
    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*SType)[i].Desc, TableText(&T, i, CName));
    (*SType)[i].Index = i;
    (*SType)[i].NLayers = SoilLayers;
    Soil->NLayers[i] = SoilLayers;

    (*SType)[i].DepthThresh = TableFloat(&T, i, CDepthThresh);
    (*SType)[i].KsAnisotropy = TableFloat(&T, i, CAnisotropy);

    if (!Options->UseKsatAnisotropy)
      (*SType)[i].MaxInfiltrationRate = TableFloat(&T, i, CMaxInfilt);

    if (Options->Infiltration == DYNAMIC)
      (*SType)[i].G_Infilt = TableFloat(&T, i, CCapillary);
    else
      (*SType)[i].G_Infilt = NOT_APPLICABLE;

    /* m/yr as read; ApplyTableParameters() adds the offset and converts it
       to m per time step */
    (*SType)[i].DeepFlux = TableFloat(&T, i, CDeepFlux);

    (*SType)[i].Albedo = TableFloat(&T, i, CAlbedo);

    if (Options->Routing)
      (*SType)[i].Manning = TableFloat(&T, i, CManning);

    /* The class values of the soil maps formed from the table */
    (*SType)[i].Depth = OptionalValue(&T, i, CDepth);
    (*SType)[i].KsLat = OptionalValue(&T, i, CKsLat);
    (*SType)[i].KsLatExp = OptionalValue(&T, i, CKsLatExp);

    /* Root zone layers and the deep layer */
    if (!((*SType)[i].WP = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].Ks = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].Porosity = (float *)calloc(NLayers, sizeof(float))) ||
        !((*SType)[i].FCap = (float *)calloc(NLayers, sizeof(float))))
      ReportError((char *)Routine, 1);

    TableFloats(&T, i, CWP, NLayers, (*SType)[i].WP);
    if (!Options->UseKsatAnisotropy)
      TableFloats(&T, i, CKs, NLayers, (*SType)[i].Ks);
    if (Options->FromTable[CM_POROSITY])
      TableFloats(&T, i, CPorosity, NLayers, (*SType)[i].Porosity);
    if (Options->FromTable[CM_FCAP])
      TableFloats(&T, i, CFCap, NLayers, (*SType)[i].FCap);
  }

  free(CWP);
  free(CKs);
  free(CPorosity);
  free(CFCap);
  FreeTable(&T);
  return NSoils;
}

/* Columns of the story parameters of a vegetation table (-1 for those the
   stomatal scheme does not read, and for the cover, LAI and height of maps
   with an [INPUTS] key) */
typedef struct {
  int Deciduous, OnMid, OnSlope, OffMid, OffSlope, MinLAI;
  int *RootFract;
  int Albedo[12];
  int RsMax, RsMin, Moist, Vpd, Rpc;
  int Vcmax25, G1, G0;
  int Cover, LAI, Height;
} STORYCOLUMNS;

/* Cover: the story's cover map in ClassMaps[] (CM_OVERCOVER or
   CM_UNDERCOVER), followed by its LAI and height maps */
static void StoryColumns(TABLE *T, OPTIONSTRUCT *Options, int Cover, int SoilLayers,
                         STORYCOLUMNS *C)
{
  int Scheme = Options->StomatalScheme;
  int i;
  int AnyDeciduous = FALSE;

  if (!(C->RootFract = calloc(SoilLayers, sizeof(int))))
    ReportError("StoryColumns", 1);

  /* Phenology columns: required when some row is deciduous, read only on
     such rows */
  C->Deciduous = TableColumn(T, "Deciduous", "-", TRUE);
  for (i = 0; i < T->NRows; i++) {
    if (TableInt(T, i, C->Deciduous) != 0)
      AnyDeciduous = TRUE;
  }
  C->OnMid = C->OnSlope = C->OffMid = C->OffSlope = C->MinLAI = -1;
  if (AnyDeciduous) {
    C->OnMid = TableColumn(T, "LeafOnMid", "day", TRUE);
    C->OnSlope = TableColumn(T, "LeafOnSlope", "logit/day", TRUE);
    C->OffMid = TableColumn(T, "LeafOffMid", "day", TRUE);
    C->OffSlope = TableColumn(T, "LeafOffSlope", "logit/day", TRUE);
    C->MinLAI = TableColumn(T, "MinLAIFraction", "fraction", TRUE);
  }
  else {
    TableKnown(T, "LeafOnMid");
    TableKnown(T, "LeafOnSlope");
    TableKnown(T, "LeafOffMid");
    TableKnown(T, "LeafOffSlope");
    TableKnown(T, "MinLAIFraction");
  }
  TableLayerColumns(T, "RootFraction", MAP_ONE, SoilLayers, "fraction", TRUE, C->RootFract);
  TableLayerColumns(T, "Albedo", MAP_MONTH, SoilLayers, "fraction", TRUE, C->Albedo);
  ClassColumns(T, Options, Cover, SoilLayers, &(C->Cover));
  ClassColumns(T, Options, Cover + 1, SoilLayers, &(C->LAI));
  ClassColumns(T, Options, Cover + 2, SoilLayers, &(C->Height));
  C->RsMax = C->RsMin = C->Vpd = C->Rpc = C->Vcmax25 = C->G1 = C->G0 = -1;
  if (Scheme == JARVIS) {
    C->RsMax = TableColumn(T, "MaxResistance", "s/m", TRUE);
    C->RsMin = TableColumn(T, "MinResistance", "s/m", TRUE);
  }
  C->Moist = TableColumn(T, "MoistureThreshold", "-", TRUE);
  if (Scheme == JARVIS) {
    C->Vpd = TableColumn(T, "VpdThreshold", "Pa", TRUE);
    C->Rpc = TableColumn(T, "Rpc", "W/m2", TRUE);
    TableKnown(T, "Vcmax25");
    TableKnown(T, "StomatalSlope");
    TableKnown(T, "StomatalIntercept");
  }
  else {
    C->Vcmax25 = TableColumn(T, "Vcmax25", "umol/m2/s", TRUE);
    C->G1 = TableColumn(T, "StomatalSlope", "kPa^0.5", TRUE);
    C->G0 = TableColumn(T, "StomatalIntercept", "mol/m2/s", TRUE);
    TableKnown(T, "MaxResistance");
    TableKnown(T, "MinResistance");
    TableKnown(T, "VpdThreshold");
    TableKnown(T, "Rpc");
  }
}

/* The story parameters of row i (NOT_APPLICABLE where a column is not
   read) */
static void ReadStory(TABLE *T, int i, STORYCOLUMNS *C, int SoilLayers, STORYTABLE *Story)
{
  if (!(Story->RootFract = (float *) calloc(SoilLayers, sizeof(float))))
    ReportError("ReadStory", 1);
  Story->Deciduous = TableInt(T, i, C->Deciduous);
  if (Story->Deciduous != 0 && Story->Deciduous != 1)
    ReportErrorCode(106, "%s: line %d: Deciduous = %d must be 0 or 1", T->Path, i + 3,
                    Story->Deciduous);
  if (Story->Deciduous) {
    Story->LeafOnMid = TableFloat(T, i, C->OnMid);
    Story->LeafOnSlope = TableFloat(T, i, C->OnSlope);
    Story->LeafOffMid = TableFloat(T, i, C->OffMid);
    Story->LeafOffSlope = TableFloat(T, i, C->OffSlope);
    Story->MinLAIFraction = TableFloat(T, i, C->MinLAI);
  }
  else {
    Story->LeafOnMid = Story->LeafOnSlope = NOT_APPLICABLE;
    Story->LeafOffMid = Story->LeafOffSlope = NOT_APPLICABLE;
    Story->MinLAIFraction = 1.0;
  }
  TableFloats(T, i, C->RootFract, SoilLayers, Story->RootFract);
  TableFloats(T, i, C->Albedo, 12, Story->AlbedoMonthly);
  Story->RsMax = OptionalValue(T, i, C->RsMax);
  Story->RsMin = OptionalValue(T, i, C->RsMin);
  Story->MoistThres = TableFloat(T, i, C->Moist);
  Story->VpdThres = OptionalValue(T, i, C->Vpd);
  Story->Rpc = OptionalValue(T, i, C->Rpc);
  Story->Vcmax25 = OptionalValue(T, i, C->Vcmax25);
  Story->G1 = OptionalValue(T, i, C->G1);
  Story->G0 = OptionalValue(T, i, C->G0);
  Story->Cover = OptionalValue(T, i, C->Cover);
  Story->LAI = OptionalValue(T, i, C->LAI);
  Story->Height = OptionalValue(T, i, C->Height);
}

/* The name of row i */
static void ReadName(TABLE *T, int i, int CName, char *Desc)
{
  if (strlen(TableText(T, i, CName)) > BUFSIZE)
    ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", T->Path,
                    i + 3, BUFSIZE);
  strcpy(Desc, TableText(T, i, CName));
}

/********************************************************************************
Function Name: InitOverstoryTable()

Purpose      : Initialize the overstory lookup table from the overstory class
               table (one row per class, Class = row number, Veg.OverstoryClass
               0 = no overstory)

Returns      : Number of overstory classes
********************************************************************************/
static int InitOverstoryTable(OVERTABLE **OType, const char *Path, OPTIONSTRUCT *Options,
                              int SoilLayers)
{
  const char *Routine = "InitOverstoryTable";
  TABLE T;
  int i;
  int NRows;
  int CName, CTrunk, CAeroAtt, CMaxSnowInt, CDrip, CSnowEff, CExtn;
  STORYCOLUMNS CStory;
  OVERTABLE *Over;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", 0, TABLE_MAXROWS);
  NRows = T.NRows;

  CName = TableColumn(&T, "Name", "-", TRUE);
  CTrunk = TableColumn(&T, "TrunkSpace", "fraction", TRUE);
  CAeroAtt = TableColumn(&T, "AerodynamicAttenuation", "-", TRUE);
  CMaxSnowInt = TableColumn(&T, "MaxSnowInterception", "m", TRUE);
  CDrip = TableColumn(&T, "MassReleaseDripRatio", "-", TRUE);
  CSnowEff = TableColumn(&T, "SnowInterceptionEff", "fraction", TRUE);
  StoryColumns(&T, Options, CM_OVERCOVER, SoilLayers, &CStory);
  CExtn = TableColumn(&T, "LightExtinction", "-", TRUE);
  TableWarnUnknown(&T);

  if (!(*OType = (OVERTABLE *) calloc((NRows > 0) ? NRows : 1, sizeof(OVERTABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NRows; i++) {
    Over = &((*OType)[i]);
    ReadName(&T, i, CName, Over->Desc);
    ReadStory(&T, i, &CStory, SoilLayers, &(Over->Story));

    Over->Trunk = TableFloat(&T, i, CTrunk);
    Over->Cn = TableFloat(&T, i, CAeroAtt);
    Over->MaxSnowInt = TableFloat(&T, i, CMaxSnowInt);
    Over->MDRatio = TableFloat(&T, i, CDrip);
    Over->SnowIntEff = TableFloat(&T, i, CSnowEff);
    Over->ExtnCoeff = TableFloat(&T, i, CExtn);
  }

  free(CStory.RootFract);
  FreeTable(&T);
  return NRows;
}

/********************************************************************************
Function Name: InitUnderstoryTable()

Purpose      : Initialize the understory lookup table from the understory
               class table (one row per class, Class = row number,
               Veg.UnderstoryClass 0 = no understory). A class holds the
               litter cover and impervious fraction of its cells whatever
               their understory cover.

Returns      : Number of understory classes
********************************************************************************/
static int InitUnderstoryTable(UNDERTABLE **UType, const char *Path, OPTIONSTRUCT *Options,
                               int SoilLayers)
{
  const char *Routine = "InitUnderstoryTable";
  TABLE T;
  int i;
  int NRows;
  int CName, CAtten, CLitter, CImperv, CDetFrac, CDetDecay;
  STORYCOLUMNS CStory;
  UNDERTABLE *Under;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Class", 0, TABLE_MAXROWS);
  NRows = T.NRows;

  CName = TableColumn(&T, "Name", "-", TRUE);
  StoryColumns(&T, Options, CM_UNDERCOVER, SoilLayers, &CStory);
  CAtten = TableColumn(&T, "RadiationAttenuation", "-", TRUE);
  CLitter = TableColumn(&T, "LitterCover", "fraction", TRUE);
  CImperv = TableColumn(&T, "ImperviousFraction", "fraction", TRUE);

  /* Some class with an impervious fraction */
  Options->Impervious = FALSE;
  for (i = 0; i < NRows; i++) {
    if (TableFloat(&T, i, CImperv) > 0.0)
      Options->Impervious = TRUE;
  }
  CDetFrac = CDetDecay = -1;
  if (Options->Impervious) {
    CDetFrac = TableColumn(&T, "DetentionFraction", "fraction", TRUE);
    CDetDecay = TableColumn(&T, "DetentionDecay", "-", TRUE);
  }
  else {
    TableKnown(&T, "DetentionFraction");
    TableKnown(&T, "DetentionDecay");
  }
  TableWarnUnknown(&T);

  if (!(*UType = (UNDERTABLE *) calloc((NRows > 0) ? NRows : 1, sizeof(UNDERTABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NRows; i++) {
    Under = &((*UType)[i]);
    ReadName(&T, i, CName, Under->Desc);
    ReadStory(&T, i, &CStory, SoilLayers, &(Under->Story));

    Under->RadAtten = TableFloat(&T, i, CAtten);
    Under->LitterCover = TableFloat(&T, i, CLitter);
    Under->ImpervFrac = TableFloat(&T, i, CImperv);
    if (Under->ImpervFrac > 0) {
      Under->DetentionFrac = TableFloat(&T, i, CDetFrac);
      Under->DetentionDecay = TableFloat(&T, i, CDetDecay);
    }
    else {
      Under->DetentionFrac = 0.;
      Under->DetentionDecay = 0.;
    }
  }

  free(CStory.RootFract);
  FreeTable(&T);
  return NRows;
}

/********************************************************************************
InitLakeTable()

The lake table: one row per lake (Lake = row number, the Lake.ID value),
its name, the segment it pours into and the scale and exponent of its
outflow power law.
********************************************************************************/
int InitLakeTable(LAKETABLE **LType, const char *Path, OPTIONSTRUCT *Options)
{
  const char *Routine = "InitLakeTable";
  TABLE T;
  int i;
  int NumLakes;
  int CName, COutlet, CScale, CExponent;

  ReadTable(Path, &T);
  TableCheckRows(&T, "Lake", 1, TABLE_MAXROWS);
  NumLakes = T.NRows;
  CName = TableColumn(&T, "Name", "-", TRUE);
  COutlet = TableColumn(&T, "OutletSegment", "-", TRUE);
  CScale = TableColumn(&T, "OutflowScale", "-", TRUE);
  CExponent = TableColumn(&T, "OutflowExponent", "-", TRUE);
  TableWarnUnknown(&T);

  if (!(*LType = (LAKETABLE *)calloc(NumLakes, sizeof(LAKETABLE))))
    ReportError((char *)Routine, 1);

  for (i = 0; i < NumLakes; i++) {
    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*LType)[i].Name, TableText(&T, i, CName));
    (*LType)[i].Index = i+1;
    (*LType)[i].Outlet = TableInt(&T, i, COutlet);
    (*LType)[i].PowLawScale = TableFloat(&T, i, CScale);
    (*LType)[i].PowLawExponent = TableFloat(&T, i, CExponent);
  }
  strcpy(Options->LakeTablePath, Path);

  FreeTable(&T);
  return NumLakes;
}

/********************************************************************************
InitCO2Table()

The CO2 file (MEDLYN): the atmospheric CO2 (umol/mol) of calendar years,
columns Year and CO2. Every year from Start to the last time step needs one
row; other years are not read.
********************************************************************************/
static void InitCO2Table(const char *Path, TIMESTRUCT *Time)
{
  TABLE T;
  int CYear, CCO2;
  int Last, NYears;
  int *Found;
  int i, y;

  ReadTable(Path, &T);
  CYear = TableColumn(&T, "Year", "-", TRUE);
  CCO2 = TableColumn(&T, "CO2", "umol/mol", TRUE);
  TableWarnUnknown(&T);

  /* The year of the last time step: End's, or the one before at its start */
  Last = Time->End.Year;
  if (Time->End.Month == 1 && Time->End.Day == 1 && Time->End.Hour == 0 && Time->End.Min == 0)
    Last--;
  CO2FirstYear = Time->Start.Year;
  NYears = Last - CO2FirstYear + 1;
  if (!(CO2OfYear = (float *) calloc(NYears, sizeof(float))) ||
      !(Found = (int *) calloc(NYears, sizeof(int))))
    ReportError("InitCO2Table", 1);

  for (i = 0; i < T.NRows; i++) {
    y = TableInt(&T, i, CYear) - CO2FirstYear;
    if (y < 0 || y >= NYears)
      continue;
    if (Found[y])
      ReportErrorCode(106, "%s: line %d: a second row for year %d", Path, i + 3,
                      y + CO2FirstYear);
    CO2OfYear[y] = TableFloat(&T, i, CCO2);
    Found[y] = TRUE;
  }
  for (y = 0; y < NYears; y++) {
    if (!Found[y])
      ReportErrorCode(106, "%s: no row for year %d (the run's time steps are in %d-%d)", Path,
                      y + CO2FirstYear, CO2FirstYear, Last);
  }

  free(Found);
  FreeTable(&T);
}

/* The atmospheric CO2 of a year of the run (umol/mol) */
float YearCO2(int Year)
{
  return CO2OfYear[Year - CO2FirstYear];
}
