/*
 * DESCRIPTION:  The output files of a run in the output directory:
 *               Basin.csv, Point_<Name>.csv, Class_<O>_<U>.csv, Streamflow.csv,
 *               Channels.csv, NetworkDynamics.csv (comma-separated, a names
 *               line, a units line, one row per time step labelled with the
 *               step's start; numbers %.9g, NA where a cell lacks a
 *               layer), the map stacks (OutputMaps.c) and the
 *               snapshots (StoreModelState.c).
 *
 *               A row stamped T holds the value of an Interval variable over
 *               the step [T, T + Dt] and the value of an Instant variable at
 *               T: the row of a step is written when the step has ended, with
 *               the Instant values held from the end of the step before (from
 *               Start for the first row).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "output.h"
#include "registry.h"

#define NETWORK_COLUMNS "Time,Segment,Channel.Storage,Channel.Inflow,Channel.LateralInflow," \
  "Channel.Outflow,Channel.Infiltration,Channel.Evaporation"
#define NETWORK_UNITS "step start,-,m3,m3/s,m3/s,m3/s,m3/s,m3/s"
#define NETWORK_NCOLUMNS 6

static const int NetworkCodes[NETWORK_NCOLUMNS] = {
  V_CHANNEL_STORAGE, V_CHANNEL_INFLOW, V_CHANNEL_LATERALINFLOW, V_CHANNEL_OUTFLOW,
  V_CHANNEL_INFILTRATION, V_CHANNEL_EVAPORATION
};

static void CreateSeries(FILES *File, char *Path, const char *Name, const char *First,
                         const char *FirstUnits, OUTLIST *List);
static HELDVALUE *HeldValues(int N);
static void WriteValue(FILE *File, double Value, int Has);
static void SeriesColumn(FILE *File, int Write, int Instant, double Value, int Has,
                         HELDVALUE *Held);
static void CloseFile(FILES *File);

/*****************************************************************************
  CreateOutputs()

  Creates (empties) every requested output file with its names and units
  lines; a file that cannot be written is error 113. Sets up the held
  Instant values of each series file.
*****************************************************************************/
void CreateOutputs(MAPSIZE *Map, TIMESTRUCT *Time, Channel *Streams, DUMPSTRUCT *Dump)
{
  char Name[2 * BUFSIZE + 1];
  Channel *Segment;
  int NSegments = 0;
  int i;

  Dump->BasinHeld = NULL;
  Dump->ChannelsHeld = NULL;
  Dump->NetworkHeld = NULL;
  for (Segment = Streams; Segment != NULL; Segment = Segment->next)
    NSegments++;

  if (Dump->List[OUT_BASIN].N > 0) {
    CreateSeries(&(Dump->Basin), Dump->Path, "Basin.csv", "Time", "step start",
                 &(Dump->List[OUT_BASIN]));
    Dump->BasinHeld = HeldValues(Dump->List[OUT_BASIN].N);
  }
  if (Dump->List[OUT_POINT].N > 0) {
    for (i = 0; i < Dump->NPoints; i++) {
      snprintf(Name, sizeof(Name), "Point_%s.csv", Dump->Points[i].Name);
      CreateSeries(&(Dump->Points[i].File), Dump->Path, Name, "Time", "step start",
                   &(Dump->List[OUT_POINT]));
      Dump->Points[i].Held = HeldValues(Dump->List[OUT_POINT].N);
    }
  }
  for (i = 0; i < Dump->NClasses; i++) {
    snprintf(Name, sizeof(Name), "Class_%d_%d.csv", Dump->Classes[i].Over,
             Dump->Classes[i].Under);
    CreateSeries(&(Dump->Classes[i].File), Dump->Path, Name, "Time", "step start",
                 &(Dump->List[OUT_CLASS]));
    Dump->Classes[i].Held = HeldValues(Dump->List[OUT_CLASS].N);
  }
  if (Dump->Streamflow) {
    CreateSeries(&(Dump->StreamflowFile), Dump->Path, "Streamflow.csv", "Time", "step start",
                 NULL);
    for (Segment = Streams; Segment != NULL; Segment = Segment->next)
      if (Segment->record)
        fprintf(Dump->StreamflowFile.FilePtr, ",%s", Segment->record_name);
    fprintf(Dump->StreamflowFile.FilePtr, "\nstep start");
    for (Segment = Streams; Segment != NULL; Segment = Segment->next)
      if (Segment->record)
        fprintf(Dump->StreamflowFile.FilePtr, ",m3/s");
    fprintf(Dump->StreamflowFile.FilePtr, "\n");
  }
  if (Dump->List[OUT_CHANNELS].N > 0) {
    CreateSeries(&(Dump->Channels), Dump->Path, "Channels.csv", "Time,Segment",
                 "step start,-", &(Dump->List[OUT_CHANNELS]));
    Dump->ChannelsHeld = HeldValues(NSegments * Dump->List[OUT_CHANNELS].N);
  }
  if (Dump->NetworkSteps > 0) {
    CreateSeries(&(Dump->Network), Dump->Path, "NetworkDynamics.csv", NULL, NULL, NULL);
    fprintf(Dump->Network.FilePtr, "%s\n%s\n", NETWORK_COLUMNS, NETWORK_UNITS);
    Dump->NetworkHeld = HeldValues(NSegments * NETWORK_NCOLUMNS);
  }
  CreateMapStacks(Map, Time, Dump);
}

/* Opens <Path><Name> for writing; with List, writes the names line (First,
   then the columns) and the units line (FirstUnits, then the units) */
static void CreateSeries(FILES *File, char *Path, const char *Name, const char *First,
                         const char *FirstUnits, OUTLIST *List)
{
  int i;

  MakePath(File->FileName, sizeof(File->FileName), "%s%s", Path, Name);
  if (!(File->FilePtr = fopen(File->FileName, "wb")))
    ReportErrorCode(113, "%s: cannot write", File->FileName);
  if (First != NULL)
    fprintf(File->FilePtr, "%s", First);
  if (List == NULL)
    return;
  for (i = 0; i < List->N; i++)
    fprintf(File->FilePtr, ",%s", List->Vars[i].Name);
  fprintf(File->FilePtr, "\n%s", FirstUnits);
  for (i = 0; i < List->N; i++)
    fprintf(File->FilePtr, ",%s", RegistryRow(List->Vars[i].Var)->Units);
  fprintf(File->FilePtr, "\n");
}

/* N held values (NULL for none) */
static HELDVALUE *HeldValues(int N)
{
  HELDVALUE *Held = NULL;

  if (N > 0 && !(Held = (HELDVALUE *) calloc(N, sizeof(HELDVALUE))))
    ReportError("CreateOutputs", 1);
  return Held;
}

/*****************************************************************************
  WriteNumber()

  A value as written in the text outputs: %.9g
*****************************************************************************/
void WriteNumber(FILE *File, double Value)
{
  fprintf(File, "%.9g", Value);
}

/* A field of a series row: the value, or NA */
static void WriteValue(FILE *File, double Value, int Has)
{
  fputc(',', File);
  if (Has)
    WriteNumber(File, Value);
  else
    fputs("NA", File);
}

/* One column of a series row with the column's value now, at the end of
   the row's step: an Interval column writes it; an Instant column writes
   the value held at the row's time and holds the value now for the next
   row. Write FALSE (at Start) only holds. */
static void SeriesColumn(FILE *File, int Write, int Instant, double Value, int Has,
                         HELDVALUE *Held)
{
  if (Instant) {
    if (Write)
      WriteValue(File, Held->Value, Held->Has);
    Held->Value = Value;
    Held->Has = Has;
  }
  else if (Write)
    WriteValue(File, Value, Has);
}

/*****************************************************************************
  WriteSeriesRows()

  The rows of the time step starting at Current, written when it has
  ended: basin means and basin variables, class means, point values.
  Write FALSE: only the Instant values are held (at Start).
*****************************************************************************/
void WriteSeriesRows(DATE *Current, MODELVIEW *Model, TOTALS *Totals, DUMPSTRUCT *Dump,
                     int Write)
{
  OUTLIST *List;
  OUTVAR *Var;
  FILE *File = NULL;
  char Time[NAMESIZE + 1];
  float Value;
  int Has;
  int i, k;

  if (Write)
    SPrintIsoDate(Current, Time);

  List = &(Dump->List[OUT_BASIN]);
  if (List->N > 0) {
    if (Write) {
      File = Dump->Basin.FilePtr;
      fputs(Time, File);
    }
    for (i = 0; i < List->N; i++) {
      Var = &(List->Vars[i]);
      if (Write || Var->Instant)
        SeriesColumn(File, Write, Var->Instant, BasinValue(Totals, Var->Code, Var->Layer), TRUE,
                     &(Dump->BasinHeld[i]));
    }
    if (Write)
      fputc('\n', File);
  }

  List = &(Dump->List[OUT_CLASS]);
  for (k = 0; k < Dump->NClasses; k++) {
    if (Write) {
      File = Dump->Classes[k].File.FilePtr;
      fputs(Time, File);
    }
    for (i = 0; i < List->N; i++) {
      Var = &(List->Vars[i]);
      if (Write || Var->Instant)
        SeriesColumn(File, Write, Var->Instant,
                     AggregateValue(&(Totals->Class[k]), Var->Code, Var->Layer),
                     Totals->Class[k].NCells > 0, &(Dump->Classes[k].Held[i]));
    }
    if (Write)
      fputc('\n', File);
  }

  List = &(Dump->List[OUT_POINT]);
  if (List->N > 0) {
    for (k = 0; k < Dump->NPoints; k++) {
      if (Write) {
        File = Dump->Points[k].File.FilePtr;
        fputs(Time, File);
      }
      for (i = 0; i < List->N; i++) {
        Var = &(List->Vars[i]);
        if (!(Write || Var->Instant))
          continue;
        Value = CellValue(Var->Code, Var->Layer, Model, Dump->Points[k].y, Dump->Points[k].x,
                          &Has);
        SeriesColumn(File, Write, Var->Instant, Value, Has, &(Dump->Points[k].Held[i]));
      }
      if (Write)
        fputc('\n', File);
    }
  }
}

/*****************************************************************************
  WriteChannelRows()

  The rows of time step Step (0 = the first), starting at Current, written
  when it has ended: outflow of the saved segments, the requested variables
  of every segment, and the network dynamics every NetworkSteps steps (its
  storage held every step). Write FALSE: only the Instant values are held
  (at Start).
*****************************************************************************/
void WriteChannelRows(DATE *Current, int Step, MODELVIEW *Model, DUMPSTRUCT *Dump,
                      int Write)
{
  OUTLIST *List = &(Dump->List[OUT_CHANNELS]);
  OUTVAR *Var;
  Channel *Segment;
  FILE *File = NULL;
  char Time[NAMESIZE + 1];
  int Instant[NETWORK_NCOLUMNS];
  int WriteNetwork;
  int i, s;

  if (Write)
    SPrintIsoDate(Current, Time);

  if (Dump->Streamflow && Write) {
    File = Dump->StreamflowFile.FilePtr;
    fputs(Time, File);
    for (Segment = Model->Streams; Segment != NULL; Segment = Segment->next) {
      if (Segment->record) {
        fputc(',', File);
        WriteNumber(File, SegmentValue(V_CHANNEL_OUTFLOW, Segment, Model->Dt));
      }
    }
    fputc('\n', File);
  }

  if (List->N > 0) {
    File = Dump->Channels.FilePtr;
    for (Segment = Model->Streams, s = 0; Segment != NULL; Segment = Segment->next, s++) {
      if (Write)
        fprintf(File, "%s,%d", Time, (int) Segment->id);
      for (i = 0; i < List->N; i++) {
        Var = &(List->Vars[i]);
        if (Write || Var->Instant)
          SeriesColumn(File, Write, Var->Instant, SegmentValue(Var->Code, Segment, Model->Dt),
                       TRUE, &(Dump->ChannelsHeld[s * List->N + i]));
      }
      if (Write)
        fputc('\n', File);
    }
  }

  if (Dump->NetworkSteps > 0) {
    File = Dump->Network.FilePtr;
    WriteNetwork = Write && (Step % Dump->NetworkSteps == 0);
    for (i = 0; i < NETWORK_NCOLUMNS; i++)
      Instant[i] = !RegistryIsInterval(CodeVar(NetworkCodes[i]));
    for (Segment = Model->Streams, s = 0; Segment != NULL; Segment = Segment->next, s++) {
      if (WriteNetwork)
        fprintf(File, "%s,%d", Time, (int) Segment->id);
      for (i = 0; i < NETWORK_NCOLUMNS; i++) {
        if (WriteNetwork || Instant[i])
          SeriesColumn(File, WriteNetwork, Instant[i],
                       SegmentValue(NetworkCodes[i], Segment, Model->Dt), TRUE,
                       &(Dump->NetworkHeld[s * NETWORK_NCOLUMNS + i]));
      }
      if (WriteNetwork)
        fputc('\n', File);
    }
  }
}

/*****************************************************************************
  HoldStartValues()

  The Instant values at Start of every series file (after the aggregation
  before the first step), held for the first rows
*****************************************************************************/
void HoldStartValues(MODELVIEW *Model, TOTALS *Totals, DUMPSTRUCT *Dump)
{
  WriteSeriesRows(NULL, Model, Totals, Dump, FALSE);
#ifndef SNOW_ONLY
  WriteChannelRows(NULL, 0, Model, Dump, FALSE);
#endif
}

/*****************************************************************************
  ExecDump()

  Everything written after time step StepsDone (1 = the first) has ended:
  the series rows of the step starting at Current, the snapshots valid at
  its end, the map layers due
*****************************************************************************/
void ExecDump(MAPSIZE *Map, DATE *Current, int StepsDone, MODELVIEW *Model, TOTALS *Totals,
              DUMPSTRUCT *Dump, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, CHANNEL *ChannelData,
              LAKETABLE *LType)
{
  WriteSeriesRows(Current, Model, Totals, Dump, TRUE);
#ifndef SNOW_ONLY
  WriteChannelRows(Current, StepsDone - 1, Model, Dump, TRUE);
#endif
  StoreSnapshots(StepsDone, Map, Dump, Model->SoilFreezing, Model->Topo, PrecipMap,
                 SnowMap, Model->Veg, Model->Soil, ChannelData, LType);
  WriteMapLayers(StepsDone, Model, Dump);
}

/*****************************************************************************
  StoreSnapshots()

  Writes the snapshots ([RUN] State Time n) valid after StepsDone time steps
*****************************************************************************/
void StoreSnapshots(int StepsDone, MAPSIZE *Map, DUMPSTRUCT *Dump, int SoilFreezing,
  TOPOPIX **TopoMap, PRECIPPIX **PrecipMap, SNOWPIX **SnowMap, VEGPIX **VegMap,
  SOILPIX **SoilMap, CHANNEL *ChannelData, LAKETABLE *LType)
{
  int i;

  for (i = 0; i < Dump->NStates; i++) {
    if (Dump->StateStep[i] == StepsDone)
      StoreModelState(Dump->Path, &(Dump->DState[i]), Map, SoilFreezing, TopoMap,
                      PrecipMap, SnowMap, VegMap, SoilMap, ChannelData, LType);
  }
}

/*****************************************************************************
  CloseOutputs()

  Closes every output file; a write error is error 113
*****************************************************************************/
void CloseOutputs(DUMPSTRUCT *Dump)
{
  int i;

  if (Dump->List[OUT_BASIN].N > 0)
    CloseFile(&(Dump->Basin));
  if (Dump->List[OUT_POINT].N > 0)
    for (i = 0; i < Dump->NPoints; i++)
      CloseFile(&(Dump->Points[i].File));
  for (i = 0; i < Dump->NClasses; i++)
    CloseFile(&(Dump->Classes[i].File));
  if (Dump->Streamflow)
    CloseFile(&(Dump->StreamflowFile));
  if (Dump->List[OUT_CHANNELS].N > 0)
    CloseFile(&(Dump->Channels));
  if (Dump->NetworkSteps > 0)
    CloseFile(&(Dump->Network));
}

static void CloseFile(FILES *File)
{
  if (ferror(File->FilePtr) || fclose(File->FilePtr) != 0)
    ReportErrorCode(113, "%s: cannot write", File->FileName);
}
