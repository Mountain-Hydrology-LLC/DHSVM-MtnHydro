
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "Calendar.h"
#include "fileio.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "getinit.h"
#include "constants.h"
#include "mapio.h"
#include "tables.h"

static void OpenStationFile(METLOCATION *Station, TIMESTRUCT *Time);
static void LastRecordTime(METLOCATION *Station, DATE *Last);
static void StepBefore(DATE *Day, int Minutes, DATE *Before);

/* Columns of a station data file and their units */
static const char *MetColumns[NMETCOLUMNS] = {
  "Time", "AirTemperature", "WindSpeed", "RelativeHumidity", "ShortwaveIn",
  "LongwaveIn", "Precipitation"
};
static const char *MetUnits[NMETCOLUMNS] = { "-", "C", "m/s", "%", "W/m2", "W/m2", "m/step" };

 /*******************************************************************************
   Function name: InitMetSources()

   Purpose      : Read the station table and open the station data files
                  ([INPUTS] Station Table, Station Data Directory). Every
                  station of the table is used, wherever it lies.

   Required     :
     LISTPTR Input           - Linked list with input strings
     OPTIONSTRUCT *Options   - Structure with different program options
     MAPSIZE *Map            - Coverage and resolution of model area
     TIMESTRUCT *Time        - Begin and end times, model timestep
     int *NStats             - Number of meteorological stations
     METLOCATION **Stat      - Station information

 *******************************************************************************/
void InitMetSources(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  TIMESTRUCT *Time, int *NStats, METLOCATION **Stat)
{
  char *Routine = "InitMetSources";
  char Path[BUFSIZE + 1];
  char Dir[BUFSIZE + 1];
  TABLE T;
  int i, k;
  int CName, CFile, CNorth, CEast, CElev;
  int CLapse[12], CNormal[12], CSnow, CMeanT;
  size_t n;
  const char *Separator;
  float East;
  float North;

  RequireConfigValue(Input, "INPUTS", "Station Table", Path);
  RequireConfigValue(Input, "INPUTS", "Station Data Directory", Dir);
  n = strlen(Dir);
  Separator = (n > 0 && (Dir[n - 1] == '/' || Dir[n - 1] == '\\')) ? "" : "/";

  ReadTable(Path, &T);
  if (T.NRows < 1)
    ReportErrorCode(106, "%s: no data rows", Path);
  CName = TableColumn(&T, "Name", "-", TRUE);
  CFile = TableColumn(&T, "File", "-", TRUE);
  CNorth = TableColumn(&T, "North", "m", TRUE);
  CEast = TableColumn(&T, "East", "m", TRUE);
  CElev = TableColumn(&T, "Elevation", "m", TRUE);
  TableLayerColumns(&T, "LapseRate", MAP_MONTH, Map->SoilLayers, "C/m", TRUE, CLapse);
  if (Options->Prism)
    TableLayerColumns(&T, "PrecipNormal", MAP_MONTH, Map->SoilLayers, "-", TRUE, CNormal);
  else
    TableKnownLayers(&T, "PrecipNormal", MAP_MONTH, Map->SoilLayers);
  CSnow = -1;
  if (Options->SnowPattern)
    CSnow = TableColumn(&T, "SnowNormal", "-", TRUE);
  else
    TableKnown(&T, "SnowNormal");
  CMeanT = -1;
  if (Options->SoilFreezing)
    CMeanT = TableColumn(&T, "MeanAirTemperature", "C", TRUE);
  else
    TableKnown(&T, "MeanAirTemperature");
  TableWarnUnknown(&T);

  *NStats = T.NRows;
  Message("\nReading %d met stations\n", *NStats);

  /* Allocate memory for the stations */
  if (!(*Stat = (METLOCATION *)calloc(*NStats, sizeof(METLOCATION))))
    ReportError(Routine, 1);

  for (i = 0; i < *NStats; i++) {
    if (strlen(TableText(&T, i, CName)) > BUFSIZE)
      ReportErrorCode(106, "%s: line %d: column Name is longer than %d characters", Path,
                      i + 3, BUFSIZE);
    strcpy((*Stat)[i].Name, TableText(&T, i, CName));

    North = TableFloat(&T, i, CNorth);
    East = TableFloat(&T, i, CEast);
    (*Stat)[i].Loc.N = Round(((Map->Yorig - 0.5 * Map->DY) - North) / Map->DY);
    (*Stat)[i].Loc.E = Round((East - (Map->Xorig + 0.5 * Map->DX)) / Map->DX);

    (*Stat)[i].Elev = TableFloat(&T, i, CElev);

    for (k = 0; k < 12; k++)
      (*Stat)[i].LapseRate[k] = TableFloat(&T, i, CLapse[k]);

    if (Options->Prism) {
      for (k = 0; k < 12; k++)
        (*Stat)[i].PrismPrecip[k] = TableFloat(&T, i, CNormal[k]);
    }
    if (Options->SnowPattern)
      (*Stat)[i].SnowPatternBase = TableFloat(&T, i, CSnow);
    if (Options->SoilFreezing)
      (*Stat)[i].MeanAirTemp = TableFloat(&T, i, CMeanT);

    MakePath((*Stat)[i].MetFile.FileName, sizeof((*Stat)[i].MetFile.FileName), "%s%s%s",
             Dir, Separator, TableText(&T, i, CFile));
    OpenStationFile(&((*Stat)[i]), Time);
  }
  Message("Forced to include all %d stations \n", *NStats);

  FreeTable(&T);
}

/* Splits the line in Station->Record into Station->Fields; FALSE if it does
   not have the fields of the header */
static int SplitRecord(METLOCATION *Station)
{
  return (SplitFields(Station->Record, Station->Fields, Station->NFields) == Station->NFields);
}

/*****************************************************************************
  OpenStationFile()

  Opens a station data file and checks its header (the seven columns and
  their units), the form of the times, and that its records cover the run
  (the first at or before Start, the last at or after the start of the last
  step); reads ahead to the first record at or after Start
*****************************************************************************/
static void OpenStationFile(METLOCATION *Station, TIMESTRUCT *Time)
{
  char *Path = Station->MetFile.FileName;
  char Text[2][BUFSIZE + 1];
  DATE First, Last, Record, LastStep;
  int j, k, n;
  long long StartMin;

  if (!(Station->Record = malloc(TABLE_LINESIZE)))
    ReportError("OpenStationFile", 1);
  if (!(Station->MetFile.FilePtr = fopen(Path, "rb")))
    ReportErrorCode(104, "%s: cannot open", Path);

  /* Column names: the seven columns by name, any other is warned */
  if (!ReadTextLine(Station->MetFile.FilePtr, Station->Record, TABLE_LINESIZE, Path, 1, 107))
    ReportErrorCode(107, "%s: line 1 (column names) is missing", Path);
  for (n = 1, j = 0; Station->Record[j] != '\0'; j++)
    if (Station->Record[j] == ',')
      n++;
  Station->NFields = n;
  if (!(Station->Fields = calloc(n, sizeof(char *))))
    ReportError("OpenStationFile", 1);
  SplitRecord(Station);
  for (k = 0; k < NMETCOLUMNS; k++)
    Station->Column[k] = -1;
  for (j = 0; j < n; j++) {
    for (k = 0; k < NMETCOLUMNS; k++) {
      if (strcmp(Station->Fields[j], MetColumns[k]) == 0) {
        if (Station->Column[k] >= 0)
          ReportErrorCode(107, "%s: column %s repeated", Path, MetColumns[k]);
        Station->Column[k] = j;
      }
    }
  }
  for (k = 0; k < NMETCOLUMNS; k++) {
    if (Station->Column[k] < 0)
      ReportErrorCode(107, "%s: column %s missing", Path, MetColumns[k]);
  }
  for (j = 0; j < n; j++) {
    for (k = 0; k < NMETCOLUMNS && Station->Column[k] != j; k++)
      ;
    if (k == NMETCOLUMNS)
      fprintf(stderr, "DHSVM WARNING: %s: column %s not used\n", Path, Station->Fields[j]);
  }

  /* Units */
  if (!ReadTextLine(Station->MetFile.FilePtr, Station->Record, TABLE_LINESIZE, Path, 2, 107))
    ReportErrorCode(107, "%s: line 2 (units) is missing", Path);
  if (!SplitRecord(Station))
    ReportErrorCode(107, "%s: line 2 has a different number of fields than line 1", Path);
  for (k = 0; k < NMETCOLUMNS; k++) {
    if (strcmp(Station->Fields[Station->Column[k]], MetUnits[k]) != 0)
      ReportErrorCode(107, "%s: column %s has units \"%s\", expected \"%s\"", Path,
                      MetColumns[k], Station->Fields[Station->Column[k]], MetUnits[k]);
  }
  Station->LineNo = 2;

  /* Records before Start are skipped; the first record must be at or
     before Start */
  StartMin = MinutesOf(&(Time->Start));
  Station->Pending = FALSE;
  n = 0;
  while (ReadTextLine(Station->MetFile.FilePtr, Station->Record, TABLE_LINESIZE, Path,
                      ++(Station->LineNo), 107)) {
    if (!SplitRecord(Station))
      ReportErrorCode(107, "%s: line %d has a different number of fields than line 1", Path,
                      Station->LineNo);
    if (!SScanIsoDate(Station->Fields[Station->Column[0]], &Record))
      ReportErrorCode(107, "%s: line %d: time \"%s\" is not YYYY-MM-DD HH:MM", Path,
                      Station->LineNo, Station->Fields[Station->Column[0]]);
    if (n++ == 0) {
      First = Record;
      if (MinutesOf(&First) > StartMin) {
        SPrintIsoDate(&First, Text[0]);
        SPrintIsoDate(&(Time->Start), Text[1]);
        ReportErrorCode(107, "%s: first record %s is after Start %s", Path, Text[0], Text[1]);
      }
    }
    if (MinutesOf(&Record) >= StartMin) {
      Station->Pending = TRUE;
      break;
    }
  }
  if (n == 0)
    ReportErrorCode(107, "%s: no records", Path);

  /* The last record must reach the start of the last time step */
  LastRecordTime(Station, &Last);
  StepBefore(&(Time->End), Time->Dt / SECPMIN, &LastStep);
  if (MinutesOf(&Last) < MinutesOf(&LastStep)) {
    SPrintIsoDate(&Last, Text[0]);
    SPrintIsoDate(&LastStep, Text[1]);
    ReportErrorCode(107, "%s: last record %s is before %s", Path, Text[0], Text[1]);
  }
}

/*****************************************************************************
  StepBefore()

  The instant Minutes (at most a day) before Day
*****************************************************************************/
static void StepBefore(DATE *Day, int Minutes, DATE *Before)
{
  int DaysPerMonth[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

  *Before = *Day;
  Before->Min -= Minutes;
  while (Before->Min < 0) {
    Before->Min += MINPHOUR;
    Before->Hour--;
  }
  if (Before->Hour < 0) {
    Before->Hour += HOURPDAY;
    Before->Day--;
    if (Before->Day < 1) {
      Before->Month--;
      if (Before->Month < 1) {
        Before->Month = 12;
        Before->Year--;
      }
      if (IsLeapYear(Before->Year))
        DaysPerMonth[1] = 29;
      Before->Day = DaysPerMonth[Before->Month - 1];
    }
  }
}

/*****************************************************************************
  LastRecordTime()

  Time of the last non-blank line of a station data file; the file position
  is kept
*****************************************************************************/
static void LastRecordTime(METLOCATION *Station, DATE *Last)
{
  FILE *File = Station->MetFile.FilePtr;
  char *Path = Station->MetFile.FileName;
  char *Tail;
  char **Fields;
  long long Size = 0, Start, Here;
  size_t n;
  char *Line;
  int i;

  if (!(Tail = malloc(TABLE_LINESIZE + 1)) ||
      !(Fields = calloc(Station->NFields, sizeof(char *))))
    ReportError("LastRecordTime", 1);
#ifdef _WIN32
  Here = _ftelli64(File);
  if (_fseeki64(File, 0, SEEK_END) != 0 || (Size = _ftelli64(File)) < 0)
#else
  Here = (long long) ftello(File);
  if (fseeko(File, 0, SEEK_END) != 0 || (Size = (long long) ftello(File)) < 0)
#endif
    ReportErrorCode(107, "%s: cannot read", Path);
  Start = (Size > TABLE_LINESIZE) ? Size - TABLE_LINESIZE : 0;
#ifdef _WIN32
  _fseeki64(File, Start, SEEK_SET);
#else
  fseeko(File, (off_t) Start, SEEK_SET);
#endif
  n = fread(Tail, 1, (size_t) (Size - Start), File);
  Tail[n] = '\0';
#ifdef _WIN32
  _fseeki64(File, Here, SEEK_SET);
#else
  fseeko(File, (off_t) Here, SEEK_SET);
#endif

  /* Drop the line ends and blanks at the end, then find the line start */
  i = (int) n - 1;
  while (i >= 0 && (Tail[i] == '\n' || Tail[i] == '\r' || Tail[i] == ' ' || Tail[i] == '\t'))
    Tail[i--] = '\0';
  while (i >= 0 && Tail[i] != '\n')
    i--;
  Line = &Tail[i + 1];
  if (SplitFields(Line, Fields, Station->NFields) != Station->NFields)
    ReportErrorCode(107, "%s: the last line has a different number of fields than line 1", Path);
  if (!SScanIsoDate(Fields[Station->Column[0]], Last))
    ReportErrorCode(107, "%s: last line: time \"%s\" is not YYYY-MM-DD HH:MM", Path,
                    Fields[Station->Column[0]]);
  free(Fields);
  free(Tail);
}
