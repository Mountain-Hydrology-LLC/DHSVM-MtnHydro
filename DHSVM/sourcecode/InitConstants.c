
#include <ctype.h>
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
#include "output.h"
#include "rad.h"

/* [RUN] time keys as read, checked once [GRID] is read */
typedef struct {
  char StepText[BUFSIZE + 1];	/* Time Step as written (hours) */
  double Seconds;		/* Time Step in seconds */
  char StartText[BUFSIZE + 1];	/* Start as written */
  char EndText[BUFSIZE + 1];	/* End as written */
  DATE Start;
  DATE End;
} RUNTIMES;

static void InitRun(LISTPTR Input, OPTIONSTRUCT *Options, RUNTIMES *Run);
static void InitRunTimes(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, RUNTIMES *Run,
                         TIMESTRUCT *Time, DUMPSTRUCT *Dump);
static void DirectoryPath(char *Buffer, size_t Size, const char *Dir);
static void InitGrid(LISTPTR Input, MAPSIZE *Map, SOLARGEOMETRY *SolarGeo);
static void InitOptions(LISTPTR Input, OPTIONSTRUCT *Options);
static void InitInputKeys(LISTPTR Input, OPTIONSTRUCT *Options);
static void InitParameters(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map);
static void CheckCombinations(OPTIONSTRUCT *Options);
static int ConfigChoice(LISTPTR Input, const char *Section, const char *Key,
                        const char **Choices, int NChoices);
static float Parameter(LISTPTR Input, const char *Key);
static float OptionalParameter(LISTPTR Input, const char *Key, int Required);

/*****************************************************************************
  Function name: InitConstants()

  Purpose      : Initialize constants and settings for DHSVM run
                 Processes the following sections in InFile:
                 [RUN]
                 [GRID]
                 [OPTIONS]
                 [INPUTS] (which keys are present; the files are read later)
                 [PARAMETERS]
                 [OUTPUT] (the form of the keys; the requests are checked
                 once the model is set up)
                 Every key is checked (present, non-empty, a valid value) in
                 the order of the sections, then the option combinations,
                 then the time step, run period, snapshot instants, map
                 times and network dynamics interval, then the output
                 directory.

  Required     :
    LISTPTR Input          - Linked list with input strings
    OPTIONSTRUCT *Options   - Structure with different program options
    MAPSIZE *Map            - Coverage and resolution of model area
    SOLARGEOMETRY *SolarGeo - Solar geometry information
    TIMESTRUCT *Time        - Begin and end times, model timestep
    DUMPSTRUCT *Dump        - Output directory, start state, snapshots and
                              output requests

*****************************************************************************/
void InitConstants(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
		   SOLARGEOMETRY *SolarGeo, TIMESTRUCT *Time, DUMPSTRUCT *Dump)
{
  RUNTIMES Run;			/* run period and time step as read */

  InitRun(Input, Options, &Run);
  InitGrid(Input, Map, SolarGeo);
  InitOptions(Input, Options);
  InitInputKeys(Input, Options);
  InitParameters(Input, Options, Map);
  InitOutputKeys(Input, Dump);
  CheckCombinations(Options);
  InitRunTimes(Input, Options, Map, &Run, Time, Dump);
}

/*****************************************************************************
  InitRun()

  [RUN]: model version, run period and time step (read here, checked by
  InitRunTimes() once [GRID] is read), start state, the form of the snapshot
  instants, groundwater spin-up, progress messages, the rebuilding of
  derived states from the snapshot
*****************************************************************************/
static void InitRun(LISTPTR Input, OPTIONSTRUCT *Options, RUNTIMES *Run)
{
  char Value[BUFSIZE + 1];
  char KeyName[BUFSIZE + 1];
  char *EndPtr;
  int n;
  DATE Instant;

  RequireConfigValue(Input, "RUN", "Model Version", Value);
  if (strcmp(Value, DHSVM_VERSION) != 0)
    ReportErrorCode(103, "[RUN] Model Version = %s: this is DHSVM %s", Value, DHSVM_VERSION);

  RequireConfigValue(Input, "RUN", "Start", Run->StartText);
  if (!SScanIsoDate(Run->StartText, &(Run->Start)))
    ReportErrorCode(103, "[RUN] Start = %s: expected YYYY-MM-DD HH:MM", Run->StartText);

  RequireConfigValue(Input, "RUN", "End", Run->EndText);
  if (!SScanIsoDate(Run->EndText, &(Run->End)))
    ReportErrorCode(103, "[RUN] End = %s: expected YYYY-MM-DD HH:MM", Run->EndText);

  /* Time step in hours */
  RequireConfigValue(Input, "RUN", "Time Step", Run->StepText);
  Run->Seconds = strtod(Run->StepText, &EndPtr) * 3600.0;
  if (EndPtr == Run->StepText || *EndPtr != '\0')
    ReportErrorCode(103, "[RUN] Time Step = %s: expected a number", Run->StepText);

  RequireConfigValue(Input, "RUN", "Output Directory", Value);

  /* Start state: a cold start (NONE) or the snapshot stamped Start */
  RequireConfigValue(Input, "RUN", "Initial State", Value);
  Options->ColdStart = IsChoice(Value, "NONE");

  /* Snapshot instants State Time 1, 2, ... (the list ends at the first
     missing number): the form here, the grid, range and order in
     InitStateDump() */
  for (n = 1; ; n++) {
    snprintf(KeyName, sizeof(KeyName), "State Time %d", n);
    if (!GetConfigValue(Input, "RUN", KeyName, Value))
      break;
    RequireConfigValue(Input, "RUN", KeyName, Value);
    if (!SScanIsoDate(Value, &Instant))
      ReportErrorCode(103, "[RUN] %s = %s: expected YYYY-MM-DD HH:MM", KeyName, Value);
  }

  /* Groundwater spin-up, only from a cold start (checked with the other
     combinations) */
  RequireConfigValue(Input, "RUN", "Spinup Years", Value);
  Options->GW_SPINUP_YRS = ConfigInt("RUN", "Spinup Years", Value);
  if (Options->GW_SPINUP_YRS < 0)
    ReportErrorCode(103, "[RUN] Spinup Years = %s: expected an integer >= 0", Value);
  Options->GW_SPINUP = (Options->GW_SPINUP_YRS > 0);
  if (Options->GW_SPINUP) {
    RequireConfigValue(Input, "RUN", "Spinup Recharge", Value);
    Options->GW_SPINUP_RECHARGE = ConfigFloat("RUN", "Spinup Recharge", Value);
  }
  else {
    GetConfigValue(Input, "RUN", "Spinup Recharge", Value);
    Options->GW_SPINUP_RECHARGE = 0.0;
  }

  RequireConfigValue(Input, "RUN", "Verbose", Value);
  Options->Verbose = ConfigBool("RUN", "Verbose", Value);
  SetVerbose(Options->Verbose);

  /* States that follow from other states rebuilt from the snapshot (absent
     = FALSE; read and ignored at a cold start) */
  Options->RecomputeDerived = FALSE;
  if (GetConfigValue(Input, "RUN", "Recompute Derived State", Value)) {
    RequireConfigValue(Input, "RUN", "Recompute Derived State", Value);
    Options->RecomputeDerived = ConfigBool("RUN", "Recompute Derived State", Value) &&
      !Options->ColdStart;
  }
}

/*****************************************************************************
  InitRunTimes()

  Checks the time step (whole minutes that divide a day), Start (on the step
  grid counted from midnight), End (the end of the last time step, a whole
  number of steps after Start), the snapshot instants, the map times and
  the network dynamics interval; sets up the run's time steps. Then checks
  that the output directory exists, before any input file is read.
*****************************************************************************/
static void InitRunTimes(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map, RUNTIMES *Run,
                         TIMESTRUCT *Time, DUMPSTRUCT *Dump)
{
  char Value[BUFSIZE + 1];
  long Dt;
  long long Minutes;

  /* Not a number, not positive or far above a day: rejected before lround()
     so that it cannot overflow; the rule itself is checked on Dt */
  if (!(Run->Seconds > 0.0 && Run->Seconds < 2.0 * SECPDAY))
    ReportErrorCode(110, "[RUN] Time Step = %s: must divide 24 hours into whole minutes",
                    Run->StepText);
  Dt = lround(Run->Seconds);
  if (!(fabs(Run->Seconds - (double) Dt) < 0.01) || Dt <= 0 || Dt % SECPMIN != 0 ||
      SECPDAY % Dt != 0)
    ReportErrorCode(110, "[RUN] Time Step = %s: must divide 24 hours into whole minutes",
                    Run->StepText);

  if (((long) (Run->Start.Hour * MINPHOUR + Run->Start.Min) * SECPMIN) % Dt != 0)
    ReportErrorCode(110, "[RUN] Start = %s: not on the %s-hour step grid from midnight",
                    Run->StartText, Run->StepText);

  Minutes = MinutesOf(&(Run->End)) - MinutesOf(&(Run->Start));
  if (Minutes <= 0 || (Minutes * SECPMIN) % Dt != 0 ||
      (Minutes * SECPMIN) / Dt > 2147483647LL)
    ReportErrorCode(110, "[RUN] End = %s: must be after Start by a whole number of steps",
                    Run->EndText);

  InitTime(Time, &(Run->Start), &(Run->End), (int) Dt);

  InitStateDump(Input, Time, Dump);
  CheckOutputTimes(Options, Map, Time, Dump);
  InitVegUpdateTimes(Input, Options, Time);

  /* Every output file is written into the output directory */
  RequireConfigValue(Input, "RUN", "Output Directory", Value);
  if (!IsDirectory(Value))
    ReportErrorCode(104, "[RUN] Output Directory %s does not exist", Value);
  DirectoryPath(Dump->Path, sizeof(Dump->Path), Value);

  /* Directory of the snapshot stamped Start, unless a cold start */
  RequireConfigValue(Input, "RUN", "Initial State", Value);
  if (Options->ColdStart)
    Dump->InitStatePath[0] = '\0';
  else
    DirectoryPath(Dump->InitStatePath, sizeof(Dump->InitStatePath), Value);
}

/*****************************************************************************
  DirectoryPath()

  A directory as a path prefix: the value, followed by a separator unless it
  ends in one
*****************************************************************************/
static void DirectoryPath(char *Buffer, size_t Size, const char *Dir)
{
  size_t n = strlen(Dir);

  if (n > 0 && (Dir[n - 1] == '/' || Dir[n - 1] == '\\'))
    MakePath(Buffer, Size, "%s", Dir);
  else
    MakePath(Buffer, Size, "%s/", Dir);
}

/*****************************************************************************
  InitGrid()

  [GRID]: grid size and position, sun geometry location, number of soil
  layers
*****************************************************************************/
static void InitGrid(LISTPTR Input, MAPSIZE *Map, SOLARGEOMETRY *SolarGeo)
{
  char Value[BUFSIZE + 1];

  RequireConfigValue(Input, "GRID", "Rows", Value);
  Map->NY = ConfigInt("GRID", "Rows", Value);
  if (Map->NY < 1 || Map->NY > 1000000)
    ReportErrorCode(103, "[GRID] Rows = %s: expected an integer from 1 to 1000000", Value);

  RequireConfigValue(Input, "GRID", "Columns", Value);
  Map->NX = ConfigInt("GRID", "Columns", Value);
  if (Map->NX < 1 || Map->NX > 1000000)
    ReportErrorCode(103, "[GRID] Columns = %s: expected an integer from 1 to 1000000", Value);

  RequireConfigValue(Input, "GRID", "Cell Size", Value);
  Map->DY = ConfigFloat("GRID", "Cell Size", Value);
  if (!(Map->DY > 0.0))
    ReportErrorCode(103, "[GRID] Cell Size = %s: expected a number > 0", Value);

  RequireConfigValue(Input, "GRID", "North", Value);
  Map->Yorig = ConfigDouble("GRID", "North", Value);

  RequireConfigValue(Input, "GRID", "West", Value);
  Map->Xorig = ConfigDouble("GRID", "West", Value);

  RequireConfigValue(Input, "GRID", "Center Latitude", Value);
  SolarGeo->Latitude = ConfigFloat("GRID", "Center Latitude", Value);
  SolarGeo->Latitude *= (float) RADPDEG;

  RequireConfigValue(Input, "GRID", "Center Longitude", Value);
  SolarGeo->Longitude = ConfigFloat("GRID", "Center Longitude", Value);
  SolarGeo->Longitude *= (float) RADPDEG;

  RequireConfigValue(Input, "GRID", "Time Zone Meridian", Value);
  SolarGeo->StandardMeridian = ConfigFloat("GRID", "Time Zone Meridian", Value);
  SolarGeo->StandardMeridian *= (float) RADPDEG;

  RequireConfigValue(Input, "GRID", "Soil Layers", Value);
  Map->SoilLayers = ConfigInt("GRID", "Soil Layers", Value);
  if (Map->SoilLayers < 1)
    ReportErrorCode(103, "[GRID] Soil Layers = %s: expected an integer >= 1", Value);

  Map->DX = Map->DY;
  Map->DXY = (float) sqrt(Map->DX * Map->DX + Map->DY * Map->DY);
  Map->X = 0;
  Map->Y = 0;
  Map->OffsetX = 0;
  Map->OffsetY = 0;
  Map->NumCells = 0;
}

/*****************************************************************************
  ConfigChoice()

  Index of the value of a required key among Choices (matched exactly,
  ignoring case)
*****************************************************************************/
static int ConfigChoice(LISTPTR Input, const char *Section, const char *Key,
                        const char **Choices, int NChoices)
{
  char Value[BUFSIZE + 1];
  int i;

  RequireConfigValue(Input, Section, Key, Value);
  for (i = 0; i < NChoices; i++) {
    if (IsChoice(Value, Choices[i]))
      return i;
  }
  ReportErrorCode(103, "[%s] %s = %s: expected %s or %s", Section, Key, Value,
                  Choices[0], Choices[1]);
  return 0;
}

/*****************************************************************************
  InitOptions()

  [OPTIONS]: the physics switches (all required)
*****************************************************************************/
static void InitOptions(LISTPTR Input, OPTIONSTRUCT *Options)
{
  const char *Gradient[] = { "TOPOGRAPHY", "WATERTABLE" };
  const char *Routing[] = { "CONVENTIONAL", "KINEMATIC" };
  const char *Ksat[] = { "ANISOTROPY", "TABLE" };
  const char *Infiltration[] = { "STATIC", "DYNAMIC" };
  const char *Stomatal[] = { "JARVIS", "MEDLYN" };
  char Value[BUFSIZE + 1];

  /* Subsurface flow gradient: land surface (TOPOGRAPHY) or water table */
  Options->FlowGradient =
    (ConfigChoice(Input, "OPTIONS", "Gradient", Gradient, 2) == 0) ? TOPOGRAPHY : WATERTABLE;

  /* Multiple flow directions or steepest descent over the 8 neighbors */
  RequireConfigValue(Input, "OPTIONS", "Multiple Flow Directions", Value);
  Options->MultiFlowDir = ConfigBool("OPTIONS", "Multiple Flow Directions", Value);

  /* Overland flow routing method */
  Options->Routing =
    (ConfigChoice(Input, "OPTIONS", "Overland Routing", Routing, 2) == 1) ? TRUE : FALSE;

  RequireConfigValue(Input, "OPTIONS", "Lake Dynamics", Value);
  Options->LakeDynamics = ConfigBool("OPTIONS", "Lake Dynamics", Value);

  /* Unsaturated lateral flow (interflow) */
  RequireConfigValue(Input, "OPTIONS", "Unsaturated Lateral Flow", Value);
  Options->UseInterflow = ConfigBool("OPTIONS", "Unsaturated Lateral Flow", Value);

  /* Vertical conductivity from lateral conductivity and anisotropy, or from
     the soil table */
  Options->UseKsatAnisotropy =
    (ConfigChoice(Input, "OPTIONS", "Vertical Ksat Source", Ksat, 2) == 0) ? TRUE : FALSE;

  /* Static or dynamic maximum infiltration rate */
  Options->Infiltration =
    (ConfigChoice(Input, "OPTIONS", "Infiltration", Infiltration, 2) == 0) ? STATIC : DYNAMIC;

  /* Monthly precipitation pattern map and station normals */
  RequireConfigValue(Input, "OPTIONS", "Precipitation Pattern", Value);
  Options->Prism = ConfigBool("OPTIONS", "Precipitation Pattern", Value);

  /* Snow pattern map to redistribute snowfall */
  RequireConfigValue(Input, "OPTIONS", "Snow Pattern", Value);
  Options->SnowPattern = ConfigBool("OPTIONS", "Snow Pattern", Value);

  /* Melt scale map multiplying the melt of every cell */
  RequireConfigValue(Input, "OPTIONS", "Melt Scale", Value);
  Options->MeltScale = ConfigBool("OPTIONS", "Melt Scale", Value);

  /* Wind scale map multiplying the station wind of every cell */
  RequireConfigValue(Input, "OPTIONS", "Wind Scale", Value);
  Options->WindScale = ConfigBool("OPTIONS", "Wind Scale", Value);

  /* Terrain shading of the beam radiation */
  RequireConfigValue(Input, "OPTIONS", "Shading", Value);
  Options->Shading = ConfigBool("OPTIONS", "Shading", Value);

  /* Canopy conductance: the resistance factors of Jarvis or the
     photosynthesis and stomata of Medlyn */
  Options->StomatalScheme =
    (ConfigChoice(Input, "OPTIONS", "Stomatal Scheme", Stomatal, 2) == 0) ? JARVIS : MEDLYN;
}

/*****************************************************************************
  InitInputKeys()

  [INPUTS]: the keys the run's options require must be present and not
  empty; the other keys are read and ignored. The files are opened and
  checked when they are read. A soil or vegetation map of ClassMaps[]
  without its key is formed from its class table (a key present must not be
  empty).
*****************************************************************************/
static void InitInputKeys(LISTPTR Input, OPTIONSTRUCT *Options)
{
  const char *Always[] = {
    "Terrain.Elevation", "Terrain.Mask", "Soil.Class", "Soil.LayerThickness",
    "Veg.OverstoryClass", "Veg.UnderstoryClass", "Soil Table", "Overstory Table",
    "Understory Table", "Channel Class Table", "Channel Segment Table", "Channel Cell Table",
    "Station Table", "Station Data Directory"
  };
  char Value[BUFSIZE + 1];
  int i;

  for (i = 0; i < (int) (sizeof(Always) / sizeof(Always[0])); i++)
    RequireConfigValue(Input, "INPUTS", Always[i], Value);

  /* The soil and vegetation maps: a map, or formed from the class table */
  for (i = 0; i < NCLASSMAPS; i++) {
    Options->FromTable[i] = !GetConfigValue(Input, "INPUTS", ClassMaps[i].Name, Value);
    if (!Options->FromTable[i] && IsEmptyStr(Value))
      ReportErrorCode(103, "[INPUTS] %s is empty", ClassMaps[i].Name);
  }

  if (Options->Prism)
    RequireConfigValue(Input, "INPUTS", "Met.PrecipPattern", Value);
  else
    GetConfigValue(Input, "INPUTS", "Met.PrecipPattern", Value);

  if (Options->SnowPattern)
    RequireConfigValue(Input, "INPUTS", "Met.SnowPattern", Value);
  else
    GetConfigValue(Input, "INPUTS", "Met.SnowPattern", Value);

  if (Options->MeltScale)
    RequireConfigValue(Input, "INPUTS", "Snow.MeltScale", Value);
  else
    GetConfigValue(Input, "INPUTS", "Snow.MeltScale", Value);

  if (Options->WindScale)
    RequireConfigValue(Input, "INPUTS", "Met.WindScale", Value);
  else
    GetConfigValue(Input, "INPUTS", "Met.WindScale", Value);

  if (Options->LakeDynamics) {
    RequireConfigValue(Input, "INPUTS", "Lake.ID", Value);
    RequireConfigValue(Input, "INPUTS", "Lake Table", Value);
  }
  else {
    GetConfigValue(Input, "INPUTS", "Lake.ID", Value);
    GetConfigValue(Input, "INPUTS", "Lake Table", Value);
  }

  /* Atmospheric CO2 of each year */
  if (Options->StomatalScheme == MEDLYN)
    RequireConfigValue(Input, "INPUTS", "CO2 File", Value);
  else
    GetConfigValue(Input, "INPUTS", "CO2 File", Value);

  /* Initial snow water equivalent: optional, used only at a cold start */
  Options->InitialSWEPath[0] = '\0';
  if (GetConfigValue(Input, "INPUTS", "Snow.InitialSWE", Value) && Options->ColdStart) {
    if (IsEmptyStr(Value))
      ReportErrorCode(103, "[INPUTS] Snow.InitialSWE is empty");
    strcpy(Options->InitialSWEPath, Value);
  }

  /* The point table: checked with the output requests */
  GetConfigValue(Input, "INPUTS", "Point Table", Value);

  /* Vegetation updates: the directory of their maps (optional; the
     instants are read by InitVegUpdateTimes()) */
  Options->VegUpdates.Directory[0] = '\0';
  if (GetConfigValue(Input, "INPUTS", "Vegetation Update Directory", Value)) {
    if (IsEmptyStr(Value))
      ReportErrorCode(103, "[INPUTS] Vegetation Update Directory is empty");
    if (!IsDirectory(Value))
      ReportErrorCode(104, "[INPUTS] Vegetation Update Directory %s does not exist", Value);
    DirectoryPath(Options->VegUpdates.Directory, sizeof(Options->VegUpdates.Directory), Value);
  }

  /* Horizon angles and sky view of the shading */
  if (Options->Shading) {
    RequireConfigValue(Input, "INPUTS", "Terrain.HorizonAngle", Value);
    RequireConfigValue(Input, "INPUTS", "Terrain.SkyView", Value);
  }
  else {
    GetConfigValue(Input, "INPUTS", "Terrain.HorizonAngle", Value);
    GetConfigValue(Input, "INPUTS", "Terrain.SkyView", Value);
  }
}

/* A required [PARAMETERS] value */
static float Parameter(LISTPTR Input, const char *Key)
{
  char Value[BUFSIZE + 1];

  RequireConfigValue(Input, "PARAMETERS", Key, Value);
  return ConfigFloat("PARAMETERS", Key, Value);
}

/* A [PARAMETERS] value required only under some options: read and ignored
   (NOT_APPLICABLE) otherwise */
static float OptionalParameter(LISTPTR Input, const char *Key, int Required)
{
  char Value[BUFSIZE + 1];

  if (Required)
    return Parameter(Input, Key);
  GetConfigValue(Input, "PARAMETERS", Key, Value);
  return NOT_APPLICABLE;
}

/*****************************************************************************
  InitParameters()

  [PARAMETERS]: every scalar value. The
  values are applied to the class and lake tables by ApplyTableParameters(),
  to the maps by ApplyParameters() and to the snow parameters by
  InitSnowParameters().
*****************************************************************************/
static void InitParameters(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map)
{
  /* Soil and vegetation maps */
  SOIL_DEPTH_ADJ = Parameter(Input, "Soil Depth Offset");
  SOIL_KSAT_ADJ = Parameter(Input, "Lateral Conductivity Scale");
  SOIL_EXPDEC_ADJ = Parameter(Input, "Exponential Decrease Scale");
  SOIL_FIELDCAP_ADJ = Parameter(Input, "Field Capacity Scale");
  SOIL_POROSITY_ADJ = Parameter(Input, "Porosity Scale");
  VEG_FC_ADJ = Parameter(Input, "Overstory Cover Scale");
  VEG_LAI_ADJ = Parameter(Input, "Overstory LAI Scale");
  VEG_UNDER_FC_ADJ = Parameter(Input, "Understory Cover Scale");
  VEG_UNDER_LAI_ADJ = Parameter(Input, "Understory LAI Scale");

  /* Class and lake tables */
  SOIL_DEEPFLUX_ADJ = Parameter(Input, "Deep Flux Offset");
  VEG_RSMIN_ADJ = OptionalParameter(Input, "Min Stomatal Resistance Scale",
                                    Options->StomatalScheme == JARVIS);
  VEG_SNOWINT_ADJ = Parameter(Input, "Snow Interception Efficiency Scale");
  VEG_MOIST_ADJ = Parameter(Input, "Moisture Threshold Scale");
  VEG_G1_ADJ = OptionalParameter(Input, "Stomatal Slope Scale", Options->StomatalScheme == MEDLYN);
  LAKE_OUTFLOW_ADJ = OptionalParameter(Input, "Lake Outflow Rate Scale", Options->LakeDynamics);
  LAKE_EXPONENT_ADJ = OptionalParameter(Input, "Lake Outflow Exponent Scale", Options->LakeDynamics);

  /* Snow */
  MIN_RAIN_TEMP = Parameter(Input, "Rain Threshold");
  MAX_SNOW_TEMP = Parameter(Input, "Snow Threshold");
  ALB_MAX = Parameter(Input, "Fresh Snow Albedo");
  ALB_ACC_LAMBDA = Parameter(Input, "Albedo Accumulation Lambda");
  ALB_MELT_LAMBDA = Parameter(Input, "Albedo Melting Lambda");
  LAMBDA_FOREST_OFFSET = Parameter(Input, "Albedo Lambda Forest Offset");
  ALB_ACC_MIN = Parameter(Input, "Albedo Accumulation Min");
  ALB_MELT_MIN = Parameter(Input, "Albedo Melting Min");

  /* Meteorology */
  PRECIP_MULTIPLIER = Parameter(Input, "Precipitation Scale");
  SNOWPAT_WEIGHT = OptionalParameter(Input, "Snow Pattern Weight", Options->SnowPattern);
  TEMPERATURE_OFFSET = Parameter(Input, "Temperature Offset");
  LAPSE_RATE_BIAS = Parameter(Input, "Lapse Rate Bias");
  LAPSE_BIAS_ELEV = Parameter(Input, "Lapse Rate Bias Elevation");

  /* Constants of the energy balance, interception and snow pack */
  Z0_GROUND = Parameter(Input, "Ground Roughness");
  Z0_SNOW = Parameter(Input, "Snow Roughness");
  LIQUID_WATER_CAPACITY = Parameter(Input, "Snow Water Capacity");
  Zref = Parameter(Input, "Reference Height");
  LAI_WATER_MULTIPLIER = Parameter(Input, "Rain LAI Multiplier");
  LAI_SNOW_MULTIPLIER = Parameter(Input, "Snow LAI Multiplier");
  MIN_INTERCEPTION_STORAGE = Parameter(Input, "Min Intercepted Snow");
  MIN_SNOW_RESET_ALBEDO = Parameter(Input, "Min Albedo Reset Snowfall");
  MAX_SURFACE_SWE = Parameter(Input, "Max Surface Snow Layer Depth");

  /* Station interpolation: maximum distance (m), used in cells */
  Options->MaxInterpDistance = Parameter(Input, "Max Interpolation Distance");
  Options->MaxInterpDist = Options->MaxInterpDistance;
  Options->MaxInterpDist /= Map->DY;

  /* Cold start: the water depth in every channel segment and the days
     since snowfall (read and ignored with a snapshot start) */
  Options->InitChannelDepth = OptionalParameter(Input, "Initial Channel Depth", Options->ColdStart);
  Options->InitDaysSinceSnowfall =
    OptionalParameter(Input, "Initial Days Since Snowfall", Options->ColdStart);
  if (!Options->ColdStart) {
    Options->InitChannelDepth = 0.0;
    Options->InitDaysSinceSnowfall = 0.0;
  }
}

/*****************************************************************************
  PrintRunSettings()

  Progress messages: routing directions, the meteorology adjustments, the
  maps formed from the class tables
*****************************************************************************/
void PrintRunSettings(OPTIONSTRUCT *Options)
{
  int m;

  if (Options->MultiFlowDir)
    Message("Using %d neighbors with multiple flow directions for surface/subsurface routing\n", NDIRS);
  else
    Message("Using %d neighbors with steepest descent for surface/subsurface routing\n", NDIRS);
  Message("\nTemperature inputs are adjusted by %.3f deg. C\n",TEMPERATURE_OFFSET);
  Message("Temperature also adjusted by %.3e deg. C / m rel. to %.0f m\n",
          LAPSE_RATE_BIAS,LAPSE_BIAS_ELEV);
  Message("Precipitation inputs are rescaled by a factor of %.3f\n",PRECIP_MULTIPLIER);
  if (Options->WindScale)
    Message("Wind speed scaled by the Met.WindScale map\n");
  Message("Stomatal scheme %s\n", (Options->StomatalScheme == MEDLYN) ? "MEDLYN" : "JARVIS");
  for (m = 0; m < NCLASSMAPS; m++) {
    if (Options->FromTable[m])
      Message("%s formed from the %s (%s%s)\n", ClassMaps[m].Name, ClassMaps[m].Table,
              ClassMaps[m].Column, (ClassMaps[m].LayerSet == MAP_SOIL) ? "_<layer>" : "");
  }
}

/*****************************************************************************
  CheckCombinations()

  Option combinations that are not allowed
*****************************************************************************/
static void CheckCombinations(OPTIONSTRUCT *Options)
{
  if (Options->SnowPattern && !Options->Prism)
    ReportErrorCode(103, "[OPTIONS] Snow Pattern = TRUE needs Precipitation Pattern = TRUE");
  if (Options->GW_SPINUP && !Options->ColdStart)
    ReportErrorCode(103, "[RUN] Spinup Years > 0 needs Initial State = NONE");
}
