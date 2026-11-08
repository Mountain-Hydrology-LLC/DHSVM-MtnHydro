/*
 * DESCRIPTION:  The variable registry: every quantity the model reads as a
 *               map, writes as an output or keeps in a snapshot, one row per
 *               variable (name, layer set, units, number type, dimension,
 *               timing, amount per step, how basin and class means are
 *               formed, the output kinds that accept it, snapshot).
 */
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "registry.h"

static const REGISTRYROW Rows[] = {
/* REGISTRY TABLE BEGIN */
  {"Terrain.Elevation", "-", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Terrain.Mask", "-", "-", "BYTE", "cell", "Static", "N", "-", "Input", "N"},
  {"Terrain.HorizonAngle", "Azimuth", "degrees", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Terrain.SkyView", "-", "fraction", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.Class", "-", "-", "BYTE", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.Depth", "-", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.LateralConductivity", "-", "m/s", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.ExponentialDecrease", "-", "1/m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.Porosity", "Soil", "m3/m3", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.FieldCapacity", "Soil", "m3/m3", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Soil.LayerThickness", "Root", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.OverstoryClass", "-", "-", "BYTE", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.UnderstoryClass", "-", "-", "BYTE", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.OverstoryCover", "-", "fraction", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.OverstoryLAI", "-", "m2/m2", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.OverstoryHeight", "-", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.UnderstoryCover", "-", "fraction", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.UnderstoryLAI", "-", "m2/m2", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Veg.UnderstoryHeight", "-", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Met.PrecipPattern", "Month", "-", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Met.SnowPattern", "-", "-", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Snow.MeltScale", "-", "-", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Met.WindScale", "-", "-", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Lake.ID", "-", "-", "BYTE", "cell", "Static", "N", "-", "Input", "N"},
  {"Snow.InitialSWE", "-", "m", "FLOAT", "cell", "Static", "N", "-", "Input", "N"},
  {"Met.Precipitation", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Met.Snowfall", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Met.Rainfall", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Met.AirTemperature", "-", "C", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Met.WindSpeed", "-", "m/s", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Met.RelativeHumidity", "-", "%", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Met.PrecipitationWaterYear", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Met.PrecipitationStart", "-", "-", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Rad.ShortwaveInObserved", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.ShortwaveIn", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.NetShortwave", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.LongwaveIn", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.CanopyNetShortwave", "Canopy", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.CanopyLongwaveIn", "Canopy", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Rad.NetRadiation", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Snow.Covered", "-", "fraction", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Snow.DaysSinceSnowfall", "-", "days", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.SWE", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Snow.Outflow", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Snow.PackLiquid", "-", "m", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.PackTemperature", "-", "C", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.SurfaceLiquid", "-", "m", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.SurfaceTemperature", "-", "C", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.ColdContent", "-", "J/m2", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.AccumulationSeason", "-", "-", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Snow.Albedo", "-", "fraction", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "Y"},
  {"Snow.Depth", "-", "m", "FLOAT", "cell", "Instant", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.Vapor", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Snow.NetShortwave", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.NetLongwave", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.SensibleHeat", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.LatentHeat", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.RainHeat", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.MeltEnergy", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "SnowMean", "Basin Point Class Map", "N"},
  {"Snow.PeakSWE", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Snow.PeakSWEDay", "-", "day", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Snow.MeltOutDay", "-", "day", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Snow.AccumulationWaterYear", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Snow.OutflowWaterYear", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Canopy.InterceptedRain", "Canopy", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Canopy.InterceptedSnow", "Canopy", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Canopy.TemporaryStorage", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Point Map", "Y"},
  {"Canopy.Temperature", "-", "C", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Canopy.Storage", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Canopy.SnowVapor", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Canopy.Resistance", "Canopy", "s/m", "FLOAT", "cell", "Interval", "N", "-", "Point Map", "N"},
  {"Canopy.Assimilation", "Canopy", "umol/m2/s", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Canopy.MoistureStress", "Canopy", "fraction", "FLOAT", "cell", "Interval", "N", "-", "Point Map", "N"},
  {"Canopy.Acclimation", "-", "fraction", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Canopy.AcclimationTemperature", "-", "C", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Canopy.GrowthTemperature", "-", "C", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Evap.Total", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Potential", "CanopySoil", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Actual", "CanopySoil", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Interception", "Canopy", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Transpiration", "CanopyRoot", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Soil", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Evap.Channel", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Soil.Moisture", "Soil", "m3/m3", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.Percolation", "Root", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.Interflow", "Root", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Soil.TableDepth", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.SatFlow", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.WaterLevel", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Soil.WaterLevelLast", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Soil.SurfaceWater", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.InfiltrationAccumulated", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.InfiltrationStartMoisture", "-", "m3/m3", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Soil.DesorptionAccumulated", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Soil.Water", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.DeepFlux", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Soil.LatentHeat", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.Ice", "Soil", "m3/m3", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.Temperature", "Soil", "C", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.NodeTemperature", "Node", "C", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Soil.NodeIce", "Node", "m3/m3", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Soil.SurfaceTemperature", "-", "C", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Soil.FrostDepth", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.FrostTop", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.GroundHeatFlux", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Soil.SensibleHeatFlux", "-", "W/m2", "FLOAT", "cell", "Interval", "N", "Mean", "Basin Point Class Map", "N"},
  {"Surface.Detention", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "Y"},
  {"Surface.Runoff", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "Y"},
  {"Surface.RunoffStart", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Stream.Interception", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Stream.Infiltration", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Lake.Outflow", "-", "m/step", "FLOAT", "cell", "Instant", "N", "-", "-", "Y"},
  {"Lake.Storage", "-", "m", "FLOAT", "cell", "Instant", "N", "Mean", "Basin Point Class Map", "N"},
  {"Lake.Level", "-", "m", "FLOAT", "cell", "Instant", "N", "-", "Point Map", "N"},
  {"Balance.WaterInput", "-", "m/step", "FLOAT", "cell", "Interval", "Y", "Mean", "Basin Point Class Map", "N"},
  {"Balance.Storage", "-", "m", "FLOAT", "basin", "Instant", "N", "Basin", "Basin", "N"},
  {"Balance.Error", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Balance.CumulativeError", "-", "m", "FLOAT", "basin", "Instant", "N", "Basin", "Basin", "N"},
  {"Balance.RainOrMelt", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Soil.SaturatedFraction", "-", "fraction", "FLOAT", "basin", "Instant", "N", "Basin", "Basin", "N"},
  {"Channel.NetworkLateralInflow", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkOutflow", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkStorageChange", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkInfiltration", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkEvaporation", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkError", "-", "m/step", "FLOAT", "basin", "Interval", "Y", "Basin", "Basin", "N"},
  {"Channel.NetworkStorage", "-", "m", "FLOAT", "basin", "Instant", "N", "Basin", "Basin", "N"},
  {"Channel.Inflow", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Channels NetworkDynamics", "N"},
  {"Channel.LateralInflow", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Channels NetworkDynamics", "N"},
  {"Channel.Outflow", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Streamflow Channels NetworkDynamics", "N"},
  {"Channel.Storage", "-", "m3", "FLOAT", "segment", "Instant", "N", "-", "Channels NetworkDynamics", "Y"},
  {"Channel.StorageChange", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Channels", "N"},
  {"Channel.Infiltration", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Channels NetworkDynamics", "N"},
  {"Channel.Evaporation", "-", "m3/s", "FLOAT", "segment", "Interval", "N", "-", "Channels NetworkDynamics", "N"},
  {"Channel.RoutingK", "-", "1/s", "FLOAT", "segment", "Instant", "N", "-", "-", "Y"},
  {"Channel.TopWaterDepth", "-", "m", "FLOAT", "segment", "Instant", "N", "-", "-", "Y"},
  {"Channel.CellTableDepth", "-", "m", "FLOAT", "channelcell", "Instant", "N", "-", "-", "Y"},
  {"Channel.CellAvailableWater", "-", "m", "FLOAT", "channelcell", "Instant", "N", "-", "-", "Y"},
/* REGISTRY TABLE END */
};

static const char *MonthTokens[12] = {
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
};

/*****************************************************************************
  RegistryCount(), RegistryRow(), RegistryFind()

  Number of variables, the row of a variable (0-based, table order), and the
  variable named Name (-1 if none)
*****************************************************************************/
int RegistryCount(void)
{
  return (int) (sizeof(Rows) / sizeof(Rows[0]));
}

const REGISTRYROW *RegistryRow(int Var)
{
  return &(Rows[Var]);
}

int RegistryFind(const char *Name)
{
  int i;

  for (i = 0; i < RegistryCount(); i++) {
    if (strcmp(Rows[i].Name, Name) == 0)
      return i;
  }
  return -1;
}

/*****************************************************************************
  RegistryLayerSet(), RegistryLayerCount()

  Layer set of a variable and its number of layers (1 for a variable
  without layers); SoilLayers is the number of root zone layers
*****************************************************************************/
int RegistryLayerSet(int Var)
{
  const char *Set = Rows[Var].Layers;

  if (strcmp(Set, "Soil") == 0)
    return SET_SOIL;
  if (strcmp(Set, "Root") == 0)
    return SET_ROOT;
  if (strcmp(Set, "Canopy") == 0)
    return SET_CANOPY;
  if (strcmp(Set, "CanopySoil") == 0)
    return SET_CANOPYSOIL;
  if (strcmp(Set, "CanopyRoot") == 0)
    return SET_CANOPYROOT;
  if (strcmp(Set, "Month") == 0)
    return SET_MONTH;
  if (strcmp(Set, "Azimuth") == 0)
    return SET_AZIMUTH;
  if (strcmp(Set, "Node") == 0)
    return SET_NODE;
  return SET_NONE;
}

int RegistryLayerCount(int Var, int SoilLayers)
{
  switch (RegistryLayerSet(Var)) {
  case SET_SOIL:
    return SoilLayers + 1;
  case SET_ROOT:
    return SoilLayers;
  case SET_CANOPY:
    return 2;
  case SET_CANOPYSOIL:
    return 3;
  case SET_CANOPYROOT:
    return 2 * SoilLayers;
  case SET_MONTH:
    return 12;
  case SET_AZIMUTH:
    return 36;
  case SET_NODE:
    return SOIL_NODES;
  default:
    return 1;
  }
}

/*****************************************************************************
  RegistryLayerName()

  Full name of layer Layer (0-based) of a variable: the name itself without
  layers, else <Name>.<token> (1..L, Deep; Overstory, Understory, Soil;
  Overstory_1..L, Understory_1..L; Jan..Dec; Az000..Az350; 1..SOIL_NODES)
*****************************************************************************/
void RegistryLayerName(int Var, int Layer, int SoilLayers, char *Buffer, size_t Size)
{
  const char *Canopy[3] = { "Overstory", "Understory", "Soil" };
  const char *Name = Rows[Var].Name;

  switch (RegistryLayerSet(Var)) {
  case SET_SOIL:
    if (Layer < SoilLayers)
      snprintf(Buffer, Size, "%s.%d", Name, Layer + 1);
    else
      snprintf(Buffer, Size, "%s.Deep", Name);
    break;
  case SET_ROOT:
  case SET_NODE:
    snprintf(Buffer, Size, "%s.%d", Name, Layer + 1);
    break;
  case SET_CANOPY:
  case SET_CANOPYSOIL:
    snprintf(Buffer, Size, "%s.%s", Name, Canopy[Layer]);
    break;
  case SET_CANOPYROOT:
    snprintf(Buffer, Size, "%s.%s_%d", Name, Canopy[Layer / SoilLayers],
             Layer % SoilLayers + 1);
    break;
  case SET_MONTH:
    snprintf(Buffer, Size, "%s.%s", Name, MonthTokens[Layer]);
    break;
  case SET_AZIMUTH:
    snprintf(Buffer, Size, "%s.Az%03d", Name, 10 * Layer);
    break;
  default:
    snprintf(Buffer, Size, "%s", Name);
    break;
  }
}

/*****************************************************************************
  RegistryFindLayer()

  The variable and layer of a full layer name (<Name>.<token>); TRUE if
  Name is one
*****************************************************************************/
int RegistryFindLayer(const char *Name, int SoilLayers, int *Var, int *Layer)
{
  char Buffer[128];
  const char *Dot;
  size_t Length;
  int i, k;

  Dot = strrchr(Name, '.');
  if (Dot == NULL)
    return 0;
  Length = (size_t) (Dot - Name);
  for (i = 0; i < RegistryCount(); i++) {
    if (RegistryLayerSet(i) == SET_NONE || strlen(Rows[i].Name) != Length ||
        strncmp(Rows[i].Name, Name, Length) != 0)
      continue;
    for (k = 0; k < RegistryLayerCount(i, SoilLayers); k++) {
      RegistryLayerName(i, k, SoilLayers, Buffer, sizeof(Buffer));
      if (strcmp(Buffer, Name) == 0) {
        *Var = i;
        *Layer = k;
        return 1;
      }
    }
  }
  return 0;
}

/*****************************************************************************
  RegistryKinds()

  The output kinds that accept a variable (KIND_* bits)
*****************************************************************************/
int RegistryKinds(int Var)
{
  const char *Words[7] = { "Basin", "Point", "Class", "Map", "Streamflow", "Channels",
                           "NetworkDynamics" };
  const char *Outputs = Rows[Var].Outputs;
  const char *p;
  size_t n;
  int Kinds = 0;
  int j;

  for (j = 0; j < 7; j++) {
    n = strlen(Words[j]);
    for (p = strstr(Outputs, Words[j]); p != NULL; p = strstr(p + 1, Words[j])) {
      if ((p == Outputs || p[-1] == ' ') && (p[n] == ' ' || p[n] == '\0')) {
        Kinds |= (1 << j);
        break;
      }
    }
  }
  return Kinds;
}

/* TRUE for a value over a time step (Timing Interval) */
int RegistryIsInterval(int Var)
{
  return strcmp(Rows[Var].Timing, "Interval") == 0;
}

/* TRUE if the variable's group (the part of its name before the dot) is the
   Length characters of Group */
int RegistryGroupMatches(int Var, const char *Group, size_t Length)
{
  return strncmp(Rows[Var].Name, Group, Length) == 0 && Rows[Var].Name[Length] == '.';
}
