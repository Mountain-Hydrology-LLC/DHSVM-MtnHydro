#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "constants.h"
#include "getinit.h"
#include "DHSVMChannel.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "settings.h"
#include "errorhandler.h"
#include "fileio.h"
#include "tables.h"

/* -----------------------------------------------------------------------------
   ReadChannelTables
   Reads the channel class, segment and cell tables ([INPUTS] Channel Class
   Table, Channel Segment Table, Channel Cell Table); their contents are
   checked by ValidateChannels()
   -------------------------------------------------------------------------- */
void ReadChannelTables(LISTPTR Input, CHANNELTABLES *Tables)
{
  const char *Routine = "ReadChannelTables";
  TABLE T;
  int i, c[7];

  memset(Tables, 0, sizeof(CHANNELTABLES));

  /* Channel classes */
  RequireConfigValue(Input, "INPUTS", "Channel Class Table", Tables->ClassPath);
  ReadTable(Tables->ClassPath, &T);
  c[0] = TableColumn(&T, "Class", "-", TRUE);
  c[1] = TableColumn(&T, "Width", "m", TRUE);
  c[2] = TableColumn(&T, "BankHeight", "m", TRUE);
  c[3] = TableColumn(&T, "ManningsN", "-", TRUE);
  TableWarnUnknown(&T);
  Tables->NClasses = T.NRows;
  if (!(Tables->ClassID = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Width = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->BankHeight = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->Friction = calloc(T.NRows + 1, sizeof(float))))
    ReportError((char *) Routine, 1);
  for (i = 0; i < T.NRows; i++) {
    Tables->ClassID[i] = TableInt(&T, i, c[0]);
    Tables->Width[i] = TableFloat(&T, i, c[1]);
    Tables->BankHeight[i] = TableFloat(&T, i, c[2]);
    Tables->Friction[i] = TableFloat(&T, i, c[3]);
  }
  FreeTable(&T);

  /* Channel segments, in routing order within each stream order */
  RequireConfigValue(Input, "INPUTS", "Channel Segment Table", Tables->SegmentPath);
  ReadTable(Tables->SegmentPath, &T);
  c[0] = TableColumn(&T, "Segment", "-", TRUE);
  c[1] = TableColumn(&T, "Order", "-", TRUE);
  c[2] = TableColumn(&T, "Slope", "m/m", TRUE);
  c[3] = TableColumn(&T, "Length", "m", TRUE);
  c[4] = TableColumn(&T, "Class", "-", TRUE);
  c[5] = TableColumn(&T, "Outlet", "-", TRUE);
  c[6] = TableColumn(&T, "SaveName", "-", TRUE);
  TableWarnUnknown(&T);
  Tables->NSegments = T.NRows;
  if (!(Tables->SegmentID = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Order = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Slope = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->Length = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->Class = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Outlet = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->SaveName = calloc(T.NRows + 1, sizeof(char *))))
    ReportError((char *) Routine, 1);
  for (i = 0; i < T.NRows; i++) {
    Tables->SegmentID[i] = TableInt(&T, i, c[0]);
    Tables->Order[i] = TableInt(&T, i, c[1]);
    Tables->Slope[i] = TableFloat(&T, i, c[2]);
    Tables->Length[i] = TableFloat(&T, i, c[3]);
    Tables->Class[i] = TableInt(&T, i, c[4]);
    Tables->Outlet[i] = TableInt(&T, i, c[5]);
    if (!(Tables->SaveName[i] = strdup(T.Fields[i][c[6]])))
      ReportError((char *) Routine, 1);
  }
  FreeTable(&T);

  /* Channel cells, in the order the records of a cell are kept */
  RequireConfigValue(Input, "INPUTS", "Channel Cell Table", Tables->CellPath);
  ReadTable(Tables->CellPath, &T);
  c[0] = TableColumn(&T, "Row", "-", TRUE);
  c[1] = TableColumn(&T, "Column", "-", TRUE);
  c[2] = TableColumn(&T, "Segment", "-", TRUE);
  c[3] = TableColumn(&T, "Length", "m", TRUE);
  c[4] = TableColumn(&T, "CutHeight", "m", TRUE);
  c[5] = TableColumn(&T, "CutWidth", "m", TRUE);
  TableWarnUnknown(&T);
  Tables->NCells = T.NRows;
  if (!(Tables->Row = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Column = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->Segment = calloc(T.NRows + 1, sizeof(int))) ||
      !(Tables->CellLength = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->CutHeight = calloc(T.NRows + 1, sizeof(float))) ||
      !(Tables->CutWidth = calloc(T.NRows + 1, sizeof(float))))
    ReportError((char *) Routine, 1);
  for (i = 0; i < T.NRows; i++) {
    Tables->Row[i] = TableInt(&T, i, c[0]);
    Tables->Column[i] = TableInt(&T, i, c[1]);
    Tables->Segment[i] = TableInt(&T, i, c[2]);
    Tables->CellLength[i] = TableFloat(&T, i, c[3]);
    Tables->CutHeight[i] = TableFloat(&T, i, c[4]);
    Tables->CutWidth[i] = TableFloat(&T, i, c[5]);
  }
  FreeTable(&T);
}

/* -----------------------------------------------------------------------------
   InitChannel
   Builds the network from the channel tables (checked by ValidateChannels())
   and binds it to the grid cells
   -------------------------------------------------------------------------- */
void
InitChannel(MAPSIZE *Map, int deltat, CHANNEL *ChannelData, CHANNELTABLES *Tables,
	    SOILTABLE *SType, SOILPIX ** SoilMap,
	    LAKETABLE *LType, TOPOPIX **TopoMap,
	    int *MaxStreamID, OPTIONSTRUCT *Options)
{
  int i, x, y;
  ChannelMapPtr cell;

  Message("\nInitializing Stream Networks\n");

  ChannelData->stream_class = NULL;
  ChannelData->streams = NULL;
  ChannelData->stream_map = NULL;

  channel_grid_init(Map->NX, Map->NY);

  Message("\tReading Stream data\n");

  ChannelData->stream_class =
    channel_create_classes(Tables->NClasses, Tables->ClassID, Tables->Width,
                           Tables->BankHeight, Tables->Friction);
  ChannelData->streams =
    channel_create_network(Tables->NSegments, Tables->SegmentID, Tables->Order, Tables->Slope,
                           Tables->Length, Tables->Class, Tables->Outlet, Tables->SaveName,
                           ChannelData->stream_class, MaxStreamID);
  ChannelData->stream_map =
    channel_grid_create_records(ChannelData->streams, Tables->NCells, Tables->Row,
                                Tables->Column, Tables->Segment, Tables->CellLength,
                                Tables->CutHeight, Tables->CutWidth, SType, SoilMap);

  /* Associate each channel segment with its constituent map cells */
  channel_combine_map_network(ChannelData->streams, ChannelData->stream_map, Map);

  if (Options->LakeDynamics) {

    /* Find channel segments that intersect lake cells */
    for (y = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++) {
        cell = ChannelData->stream_map[x][y];
        while (cell != NULL) {
          if (TopoMap[y][x].LakeID != 0) {
            cell->channel->IntersectsLake = TRUE;
            cell->channel->lake = &(LType[TopoMap[y][x].LakeID - 1]);
          }
          cell = cell->next;
        }
      }
    }

    /* Find outlets from lakes */
    for (i = 0; i < Map->NumLakes; i++) {
      LType[i].OutletID = (SegmentID) LType[i].Outlet;
      LType[i].outlet = channel_find_segment(ChannelData->streams, LType[i].OutletID);
      LType[i].outlet->IntersectsLake = FALSE;
      LType[i].outlet->lake = NULL;
    }
  } /* End of lake dynamics */

  error_handler(ERRHDL_STATUS,
		"InitChannel: computing stream network routing coefficients");
  channel_routing_parameters(ChannelData->streams, (double) deltat);
}

/* -------------------------------------------------------------
   RouteChannel
   ------------------------------------------------------------- */
void
RouteChannel(CHANNEL *ChannelData, TIMESTRUCT *Time, MAPSIZE *Map,
	    TOPOPIX **TopoMap, SOILPIX **SoilMap,
	     OPTIONSTRUCT *Options, NETSTRUCT **Network, SOILTABLE *SType,
	     VEGPIX **VegMap, EVAPPIX **Evap, LAKETABLE *LType)
{
  int k, x, y;
  float max_bank_height = 0.0;
  float AdjTableDepth;
  float Transmissivity;
  float ChannelTableDepth;
  float Depth;
  float EffThickness;
  float MaxInfiltrationCap;
  int i;
  float StreamInfiltration, StreamEvap;
  
  /* Add IExcess and ChannelInt to stream channels (by descending elevation) */
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      channel_grid_inc_inflow(ChannelData->stream_map, x, y,
                              SoilMap[y][x].IExcess * Map->DX * Map->DY);
      
      channel_grid_satflow(ChannelData->stream_map, x, y);
      
      SoilMap[y][x].ChannelInt += SoilMap[y][x].IExcess;
      SoilMap[y][x].IExcess = 0.0f;
    }
  }
  
  /* Route stream channels */
  /* Account for infiltration out of streams that are above the water table */
  /* Loop thru all of the cells in descending order of elevation */
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      
      /* Only allow stream re-infiltration if local water table is 1 mm below deepest channel */
      /* Update water table depth directly below channel based on lateral diffusion on previous timestep */
      AdjTableDepth = TopoMap[y][x].Dem - SoilMap[y][x].WaterLevel;
      max_bank_height = channel_grid_cell_maxbankht(ChannelData->stream_map, x, y);
      Transmissivity = CalcTransmissivity(AdjTableDepth, max_bank_height,
                                          SoilMap[y][x].KsLat,
                                          SoilMap[y][x].KsLatExp,
                                          SType[SoilMap[y][x].Soil - 1].DepthThresh);
      ChannelTableDepth = channel_grid_table_depth(ChannelData->stream_map, x, y, Time->Dt,
                                                   AdjTableDepth, Transmissivity,
                                                   (SoilMap[y][x].Porosity[Network[y][x].CutBankZone] -
                                                     SoilMap[y][x].FCap[Network[y][x].CutBankZone]),
                                                     Map->DX);
      
      if (ChannelTableDepth > max_bank_height) {
        
        /* Find maximum amount of water that can be added to subsurface */
        /* So that water table immediately below channel just reaches bottom of lowest channel cut */
        /* With soil freezing a layer holding ice takes at most its free pore space */
        MaxInfiltrationCap = 0.0;
        Depth = 0.0;
        for (i = 0; i < SType[SoilMap[y][x].Soil - 1].NLayers && Depth < ChannelTableDepth; i++) {
          if (SoilMap[y][x].LayerThickness[i] < (SoilMap[y][x].Depth - Depth))
            Depth += SoilMap[y][x].LayerThickness[i];
          else
            Depth = SoilMap[y][x].Depth;
          
          if (Depth > max_bank_height) {
            EffThickness = ((Depth - max_bank_height) < SoilMap[y][x].LayerThickness[i] ?
                              (Depth - max_bank_height) : SoilMap[y][x].LayerThickness[i]);
            if (Depth < ChannelTableDepth || (Options->SoilFreezing && SoilMap[y][x].Ice[i] > 0.0))
              MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * EffThickness;
            else {
              EffThickness -= (Depth - ChannelTableDepth);
              MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) * EffThickness;
            }
          }
        }
        /* Also add deep layer water capacity if water table is below root zone layers */
        if (ChannelTableDepth > Depth) {
          i = SType[SoilMap[y][x].Soil - 1].NLayers;
          if (Options->SoilFreezing && SoilMap[y][x].Ice[i] > 0.0)
            MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].Moist[i]) * (ChannelTableDepth - Depth);
          else
            MaxInfiltrationCap += (SoilMap[y][x].Porosity[i] - SoilMap[y][x].FCap[i]) * (ChannelTableDepth - Depth);
        }
      } else
        MaxInfiltrationCap = 0.0;
      
      channel_grid_calc_infiltration(ChannelData->stream_map, x, y, Time->Dt,
                                     AdjTableDepth, MaxInfiltrationCap, Map->DX);
    }
  }
  
  /* Route lakes, then route streams */
  if (Options->LakeDynamics) {
    for (i = 0; i < Map->NumLakes; i++) {
      LType[i].outlet->lake_inflow += LType[i].Outflow * LType[i].Area;
    }
  }
  
  channel_route_network(ChannelData->streams, Time->Dt);
  
  for (k = (Map->NumCells - 1); k > -1;  k--) {
    y = Map->OrderedCells[k].y;
    x = Map->OrderedCells[k].x;
    if (channel_grid_has_channel(ChannelData->stream_map, x, y) && TopoMap[y][x].LakeID == 0) {
      
      StreamInfiltration = channel_grid_infiltration(ChannelData->stream_map, x, y);
      StreamInfiltration /= Map->DX * Map->DY;
      SoilMap[y][x].SatFlow += StreamInfiltration;
      SoilMap[y][x].ChannelInfiltration += StreamInfiltration;
      
      StreamEvap = channel_grid_evaporation(ChannelData->stream_map, x, y);
      StreamEvap /= Map->DX * Map->DY;
      VegMap[y][x].MoistureFlux += StreamEvap;
      Evap[y][x].ETot += StreamEvap;
      Evap[y][x].EvapChannel = StreamEvap;
    }
  }
}

/* -------------------------------------------------------------
   ChannelCut
   computes necessary parameters for cell storage adjustment from
   channel dimensions
   ------------------------------------------------------------- */
void ChannelCut(int y, int x, CHANNEL * ChannelData, NETSTRUCT * Network)
{
  float bank_height = 0.0;
  float cut_area = 0.0;

  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    bank_height = channel_grid_cell_bankht(ChannelData->stream_map, x, y);
    cut_area = channel_grid_cell_width(ChannelData->stream_map, x, y) * 
      channel_grid_cell_length(ChannelData->stream_map, x, y);
  }
  Network->Area = cut_area;
  Network->BankHeight = bank_height;
}

