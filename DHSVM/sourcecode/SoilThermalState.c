/*
 * DESCRIPTION:  State of the soil temperature nodes (SoilFreezing): the
 *               bottom boundary temperature of each cell, the cold start of
 *               the nodes, and the layer ice, layer temperatures and frost
 *               depth that follow from the node state.
 */
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "soiltemp.h"

/*****************************************************************************
  InitSoilThermal()

  The node grid, and the bottom boundary temperature of every mask cell:
  the stations' mean air temperature, weighted and lapsed to the cell like
  the air temperature (each station's lapse rate the mean of its monthly
  rates), with the temperature offset and the lapse rate bias
*****************************************************************************/
void InitSoilThermal(MAPSIZE *Map, TOPOPIX **TopoMap, SOILPIX **SoilMap,
                     METLOCATION *Stat, int NStats, uchar ***MetWeights)
{
  float *Lapse;
  float WeightSum, TBot;
  int i, k, x, y;

  InitSoilNodes();

  if (!(Lapse = (float *) calloc(NStats, sizeof(float))))
    ReportError("InitSoilThermal", 1);
  for (i = 0; i < NStats; i++) {
    for (k = 0; k < 12; k++)
      Lapse[i] += Stat[i].LapseRate[k];
    Lapse[i] /= 12.;
  }

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (!INBASIN(TopoMap[y][x].Mask))
        continue;
      WeightSum = 0.0;
      for (i = 0; i < NStats; i++)
        WeightSum += (float) MetWeights[y][x][i];
      TBot = 0.0;
      for (i = 0; i < NStats; i++)
        TBot += ((float) MetWeights[y][x][i]) / WeightSum *
          LapseT(Stat[i].MeanAirTemp, Stat[i].Elev, TopoMap[y][x].Dem, Lapse[i]);
      TBot += TEMPERATURE_OFFSET;
      TBot += LAPSE_RATE_BIAS * (TopoMap[y][x].Dem - LAPSE_BIAS_ELEV);
      SoilMap[y][x].TBot = TBot;
    }
  }
  free(Lapse);
}

/*****************************************************************************
  ColdStartSoilThermal()

  Cold start of a cell's nodes: every node and the soil surface at the
  bottom temperature, each node's ice in equilibrium with it
*****************************************************************************/
void ColdStartSoilThermal(SOILPIX *Soil, SOILTABLE *SType)
{
  SOILNODES Nodes;
  int k;

  for (k = 0; k < SOIL_NODES; k++)
    Soil->IceNode[k] = 0.0;
  NodeMoisture(SType->NLayers, Soil->LayerThickness, Soil->Depth, Soil->Moist,
               Soil->Porosity, Soil->Press, Soil->PoreDist, SType->Quartz,
               SType->Organic, Soil->IceNode, &Nodes);
  for (k = 0; k < SOIL_NODES; k++) {
    Soil->TNode[k] = Soil->TBot;
    Soil->IceNode[k] = NodeIce(Soil->TBot, Nodes.Moist[k], Nodes.Porosity[k],
                               Nodes.Press[k], Nodes.PoreDist[k]);
  }
  Soil->TSurf = Soil->TBot;
}

/*****************************************************************************
  SoilThermalDerived()

  The layer temperatures, layer ice (at most the layer's moisture, at
  least 0), frozen zones and frost depth (the bottom of the deepest zone,
  0 without frost) of a cell from its node temperatures and ice
*****************************************************************************/
void SoilThermalDerived(SOILPIX *Soil, int NLayers)
{
  double Water[SOIL_NODES];
  int l;

  LayerAverages(Soil->TNode, Soil->IceNode, NLayers, Soil->LayerThickness,
                Soil->Depth, Soil->Moist, Soil->TLayer, Soil->Ice);
  for (l = 0; l <= NLayers; l++) {
    if (Soil->Ice[l] < 0.)
      Soil->Ice[l] = 0.;
  }
  NodeWater(NLayers, Soil->LayerThickness, Soil->Depth, Soil->Moist, Soil->Porosity, Water);
  Soil->NFrost = FrostZones(Soil->TNode, Soil->IceNode, Water, Soil->Depth,
                            Soil->FrostTop, Soil->FrostBottom);
  Soil->FrostDepth = (Soil->NFrost > 0) ? Soil->FrostBottom[Soil->NFrost - 1] : 0.0;
}

/*****************************************************************************
  UpdateSoilThermalDerived()

  SoilThermalDerived() for every mask cell: at the end of every time step,
  so the outputs and a run restarted from a snapshot see the same layer
  values
*****************************************************************************/
void UpdateSoilThermalDerived(MAPSIZE *Map, TOPOPIX **TopoMap, SOILPIX **SoilMap,
                              SOILTABLE *SType)
{
  int x, y;

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask))
        SoilThermalDerived(&(SoilMap[y][x]), SType[SoilMap[y][x].Soil - 1].NLayers);
    }
  }
}
