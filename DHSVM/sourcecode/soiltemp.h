#ifndef SOILTEMP_H
#define SOILTEMP_H

#include "constants.h"

/* Soil and thermal properties of the soil temperature nodes for one cell
   and time step (SoilFreezing). Nodes below the soil depth are saturated
   and take the deep layer's properties. */
typedef struct {
  double Moist[SOIL_NODES];     /* Total water, liquid and ice (m3/m3) */
  double Porosity[SOIL_NODES];  /* Porosity (m3/m3) */
  double Press[SOIL_NODES];     /* Brooks-Corey bubbling pressure (m) */
  double PoreDist[SOIL_NODES];  /* Brooks-Corey pore size distribution index */
  double Quartz[SOIL_NODES];    /* Quartz fraction of the mineral solids */
  double Organic[SOIL_NODES];   /* Organic fraction of the solids by volume */
  double Kappa[SOIL_NODES];     /* Thermal conductivity at the start of the step (W/(m C)) */
  double Cs[SOIL_NODES];        /* Volumetric heat capacity at the start of the step (J/(m3 C)) */
} SOILNODES;

extern double SoilNodeDepth[SOIL_NODES];  /* Depth of each node (m), set by InitSoilNodes() */

/* SoilProperties.c */
double UnfrozenWater(double T, double Porosity, double Press, double PoreDist);
double SoilConductivity(double Moist, double Liquid, double Porosity,
                        double Quartz, double Organic);
double SoilHeatCapacity(double Porosity, double Liquid, double Ice,
                        double Organic);
double NodeIce(double T, double Moist, double Porosity, double Press,
               double PoreDist);
void NodeWater(int NLayers, float *LayerThickness, float SoilDepth,
               float *Moist, float *Porosity, double *Water);
void NodeMoisture(int NLayers, float *LayerThickness, float SoilDepth,
                  float *Moist, float *Porosity, float *Press, float *PoreDist,
                  float *Quartz, float *Organic, float *IceNode,
                  SOILNODES *Nodes);

/* SoilProfile.c */
void InitSoilNodes(void);
double SolveProfile(float *TNode, float *IceNode, float TBot, SOILNODES *Nodes,
                    double Dt, double H0, double Beta, double *Storage,
                    double *Rest);
int FrostZones(float *TNode, float *IceNode, double *Water, float SoilDepth,
               float *TopZ, float *BottomZ);
void LayerAverages(float *TNode, float *IceNode, int NLayers,
                   float *LayerThickness, float SoilDepth, float *Moist,
                   float *TLayer, float *Ice);

#endif
