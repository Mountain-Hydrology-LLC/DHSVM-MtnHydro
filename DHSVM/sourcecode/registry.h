#ifndef REGISTRY_H
#define REGISTRY_H

#include <stddef.h>

/* Layer sets of the variables */
enum LAYERSET {
  SET_NONE, SET_SOIL, SET_ROOT, SET_CANOPY, SET_CANOPYSOIL, SET_CANOPYROOT,
  SET_MONTH, SET_AZIMUTH, SET_NODE
};

/* Output kinds that accept a variable (bits) */
#define KIND_BASIN       0x01
#define KIND_POINT       0x02
#define KIND_CLASS       0x04
#define KIND_MAP         0x08
#define KIND_STREAMFLOW  0x10
#define KIND_CHANNELS    0x20
#define KIND_NETWORK     0x40

/* One variable: name and nine attributes, as text */
typedef struct {
  const char *Name;
  const char *Layers;
  const char *Units;
  const char *NumberType;
  const char *Dimension;
  const char *Timing;
  const char *SeriesFlux;
  const char *Aggregation;
  const char *Outputs;
  const char *Snapshot;
} REGISTRYROW;

int RegistryCount(void);
const REGISTRYROW *RegistryRow(int Var);
int RegistryFind(const char *Name);
int RegistryLayerSet(int Var);
int RegistryLayerCount(int Var, int SoilLayers);
void RegistryLayerName(int Var, int Layer, int SoilLayers, char *Buffer, size_t Size);
int RegistryFindLayer(const char *Name, int SoilLayers, int *Var, int *Layer);
int RegistryKinds(int Var);
int RegistryIsInterval(int Var);
int RegistryGroupMatches(int Var, const char *Group, size_t Length);

#endif
