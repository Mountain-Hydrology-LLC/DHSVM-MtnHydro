#ifndef MAPIO_H
#define MAPIO_H

#include <stdio.h>
#include "data.h"
#include "getinit.h"

#define MAP_NAMESIZE 64		/* Longest layer name + 1 */
#define MAP_NODATA_BYTE 255		/* No-data value of a BYTE map */
#define MAP_NODATA_FLOAT (-9999.0f)	/* No-data value of a FLOAT map */

/* Layer sets of the input maps */
#define MAP_ONE   0			/* one layer, named like the map */
#define MAP_SOIL  1			/* soil layers 1..L and Deep */
#define MAP_MONTH 2			/* Jan..Dec */
#define MAP_AZIMUTH 3			/* Az000..Az350 */
#define MAP_ROOT  4			/* soil layers 1..L (no Deep) */

#define MAP_NAZIMUTHS 36		/* Azimuths of a horizon angle map */

/* The soil and vegetation maps a class table can give (CM_* order) */
extern const CLASSMAP ClassMaps[NCLASSMAPS];

void CheckMap(MAPFILE *File, const char *Path, const char *Name, int NumberType,
              const char *Units, int NLayers, char (*LayerNames)[MAP_NAMESIZE],
              const char *Time, MAPSIZE *Map, int Code);
void CloseMap(MAPFILE *File);
void CloseMapFile(FILE *File, const char *Path);
FILE *CreateMap(const char *Path, const char *Name, const char *Units,
                int NLayers, char (*LayerNames)[MAP_NAMESIZE],
                const char *Time, MAPSIZE *Map);
int MapLayerCount(int LayerSet, int SoilLayers);
void MapLayerNames(const char *Name, int LayerSet, int SoilLayers,
                   char (*LayerNames)[MAP_NAMESIZE]);
void OpenInputMap(MAPFILE *File, const char *Path, const char *Name,
                  int NumberType, const char *Units, int LayerSet, MAPSIZE *Map);
void ReadInputMap(LISTPTR Input, const char *Name, int NumberType,
                  const char *Units, int LayerSet, MAPSIZE *Map, void *Values);
void ReadMapLayer(MAPFILE *File, int Layer, void *Values, MAPSIZE *Map);
void WriteMapLayer(FILE *File, float *Values, MAPSIZE *Map, const char *Path);

#endif
