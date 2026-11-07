
#ifndef DATA_H
#define DATA_H

#include "settings.h"
#include "Calendar.h"
#include "channel.h"

typedef struct {
  int N;			/* Northing */
  int E;			/* Easting */
} COORD;

typedef struct {
  char FileName[BUFSIZE * 2 + 1];
  FILE *FilePtr;
} FILES;

#define OUT_NAMESIZE 64		/* Longest variable layer name + 1 */

/* One column of an output: a layer of a variable */
typedef struct {
  int Var;			/* Row of the variable registry */
  int Layer;			/* Layer within the variable's layer set (0 without layers) */
  int Code;			/* What the value is (VALUECODE, output.h) */
  int Instant;			/* TRUE: a value at an instant, FALSE: over a time step */
  char Name[OUT_NAMESIZE];	/* Name of the layer */
} OUTVAR;

/* The value of an Instant column at the time of the next series row */
typedef struct {
  double Value;
  int Has;			/* FALSE: NA (the cell lacks the layer) */
} HELDVALUE;

/* A variable list of [OUTPUT]: its values as read (the keys Basin, Basin 2,
   ..., then the request table's rows) and the columns they expand to */
typedef struct {
  int NKeys;
  char (*Key)[2 * BUFSIZE + 1];	/* Origin of each value, as messages name it */
  char (*Value)[BUFSIZE + 1];
  int N;
  OUTVAR *Vars;
} OUTLIST;

/* The variable lists of [OUTPUT] */
enum OUTLISTS { OUT_BASIN, OUT_POINT, OUT_CLASS, OUT_CHANNELS, OUT_MAPS, NOUTLISTS };

typedef struct {
  char Name[BUFSIZE + 1];	/* Point name (file Point_<Name>.csv) */
  int Line;			/* Line of the point table */
  double North;			/* Coordinates as in the point table */
  double East;
  int y;			/* Row of the point's cell */
  int x;			/* Column of the point's cell */
  FILES File;
  HELDVALUE *Held;		/* Its Instant columns at the next row's time */
} POINTOUT;

/* A combination of an overstory and an understory class with means: every
   mask cell with these two classes */
typedef struct {
  int Over;			/* Overstory class (0 = none) */
  int Under;			/* Understory class (0 = none) */
  int Present;			/* TRUE if a vegetation update gives the combination to a mask cell */
  FILES File;			/* Class_<Over>_<Under>.csv */
  HELDVALUE *Held;		/* Its Instant columns at the next row's time */
} CLASSOUT;

typedef struct {
  char FileName[BUFSIZE * 2 + 1];	/* <Output Directory>/<layer name>.map */
  int *Step;			/* Time steps done when each layer is written */
} MAPOUT;

typedef struct {
  char Path[BUFSIZE * 2 + 1];		/* Output directory, ending in a separator */
  char InitStatePath[BUFSIZE * 2 + 1];	/* Directory of the initial snapshot, ending in a separator */
  int NStates;				/* Number of snapshots */
  DATE *DState;				/* Instants of the snapshots */
  int *StateStep;			/* Number of time steps done at each snapshot */
  OUTLIST List[NOUTLISTS];		/* Basin, Point, Class, Channels and Maps lists */
  char ClassesText[BUFSIZE + 1];	/* [OUTPUT] Classes as read */
  int NClasses;				/* Class combinations with a Class_<O>_<U>.csv */
  CLASSOUT *Classes;
  int NPoints;				/* Points of the point table */
  POINTOUT *Points;
  int Streamflow;			/* TRUE: Streamflow.csv */
  char NetworkText[BUFSIZE + 1];	/* [OUTPUT] Network Dynamics Interval as read, empty if absent */
  int NetworkSteps;			/* Network dynamics interval in time steps, 0 = none */
  int NMapTimes;			/* Layer times of the map stacks */
  DATE *MapTime;
  MAPOUT *Maps;				/* One per column of the Maps list */
  FILES Basin;				/* Basin.csv */
  FILES StreamflowFile;			/* Streamflow.csv */
  FILES Channels;			/* Channels.csv */
  FILES Network;			/* NetworkDynamics.csv */
  HELDVALUE *BasinHeld;		/* Instant columns at the next row's time: Basin.csv, */
  HELDVALUE *ChannelsHeld;		/* Channels.csv (by segment, then column) and */
  HELDVALUE *NetworkHeld;		/* NetworkDynamics.csv (by segment, then column) */
} DUMPSTRUCT;

typedef struct {
  float ETot;			/* Total amount of evapotranspiration */
  float *EPot;			/* Potential transpiration from each vegetation/soil layer */
  float *EAct;			/* Actual transpiration from each vegetation soil layer */
  float *EInt;			/* Evaporation from interception for each vegetation layer */
  float **ESoil;		/* Transpiration for each vegetation layer from each soil zone */
  float EvapSoil;		/* Evaporation from the upper soil layer */
  float EvapChannel; /* Evaporation from the stream channel */
  float Rc[MAXVEGLAYERS];	/* Canopy resistance of each vegetation layer (s/m), NA if none */
  float An[MAXVEGLAYERS];	/* Net assimilation of each vegetation layer per cell area
				   (umol/m2/s) (MEDLYN) */
  float Beta[MAXVEGLAYERS];	/* Soil moisture factor of each vegetation layer, NA if none */
} EVAPPIX;

typedef struct {
  int TimeStep;
  float Fraction;
} UNITHYDR;

typedef struct {
  int MaxTravelTime;
  int TotalWaveLength;
  int *WaveLength;
} UNITHYDRINFO;

typedef struct {
  char Const[BUFSIZE + 1];	/* Filename for main input file  */
} INPUTFILES;

typedef struct {
  int NTypes;
  int *NLayers;
  int MaxLayers;
} LAYER;

typedef struct {
  float Tair;			/* Air temperature (C) */
  float Rh;			    /* Relative humidity (%) */
  float Wind;			/* Wind (m/s) */
  float VICSin;         /* Observed incoming shortwave radiation 
                           without topographic or canopy shading (W/m2) */
  float Sin;			/* Incoming shortwave adjusted for topographic shading (W/m^2) */
  float SinBeam;		/* Incoming beam radiation (W/m^2) */
  float SinDiffuse;		/* Incoming diffuse radiation (W/m^2) */
  float Lin;			/* Incoming longwave (W/m^2) */
  float AirDens;		/* Air density on kg/m^3 */
  float Lv;			    /* Latent heat of vaporization (J/kg) */
  float Press;			/* Atmospheric pressure (Pa) */
  float Gamma;			/* Psychrometric constant (Pa/C) */
  float Es;			    /* Saturated vapor pressure (Pa) */
  float Eact;			/* Actual vapor pressure (Pa) */
  float Slope;			/* Slope of vapor pressure curve (Pa/C) */
  float Vpd;			/* Vapor pressure deficit (Pa) */
} PIXMET;

typedef struct {
  float Rank;
  int   x;
  int   y;
} ITEM;

typedef struct {
  double Xorig;					 /* X coordinate of Northwest corner */
  double Yorig;					 /* Y coordinate of Northwest corner */
  int X;						 /* Current x position */
  int Y;						 /* Current y position */
  int NX;						 /* Number of pixels in x direction */
  int NY;						 /* Number of pixels in y direction */
  float DX;						 /* Pixel spacing in x-direction */
  float DY;						 /* Pixel spacing in y-direction */
  float DXY;					 /* Pixel spacing in diagonal */
  int OffsetX;					 /* Offset in x-direction compared to basemap */
  int OffsetY;					 /* Offset in y-direction compared to basemap */
  int NumCells;                  /* Number of cells within the basin */
  int NumLakes;                  /* Number of lakes within the basin */
  int SoilLayers;                /* Number of root zone layers */
  ITEM *OrderedCells;            /* Structure array to hold the ranked elevations; NumCells in size */
} MAPSIZE;

typedef struct {
  char Path[BUFSIZE * 2 + 1];	 /* Map file */
  FILE *Ptr;					 /* Open file, NULL when closed */
  long long HeaderBytes;		 /* Offset of the first value */
  int NumberType;				 /* NC_BYTE or NC_FLOAT */
  int Layers;					 /* Number of layers */
  int ErrorCode;				 /* Exit code of a read error */
} MAPFILE;

typedef struct {
  float Tair;					/* Air temperature (C) */
  float TempLapse;				/* Temperature lapse rate (C/m) */
  float Rh;						/* Relative humidity (%) */
  float Wind;					/* Wind (m/s) */
  float Sin;					/* Incoming shortwave (W/m^2) */
  float SinBeamObs;				/* Observed incoming beam radiation (W/m^2) */
  float SinDiffuseObs;			/* Observed incoming diffuse radiation (W/m^2) */
  float SinBeamMod;				/* Modeled incoming beam radiation (W/m^2) */
  float SinDiffuseMod;			/* Modeled incoming diffuse radiation (W/m^2) */
  float BeamRatio;				/* Ratio of observed beam to modeled beam */
  float DiffuseRatio;			/* Ratio of observed diffuse to modeled diffuse */
  float Lin;					/* Incoming longwave (W/m^2) */
  float ClearIndex;				/* Cloudiness index */
  float Precip;					/* Precipitation if available (m) */
  float Q;					/* Specific humidity at the station (kg/kg) */
} MET;

#define NMETCOLUMNS 7			/* Time and the six values of a station data file */

typedef struct {
  char Name[BUFSIZE + 1];		  /* Station name */
  COORD Loc;					        /* Station locations */
  float Elev;					        /* Station elevations */
  float LapseRate[12];		    /* Temperature lapse rate of each month (C/m) */
  float PrismPrecip[12];		  /* Precipitation normals: precipitation pattern map value at the station for each month */
  float SnowPattern;		      /* Snow pattern average for each station - after re-weighting with precip pattern */
  float SnowPatternBase;		  /* Snow normal: snow pattern map value at the station */
  FILES MetFile;				      /* File with observations */
  int Column[NMETCOLUMNS];	  /* Field of Time and of each value in a record */
  int NFields;				        /* Fields of a record */
  int LineNo;				          /* Line of the record read last */
  int Pending;				        /* TRUE if Fields hold the next record, read ahead */
  char *Record;				        /* Line buffer */
  char **Fields;			        /* Fields of the line in Record */
  MET Data;
} METLOCATION;

typedef struct {
  int NGrids;            
  int Decimal;  
  float LatNorth;               /* extreme north latitude */
  float LonEast;                /* extreme east longitude */
  float LatSouth;               /* extreme south latitude */
  float LonWest;                /* extreme west longitude */
  char filepath[BUFSIZE + 1];   /* file path */
  char fileprefix[BUFSIZE + 1]; /* file path */
} GRID;

/* Vegetation updates: at each instant the eight Veg.* maps of the cells are
   replaced by the maps <Directory>Veg.<Variable>.<stamp>.map (the two class
   maps and the maps whose [INPUTS] keys are given; the others are formed from
   the class tables) */
typedef struct {
  int N;                        /* Updates in the config */
  DATE *Time;                   /* Instant of each update */
  long long *Step;              /* Time steps from Start to each update (negative before Start) */
  char Directory[BUFSIZE + 1];  /* Directory of the update maps, ending in a separator */
  int Next;                     /* First update applied during the run */
} VEGUPDATES;

/* The eight vegetation maps of the cells, as read or formed from the class
   tables (each array NY x NX, row after row): the input maps or the maps of
   one vegetation update */
typedef struct {
  unsigned char *OverstoryClass;	/* Row of the overstory table, 0 = none */
  unsigned char *UnderstoryClass;	/* Row of the understory table, 0 = none */
  float *OverstoryCover;
  float *OverstoryLAI;
  float *OverstoryHeight;
  float *UnderstoryCover;
  float *UnderstoryLAI;			/* LAI per covered area */
  float *UnderstoryHeight;
} VEGMAPS;

/* The soil and vegetation maps a class table can give (ClassMaps[] and
   OPTIONSTRUCT.FromTable[] order): without its [INPUTS] key a map is formed
   from a column of the class table through the class map */
#define CM_DEPTH        0
#define CM_KSLAT        1
#define CM_EXPDEC       2
#define CM_POROSITY     3
#define CM_FCAP         4
#define CM_OVERCOVER    5	/* the three overstory maps in this order */
#define CM_OVERLAI      6
#define CM_OVERHEIGHT   7
#define CM_UNDERCOVER   8	/* the three understory maps in this order */
#define CM_UNDERLAI     9
#define CM_UNDERHEIGHT 10
#define NCLASSMAPS     11

typedef struct {
  const char *Name;		/* Map name and [INPUTS] key */
  const char *Table;		/* [INPUTS] key of the class table */
  const char *Column;		/* Class column (base of the soil layer columns) */
  const char *Units;
  int LayerSet;			/* MAP_ONE or MAP_SOIL */
} CLASSMAP;

typedef struct {
  int FlowGradient;			/* Specifies whether the flow gradient is based
                           on the terrain elevation (TOPOGRAPHY) or the 
                           water table elevation (WATERTABLE).  The 
                           TOPOGRAPHY method is much faster, since the 
                           flow direction and gradient do not have to 
                           be recalculated every timestep */
  int MultiFlowDir;     /* TRUE to use multiple flow directions or FALSE for steepest descent */
  int Routing;          /* Overland flow routing indicator, either CONVENTIONAL (FALSE) or KINEMATIC (TRUE) */
  int LakeDynamics;		  /* If TRUE, lake dynamics will be simulated using power law storage relationships */
  int UseInterflow;		  /* If TRUE, unsaturated lateral flow is simulated */
  int UseKsatAnisotropy;/* Vertical Ksat from lateral Ksat and anisotropy (TRUE) or default table (FALSE) */
  int Infiltration;     /* Specifies static or dynamic maximum infiltration rate */
  int StomatalScheme;   /* Canopy conductance: JARVIS or MEDLYN */
  int Prism;					  /* Precipitation Pattern: if TRUE, the monthly precipitation pattern map and the
                           station normals are used to interpolate precipitation */
  int SnowPattern;		  /* If TRUE, user supplied snow map will be used to redistribute precipitation */
  int MeltScale;        /* TRUE: the melt of every cell times its Snow.MeltScale value */
  int WindScale;        /* TRUE: the wind of every cell times its Met.WindScale value */
  int Shading;					/* if TRUE then terrain shading for solar is on */
  float MaxInterpDistance; /* Maximum station distance of the inverse distance interpolation (m) */
  float MaxInterpDist;  /* Maximum station distance of the inverse distance interpolation (cells) */
  int Impervious;       /* TRUE if some understory class has an impervious fraction */
  int GW_SPINUP;        /* Whether to spinup groundwater state prior to launching run */
  int GW_SPINUP_YRS;    /* Number of years in groundwater spinup */
  float GW_SPINUP_RECHARGE; /* Yearly groundwater recharge rate during spinup (m/yr) */
  int Verbose;          /* TRUE for progress messages */
  int ColdStart;        /* TRUE: Initial State = NONE; FALSE: start from a snapshot */
  int RecomputeDerived; /* TRUE: the states that follow from others are rebuilt from a snapshot */
  float InitChannelDepth;      /* Water depth in every channel segment at a cold start (m) */
  float InitDaysSinceSnowfall; /* Days since snowfall in every cell at a cold start */
  char InitialSWEPath[BUFSIZE + 1];  /* Snow water equivalent map of a cold start, empty if none */
  MAPFILE PrecipPattern;       /* Monthly precipitation pattern map, open during the run */
  MAPFILE HorizonAngle;        /* Horizon angle map, open until the shading is set up */
  char LakeTablePath[BUFSIZE + 1];    /* Lake table, for messages */
  VEGUPDATES VegUpdates;       /* Vegetation maps replaced during the run */
  int FromTable[NCLASSMAPS];   /* TRUE: the map (and each update's) is formed from its class table */
} OPTIONSTRUCT;

/* The input maps as read or formed from the class tables, before the
   parameters apply (each array NY x NX per layer, row after row; NULL when
   the run does not use the map) */
typedef struct {
  float *Elevation;
  unsigned char *Mask;
  unsigned char *SoilClass;
  unsigned char *LakeID;
  float *Depth;
  float *LateralConductivity;
  float *ExponentialDecrease;
  float *Porosity;		/* L + 1 layers */
  float *FieldCapacity;		/* L + 1 layers */
  float *LayerThickness;	/* L layers */
  VEGMAPS Veg;			/* The eight vegetation maps */
  float *SnowPattern;
  float *MeltScale;
  float *WindScale;
  float *InitialSWE;
  float *SkyView;
} INPUTMAPS;

typedef struct {
  float Precip;					    /* Total amount of precipitation at pixel (m) */
  float SumPrecip;          /* Accumulated precipitation at pixel (m) */
  float SnowAccum;          /* Cumulative SWE accumulation for the water year */
  float SnowMelt;           /* Cumulative SWE melt for the water year */
  float RainFall;		        /* Amount of rainfall (m) */
  float SnowFall;		        /* Amount of snowfall determined by air temperature (m) */
  float *IntRain;		        /* Rain interception by each vegetation layer (m) */
  float *IntSnow;		        /* Snow interception by each vegetation layer (m) */
  float TempIntStorage;			/* Temporary snow and rain interception storage, used by MassRelease() */
  int PrecipStart;          /* TRUE if there was surface water in the last time step */
 } PRECIPPIX;

typedef struct {
  float NetShort[2];        /* Shortwave radiation for vegetation surfaces and ground/snow surface W/m2;
                               the lower surface's entry holds only the understory's part */
  float NetShortGround;     /* Net shortwave of the ground (soil or snow) under and outside the
                               understory, per cell area W/m2 */
  float ShortInLower;       /* Shortwave reaching the lower surface (the top of the understory),
                               per cell area W/m2 */
  float LongIn[2];		      /* Incoming longwave radiation for vegetation surfaces and ground/snow surface W/m2 */
  float LongOut[2];		      /* Outgoing longwave radiation for vegetation surfaces and ground/snow surface W/m2;
                               the lower surface's entry is the upward longwave above the understory */
  float LongInGround;       /* Incoming longwave at the ground (soil or snow) W/m2 */
  float LongOutGround;      /* Longwave emitted by the ground (soil or snow) W/m2 */
  float PixelNetShort;	    /* Net shortwave for the entire pixel W/m2 */
  float NetRadiation[2];    /* Net radiation received by the entire pixel W/m2 */
  float PixelLongIn;	      /* Incoming longwave for entire pixel W/m2 */
  float PixelLongOut;	      /* Outgoing longwave for entire pixel W/m2 */
  float ObsShortIn;         /* Incoming shortwave radiation straight from the weather file without topographic or canopy shading */
  float BeamIn;             /* Incoming beam radiation */
  float DiffuseIn;          /* Incomning diffuse radiation */
  float Tair;               /* Air temperature */
  float Wind;               /* Wind speed at the reference height (m/s) */
  float Rh;                 /* Relative humidity (%) */
} PIXRAD;

typedef struct {
  float Area;			    /* Area of channel cut (m) */
  float BankHeight;		/* Height of channel cut (m) */
  int   CutBankZone;	/* Number of the soil layer that contains the bottom of the channel cut */
  float *PercArea;		/* Area of percolation zone for each soil layer, corrected for the channel cut,
                         divided by the grid cell area (0-1)  */
  float *Adjust;		  /* Array with coefficients to correct for loss of soil storage due to channel for each soil layer.
                         Multiplied with RootDepth to give the zone thickness for use in calculating soil moisture */
} NETSTRUCT;

typedef struct {
  float SolarAzimuth;		/* solar azimuth */
  float Latitude;		    /* Latitude of center of study area */
  float Longitude;			/* Longitude of center of study area */
  float StandardMeridian;	/* Standard meridian for current time zone */
  float NoonHour;		    /* Time at which solar noon occurs for current location */

  float Declination;		/* Solar declination */
  float HalfDayLength;		/* Length of half day in hours */
  float Sunrise;		    /* Hour of sunrise */
  float Sunset;				/* Hour of sunset */
  float TimeAdjustment;		/* Time adjustment to be made between center of study area and standard meridian */

  float SunEarthDistance;	/* Distance from Sun to Earth */
  float SineSolarAltitude;	/* Sine of sun's SolarAltitude  */
  int DayLight;				/* FALSE: measured solar radiation and the sun is below the horizon.
                         TRUE: sun is above the horizon */
  float SolarTimeStep;		/* Fraction of the timestep the sun is above the horizon  */

  float SunMax;				/* Calculated solar radiation at the top of the atmosphere (W/m^2) */
} SOLARGEOMETRY;

/* Beam shading computed each time step from the horizon angles: the sun
   positions of the step (shared by all cells) and the horizon of every cell */
typedef struct {
  int NX;                   /* Number of grid columns */
  unsigned short *Horizon;  /* Horizon angles (hundredths of a degree), NAZIMUTHS per
                               cell, cells row after row; set in mask cells only */
  float *SkyView;           /* Sky view factor of each cell (scales the diffuse shortwave),
                               cells row after row */
  int NSun;                 /* Sun positions per time step */
  int NUp;                  /* Positions of the current step with the sun above the horizon */
  double SumSinAlt;         /* Sum of the sine of the sun altitude over those positions */
  double *SinAlt;           /* Sine of the sun altitude */
  double *CosAlt;           /* Cosine of the sun altitude */
  double *AltDeg;           /* Sun altitude (degrees) */
  double *CosAz;            /* Cosine of the sun azimuth */
  double *SinAz;            /* Sine of the sun azimuth */
  int *Az0;                 /* Horizon azimuth just counterclockwise of the sun azimuth */
  int *Az1;                 /* Horizon azimuth just clockwise of it */
  double *AzWeight;         /* Interpolation weight of Az1 */
} SHADE;


typedef struct {
  uchar HasSnow;			/* Snow cover flag determined by SWE */
  uchar SnowCoverOver;		/* Flag overstory can be covered */
  uchar AccumSeason;      /* Flag for whether Tsurf > 0 since last snowfall */
  float LastSnow;			/* Days since last snowfall */
  float Swq;				/* Snow water equivalent */
  float OldSwq;             /* Snow water equivalent from previous time step */
  float Melt;				/* Snow Melt */
  float Outflow;		    /* Snow pack outflow (m) */
  float PackWater;			/* Liquid water content of snow pack */
  float TPack;				/* Temperature of snow pack */
  float SurfWater;			/* Liquid water content of surface layer */
  float TSurf;				/* Temperature of snow pack surface layer */
  float Albedo;				/* Albedo of snow pack */
  float AlbedoGround; /* Albedo of the ground or understory under the snowpack */
  float Depth;				/* Snow depth; Does not appear to be calculated
							   or used anywhere */
  float VaporMassFlux;		/* Vapor mass flux to/from snow pack,(m/timestep). 
                               A negative value indicates flux from snow -- sublimiation */
  float CanopyVaporMassFlux;/* Vapor mass flux to/from intercepted snow in the canopy (m/timestep) */
  float Qsw;                /* Net shortwave radiation exchange at surface */
  float Qlw;                /* Net longwave radiation exchange at surface */
  float Qs;				    /* Sensible heat exchange */
  float Qe;				    /* Latent heat exchange */
  float Qp;                 /* advected heat from rain input */
  float MeltEnergy;			/* Energy used to melt snow and change of cold content of snow pack */
  
  // spatial parameters
  float Ts;	   				    /* snow temperature threshold */
  float Tr;	    				/* rain temperature threshold */
  float amax;					/* fresh snow albedo */
  float LamdaAcc;               /* lambda for accumulation period */
  float LamdaMelt;              /* lambda for melt period */
  float AccMin;                 /* minimum albedo for accumulation period*/
  float MeltMin;                /* minimum albedo for melt period*/
  float MeltScale;              /* multiplier on the melt of the pack (Snow.MeltScale map, else 1) */

  // SWE stats 
  float MaxSwe;         /* Peak SWE of the water year*/
  unint MaxSweDate;       /* Peak SWE date/timestep of the water year*/
  unint MeltOutDate;    /* Last day of SWE of the water year */
} SNOWPIX;

typedef struct {
  int   Soil;			/* Soil type */
  float Depth;			/* Depth of total soil zone, including all root
						zone layers, and the saturated zone */
  float *Moist;			/* Soil moisture content in layers (0-1) */
  float *Perc;			/* Percolation from layers */
  float *Temp;			/* Temperature in each layer (C) */
  float TableDepth;		/* Depth of water table below ground surface (m) */
  float WaterLevel;		/* Absolute height of the watertable above datum (m), 
						          i.e. corrected for terrain elevation */
  float WaterLevelLast;
  float SatFlow;		/* Amount of saturated flow generated */
  float *InterFlow;		/* Amount of unsaturated lateral interflow entering each layer */
  float MaxDeepFlux; /* Theoretical deep flux before accounting for available water */
  float DeepFlux; /* Amount of water entering (+) or leaving (-) bottom of soil column (m) */
  float IExcess;		/* Amount of surface runoff (m) generated from HOF and Return flow */
  float Runoff;         /* Surface water flux (m) from the grid cell. */
  float ChannelInt;		/* Amount of subsurface flow intercepted by the channel */
  float ChannelInfiltration;		/* Amount of channel storage returned to the subsurface */
  float TSurf;			/* Soil surface temperature */
  float Qnet;			/* Net radiation exchange at surface */
  float Qrest;			/* Rest term for energy balance (should be 0) */
  float Qs;				/* Sensible heat exchange */
  float Qe;				/* Latent heat exchange */
  float Qg;				/* Ground heat exchange */
  float Qst;			/* Ground heat storage */
  float Ra;				/* Soil surface aerodynamic resistance (s/m) */
  float InfiltAcc;               /* Accumulated water in the top layer (m) */
  float MoistInit;               /* Initial moisture content when ponding begins (0-1) */
  float DesorptionAcc;           /* Soil evaporation per unit ground area since the top layer
                                    was last wetted (m) */
  float startRunoff;             /* Surface water flux from the previous (sub) time step. Used for kinematic wave routing */
  float DetentionStorage;        /* amount of water kept in detention storage when impervious fraction > 0 */
  float DetentionIn;			 /* detention storage change in current time step */
  float DetentionOut;            /* water flow out of detention storage */
  float MaxInfiltrationRate; /* Maximum surface infiltration rate (m/s) */
  float KsLat;      		 /* Soil Lateral Conductivity */
  float KsLatExp;        /* Exponent for vertical change of KsLat */
  float *KsVert;         /* Soil Vertical Conductivity - each layer */
  float *Porosity;       /* Soil Porosity - each layer */
  float *FCap;			     /* Soil Field Capacity - each layer */
  float *PoreDist;       /* Brooks-Corey pore size distribution index - each layer */
  float *Press;          /* Brooks-Corey bubbling pressure (m) - each layer */
  float *LayerThickness; /* Thickness of each root zone layer (m) */
  float TotalThickness;  /* Sum of LayerThickness from 0.0 in layer order (m) */
} SOILPIX;

typedef struct {
  char Desc[BUFSIZE + 1];	/* Soil type */
  int Index;
  int NLayers;				/* Number of root zone layers (the deep layer is index NLayers) */
  float Albedo;				/* Albedo of the soil surface */
  float Manning;      /* Manning's roughness of the soil surface */
  float *WP;				/* Wilting point for each layer and the deep layer */
  float *Ks;				/* Saturated hydraulic conductivity (vertical) for each layer and the deep layer */
  float *Porosity;			/* Porosity of each layer and the deep layer (Soil.Porosity from the table) */
  float *FCap;				/* Field capacity of each layer and the deep layer (Soil.FieldCapacity from the table) */
  float Depth;				/* Soil depth (m), Soil.Depth from the table */
  float KsLat;				/* Lateral conductivity (m/s), Soil.LateralConductivity from the table */
  float KsLatExp;			/* Exponential decrease (1/m), Soil.ExponentialDecrease from the table */
  float DepthThresh;    /* Threshold water table depth, beyond which transmissivity decays linearly with water table depth */
  float KsAnisotropy; /* KsVert = KsLat / KsAnisotropy */
  float MaxInfiltrationRate; /* Maximum infiltration rate for upper layer (m/s) */
  float DeepFlux; /* Net groundwater flux upwelling (+) or downwelling (-) (m/yr as read; m per time step after ApplyTableParameters) */
  float G_Infilt;                /* Mean capillary drive for dynamic maximum infiltration rate (m)   */
} SOILTABLE;

struct LAKETABLE {
  char Name[BUFSIZE + 1];	/* Lake name */
  int Index;
  int Outlet;             /* Segment into which the lake pours, as in the lake table */
  SegmentID OutletID;     /* ID of channel into which the lake pours */
  Channel *outlet;		    /* Pointer to outlet segment record */
  float PowLawScale;      /* Scale of the power law of the lake outflow */
  float PowLawExponent;   /* Exponent of the power law of the lake outflow */
  float Area;             /* Lake area (m^2) */
  float Storage;          /* Lake storage area-normalized (m) */
  float Inflow;           /* Lake inflow area-normalized (m/timestep) */
  float Outflow;          /* Lake outflow area-normalized (m/timestep) */
};

typedef struct {
  char Distribution[BUFSIZE+1];	/* Distribution type */
  float mean;
  float stdev;
  float min;
  float max;
  float mode;
} STATSTABLE;


typedef struct {
  float Dem;					/* Elevations */
  uchar Mask;					/* Mask for modeled area */
  unshort Travel;				/* Travel time */
  float Grad;					/* Sum of downslope slope-width products */
  float Slope;					/* Land surface slope */
  float Aspect;					/* Land surface slope direction */
  float FlowGrad;				/* Magnitude of subsurface flow gradient slope * width */
  unsigned char Dir[MAXDIRS];		/* Fraction of surface flux moving in each direction*/
  unsigned int TotalDir;	    /* Sum of Dir array */
  int drains_x;					/* x-loc of cell to which this impervious cell drains */
  int drains_y;					/* y-loc of cell to which this impervious cell drains */
  int LakeID;           /* Unique ID of each discrete lake */
  ITEM *OrderedTopoIndex;       /* Structure array to hold the ranked topoindex for fine pixels in a coarse pixel */
  uchar LateralDir;    /* Direction of downhill grid cell for unsaturated flow */
  float CosSlope;    /* Cosine of slope for unsaturated flow partitioned downslope */
  float SinSlope;    /* Sine of slope for unsaturated flow partitioned downslope */
  float WindScale;   /* Multiplier on the station wind (Met.WindScale map, else 1) */
} TOPOPIX;

/* The vegetation of a cell. The canopy arrays hold the layers present in
   position order, overstory first (CanopyIndex()) */
typedef struct {
  int OverClass;			    /* Row of the overstory table, 0 = none */
  int UnderClass;			    /* Row of the understory table, 0 = none */
  unsigned char OverStory;	    /* TRUE if the cell has an overstory layer (OverClass > 0) */
  unsigned char UnderStory;	    /* TRUE if the cell has an understory layer (UnderClass > 0
                                   and understory cover > 0) */
  int NVegLayers;		        /* Vegetation layers of the cell (0-2) */
  float UnderstoryCover;	    /* Understory cover: its map times its scale, at most 1 (0-1),
                                   0 without an understory layer */
  float Tcanopy;		        /* Canopy temperature (C) */
  float TAcclim;		        /* Cold-acclimation temperature of the photosynthesis (C) */
  float TGrowth;		        /* Growth temperature: the 10-day mean air temperature (C) */
  float MoistureFlux;		    /* Amount of water transported from the pixel
                                   to the atmosphere (m/timestep) */
  float MeltEnergy;			    /* Energy used to melt snow and change of cold content
                                of snow pack */
  float *Fract;                  /* Fractional Coverage*/
  float *LAI;                  /* LAI of the day (SetCellLAI) */
  float *Height;                /* Height of each vegetation layer (m), from the height maps */
  float *LAIFull;              /* Full-leaf LAI of each layer: its LAI map times its scale */
  float *MaxInt;                  /* Max Interception*/   
  float Taud;                     /* Transmittance of the overstory for diffuse radiation in the month */
  float TauU;                     /* Shortwave transmittance of the understory in the month
                                     (1 without an understory layer) */
  float U[2];                     /* Wind profile of the cell for a reference wind of 1 m/s (CalcAerodynamic) */
  float USnow;                    /* Wind 2 m above the snow surface, same normalization */
  float Ra[2];                    /* Aerodynamic resistances in the absence of snow, same normalization */
  float RaSnow;                   /* Aerodynamic resistance above the snow surface, same normalization */
} VEGPIX;

/* The parameters of one story (overstory or understory) of a vegetation
   table row */
typedef struct {
  float RsMax;			/* Maximum stomatal resistance (s/m), JARVIS */
  float RsMin;			/* Minimum stomatal resistance (s/m), JARVIS */
  float MoistThres;		/* Soil moisture threshold above which soil
				   moisture does not restrict transpiration,
				   relative: 0 at the wilting point, 1 at
				   field capacity */
  float VpdThres;		/* Vapor pressure deficit threshold above which
				   stomatal closure occurs (Pa), JARVIS */
  float Rpc;			/* Visible radiation of half the maximum light factor (W/m2), JARVIS */
  float Vcmax25;		/* Maximum carboxylation rate at 25 C at the top of the
				   canopy (umol/m2/s), MEDLYN */
  float G1;			/* Stomatal slope (kPa^0.5), MEDLYN */
  float G0;			/* Stomatal intercept (mol/m2/s), MEDLYN */
  float *RootFract;		/* Fraction of roots in each root zone layer */
  float AlbedoMonthly[12];	/* Albedo of each month */
  float Albedo;			/* Albedo of the current month */
  int Deciduous;		/* 1: the LAI follows the leaf-on and leaf-off logistics, 0: full LAI */
  float LeafOnMid;		/* Day of year of the leaf-on logistic's midpoint (January 1 = 1), deciduous */
  float LeafOnSlope;		/* Slope of the leaf-on logistic (logit/day), deciduous */
  float LeafOffMid;		/* Day of year of the leaf-off logistic's midpoint, deciduous */
  float LeafOffSlope;		/* Slope of the leaf-off logistic (logit/day, negative), deciduous */
  float MinLAIFraction;		/* LAI without leaves as a fraction of the full LAI, deciduous */
  float Cover;			/* Cover (fraction): the story's cover map from the table */
  float LAI;			/* LAI (m2/m2): the story's LAI map from the table */
  float Height;			/* Height (m): the story's height map from the table */
} STORYTABLE;

/* One row of the overstory table */
typedef struct {
  char Desc[BUFSIZE + 1];	/* Name */
  STORYTABLE Story;		/* Stomatal, root, albedo and LAI parameters */
  float ExtnCoeff;		/* Shortwave extinction coefficient per unit LAI (LightExtinction) */
  float Cn;			/* Canopy attenuation coefficient for wind profile */
  float MaxSnowInt;		/* Maximum snow interception capacity */
  float MDRatio;		/* Ratio of Mass Release to Meltwater drip from int snow */
  float SnowIntEff;		/* Efficiency of snow interception process */
  float Trunk;			/* Fraction of overstory height that identifies the top of the trunk space */
} OVERTABLE;

/* One row of the understory table */
typedef struct {
  char Desc[BUFSIZE + 1];	/* Name */
  STORYTABLE Story;		/* Stomatal, root, albedo and LAI parameters */
  float RadAtten;		/* Shortwave attenuation coefficient of the understory per unit LAI */
  float LitterCover;		/* Fraction of the cell under litter, where the soil does not evaporate */
  float ImpervFrac;		/* fraction of pixel that is impervious */
  float DetentionFrac;		/* fraction of flow on impervious area goes to flood detention storage */
  float DetentionDecay;		/* Decay coefficient of linear reservoir storage for impervious surface detention facility */
} UNDERTABLE;

/* The overstory and understory tables (class = row number) */
typedef struct {
  int NOver;			/* Rows of the overstory table */
  OVERTABLE *Over;
  int NUnder;			/* Rows of the understory table */
  UNDERTABLE *Under;
} VEGTABLES;

/* Means over mask cells (the basin, or the cells of one class combination) */
typedef struct {
  int NCells;			/* Mask cells averaged over */
  int NSnow;			/* Of these, cells with snow at the end of the step */
  float *Sum;			/* Running totals, by slot (output.h) */
  float *Mean;			/* Means, by slot */
  float OldSwq;			/* Mean snow water equivalent at the start of the step */
} AGGREGATE;

/* Basin and class means of a time step, the basin water balance and the
   sums of the run */
typedef struct {
  AGGREGATE Basin;
  int NClasses;			/* Class combinations with means */
  AGGREGATE *Class;
  int *ClassOver;		/* Overstory class of each combination */
  int *ClassUnder;		/* Understory class of each combination */
  int NSaturated;		/* Mask cells with the water table above 90% of the soil depth */
  float SaturatedFraction;
  float Storage;		/* Water in the basin at the end of the step (m) */
  float OldStorage;		/* ... at the start of the step */
  float Error;			/* Water balance error of the step (m) */
  float CumulativeError;	/* Sum of the errors of the run (m) */
  float RainOrMelt;		/* Rain, or snowpack outflow when there is snow (m) */
  double Network[7];		/* Channel network sums over the basin area (output.h) */
  float StartStorage;		/* Water in the basin at Start (m) */
  float CumPrecip;		/* Sums of the run (m) */
  float CumSnowVapor;
  float CumDeepFlux;
  float CumET;
  float CumChannelInt;
  float CumChannelInfiltration;
  float CumChannelEvap;
  double CumOutlet;		/* Outflow of the basin outlet (m) */
  double CumNetworkError;	/* Sum of the channel network errors of the run (m) */
} TOTALS;

typedef struct node node;
struct node {
  node *next;
  int x;
  int y;
};

/* The model state the outputs read */
typedef struct {
  MAPSIZE *Map;
  TOPOPIX **Topo;
  SOILPIX **Soil;
  SOILTABLE *SType;
  VEGPIX **Veg;
  SNOWPIX **Snow;
  PRECIPPIX **Precip;
  EVAPPIX **Evap;
  PIXRAD **Rad;
  NETSTRUCT **Network;
  Channel *Streams;		/* Channel segments in segment table order */
  int Dt;			/* Time step (s) */
  int NumLakes;			/* Lakes of the lake table (0 without Lake Dynamics) */
  float *LakeLevel;		/* Level of each lake: the mean surface water of its cells (m), set by Aggregate */
  int *LakeCells;		/* Mask cells of each lake, counted by Aggregate */
} MODELVIEW;

#endif
