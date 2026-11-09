
#ifndef CONSTANTS_H
#define CONSTANTS_H

#define CH_ICE     (2100.0e3)	/* Volumetric heat capacity (J/(m3*C) of ice (0C) */
#define CH_WATER   (4186.8e3)	/* Volumetric heat capacity (J/(m3*C) of water */
#define CP         1013.0		/* Specific heat of moist air at constant pressure (J/(kg*C)) */
#define DELTAT     50.		    /* Used in SnowMelt to bracket the effective surface temperature (C) */
#define DEGPRAD     57.29578	/* degree per radian */
#define D0_MULTIPLIER  0.63		/* Multiplier for vegetation height to get displacement height (m) */
#define EPS         0.622		/* ratio of molecular weight of water vapor to that for dry air */
#define G           9.81		/* gravitational accelleration (m/(s^2)) */
#define GRAMSPKG    1000.	    /* grams per kilogram */
#define JOULESPCAL  4.1868	    /* Joules per calorie */
#define KhH2O       0.58		/* Thermal conductivity of water (W/(mk)) */
#define LF            (333.7e3) /* latent heat of fusion (J/kg) */
#define MINPDEG        4.		/* minutes per degree longitude */
#undef PI
#define PI   3.14159265358979323846
#define RADPHOUR    0.2617994   /* radians per hour: Earth's Rotation (2 PI rad/day) * (1 day/24 h) */
#define RADPDEG     (PI/180.0)	/* radians per degree */
#define SOIL_NODES 10			/* Soil temperature nodes, the surface to SOIL_DAMPING_DEPTH */
#define SOIL_DAMPING_DEPTH 4.		/* Depth of the bottom soil temperature node (m) */
#define SOIL_FROST_ZONES ((SOIL_NODES + 1) / 2)	/* Most separate frozen zones of the nodes */
#define SOIL_PARTICLE_DENSITY 2685.	/* Density of the mineral soil solids (kg/m3) */
#define SOIL_TKTRIP 273.16		/* Triple point of water (K), freezing point depression */
#define SOIL_K_ICE 2.2			/* Thermal conductivity of ice in soil (W/(m C)) */
#define SOIL_K_WATER 0.57		/* Thermal conductivity of water in soil (W/(m C)) */
#define SOIL_KDRY_ORG 0.05		/* Dry thermal conductivity of organic solids (W/(m C)) */
#define SOIL_KS_ORG 0.25		/* Thermal conductivity of organic solids (W/(m C)) */
#define SOIL_K_LITTER 0.1		/* Thermal conductivity of surface litter (W/(m C)) */
#define SOIL_CS_MIN 2.0e6		/* Volumetric heat capacity of mineral solids (J/(m3 C)) */
#define SOIL_CS_ORG 2.7e6		/* Volumetric heat capacity of organic solids (J/(m3 C)) */
#define SOIL_CS_WATER 4.2e6		/* Volumetric heat capacity of water in soil (J/(m3 C)) */
#define SOIL_CS_ICE 1.9e6		/* Volumetric heat capacity of ice in soil (J/(m3 C)) */
#define SOIL_CS_AIR 1.3e3		/* Volumetric heat capacity of air in soil (J/(m3 C)) */
#define SOLARCON   1360.	    /* Solar constant (W/m^2) */
#define STEFAN    (5.6696e-8)	/* Stefan-Boltzmann constant (W/(M^2*C^4) */
#define VISFRACT    0.5		    /* part of shortwave that is in the visible range of the spectrum */
#define VON_KARMAN  0.4		    /* Von Karman's constant */
#define WATER_DENSITY 1000.		/* Density of water in kg/m3 */
#define Z0_MULTIPLIER 0.13		/* Multiplier for vegetation height to get roughness length (m) */
#define MTHRESH        0.9   /* Saturation extent is based on the number of pixels with a water table 
                                that is at least MTHRESH of soil depth. */
#define NDIRS          8     /* Neighbors used in surface/subsurface routing */

/**************** extern constants - see globals.c ****************/

extern int xdirection[];        /* Column offsets of the NDIRS neighbors */
extern int ydirection[];        /* Row offsets of the NDIRS neighbors */

extern float LAI_SNOW_MULTIPLIER;		/* Multiplier to calculate the amount of available snow interception as a function of LAI */
extern float LAI_WATER_MULTIPLIER;		/* Multiplier to determine maximum interception storage as a function of LAI  */
extern float LIQUID_WATER_CAPACITY;		/* Water holding capacity of snow as a fraction of snow-water-equivalent */
extern float MAX_SNOW_TEMP;				/* Maximum temperature at which snow can occur (C) */
extern float MIN_INTERCEPTION_STORAGE;	/* The amount of snow on the canopy that can only be melted off (m) */
extern float MIN_SNOW_RESET_ALBEDO; /* Snowfall rate (m/d) above which snowfall resets the snow albedo (a smaller rate r scales the days since snowfall by 1 - r/R per day); also the SWE depth (m) below which the ground albedo is blended in */
extern float LAMBDA_FOREST_OFFSET; /* Additional correction applied to albedo decay rate in forest */
extern float MIN_RAIN_TEMP;				/* Minimum temperature at which rain can occur (C) */
extern unsigned char OUTSIDEBASIN;		/* Mask value indicating outside the basin (0) */
extern float SNOWPAT_WEIGHT; /* Fractional weight between snow pattern and precip pattern for snowfall */
extern float TEMPERATURE_OFFSET; /* Uniform offset added to input air temperature (C) */
extern float LAPSE_RATE_BIAS; /* Additional lapse rate relative to LAPSE_BIAS_ELEV (C/m) */
extern float LAPSE_BIAS_ELEV; /* Elevation where temperature bias is zero (m) */
extern float SOIL_DEPTH_ADJ; /* Soil Depth Offset (m), added to the soil depth map */
extern float SOIL_KSAT_ADJ; /* Lateral Conductivity Scale */
extern float SOIL_EXPDEC_ADJ; /* Exponential Decrease Scale */
extern float SOIL_POROSITY_ADJ; /* Porosity Scale */
extern float SOIL_FIELDCAP_ADJ; /* Field Capacity Scale */
extern float VEG_LAI_ADJ; /* Overstory LAI Scale */
extern float VEG_FC_ADJ; /* Overstory Cover Scale */
extern float VEG_UNDER_LAI_ADJ; /* Understory LAI Scale */
extern float VEG_UNDER_FC_ADJ; /* Understory Cover Scale */
extern float SOIL_DEEPFLUX_ADJ; /* Deep Flux Offset (m/yr), added to the deep flux of every soil class */
extern float VEG_RSMIN_ADJ; /* Min Stomatal Resistance Scale (JARVIS) */
extern float VEG_MOIST_ADJ; /* Moisture Threshold Scale */
extern float VEG_G1_ADJ; /* Stomatal Slope Scale (MEDLYN) */
extern float ATMOS_CO2; /* Atmospheric CO2 of the step (umol/mol), MEDLYN */
extern float VEG_SNOWINT_ADJ; /* Snow Interception Efficiency Scale */
extern float LAKE_OUTFLOW_ADJ; /* Lake Outflow Rate Scale */
extern float LAKE_EXPONENT_ADJ; /* Lake Outflow Exponent Scale */
extern float Z0_GROUND;					/* Roughness length for bare soil (m) */
extern float Z0_SNOW;					/* Roughness length for snow (m) */
extern float Zref;						/* Reference height (m) */

/* Snow albedo decay curve */
extern float ALB_MAX;                   /* Fresh snow albedo */                                                               
extern float ALB_ACC_LAMBDA;            /* Snow freeze albedo curve control parameters */
extern float ALB_MELT_LAMBDA;           /* Snow thaw albedo curve control parameters */
extern float ALB_ACC_MIN;
extern float ALB_MELT_MIN;
extern float PRECIP_MULTIPLIER;        /* Precipitation multiplier */
extern float MAX_SURFACE_SWE; 	       /* Maximum depth of the surface layer in snow water equivalent (m) */
extern float SNOW_DENSITY;             /* Density of all snow (kg/m3) */
#endif
