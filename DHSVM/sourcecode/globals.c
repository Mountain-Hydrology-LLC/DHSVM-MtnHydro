
/* These indices are so neighbors can be looked up quickly */
int xdirection[] = {-1,  0,  1,  1,  1,  0, -1, -1};
int ydirection[] = {1,  1,  1,  0, -1, -1, -1,  0};

float LAI_SNOW_MULTIPLIER;		/* Multiplier to calculate the amount of available snow interception as a function of LAI */
float LAI_WATER_MULTIPLIER;		/* Multiplier to determine maximum interception storage as a function of LAI  */
float LIQUID_WATER_CAPACITY;		/* Water holding capacity of snow as a fraction of snow-water-equivalent */
float MAX_SNOW_TEMP;				/* Maximum temperature at which snow can occur (C) */
float MIN_INTERCEPTION_STORAGE;	/* The amount of snow on the canopy that can only be melted off (m) */
float MIN_SNOW_RESET_ALBEDO; /* Snowfall rate (m/d) above which snowfall resets the snow albedo (a smaller rate r scales the days since snowfall by 1 - r/R per day); also the SWE depth (m) below which the ground albedo is blended in */
float LAMBDA_FOREST_OFFSET; /* Additional correction applied to albedo decay rate in forest */
float MIN_RAIN_TEMP;				/* Minimum temperature at which rain can occur (C) */
unsigned char OUTSIDEBASIN = 0;		/* Mask value indicating outside the basin (0) */
float SNOWPAT_WEIGHT; /* Fractional weight between snow pattern and precip pattern for snowfall */
float TEMPERATURE_OFFSET; /* Uniform offset added to input air temperature (C) */
float LAPSE_RATE_BIAS; /* Additional lapse rate relative to LAPSE_BIAS_ELEV (C/m) */
float LAPSE_BIAS_ELEV; /* Elevation where temperature bias is zero (m) */
float SOIL_DEPTH_ADJ; /* Soil Depth Offset (m), added to the soil depth map */
float SOIL_KSAT_ADJ; /* Lateral Conductivity Scale */
float SOIL_EXPDEC_ADJ; /* Exponential Decrease Scale */
float SOIL_POROSITY_ADJ; /* Porosity Scale */
float SOIL_FIELDCAP_ADJ; /* Field Capacity Scale */
float VEG_LAI_ADJ; /* Overstory LAI Scale */
float VEG_FC_ADJ; /* Overstory Cover Scale */
float VEG_UNDER_LAI_ADJ; /* Understory LAI Scale */
float VEG_UNDER_FC_ADJ; /* Understory Cover Scale */
float SOIL_DEEPFLUX_ADJ; /* Deep Flux Offset (m/yr), added to the deep flux of every soil class */
float VEG_RSMIN_ADJ; /* Min Stomatal Resistance Scale (JARVIS) */
float VEG_MOIST_ADJ; /* Moisture Threshold Scale */
float VEG_G1_ADJ; /* Stomatal Slope Scale (MEDLYN) */
float ATMOS_CO2; /* Atmospheric CO2 of the step (umol/mol), MEDLYN */
float VEG_SNOWINT_ADJ; /* Snow Interception Efficiency Scale */
float LAKE_OUTFLOW_ADJ; /* Lake Outflow Rate Scale */
float LAKE_EXPONENT_ADJ; /* Lake Outflow Exponent Scale */
float Z0_GROUND;					/* Roughness length for bare soil (m) */
float Z0_SNOW;					/* Roughness length for snow (m) */
float Zref;						/* Reference height (m) */

/* Snow albedo decay curve */
float ALB_MAX;                   /* Fresh snow albedo */                                                               
float ALB_ACC_LAMBDA;            /* Snow freeze albedo curve control parameters */
float ALB_MELT_LAMBDA;           /* Snow thaw albedo curve control parameters */
float ALB_ACC_MIN;
float ALB_MELT_MIN;
float PRECIP_MULTIPLIER;        /* Precipitation multiplier */
float MAX_SURFACE_SWE; 	       /* Maximum depth of the surface layer in snow water equivalent (m) */
float SNOW_DENSITY;             /* Density of all snow (kg/m3) */
