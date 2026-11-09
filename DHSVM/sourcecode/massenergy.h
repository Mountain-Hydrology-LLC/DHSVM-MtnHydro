
#ifndef MASSENERGY_H
#define MASSENERGY_H

#include "data.h"
#include <stdarg.h>
#include "DHSVMChannel.h"
#include "photosynthesis.h"

/* The stomatal scheme of a vegetation layer's transpiration: what MEDLYN
   needs beside the meteorology and the story's parameters, and what both
   schemes report */
typedef struct {
  int Scheme;			/* JARVIS or MEDLYN */
  PHOTOKIN Kin;			/* Kinetics of the cell and step (MEDLYN) */
  float Ca;			/* Atmospheric CO2 (umol/mol) (MEDLYN) */
  float Acclim;			/* Cold-acclimation factor (MEDLYN) */
  int TwoLeaf;			/* TRUE: sunlit and shaded leaves (MEDLYN overstory) */
  float SinAlt;			/* Sine of the solar altitude (TwoLeaf) */
  float VisBeam;		/* Incident visible beam above the canopy (W/m2) (TwoLeaf) */
  float VisDiffuse;		/* Incident visible diffuse above the canopy (W/m2) (TwoLeaf) */
  float Rc;			/* Out: canopy resistance (s/m) */
  float An;			/* Out: net assimilation per covered area (umol/m2/s) (MEDLYN) */
  float Beta;			/* Out: soil moisture factor */
} CANOPYSCHEME;

float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
		       float VpdThres, float MoistThres, float WP, float FCap,
		       float TSoil, float SoilMoisture, float Vpd, float Rp);

void CanopyResistanceLayers(int Layer, PIXMET *Met, float Rp, STORYTABLE *Story,
                            SOILTABLE *SType, float *FCap, float *Moist, float *SoilTemp,
                            VEGPIX *LocalVeg, CANOPYSCHEME *Canopy, float *Rc,
                            float *Share);

float Desorption(int Dt, float Moisture, float Porosity, float Ks, 
			   float Press, float m);

void EvapoTranspiration(int Layer, int Dt, PIXMET *Met,
              float NetRad, float Rp, STORYTABLE *Story, SOILTABLE *SType,
              float MoistureFlux, float *Moist, float *Temp, float *Int,
              float *EPot, float *EInt, float **ESoil, float *EAct,
              float *ETot, float *Adjust, float *LayerThickness, float *FCap, float Ra,
              VEGPIX *LocalVeg, CANOPYSCHEME *Canopy, float Overflow, float *Throughfall);

void InterceptionStorage(int NAct, float *MaxInt, float *Fract, float *Int,
               float *Precip, float *Overflow);

void LongwaveBalance(unsigned char OverStory,
			   float F, float Ld, float Tcanopy, float Tsurf,
               float Cu, float Tair, PIXRAD *LocalRad);

float OverstoryViewFactor(VEGPIX *LocalVeg);

void NoSensibleHeatFlux(int Dt, PIXMET *LocalMet, float ETot, SOILPIX *LocalSoil);

float GroundEnergyBalance(int y, int x, int Dt, PIXMET *LocalMet, SOILTABLE *SType,
                          SOILPIX *Soil, int Snow, float SnowTemp, float SnowCond,
                          float NetRadGround, float LongOutGround, float LowerRa,
                          float LowerWind, float GroundEvap, float LitterFrac, float LitterThick);

void RadiationBalance(OPTIONSTRUCT *Options,
              int Overstory, float FuR, float Cu,
              float SineSolarAltitude, float VICRs,
              float Rs, float Rsb, float Rsd, float Ld, float Tair, float TGround, float Tcanopy,
              float SoilAlbedo, OVERTABLE *Over, UNDERTABLE *Under, SNOWPIX *LocalSnow,
              PIXRAD *LocalRad, VEGPIX *LocalVeg);

void ShortwaveBalance(unsigned char OverStory,
              float F, float Rs, float Rsb, float Rsd, float Tau,
              float Taud, float *Albedo, float FuR, float AlbedoU, float AlbedoG,
              float TauU, PIXRAD * LocalRad);

float PondEvaporation(int Dt, float DXDY, float ChannelArea,
                      float Temp, float Slope, float Gamma,
                      float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                      float Evapotranspiration, float *IExcess);

void ChannelEvaporation(int Dt, float DXDY,
                        float Temp, float Slope, float Gamma,
                        float Lv, float AirDens, float Vpd, float NetRad, float LowerRa,
                        float Evapotranspiration, int x, int y, CHANNEL *ChannelData);

float ChannelSoilEvaporation(int Dt, float DXDY,
                             float Temp, float Slope, float Gamma, float Lv,
                             float AirDens, float Vpd, float NetRad, float RaSoil,
                             float Evapotranspiration, float *Porosity, float *FCap, float *Ks,
                             float *Press, float *m, float LayerThickness,
                             float *MoistContent, float *Adjust,
                             int x, int y, CHANNEL *ChannelData, int CutBankZone);

float SoilEvaporation(int Dt,
                      float Temp, float Slope, float Gamma, float Lv,
                      float AirDens, float Vpd, float NetRad, float RaSoil,
                      float Evapotranspiration, float Porosity, float FCap, float Ks,
                      float Press, float m, float WP, float RootDepth,
                      float *MoistContent, float Adjust, float BareFract, float *PotEvap,
                      float *DesorptionAcc);

float StabilityCorrection(float Z, float d, float Tsurf, float Tair,
			  float Wind, float Z0);

#endif
