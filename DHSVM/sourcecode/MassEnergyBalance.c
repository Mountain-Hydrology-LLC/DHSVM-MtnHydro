
#ifdef SNOW_ONLY
  #define NO_SOIL
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "snow.h"
#include "constants.h"
#include "soilmoisture.h"
#include "Calendar.h"
#include "photosynthesis.h"

 /*****************************************************************************
   Function name: MassEnergyBalance()

   Reference    :
     Epema, G.F. and H.T. Riezbos, 1983, Fall Velocity of waterdrops at different
 heights as a factor influencing erosivity of simulated rain. Rainfall simulation,
 Runoff and Soil Erosion. Catena suppl. 4, Braunschweig. Jan de Ploey (Ed), 1-17.

   Laws, J.o., and D.A. Parsons, 1943, the relation of raindrop size to intensity.
 Trans. Am. Geophys. Union, 24: 452-460.

 *****************************************************************************/
void MassEnergyBalance(OPTIONSTRUCT *Options, int y, int x,
  float SineSolarAltitude, float DX, float DY,
  int Dt, int InfiltOption, int MaxSoilLayers, int MaxVegLayers, PIXMET *LocalMet,
  NETSTRUCT *LocalNetwork, PRECIPPIX *LocalPrecip,
  OVERTABLE *Over, UNDERTABLE *Under, VEGPIX *LocalVeg, SOILTABLE *SType,
  SOILPIX *LocalSoil, SNOWPIX *LocalSnow, PIXRAD *LocalRad,
  EVAPPIX *LocalEvap, CHANNEL *ChannelData,
  VEGPIX **VegMap, NETSTRUCT **Network, SOILPIX **SoilMap,
  TOPOPIX *LocalTopo, TOPOPIX **TopoMap, MAPSIZE *Map)
{
#ifndef NO_SOIL
  float SurfaceWater;		/* Pixel average depth of water before infiltration is calculated (m) */
  float ChannelWater;       /* Precip that hits the channel */
  float Infiltration;		/* Infiltration into the top soil layer (m) */
  float Infiltrability;     /* Dynamic infiltration capacity (m/s)*/
  float B;                  /* Capillary drive and soil saturation deficit used
                               in dynamic infiltration calculation*/
  float MaxInfiltration;	/* Maximum infiltration into the top soil layer (m) */
  float PercArea;           /* Surface area of percolation corrected for
                               channel area, divided by the grid cell area (0-1)  */
#endif
  float LowerRa;		    /* Aerodynamic resistance for lower layer (s/m) */
  float NetRadiation;		/* Total Net long- and shortwave radiation for each veg layer (W/m2) */
  float Rp;					/* radiation flux in visible part of the spectrum (W/m^2) */
  float UpperRa;		    /* Aerodynamic resistance for upper layer (s/m) */
  float UpperWind;			/* Wind for upper layer (m/s) */
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float Tsurf;				/* Surface temperature used in LongwaveBalance (C) */
  int   NVegLActual;		/* Number of vegetation layers above snow */
  float RainOverflow[MAXVEGLAYERS] = {0.0, 0.0}; /* Rain each layer's full interception storage passed on (m) */
  float DryChannelEvap;
  int   NVeg;               /* Vegetation layers of the cell (the soil's entry of EPot and EAct) */

  /* Understory cover and burial by the snow */
  int   Lower;              /* Radiation index of the lower surface (understory, soil or snow) */
  float Fu;                 /* Understory cover (0 without an understory layer) */
  float FuR;                /* Fraction of the lower surface covered by the understory for the
                               shortwave of the step (0 with the understory buried) */
  float Cu;                 /* Fraction of the lower surface opaque to longwave by the
                               understory (0 with the understory buried) */
  float BurySwe;            /* Snow water equivalent above which the understory is buried (m) */
  int   Exposed;            /* TRUE if the understory is above the snow at the start of the step,
                               with the step's snowfall */
  float NetShortL;          /* Net shortwave of the lower surface (W/m2) */
  float NetRadGround;       /* Net radiation of the ground (soil or snow) (W/m2) */
  float LongOutU;           /* Longwave emitted by the understory (W/m2) */
  float SoilEvap;           /* Evaporation from the soil (m) */
  float SoilEPot;           /* Potential evaporation of the soil (m) */
  CANOPYSCHEME Canopy;      /* The stomatal scheme of the vegetation layers */

  /* Used in snow surface energy balance */
  float OldSnowTSurf;       /* Effective surface temperature at the end of the last timestep (C) */
  float SnowTmean;          /* Average snow surface temperature*/
  double Tmp;			    /* Temporary value */
  float Ls;			        /* Latent heat of sublimation (J/kg) */

  /* Edited by Zhuoran Duan zhuoran.duan@pnnl.gov 06/21/2006*/
  /*Add a function to modify soil moisture by add/extract SatFlow from previous time step*/
  DistributeSatflow(Dt, DX, DY, LocalSoil->SatFlow,
    SType->NLayers, LocalSoil->Depth, LocalSoil->LayerThickness,
    LocalSoil->Porosity, LocalSoil->FCap,
    LocalNetwork->Adjust, &(LocalSoil->TableDepth),
    &(LocalSoil->IExcess), LocalSoil->Moist);
  
  /* Temperature memories of the canopy photosynthesis */
  PhotoTemperatureMemory(LocalMet->Tair, (float) Dt, &(LocalVeg->TAcclim), &(LocalVeg->TGrowth));

  /* The understory covers the fraction Fu of the cell and is buried when the
     snow water equivalent exceeds SNOW_DENSITY times its height. The lower
     surface at the start of the step (with the step's snowfall): the
     understory over FuR above the ground (soil or snow), transmitting the
     shortwave with TauU; it is opaque to longwave over Cu */
  NVeg = LocalVeg->NVegLayers;
  Lower = (LocalVeg->OverStory == TRUE) ? 1 : 0;
  Fu = 0.0;
  BurySwe = 0.0;
  Exposed = FALSE;
  if (LocalVeg->UnderStory == TRUE) {
    Fu = LocalVeg->UnderstoryCover;
    BurySwe = SNOW_DENSITY * LocalVeg->Height[Lower];
    Exposed = (LocalSnow->Swq + LocalPrecip->SnowFall <= BurySwe) ? TRUE : FALSE;
  }
  FuR = (LocalSnow->HasSnow == TRUE && Exposed != TRUE) ? 0.0 : Fu;
  Cu = FuR * (1 - LocalVeg->TauU);

  /* Calculate the number of vegetation layers above the snow */
  NVegLActual = LocalVeg->NVegLayers;
  if (LocalVeg->UnderStory == TRUE && FuR == 0.0)
    --NVegLActual;

  /* initialize the total amount of evapotranspiration, and MeltEnergy */
  LocalEvap->ETot = 0.0;
  LocalVeg->MeltEnergy = 0.0;
  LocalVeg->MoistureFlux = 0.0;

  /* calculate the radiation balance for pixels */
  RadiationBalance(Options,
    LocalVeg->OverStory, FuR, Cu, SineSolarAltitude,
    LocalMet->VICSin, LocalMet->Sin, LocalMet->SinBeam,
    LocalMet->SinDiffuse, LocalMet->Lin, LocalMet->Tair, LocalVeg->Tcanopy,
    SType->Albedo, Over, Under, LocalSnow, LocalRad, LocalVeg);

  /* calculate the actual aerodynamic resistances and wind speeds */
  UpperWind = LocalVeg->U[0] * LocalMet->Wind;
  UpperRa = LocalVeg->Ra[0] / LocalMet->Wind;
  if (LocalVeg->OverStory == TRUE)
    LowerRa = LocalVeg->Ra[1] / LocalMet->Wind;
  else
    LowerRa = UpperRa;

  /* calculate the amount of interception storage, and the amount of
     throughfall. Of course this only needs to be done if there is
     vegetation present. */
#ifndef NO_SNOW

  if (LocalVeg->OverStory == TRUE &&
    (LocalPrecip->IntSnow[0] || LocalPrecip->SnowFall > 0.0)) {
    SnowInterception(Options, y, x, Dt, LocalVeg->Fract[0],
      LocalVeg->LAI[0], LocalVeg->MaxInt[0], Over->MaxSnowInt, Over->MDRatio,
      Over->SnowIntEff, UpperRa, LocalMet->AirDens,
      LocalMet->Eact, LocalMet->Lv, LocalRad, LocalMet->Press,
      LocalMet->Tair, LocalMet->Vpd, UpperWind,
      &(LocalPrecip->RainFall), &(LocalPrecip->SnowFall),
      &(LocalPrecip->IntRain[0]), &(LocalPrecip->IntSnow[0]),
      &(LocalPrecip->TempIntStorage),
      &(LocalSnow->CanopyVaporMassFlux), &(LocalVeg->Tcanopy),
      &(LocalVeg->MeltEnergy), LocalVeg->Height, LocalVeg->UnderStory);
    LocalVeg->MoistureFlux -= LocalSnow->CanopyVaporMassFlux;
    /* Because we now have a new estimate of the canopy temperature we can
       recalculate the longwave balance */
    if (LocalSnow->HasSnow == TRUE)
      Tsurf = LocalSnow->TSurf;
    else
      Tsurf = LocalMet->Tair;
    LongwaveBalance(LocalVeg->OverStory, LocalVeg->Fract[0],
      LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, Cu, LocalMet->Tair, LocalRad);

    /* An exposed understory intercepts the rain passing the overstory */
    if (NVegLActual > 1)
      InterceptionStorage(1, &(LocalVeg->MaxInt[1]), &(LocalVeg->Fract[1]),
        &(LocalPrecip->IntRain[1]), &(LocalPrecip->RainFall), &RainOverflow[1]);
  }
  /* if no snow */
  else if (LocalVeg->NVegLayers > 0) {
    LocalVeg->Tcanopy = LocalMet->Tair;
    LocalSnow->CanopyVaporMassFlux = 0.0;
    LocalPrecip->TempIntStorage = 0.0;
    InterceptionStorage(NVegLActual, LocalVeg->MaxInt, LocalVeg->Fract, LocalPrecip->IntRain,
      &(LocalPrecip->RainFall), RainOverflow);
  }

  /* if snow is present, simulate the snow pack dynamics */
  if (LocalSnow->HasSnow || LocalPrecip->SnowFall > 0.0) {
    /* Under an exposed understory the snow gets the ground's shortwave and
       longwave; otherwise all of the lower surface's */
    if (Exposed == TRUE) {
      SnowNetShort = LocalRad->NetShortGround;
      SnowLongIn = LocalRad->LongInGround;
    }
    else {
      SnowNetShort = LocalRad->NetShort[Lower] + LocalRad->NetShortGround;
      SnowLongIn = LocalRad->LongIn[Lower];
    }

    SnowWind = LocalVeg->USnow * LocalMet->Wind;
    SnowRa = LocalVeg->RaSnow / LocalMet->Wind;

    OldSnowTSurf = LocalSnow->TSurf;
    LocalSnow->Outflow =
      SnowMelt(y, x, Dt, 2. + Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, LocalPrecip->RainFall, LocalPrecip->SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &(LocalSnow->PackWater), &(LocalSnow->SurfWater),
        &(LocalSnow->Swq), &(LocalSnow->VaporMassFlux),
        &(LocalSnow->TPack), &(LocalSnow->TSurf), &(LocalVeg->MeltEnergy),
        LocalSnow->MeltScale);

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    SnowTmean = 0.5 * (OldSnowTSurf + LocalSnow->TSurf);
    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, SnowTmean, LocalMet->Tair,
        SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;
    /* convert snow surface temperature from C to K */
    Tmp = SnowTmean + 273.15;
    /* net shortwave radiation */
    LocalSnow->Qsw = SnowNetShort;
    /* net longwave radiation */
    LocalSnow->Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    LocalSnow->Qs = LocalMet->AirDens * CP * (LocalMet->Tair - SnowTmean) / SnowRa;

    /* Calculate latent heat flux */
    if (SnowTmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      LocalSnow->Qe = LocalMet->Lv * LocalSnow->VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * SnowTmean) * JOULESPCAL * GRAMSPKG;
      LocalSnow->Qe = Ls * LocalSnow->VaporMassFlux * WATER_DENSITY;
    }
    LocalSnow->Qe /= Dt;

    /* Calculate advected heat flux from rain */
    LocalSnow->Qp = (CH_WATER * LocalMet->Tair * LocalPrecip->RainFall) / Dt;

    /* Calculate change in cold content */
    LocalSnow->MeltEnergy = LocalVeg->MeltEnergy;
    /* end of snow energy tems */
    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    LocalPrecip->RainFall = 0.0;
    LocalVeg->MoistureFlux -= LocalSnow->VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
       can recalculate the longwave balance (with the understory if it is
       above the remaining snow) */
    Tsurf = LocalSnow->TSurf;
    Cu = (LocalSnow->Swq <= BurySwe) ? Fu * (1 - LocalVeg->TauU) : 0.0;
    LongwaveBalance(LocalVeg->OverStory, LocalVeg->Fract[0],
      LocalMet->Lin, LocalVeg->Tcanopy, Tsurf, Cu, LocalMet->Tair,
      LocalRad);
  }
  else {
    LocalSnow->Outflow = 0.0;
    LocalSnow->VaporMassFlux = 0.0;
    LocalSnow->Qe = 0.;
    LocalSnow->Qs = 0.;
    LocalSnow->Qp = 0.;
    LocalSnow->Qsw = 0.;
    LocalSnow->Qlw = 0.;
    LocalSnow->MeltEnergy = 0.;							   
  }

  /* Determine whether a snow pack is still present, or whether everything
  has melted */
  if (LocalSnow->Swq > 0.0)
    LocalSnow->HasSnow = TRUE;
  else
    LocalSnow->HasSnow = FALSE;

#endif

#ifndef NO_ET
  /* The stomatal scheme: under MEDLYN the kinetics of the cell and step
     (only the dark respiration without shortwave), the CO2 of the year and
     the cold acclimation */
  Canopy.Scheme = Options->StomatalScheme;
  Canopy.TwoLeaf = FALSE;
  if (Canopy.Scheme == MEDLYN && LocalVeg->NVegLayers > 0) {
    PhotoKinetics(LocalMet->Tair, LocalVeg->TGrowth, LocalMet->Press, LocalMet->Sin > 0.0,
                  &(Canopy.Kin));
    Canopy.Ca = ATMOS_CO2;
    Canopy.Acclim = PhotoAcclimation(LocalVeg->TAcclim);
  }

  /* calculate the amount of evapotranspiration from each vegetation layer
     above the ground/soil surface.  Also calculate the total amount of
     evapotranspiration from the vegetation. JARVIS: the visible radiation
     incident at the top of the layer; MEDLYN: the incident beam and diffuse
     visible radiation for the sunlit and shaded leaves */
  if (LocalVeg->OverStory == TRUE) {
    Rp = VISFRACT * LocalMet->Sin;
    NetRadiation = LocalRad->NetShort[0] +
      LocalRad->LongIn[0] - 2 * LocalVeg->Fract[0] * LocalRad->LongOut[0];
    LocalRad->NetRadiation[0] = NetRadiation;
    Canopy.TwoLeaf = TRUE;
    Canopy.SinAlt = SineSolarAltitude;
    Canopy.VisBeam = VISFRACT * LocalMet->SinBeam;
    Canopy.VisDiffuse = VISFRACT * LocalMet->SinDiffuse;
    EvapoTranspiration(0, Dt, LocalMet, NetRadiation,
      Rp, &(Over->Story), SType, LocalVeg->MoistureFlux, LocalSoil->Moist, LocalSoil->Temp,
      &(LocalPrecip->IntRain[0]), LocalEvap->EPot, LocalEvap->EInt, LocalEvap->ESoil,
      LocalEvap->EAct, &(LocalEvap->ETot), LocalNetwork->Adjust, LocalSoil->LayerThickness,
      LocalSoil->FCap, UpperRa, LocalVeg, &Canopy, RainOverflow[0], &(LocalPrecip->RainFall));
    LocalVeg->MoistureFlux += LocalEvap->EAct[0] + LocalEvap->EInt[0];
    LocalEvap->Rc[0] = Canopy.Rc;
    LocalEvap->An[0] = Canopy.An;
    LocalEvap->Beta[0] = Canopy.Beta;
  }				/* end if(LocalVeg->OverStory == TRUE) */

  /* The net shortwave of the whole lower surface, and the net radiation of
     the ground (soil or snow) under and outside the understory */
  NetShortL = LocalRad->NetShort[Lower] + LocalRad->NetShortGround;
  NetRadGround = LocalRad->NetShortGround + LocalRad->LongInGround - LocalRad->LongOutGround;

  /* Understory above the snow at the end of the step. Net radiation from
     its part of the shortwave, the longwave it absorbs from above and from
     the ground, and its emission at the air temperature up and down over
     its opaque fraction; if it was buried at the start of the step, its opaque fraction's share of
     the whole lower surface's shortwave and incoming longwave (the ground
     takes the rest). The visible radiation: JARVIS the incident at the top
     of the understory, MEDLYN the absorbed per covered area; the demand
     above it per covered area */
  if (LocalVeg->UnderStory == TRUE && LocalSnow->Swq <= BurySwe) {
    Tmp = LocalMet->Tair + 273.15;
    LongOutU = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    if (FuR > 0.0) {
      Rp = VISFRACT * LocalRad->NetShort[Lower] / FuR;
      NetRadiation =
        LocalRad->NetShort[Lower] +
        Cu * (LocalRad->LongIn[Lower] + LocalRad->LongOutGround) -
        2 * LocalVeg->Fract[Lower] * (1 - LocalVeg->TauU) * LongOutU;
    }
    else {
      /* The shortwave and incoming longwave of the lower surface shared
         between the understory (its opaque fraction, emitting at the air
         temperature up and down) and the ground below and outside it */
      float CuE = Fu * (1 - LocalVeg->TauU);
      Rp = VISFRACT * (1 - LocalVeg->TauU) * NetShortL;
      NetRadiation =
        CuE * (NetShortL + LocalRad->LongIn[Lower] + LocalRad->LongOutGround) -
        2 * LocalVeg->Fract[Lower] * (1 - LocalVeg->TauU) * LongOutU;
      NetRadGround = (1 - CuE) * (NetShortL + LocalRad->LongIn[Lower]) + CuE * LongOutU -
        LocalRad->LongOutGround;
    }
    if (Canopy.Scheme == JARVIS)
      Rp = VISFRACT * LocalRad->ShortInLower;
    LocalRad->NetRadiation[Lower] = NetRadiation;
    Canopy.TwoLeaf = FALSE;
    EvapoTranspiration(Lower, Dt, LocalMet, NetRadiation,
      Rp, &(Under->Story), SType, LocalVeg->MoistureFlux * LocalVeg->Fract[Lower],
      LocalSoil->Moist, LocalSoil->Temp,
      &(LocalPrecip->IntRain[Lower]), LocalEvap->EPot, LocalEvap->EInt, LocalEvap->ESoil,
      LocalEvap->EAct, &(LocalEvap->ETot), LocalNetwork->Adjust, LocalSoil->LayerThickness,
      LocalSoil->FCap, LowerRa, LocalVeg, &Canopy, RainOverflow[Lower], &(LocalPrecip->RainFall));
    LocalVeg->MoistureFlux += LocalEvap->EAct[Lower] + LocalEvap->EInt[Lower];
    LocalEvap->Rc[Lower] = Canopy.Rc;
    LocalEvap->An[Lower] = Canopy.An;
    LocalEvap->Beta[Lower] = Canopy.Beta;
  }
  else if (LocalVeg->UnderStory == TRUE) {
    int k;
    LocalEvap->EPot[Lower] = 0.;
    for (k = 0; k < MaxSoilLayers; k++)
      LocalEvap->ESoil[Lower][k] = 0.;
    LocalEvap->EAct[Lower] = 0.;
    LocalEvap->EInt[Lower] = 0.;
    LocalRad->NetRadiation[Lower] = 0.;
    LocalEvap->Rc[Lower] = NA;
    LocalEvap->An[Lower] = 0.;
    LocalEvap->Beta[Lower] = NA;
  }

  /* Calculate open water evaporation from surface ponding if present, with
     the ground's net radiation */
  if (LocalSnow->HasSnow != TRUE && LocalSoil->IExcess > 0.) {
    LocalEvap->EvapSoil =
      PondEvaporation(Dt, (DX*DY), LocalNetwork->Area,
                      LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadGround, LowerRa,
                      LocalVeg->MoistureFlux, &(LocalSoil->IExcess));
  }
  else
    LocalEvap->EvapSoil = 0.0;

  /* Calculate soil evaporation from the upper soil layer of the cell outside
     its understory class's litter cover (under and outside the understory)
     if no snow is present: the ground's net radiation, the demand above it
     (and the ponding evaporation) */
  if (LocalSnow->HasSnow != TRUE) {
    if (LocalVeg->UnderStory == TRUE)
      LocalRad->NetRadiation[Lower] += NetRadGround;
    else
      LocalRad->NetRadiation[Lower] = NetRadGround;
    /* New: SoilEvaporation is added to the pre-existing value of EvapSoil,
       which is either equal to SurfaceEvaporation or 0. */
    SoilEvap =
      SoilEvaporation(Dt,
                      LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                      LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadGround, LowerRa,
                      (LocalVeg->MoistureFlux + LocalEvap->EvapSoil),
                      LocalSoil->Porosity[0], LocalSoil->FCap[0], LocalSoil->KsVert[0],
                      LocalSoil->Press[0], LocalSoil->PoreDist[0], SType->WP[0],
                      LocalSoil->LayerThickness[0],
                      &(LocalSoil->Moist[0]), LocalNetwork->Adjust[0],
                      (Under != NULL) ? 1. - Under->LitterCover : 1., &SoilEPot,
                      &(LocalSoil->DesorptionAcc));
    LocalEvap->EvapSoil += SoilEvap;
    LocalEvap->EAct[NVeg] = SoilEvap;
    LocalEvap->EPot[NVeg] = SoilEPot / Dt;
  }
  else {
    if (LocalVeg->UnderStory != TRUE)
      LocalRad->NetRadiation[Lower] = 0.;
    LocalEvap->EAct[NVeg] = 0.;
    LocalEvap->EPot[NVeg] = 0.;
  }
  if (LocalVeg->OverStory != TRUE)
    LocalRad->NetRadiation[1] = 0.;

  LocalVeg->MoistureFlux += LocalEvap->EvapSoil;
  LocalEvap->ETot += LocalEvap->EvapSoil;

  /* Calculate open water and/or dry soil evaporation from stream channels */
  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    NetRadiation = NetRadGround;

    /* Calculate soil evaporation from dry channels if not snow-covered */
    if (LocalSnow->HasSnow != TRUE) {
      DryChannelEvap =
        ChannelSoilEvaporation(Dt, (DX*DY),
                               LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                               LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                               LocalVeg->MoistureFlux,
                               LocalSoil->Porosity, LocalSoil->FCap, LocalSoil->KsVert,
                               LocalSoil->Press, LocalSoil->PoreDist,
                               ((SType->NLayers == LocalNetwork->CutBankZone) ?
                                  LocalSoil->Depth - LocalSoil->TotalThickness :
                                  LocalSoil->LayerThickness[LocalNetwork->CutBankZone]),
                               LocalSoil->Moist, LocalNetwork->Adjust,
                               x, y, ChannelData, LocalNetwork->CutBankZone);
      LocalEvap->EvapSoil += DryChannelEvap;
      LocalVeg->MoistureFlux += DryChannelEvap;
      LocalEvap->ETot += DryChannelEvap;
    }
    
    /* Calculate amount of potential evaporation from flowing channels
       Note that actual subtraction and water balance is handled during routing
       Note that it is assumed that flowing channels are always snow-free */
    ChannelEvaporation(Dt, (DX*DY),
                       LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
                       LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd, NetRadiation, LowerRa,
                       LocalVeg->MoistureFlux, x, y, ChannelData);
  }
  else
    LocalEvap->EvapChannel = 0.0;
#endif

  /* add the water that was not intercepted to the upper soil layer */

#ifndef NO_SOIL

  /* Impervious fraction of the cell: its understory class's (0 without one) */
  float ImpervFrac = (Under != NULL) ? Under->ImpervFrac : 0.0;

  /* This has been modified so that PercArea for infiltration is calculated
     locally to account for the fact that some cells have streams.
     I am not sure if the old PercArea (which does not account for the fact
     that some cells have streams) needs to remain the same.
     Currently, the old PercArea is passed to UnsaturatedFlow */

  MaxInfiltration = 0.;
  ChannelWater = 0.;
  SurfaceWater = 0.;
  PercArea = 1.;
  
  /* ChannelWater is precipitation falling on the channel */
  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
    PercArea = 1. - LocalNetwork->Area / (DX*DY);
    ChannelWater = (LocalNetwork->Area / (DX*DY)) * LocalPrecip->RainFall;
  }

  /* SurfaceWater is rain falling on the hillslope +
     snowmelt on the hillslope (there is no snowmelt on the channel) +
     existing IExcess */
  SurfaceWater = (PercArea * LocalPrecip->RainFall) + LocalSnow->Outflow + LocalSoil->IExcess;
  
  if (InfiltOption == STATIC)
    MaxInfiltration = (1. - ImpervFrac) * PercArea * LocalSoil->MaxInfiltrationRate * Dt;
  else { /* InfiltOption == DYNAMIC
        Dynamic Infiltration Capacity after Parlange and Smith 1978,
        as used in KINEROS and THALES */
    Infiltration = 0.0;
    if (SurfaceWater > 0.) {
      /* Infiltration is a function of the amount of water infiltrated since
     the storm started */
      if (LocalPrecip->PrecipStart) {
        LocalSoil->MoistInit = LocalSoil->Moist[0];
        LocalSoil->InfiltAcc = 0.0;
      }
      /* Check that the B parameter > 0 */
      if ((LocalSoil->InfiltAcc > 0.) && (LocalSoil->Porosity[0] > LocalSoil->MoistInit)) {
        B = (LocalSoil->Porosity[0] - LocalSoil->MoistInit) * (SType->G_Infilt + SurfaceWater);
        Infiltrability = LocalSoil->KsVert[0] * exp((LocalSoil->InfiltAcc) / B) /
          (exp((LocalSoil->InfiltAcc) / B) - 1.);
      }
      else
        Infiltrability = SurfaceWater / Dt;

      MaxInfiltration = Infiltrability * PercArea * (1. - ImpervFrac) * Dt;
      LocalPrecip->PrecipStart = FALSE;
    }/* end  if (SurfaceWater > 0.) */
    else
      LocalPrecip->PrecipStart = TRUE;

  } /* end Dynamic MaxInfiltration calculation */

  Infiltration = (1. - ImpervFrac) * SurfaceWater;

  if (Infiltration > MaxInfiltration)
    Infiltration = MaxInfiltration;
  
  LocalSoil->IExcess = SurfaceWater - Infiltration;

  /* The infiltration rewets the top layer: it offsets the soil evaporation
     accumulated since the last wetting */
  LocalSoil->DesorptionAcc -= Infiltration;
  if (LocalSoil->DesorptionAcc < 0.0)
    LocalSoil->DesorptionAcc = 0.0;

  /*Add water that hits the channel network to the channel network */
  if (ChannelWater > 0.) {
    channel_grid_inc_inflow(ChannelData->stream_map, x, y, ChannelWater * DX * DY);
    LocalSoil->ChannelInt += ChannelWater;
  }

  /* Calculate unsaturated soil water movement, and adjust soil water table depth */
  UnsaturatedFlow(Options, Dt, DX, DY, Infiltration,
    SType->NLayers, LocalSoil->Depth,
    LocalNetwork->Area, LocalSoil->LayerThickness, LocalSoil->KsVert, SType->KsAnisotropy,
    LocalSoil->PoreDist, LocalSoil->Porosity, LocalSoil->FCap, LocalSoil->Perc,
    LocalNetwork->PercArea, LocalNetwork->Adjust, LocalNetwork->CutBankZone,
    LocalNetwork->BankHeight, &(LocalSoil->TableDepth), &(LocalSoil->IExcess),
    LocalSoil->Moist, InfiltOption,
    VegMap, Network, SoilMap,
    LocalTopo->CosSlope, LocalTopo->SinSlope, TopoMap,
    Map, y, x);

  /* Infiltration is updated in UnsaturatedFlow and accumulated
     below */
  if ((InfiltOption == DYNAMIC) && (SurfaceWater > 0.))
    LocalSoil->InfiltAcc += Infiltration;
#endif

  NoSensibleHeatFlux(Dt, LocalMet, LocalVeg->MoistureFlux, LocalSoil);
}
