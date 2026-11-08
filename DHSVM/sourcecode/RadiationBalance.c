
/*
Reference:

   Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed
   hydrology-vegetation model for complex terrain, Water Resour. Res.,
   30(6), 1665-1679, 1994.

   Nijssen and Lettenmaier, A simplified approach for predicting shortwave
   radiation transfer through boreal forest canopies, JGR, 1999.

   **references of improved radiation scheme **

   Reifsnyder, W. E., and H. W. Lull (1965), Radiant energy in relation to forests,
   USDA For. Serv. Tech. Bull., 1344, 111.

   Thyer M. et al. (2004), Diagnosing a distributed hydrologic model for two
   high-elevation forested catchments based on detailed stand- and basin-scale data.
   Water Resources Research. DOI:10.1029/2003WR002414.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

 /*****************************************************************************
   Function name: RadiationBalance()

   Purpose      : Calculate the radiation balance for the individual canopy layers

   Required     :
     float Rs           - Incoming shortwave radiation adjusted by topo shading
                          if shading option is on (W/m2)
     float Rsb          - Incoming direct shortwave radiation separated from Rs
     float Rsd          - Incoming diffuse shortwave radiation separated from Rs
     float VIC_Rs       - Input 'observed' shortwave radiation in the forcing file
     float Ld           - Incoming longwave radiation (W/m2)
     float Tair         - Ambient air temperature (C)
     float TGround      - Temperature of the snow-free ground (C): the air
                          temperature, or with SoilFreezing the soil surface
                          temperature of the last step
     float Tcanopy      - Canopy temperature from the previous timestep (C)
     float FuR          - Fraction of the lower surface covered by the
                          understory for the shortwave (0 without an
                          understory or with the understory buried)
     float Cu           - Fraction of the lower surface opaque to longwave
                          by the understory at the air temperature (0
                          without an understory or with the understory
                          buried)
     OVERTABLE *Over    - Overstory class of the cell (NULL if none)
     UNDERTABLE *Under  - Understory class of the cell (NULL if none)
     SNOWPIX LocalSnow  - Information about snow conditions at current pixel
     PIXRAD *LocalRad   - Components of radiation balance for current pixel

   Comments     :
     This routine is implemented according to Wigmosta et al. (1994), with a
     small change for the canopy surface temperature.
     The following assumptions are made:
       - No snow:
         the soil surface is at TGround
       - Snow:
         soil temperature is the snow surface temperature
       - Canopy temperature is equal to the air temperature is there is no snow
         interception
       - The understory is at the air temperature and transmits the
         shortwave with its transmittance of the month; it is opaque to
         longwave over the fraction Cu
       - The overstory is opaque to longwave over F (1 - Taud), its cover
         times one minus its diffuse transmittance (the longwave view
         factor, OverstoryViewFactor())
       - There are at most two vegetation layers

 *****************************************************************************/
void RadiationBalance(OPTIONSTRUCT *Options,
  int OverStory, float FuR, float Cu,
  float SineSolarAltitude, float VIC_Rs, float Rs,
  float Rsb, float Rsd, float Ld,
  float Tair, float TGround, float Tcanopy, float SoilAlbedo,
  OVERTABLE *Over, UNDERTABLE *Under, SNOWPIX *LocalSnow, PIXRAD *LocalRad, VEGPIX *LocalVeg)
{
  float F;			    /* Fraction of pixel covered by top canopy layer [0-1] */
  float Albedo[2];		/* Albedo of each layer */
  float AlbedoG;		/* Albedo of the lower surface outside the understory (snow or soil) */
  float AlbedoL;		/* Albedo of the lower surface */
  float AlbedoU;		/* Albedo of the understory */
  float Tau = 0.0;			/* Transmittance of the overstory for the beam */
  float Taud;			/* Transmittance of the overstory for the diffuse */
  float Tsurf;			/* Surface temperature (C) */

  F = LocalVeg->Fract[0];

  /* Determine Albedo: the lower surface is the understory over the
     fraction FuR, the snow (or the soil without snow) over the rest */
  AlbedoG = (LocalSnow->HasSnow == TRUE) ? LocalSnow->Albedo : SoilAlbedo;
  AlbedoU = (FuR > 0.0) ? Under->Story.Albedo : NOT_APPLICABLE;
  if (FuR > 0.0)
    AlbedoL = FuR * AlbedoU + (1 - FuR) * AlbedoG;
  else
    AlbedoL = AlbedoG;
  Albedo[1] = NOT_APPLICABLE;
  if (OverStory == TRUE) {
    Albedo[0] = Over->Story.Albedo;
    Albedo[1] = AlbedoL;
  }
  else
    Albedo[0] = AlbedoL;


  /* Beam transmittance of the overstory from the cell's LAI of the month
     along the sun's path (Beer's law); diffuse transmittance of the cell
     (SetCellMonth) */
  Taud = NOT_APPLICABLE;
  if (OverStory == TRUE) {
    if (SineSolarAltitude > 0. && Rs > 0.)
      Tau = exp(-Over->ExtnCoeff * LocalVeg->LAI[0] / SineSolarAltitude);
    else
      Tau = 0.;
    Taud = LocalVeg->Taud;
  }
  ShortwaveBalance(OverStory, F, Rs, Rsb, Rsd, Tau, Taud, Albedo,
    FuR, AlbedoU, AlbedoG, LocalVeg->TauU, LocalRad);

  if (LocalSnow->HasSnow == TRUE)
    Tsurf = LocalSnow->TSurf;
  else
    Tsurf = TGround;

  /* Calculate longwave radiation balance */
  LongwaveBalance(OverStory, OverstoryViewFactor(LocalVeg), Ld, Tcanopy, Tsurf, Cu,
    Tair, LocalRad);
  
  LocalRad->PixelLongIn = Ld;

  // Input raw (downward) shortwave radiation without topo shading 
  LocalRad->ObsShortIn = VIC_Rs;

}

/*****************************************************************************
  OverstoryViewFactor()

  Fraction of the cell opaque to longwave by the overstory: its cover F
  times one minus its diffuse transmittance Taud of the day (SetCellLAI),
  the hemispherical integral of the beam transmittance also used for the
  diffuse shortwave. The sky reaches the lower surface through the gaps of
  the crowns. 0 without an overstory.
*****************************************************************************/
float OverstoryViewFactor(VEGPIX *LocalVeg)
{
  if (LocalVeg->OverStory != TRUE)
    return 0.0;
  return LocalVeg->Fract[0] * (1 - LocalVeg->Taud);
}

/************************************************************************************************
  Function name: LongwaveBalance()

  Purpose      : Calculate the longwave radiation balance for the individual
                 canopy layers

  Required     :
    unsigned char OverStory
                     - Flag to indicate presence/absence of overstory
    float F          - Fraction of pixel opaque to longwave by the top
                       vegetation layer (OverstoryViewFactor())
    float Ld         - Downward (incoming) longwave radiation (W/m2)
    float Tcanopy    - Canopy temperature (C)
    float Tsurf      - Surface temperature of the ground (soil or snow) (C)
    float Cu         - Fraction of the lower surface opaque to longwave by
                       the understory at the air temperature, 0 if none
    float Tair       - Air temperature (C)
    PIXRAD &LocalRad - Components of radiation balance for current pixel

  Comments     : This function is used to update the longwave radiation
                 balance when new surface temperatures are calculated.
                 The understory absorbs Cu of the longwave reaching it
                 from above and from the ground and emits Cu of its
                 blackbody longwave up and down; the lower surface's
                 LongOut is the upward longwave above the understory.
************************************************************************************************/
void LongwaveBalance(unsigned char OverStory,
  float F, float Ld, float Tcanopy, float Tsurf, float Cu, float Tair,
  PIXRAD *LocalRad)
{
  double Tmp;
  float LongOutU = 0.0;		/* Blackbody longwave of the understory (W/m2) */
  int Lower = (OverStory == TRUE) ? 1 : 0;

  /* calculate emitted longwave for each layer */
  if (OverStory == TRUE) {
    Tmp = Tcanopy + 273.15;
    LocalRad->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    Tmp = Tsurf + 273.15;
    LocalRad->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
  }
  else {
    Tmp = Tsurf + 273.15;
    LocalRad->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    LocalRad->LongOut[1] = 0.;
  }

  /* lower surface: the ground's emission passing the understory, and the
     understory's at the air temperature over Cu */
  LocalRad->LongOutGround = LocalRad->LongOut[Lower];
  if (Cu > 0.0) {
    Tmp = Tair + 273.15;
    LongOutU = STEFAN * (Tmp * Tmp * Tmp * Tmp);
    LocalRad->LongOut[Lower] = Cu * LongOutU + (1 - Cu) * LocalRad->LongOutGround;
  }

  /* Calculate the incoming longwave for each layer */
  if (OverStory == TRUE) {
    LocalRad->LongIn[0] = (Ld + LocalRad->LongOut[1]) * F;
    LocalRad->LongIn[1] = Ld * (1-F) + LocalRad->LongOut[0] * F;
  }
  else {
    LocalRad->LongIn[0] = Ld;
    LocalRad->LongIn[1] = 0.;
  }

  /* the ground receives the longwave passing the understory and the
     understory's over Cu */
  LocalRad->LongInGround = LocalRad->LongIn[Lower];
  if (Cu > 0.0)
    LocalRad->LongInGround = (1 - Cu) * LocalRad->LongIn[Lower] + Cu * LongOutU;

  /* Calculate the radiative components for the entire pixel.  Use the
     snow/soil surface temperature as an estimate for the pixel temperature.
     LocalRad->PixelLongOut is calculated in the sensible heat flux routine,
     and is not needed otherwise.  Here it is initialized anyway, as if the
     surface temperature is already known */
  LocalRad->PixelLongIn = Ld;
  if (OverStory == TRUE) {
    LocalRad->PixelLongOut = LocalRad->LongOut[0] * F +
      LocalRad->LongOut[1] * (1 - F);
  }
  else {
    LocalRad->PixelLongOut = LocalRad->LongOut[0];
  }
}

/*****************************************************************************
  Function name: ShortwaveBalance()

  Purpose      : Calculate the shortwave radiation balance for the individual
                 pixels

  Required     :
    unsigned char OverStory
                     - Flag to indicate presence of overstory
    float F          - Fraction of pixel covered by top vegetation layer
    float Rs         - Incoming shortwave radiation (W/m2)
    float Rsb, Rsd   - Beam and diffuse parts of Rs (W/m2)
    float Tau        - Beam transmittance of the overstory
    float Taud       - Diffuse transmittance of the overstory
    float *Albedo    - Albedo of each layer
    float FuR        - Fraction of the lower surface covered by the understory
    float AlbedoU    - Albedo of the understory
    float AlbedoG    - Albedo of the ground (soil or snow)
    float TauU       - Shortwave transmittance of the understory
    PIXRAD *LocalRad - Components of radiation balance for current pixel

  Comments     :
    This function needs to be called only once for each pixel for each
    timestep, since the shortwave radiation balance is independent of
    the surface temperatures. The net shortwave of the lower surface is
    split into the understory's part (NetShort of the lower layer) and the
    ground's under and outside the understory (NetShortGround).
*****************************************************************************/
void ShortwaveBalance(unsigned char OverStory,
  float F, float Rs, float Rsb, float Rsd, float Tau, float Taud,
  float *Albedo, float FuR, float AlbedoU, float AlbedoG, float TauU,
  PIXRAD *LocalRad)
{
  int Lower = (OverStory == TRUE) ? 1 : 0;
  float NetShortL;		/* Net shortwave of the lower surface (W/m2) */
  float ShortL;			/* Shortwave reaching the lower surface (W/m2) */

  /* Calculate the net shortwave for each layer */
  /* Overstory present, i.e. two layers */
  if (OverStory == TRUE) {
    /* The beam is transmitted with Tau and the diffuse with Taud, the same
       transmitted shortwave in both layers */
    ShortL = Rs*(1-F) + F*(Rsb*Tau + Rsd*Taud);
    LocalRad->NetShort[0] = F * ((1-Albedo[0])*Rs - (1-Albedo[1])*(Rsb*Tau + Rsd*Taud));
    LocalRad->NetShort[1] =
      (1-Albedo[1])*(Rs*(1-F) + F*(Rsb*Tau + Rsd*Taud));
  }
  else {
    ShortL = Rs;
    LocalRad->NetShort[0] = Rs * (1-Albedo[0]);
    LocalRad->NetShort[1] = 0.;
  }

  LocalRad->ShortInLower = ShortL;

  /* Calculate the net shortwave for the entire pixel */
  if (OverStory == TRUE)
    LocalRad->PixelNetShort = LocalRad->NetShort[0] + LocalRad->NetShort[1];
  else
    LocalRad->PixelNetShort = LocalRad->NetShort[0];

  /* The understory absorbs over FuR what it does not transmit to the
     ground; the ground (soil or snow) the rest of the lower surface's */
  NetShortL = LocalRad->NetShort[Lower];
  if (FuR == 0.0) {
    LocalRad->NetShort[Lower] = 0.0;
    LocalRad->NetShortGround = NetShortL;
  }
  else {
    LocalRad->NetShort[Lower] = FuR * ShortL * ((1 - AlbedoU) - TauU * (1 - AlbedoG));
    LocalRad->NetShortGround = NetShortL - LocalRad->NetShort[Lower];
  }
}
