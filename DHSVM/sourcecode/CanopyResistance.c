
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "massenergy.h"
#include "constants.h"
#include "photosynthesis.h"

 /*****************************************************************************
   CanopyResistance()
 *****************************************************************************/
float CanopyResistance(float LAI, float RsMin, float RsMax, float Rpc,
  float VpdThres, float MoistThres, float WP, float FCap,
  float TSoil, float SoilMoisture, float Vpd, float Rp)
{
  float MoistFactor;	/* multiplier for resistance due to soil moisture feed-back */
  float Resistance;		/* Canopy resistance (s/m) */
  float RpFactor;		/* multiplier for resistance due to light level feed-back */
  float TFactor;		/* multiplier for resistance due to soil temperaure feed-back */
  float VpdFactor;		/* multiplier for resistance due to vapor pressure deficit feed-back */

  /* for OBS; below 0 C the factor keeps its value at 0 C */
  if (TSoil < 0)
    TSoil = 0;
  TFactor = 1.0 / (0.176 + 0.0770 * TSoil - 0.0018 * TSoil * TSoil);

  if (TFactor <= 0) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }

  /* equation 14, Wigmosta et al [1994] */

  if (Vpd >= VpdThres) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else
    VpdFactor = 1.0 / (1 - Vpd / VpdThres);

  /* equation 15, Wigmosta et al [1994 */

  RpFactor = 1.0 / ((RsMin / RsMax + Rp / Rpc) / (1 + Rp / Rpc));

  /* equation 16, Wigmosta et al [1994]; the threshold MoistThres is relative,
     0 at the wilting point and 1 at field capacity */

  MoistThres = WP + MoistThres * (FCap - WP);
  if (SoilMoisture <= WP) {
    Resistance = DHSVM_HUGE;
    return Resistance;
  }
  else if (SoilMoisture < MoistThres)
    MoistFactor = (MoistThres - WP) / (SoilMoisture - WP);
  else
    MoistFactor = 1.0;

  Resistance = TFactor * VpdFactor * RpFactor * MoistFactor * RsMin / LAI;

  return Resistance;

}

/*****************************************************************************
  CanopyResistanceLayers()

  The canopy resistance associated with each soil layer Rc (s/m) and each
  layer's share of the transpiration of vegetation layer Layer (FCap: the
  cell's field capacity of each layer), and what the scheme reports (Canopy: the canopy resistance, the net assimilation per
  covered area, the soil moisture factor).
  JARVIS: Rc of each layer from its moisture (Wigmosta et al. 1994), the
  shares the root fractions; the reported resistance is their root-weighted
  harmonic mean.
  MEDLYN: one Rc from the canopy conductance to water vapour of the
  photosynthesis (sunlit and shaded leaves of the overstory, one big leaf of
  the understory absorbing Rp) with the soil moisture factor, the shares
  root fraction x layer factor.
*****************************************************************************/
void CanopyResistanceLayers(int Layer, PIXMET *Met, float Rp, STORYTABLE *Story,
                            SOILTABLE *SType, float *FCap, float *Moist, float *SoilTemp,
                            VEGPIX *LocalVeg, CANOPYSCHEME *Canopy, float *Rc,
                            float *Share)
{
  PHOTOSTOMATA Stomata;
  double Gw;			/* Canopy conductance to water vapour (mol/m2/s) */
  double An;			/* Canopy net assimilation (umol/m2/s) */
  double Vcmax25;		/* Top of the canopy, acclimated (umol/m2/s) */
  float Lai = LocalVeg->LAI[Layer];
  float Sum;
  int i;

  if (Canopy->Scheme == JARVIS) {
    Sum = 0.0;
    for (i = 0; i < SType->NLayers; i++) {
      Rc[i] = CanopyResistance(Lai, Story->RsMin, Story->RsMax, Story->Rpc,
        Story->VpdThres, Story->MoistThres, SType->WP[i], FCap[i], SoilTemp[i], Moist[i],
        Met->Vpd, Rp);
      Sum += Story->RootFract[i] / Rc[i];
    }
    Canopy->Rc = 1.0 / Sum;
    Canopy->An = 0.0;
    Canopy->Beta = PhotoMoistureStress(SType->NLayers, Story->RootFract, Moist, SType->WP,
                                       FCap, Story->MoistThres, Share);
    for (i = 0; i < SType->NLayers; i++)
      Share[i] = Story->RootFract[i];
    return;
  }

  Canopy->Beta = PhotoMoistureStress(SType->NLayers, Story->RootFract, Moist, SType->WP,
                                     FCap, Story->MoistThres, Share);
  PhotoStomata(Canopy->Ca, 0.001 * Met->Vpd, Canopy->Beta, Story->G1, Story->G0, &Stomata);
  Vcmax25 = Story->Vcmax25 * Canopy->Acclim;
  if (Canopy->TwoLeaf)
    Gw = PhotoTwoLeafCanopy(&(Canopy->Kin), &Stomata, Vcmax25, Lai, Canopy->SinAlt,
                            Canopy->VisBeam, Canopy->VisDiffuse, &An);
  else
    Gw = PhotoBigLeafCanopy(&(Canopy->Kin), &Stomata, Vcmax25, Lai, Rp, &An);

  /* mol/m2/s to m/s */
  Gw *= PhotoMolarVolume(Met->Tair, Met->Press);
  Canopy->Rc = (Gw > 1.0 / DHSVM_HUGE) ? 1.0 / Gw : DHSVM_HUGE;
  Canopy->An = An;
  for (i = 0; i < SType->NLayers; i++)
    Rc[i] = Canopy->Rc;
}
