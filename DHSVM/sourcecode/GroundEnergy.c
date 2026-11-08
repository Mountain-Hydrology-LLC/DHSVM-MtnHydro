/*
 * DESCRIPTION:  Energy balance of the ground (SoilFreezing), as in CLM 5
 *               (Lawrence et al. 2018, "Soil and Snow Temperatures"): the
 *               top temperature node is the ground surface. Without snow it
 *               takes the net radiation, the sensible heat and the latent
 *               heat of the ground evaporation (fixed), linearised about its
 *               temperature at the start of the step; under snow the heat
 *               conducted from the snow, implicit in the ground. One solve
 *               of the node temperatures with phase change follows
 *               (SolveProfile()).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "constants.h"
#include "soiltemp.h"

#define GROUND_Z 2.0            /* Reference height of the stability correction (m) */

/*****************************************************************************
  GroundEnergyBalance()

  Node temperatures and ice of a cell at the end of the step, with the
  surface temperature, layer ice, layer temperatures and frost depth that
  follow. Snow: TRUE under snow, the ground then exchanges heat only with
  the snow at SnowTemp (C) through the conductance SnowCond (W/m2/K).
  Otherwise NetRadGround: the net radiation of the ground with the emission
  at the last step's surface temperature (LongOutGround); GroundEvap: the
  evaporation from the ground this step (m). Writes the surface terms into
  the cell's TSurf, Qnet, Qs, Qg, Qst, Qrest and Ra. Returns the heat flux
  from the ground into the snow (W/m2, 0 without snow).
*****************************************************************************/
float GroundEnergyBalance(int y, int x, int Dt, PIXMET *LocalMet, SOILTABLE *SType,
                          SOILPIX *Soil, int Snow, float SnowTemp, float SnowCond,
                          float NetRadGround, float LongOutGround, float LowerRa,
                          float LowerWind, float GroundEvap)
{
  SOILNODES Nodes;
  double T0 = Soil->TNode[0], Tk = T0 + 273.15;
  double Rn = 0., H = 0., LatentHeat = 0., Ra = DHSVM_HUGE;
  double DRn = 0., DH = 0., H0, Beta, Into, Storage, Rest, Change;

  /* Node water from the layers now; the node ice at most the node water */
  NodeMoisture(SType->NLayers, Soil->LayerThickness, Soil->Depth, Soil->Moist,
               Soil->Porosity, Soil->Press, Soil->PoreDist, SType->Quartz,
               SType->Organic, Soil->IceNode, &Nodes);

  /* The flux into the ground at the start-of-step surface temperature and
     its derivative */
  if (Snow) {
    H0 = SnowCond * (SnowTemp - T0);
    Beta = -SnowCond;
  }
  else {
    Rn = NetRadGround + LongOutGround - STEFAN * Tk * Tk * Tk * Tk;
    DRn = -4. * STEFAN * Tk * Tk * Tk;
    if (LowerWind > 0.)
      Ra = LowerRa / StabilityCorrection(GROUND_Z, 0.f, (float) T0, LocalMet->Tair,
                                         LowerWind, Z0_GROUND);
    H = LocalMet->AirDens * CP * (LocalMet->Tair - T0) / Ra;
    DH = -LocalMet->AirDens * CP / Ra;
    LatentHeat = -LocalMet->Lv * WATER_DENSITY * GroundEvap / Dt;
    H0 = Rn + H + LatentHeat;
    Beta = DRn + DH;
  }

  Into = SolveProfile(Soil->TNode, Soil->IceNode, Soil->TBot, &Nodes, (double) Dt,
                      H0, Beta, &Storage, &Rest);
  SoilThermalDerived(Soil, SType->NLayers);

#ifdef SOIL_ENERGY_DEBUG
  if (fabs(Rest) > 1e-6)
    fprintf(stderr, "GroundEnergyBalance: y %d x %d energy balance error %g W/m2\n",
            y, x, Rest);
#else
  (void) y;
  (void) x;
#endif

  /* The fluxes at the surface temperature of the conduction solve, as they
     entered the ground */
  Change = (Beta != 0.) ? (Into - H0) / Beta : 0.;
  Soil->TSurf = Soil->TNode[0];
  Soil->Ra = (float) Ra;
  Soil->Qnet = (float) (Rn + DRn * Change);
  Soil->Qs = (float) (H + DH * Change);
  Soil->Qg = (float) -Into;
  Soil->Qst = (float) Storage;
  Soil->Qrest = (float) Rest;
  return Snow ? (float) -Into : 0.0f;
}
