/*
 * DESCRIPTION:  Snow depth from the snow water equivalent at the one snow
 *               density of the run ([PARAMETERS] Snow Density), and the
 *               thermal conductivity of snow
 */
#include <stdio.h>
#include "settings.h"
#include "data.h"
#include "constants.h"
#include "snow.h"

/*****************************************************************************
  SnowDepth()

  Depth (m) of Swq (m) of snow at SNOW_DENSITY (kg/m3)
*****************************************************************************/
float SnowDepth(float Swq)
{
  return (float) (WATER_DENSITY * Swq / SNOW_DENSITY);
}

/*****************************************************************************
  SnowConductivity()

  Effective thermal conductivity (W/m/K) of snow of Density (kg/m3): the
  regression of Sturm et al. (1997, J. Glaciol. 43(143), 26-41) on
  488 measurements, rho in g/cm3: k = 0.138 - 1.01 rho + 3.233 rho^2 for
  0.156 <= rho <= 0.6, k = 0.023 + 0.234 rho below 0.156 (0.0591 and
  0.0595 at 0.156)
*****************************************************************************/
float SnowConductivity(float Density)
{
  double Rho = Density / 1000.;

  if (Rho < 0.156)
    return (float) (0.023 + 0.234 * Rho);
  return (float) (0.138 - 1.01 * Rho + 3.233 * Rho * Rho);
}
