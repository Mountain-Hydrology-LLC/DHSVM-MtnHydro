/*
 * DESCRIPTION:  C3 photosynthesis and stomatal conductance of leaves and
 *               canopies: temperature kinetics, electron transport, the
 *               coupled FvCB and Medlyn solution of a leaf in closed form,
 *               the sunlit and shaded leaves of the overstory, the big leaf
 *               of the understory, the soil moisture factor and the
 *               temperature memories. Uses only math.h.
 */
#include <math.h>
#include "photosynthesis.h"

static double LargerRoot(double A, double B, double C);

/*****************************************************************************
  PhotoArrhenius(), PhotoHighTemperature()

  Temperature factors, 1 at 25 C (CLM5): exp(Ha/(R Tref) (1 - Tref/Tk)) and
  [1 + exp((Tref dS - Hd)/(R Tref))] / [1 + exp((dS Tk - Hd)/(R Tk))]
*****************************************************************************/
double PhotoArrhenius(double Tk, double Ha)
{
  return exp(Ha / (PHOTO_R * PHOTO_TREF) * (1.0 - PHOTO_TREF / Tk));
}

double PhotoHighTemperature(double Tk, double Hd, double Ds)
{
  return (1.0 + exp((PHOTO_TREF * Ds - Hd) / (PHOTO_R * PHOTO_TREF))) /
    (1.0 + exp((Ds * Tk - Hd) / (PHOTO_R * Tk)));
}

/*****************************************************************************
  PhotoKinetics()

  Kinetics at the air temperature Tair (C) with the growth temperature
  TGrowth (C) and the air pressure Press (Pa); without Light only the dark
  respiration (the rest is not used). Kc and Ko are constant in Pa and
  Gamma* as a mole fraction, so Km in mole fraction rises with altitude
  (rpmodel). The entropy terms and the Jmax25:Vcmax25 ratio follow the growth
  temperature within 11-35 C (Kattge and Knorr 2007, CLM5).
*****************************************************************************/
void PhotoKinetics(float Tair, float TGrowth, float Press, int Light, PHOTOKIN *K)
{
  double Tk = (double) Tair + PHOTO_TFRZ;
  double Tg = TGrowth;
  double Kc, Ko;

  K->Rd = PHOTO_RD_RATIO * PhotoArrhenius(Tk, PHOTO_HA_RD) *
    PhotoHighTemperature(Tk, PHOTO_HD_RD, PHOTO_DS_RD);
  if (!Light) {
    K->Vcmax = K->Jmax = K->Km = K->GammaStar = 0.0;
    return;
  }

  if (Tg < PHOTO_TGROWTH_MIN)
    Tg = PHOTO_TGROWTH_MIN;
  if (Tg > PHOTO_TGROWTH_MAX)
    Tg = PHOTO_TGROWTH_MAX;

  Kc = PHOTO_KC25 * PhotoArrhenius(Tk, PHOTO_HA_KC);
  Ko = PHOTO_KO25 * PhotoArrhenius(Tk, PHOTO_HA_KO);
  K->Km = Kc * (1.0 + PHOTO_O2 * Press / Ko) / Press * 1.0e6;
  K->GammaStar = PHOTO_GAMMASTAR25 * PhotoArrhenius(Tk, PHOTO_HA_GAMMASTAR);

  K->Vcmax = PhotoArrhenius(Tk, PHOTO_HA_VCMAX) *
    PhotoHighTemperature(Tk, PHOTO_HD_VCMAX, PHOTO_DS_VCMAX - PHOTO_DS_VCMAX_T * Tg);
  K->Jmax = (PHOTO_JV_RATIO - PHOTO_JV_RATIO_T * Tg) * PhotoArrhenius(Tk, PHOTO_HA_JMAX) *
    PhotoHighTemperature(Tk, PHOTO_HD_JMAX, PHOTO_DS_JMAX - PHOTO_DS_JMAX_T * Tg);
}

/*****************************************************************************
  PhotoStomata()

  The stomatal terms of a canopy: D = max(VPD, 0.05 kPa); the moisture
  factor Beta multiplies both g1 and g0 (CABLE)
*****************************************************************************/
void PhotoStomata(double Ca, double VpdKpa, double Beta, double G1, double G0,
                  PHOTOSTOMATA *S)
{
  double D = (VpdKpa > PHOTO_VPD_MIN) ? VpdKpa : PHOTO_VPD_MIN;

  S->Ca = Ca;
  S->Beta = Beta;
  S->G0 = Beta * G0;
  S->Slope = 1.0 + Beta * G1 / sqrt(D);
}

/*****************************************************************************
  PhotoElectronTransport()

  J (umol/m2/s) from Jmax and the absorbed PAR (umol/m2/s): the smaller root
  of theta J^2 - (I2 + Jmax) J + I2 Jmax = 0, I2 = 0.5 (1 - f) PAR (de Pury
  and Farquhar 1997, CLM5), in the form without cancellation
*****************************************************************************/
double PhotoElectronTransport(double Jmax, double ParAbs)
{
  double I2 = 0.5 * (1.0 - PHOTO_FNPS) * ParAbs;
  double Sum = I2 + Jmax;
  double Disc;

  if (I2 <= 0.0)
    return 0.0;
  Disc = Sum * Sum - 4.0 * PHOTO_THETA_J * I2 * Jmax;
  if (Disc < 0.0)
    Disc = 0.0;
  return 2.0 * I2 * Jmax / (Sum + sqrt(Disc));
}

/* Larger root of A x^2 + B x + C = 0 (A > 0) */
static double LargerRoot(double A, double B, double C)
{
  double Disc = B * B - 4.0 * A * C;
  double Q, R1, R2;

  if (Disc < 0.0)
    Disc = 0.0;
  Q = -0.5 * (B + ((B >= 0.0) ? sqrt(Disc) : -sqrt(Disc)));
  if (Q == 0.0)
    return -B / (2.0 * A);
  R1 = Q / A;
  R2 = C / Q;
  return (R1 > R2) ? R1 : R2;
}

/*****************************************************************************
  PhotoLeaf()

  Net assimilation (umol/m2/s, returned), stomatal conductance to water
  vapour Gs (mol/m2/s) and intercellular CO2 Ci (umol/mol) of a leaf of
  capacity Vcmax25 (umol/m2/s) absorbing ParAbs (umol/m2/s); Light: Rd in the
  light. An = min(Ac, Aj) - Rd (Farquhar et al. 1980) with
  gs = g0 + 1.6 (1 + g1/sqrt(D)) An/Ca for An > 0, else g0 (Medlyn et al.
  2011), Ci of each limitation in closed form (MAESTRA); where the gross
  RuBP-limited rate is at most Rd (below the light compensation point) that
  limitation is taken at Ci = Ca (plantecophys)
*****************************************************************************/
double PhotoLeaf(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                 double ParAbs, int Light, double *Gs, double *Ci)
{
  double Vcmax = Vcmax25 * K->Vcmax;
  double Rd = Vcmax25 * K->Rd * (Light ? PHOTO_RD_LIGHT : 1.0);
  double VJ = 0.25 * PhotoElectronTransport(Vcmax25 * K->Jmax, ParAbs);
  double Ca = S->Ca;
  double Km = K->Km;
  double Gam = K->GammaStar;
  double Gq = S->Slope / Ca;		/* gs/An on the CO2 basis, An > 0 */
  double G0c = S->G0 / PHOTO_H2O_CO2;
  double Diff = 1.0 - Ca * Gq;
  double A2, B2, C2, CiC, CiJ, AnC, AnJ, An, CiLeaf;

  /* Rubisco limited */
  A2 = G0c + Gq * (Vcmax - Rd);
  B2 = Diff * (Vcmax - Rd) + G0c * (Km - Ca) - Gq * (Vcmax * Gam + Km * Rd);
  C2 = -Diff * (Vcmax * Gam + Km * Rd) - G0c * Km * Ca;
  CiC = LargerRoot(A2, B2, C2);
  AnC = Vcmax * (CiC - Gam) / (CiC + Km) - Rd;

  /* RuBP regeneration limited; taken at Ci = Ca where its gross rate at the
     closed-form Ci is at most Rd, which is so for every Ci < Ca when its
     gross rate at Ca is (the closed form is then skipped) */
  CiJ = Ca;
  AnJ = VJ * (Ca - Gam) / (Ca + 2.0 * Gam) - Rd;
  if (AnJ > 0.0) {
    A2 = G0c + Gq * (VJ - Rd);
    B2 = Diff * (VJ - Rd) + G0c * (2.0 * Gam - Ca) - Gq * (VJ * Gam + 2.0 * Gam * Rd);
    C2 = -Diff * Gam * (VJ + 2.0 * Rd) - G0c * 2.0 * Gam * Ca;
    CiJ = LargerRoot(A2, B2, C2);
    AnJ = VJ * (CiJ - Gam) / (CiJ + 2.0 * Gam) - Rd;
    if (AnJ <= 0.0) {
      CiJ = Ca;
      AnJ = VJ * (Ca - Gam) / (Ca + 2.0 * Gam) - Rd;
    }
  }

  if (AnC < AnJ) {
    An = AnC;
    CiLeaf = CiC;
  }
  else {
    An = AnJ;
    CiLeaf = CiJ;
  }
  *Gs = (An > 0.0) ? S->G0 + PHOTO_H2O_CO2 * Gq * An : S->G0;
  if (Ci != 0)
    *Ci = CiLeaf;
  return An;
}

/*****************************************************************************
  PhotoSunShade()

  Absorbed PAR of the sunlit and shaded leaves of a canopy (de Pury and
  Farquhar 1997 eqs. 13, 18, 20, 21 and A19-A21, spherical leaf angles): from
  the incident beam Ib and diffuse Id at the canopy top (umol/m2/s), the beam
  extinction Kb, the sunlit leaf area LSun and the PAR absorbed by the sunlit
  and the shaded leaves (umol/m2 ground/s)
*****************************************************************************/
void PhotoSunShade(double Lai, double SinAlt, double Ib, double Id,
                   double *Kb, double *LSun, double *IcSun, double *IcSha)
{
  double Sq = sqrt(1.0 - PHOTO_SIGMA);
  double Kbe = PHOTO_G / SinAlt;		/* beam, black leaves */
  double Kbs = Kbe * Sq;			/* beam and its scattering */
  double Kds = PHOTO_KD;			/* diffuse and its scattering */
  double RhoH = (1.0 - Sq) / (1.0 + Sq);
  double RhoCb = 1.0 - exp(-2.0 * RhoH * Kbe / (1.0 + Kbe));
  double EbL = exp(-Kbe * Lai);
  double EbsL = exp(-Kbs * Lai);
  double EdsL = exp(-Kds * Lai);
  double Ic;

  Ic = (1.0 - RhoCb) * Ib * (1.0 - EbsL) + (1.0 - PHOTO_RHO_CD) * Id * (1.0 - EdsL);
  *IcSun = Ib * (1.0 - PHOTO_SIGMA) * (1.0 - EbL) +
    Id * (1.0 - PHOTO_RHO_CD) * (1.0 - EdsL * EbL) * Kds / (Kds + Kbe) +
    Ib * ((1.0 - RhoCb) * (1.0 - EbsL * EbL) * Kbs / (Kbs + Kbe) -
          (1.0 - PHOTO_SIGMA) * (1.0 - EbL * EbL) / 2.0);
  *IcSha = Ic - *IcSun;
  *LSun = (1.0 - EbL) / Kbe;
  *Kb = Kbe;
}

/*****************************************************************************
  PhotoTwoLeafCanopy()

  Canopy conductance to water vapour (mol/m2 ground/s, returned) and net
  assimilation An (umol/m2 ground/s) of a canopy of sunlit and shaded leaves
  from the incident visible beam and diffuse radiation (W/m2) and the sine
  of the solar altitude. Vcmax25 at the top declines with exp(-kn LAI)
  (CLM5): each leaf has the mean of its share. Sun down or no light: every
  leaf shaded, An = -Rd in the dark, gs = g0.
*****************************************************************************/
double PhotoTwoLeafCanopy(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                          double Lai, double SinAlt, double VisBeam, double VisDiffuse,
                          double *An)
{
  double Ib = PHOTO_UMOL_PER_J * VisBeam;
  double Id = PHOTO_UMOL_PER_J * VisDiffuse;
  double EnL = exp(-PHOTO_KN * Lai);
  double VcTotal = Vcmax25 * (1.0 - EnL) / PHOTO_KN;
  double Kb, LSun, LSha, IcSun, IcSha, VcSun, VcSha, Gs, AnLeaf, Gw;

  if (SinAlt <= 0.0 || Ib + Id <= 0.0) {
    *An = -VcTotal * K->Rd;
    return S->G0 * Lai;
  }

  PhotoSunShade(Lai, SinAlt, Ib, Id, &Kb, &LSun, &IcSun, &IcSha);
  LSha = Lai - LSun;
  /* exp(-(kn + kb) L), with exp(-kb L) = 1 - kb LSun */
  VcSun = Vcmax25 * (1.0 - EnL * (1.0 - Kb * LSun)) / (PHOTO_KN + Kb);
  VcSha = VcTotal - VcSun;

  AnLeaf = PhotoLeaf(K, S, VcSun / LSun, IcSun / LSun, 1, &Gs, 0);
  *An = AnLeaf * LSun;
  Gw = Gs * LSun;
  if (LSha > 1.0e-6 * Lai && VcSha > 0.0) {
    AnLeaf = PhotoLeaf(K, S, VcSha / LSha, IcSha / LSha, 1, &Gs, 0);
    *An += AnLeaf * LSha;
    Gw += Gs * LSha;
  }
  return Gw;
}

/*****************************************************************************
  PhotoBigLeafCanopy()

  The same for a canopy of one big leaf absorbing the visible radiation
  VisAbs (W/m2 ground), with the mean Vcmax25 of the canopy
*****************************************************************************/
double PhotoBigLeafCanopy(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                          double Lai, double VisAbs, double *An)
{
  double VcLeaf = Vcmax25 * (1.0 - exp(-PHOTO_KN * Lai)) / (PHOTO_KN * Lai);
  double Gs, AnLeaf;

  if (VisAbs <= 0.0) {
    *An = -VcLeaf * K->Rd * Lai;
    return S->G0 * Lai;
  }
  AnLeaf = PhotoLeaf(K, S, VcLeaf, PHOTO_UMOL_PER_J * VisAbs / Lai, 1, &Gs, 0);
  *An = AnLeaf * Lai;
  return Gs * Lai;
}

/*****************************************************************************
  PhotoMoistureStress()

  Soil moisture factor of a canopy (CABLE): per layer 0 at or below the
  wilting point, otherwise 1 at or above the threshold and linear between
  (tested in this order, as the Jarvis moisture factor, so a threshold at or
  below the wilting point gives a step at the wilting point). Threshold is
  relative: the layer's threshold is WP + Threshold x (FCap - WP); their sum
  weighted by the root fractions, at least 1e-4. Share: each layer's part of
  the transpiration, root fraction x layer factor over their sum (the root
  fractions where the sum is 0)
*****************************************************************************/
float PhotoMoistureStress(int NLayers, const float *RootFract, const float *Moist,
                          const float *WP, const float *FCap, float Threshold,
                          float *Share)
{
  float Sum = 0.0;
  float Layer;
  float Thres;
  int i;

  for (i = 0; i < NLayers; i++) {
    Thres = WP[i] + Threshold * (FCap[i] - WP[i]);
    if (Moist[i] <= WP[i])
      Layer = 0.0;
    else if (Moist[i] >= Thres)
      Layer = 1.0;
    else
      Layer = (Moist[i] - WP[i]) / (Thres - WP[i]);
    Share[i] = RootFract[i] * Layer;
    Sum += Share[i];
  }
  for (i = 0; i < NLayers; i++)
    Share[i] = (Sum > 0.0) ? Share[i] / Sum : RootFract[i];
  return (Sum > PHOTO_BETA_MIN) ? Sum : PHOTO_BETA_MIN;
}

/*****************************************************************************
  PhotoTemperatureMemory()

  The cold-acclimation temperature (time constant 200 h, Kolari et al. 2007)
  and the growth temperature (10 days, CLM5) relax to the air temperature
  over the step of Dt seconds; a memory at or below PHOTO_UNSET (a cold
  start) takes the air temperature
*****************************************************************************/
void PhotoTemperatureMemory(float Tair, float Dt, float *TAcclim, float *TGrowth)
{
  if (*TAcclim <= (float) PHOTO_UNSET)
    *TAcclim = Tair;
  else
    *TAcclim += (Tair - *TAcclim) * Dt / (float) PHOTO_TAU_ACCLIM;
  if (*TGrowth <= (float) PHOTO_UNSET)
    *TGrowth = Tair;
  else
    *TGrowth += (Tair - *TGrowth) * Dt / (float) PHOTO_TAU_GROWTH;
}

/*****************************************************************************
  PhotoAcclimation()

  Photosynthetic capacity relative to summer from the cold-acclimation
  temperature (Kolari et al. 2007, pyAPES): (X + 4.67)/(15 + 4.67) within
  [0.1, 1]
*****************************************************************************/
float PhotoAcclimation(float TAcclim)
{
  float F = (TAcclim - PHOTO_ACCLIM_TBASE) / (PHOTO_ACCLIM_TMAX - PHOTO_ACCLIM_TBASE);

  if (F > 1.0)
    F = 1.0;
  if (F < PHOTO_ACCLIM_FMIN)
    F = PHOTO_ACCLIM_FMIN;
  return F;
}

/*****************************************************************************
  PhotoMolarVolume()

  R T / P (m3/mol): a conductance in mol/m2/s times it is in m/s
*****************************************************************************/
double PhotoMolarVolume(float Tair, float Press)
{
  return PHOTO_R * ((double) Tair + PHOTO_TFRZ) / Press;
}
