#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H

/* C3 leaf photosynthesis (Farquhar, von Caemmerer and Berry 1980) coupled to
   the stomatal conductance of Medlyn et al. (2011), the canopy integration
   (sunlit and shaded leaves, de Pury and Farquhar 1997; one big leaf), the
   soil moisture factor and the temperature memories of the canopy */

/* Units, gas constant */
#define PHOTO_R             8.314     /* J/mol/K */
#define PHOTO_TFRZ          273.15    /* K */
#define PHOTO_TREF          298.15    /* K */
#define PHOTO_UMOL_PER_J    4.6       /* visible radiation, umol/J (CLM5) */
#define PHOTO_H2O_CO2       1.6       /* diffusivity ratio of the stomatal pore */

/* Rubisco kinetics at 25 C (Bernacchi et al. 2001; Kc, Ko in Pa as in
   rpmodel, Gamma* as a mole fraction) and activation energies (J/mol) */
#define PHOTO_KC25          39.97     /* Pa */
#define PHOTO_KO25          27480.0   /* Pa */
#define PHOTO_GAMMASTAR25   42.75     /* umol/mol */
#define PHOTO_O2            0.209476  /* mol/mol */
#define PHOTO_HA_KC         79430.0
#define PHOTO_HA_KO         36380.0
#define PHOTO_HA_GAMMASTAR  37830.0

/* Vcmax, Jmax and Rd temperature responses (CLM5; entropy terms and the
   Jmax25:Vcmax25 ratio of Kattge and Knorr 2007 at the growth temperature) */
#define PHOTO_HA_VCMAX      72000.0
#define PHOTO_HD_VCMAX      200000.0
#define PHOTO_DS_VCMAX      668.39    /* J/mol/K, minus PHOTO_DS_VCMAX_T x growth C */
#define PHOTO_DS_VCMAX_T    1.07
#define PHOTO_HA_JMAX       50000.0
#define PHOTO_HD_JMAX       200000.0
#define PHOTO_DS_JMAX       659.70
#define PHOTO_DS_JMAX_T     0.75
#define PHOTO_JV_RATIO      2.59      /* Jmax25/Vcmax25, minus PHOTO_JV_RATIO_T x growth C */
#define PHOTO_JV_RATIO_T    0.035
#define PHOTO_TGROWTH_MIN   11.0      /* C, growth temperature limits */
#define PHOTO_TGROWTH_MAX   35.0
#define PHOTO_RD_RATIO      0.015     /* Rd25/Vcmax25 */
#define PHOTO_HA_RD         46390.0
#define PHOTO_HD_RD         150650.0
#define PHOTO_DS_RD         490.0
#define PHOTO_RD_LIGHT      0.67      /* Rd factor in the light (CLM5) */

/* Electron transport (de Pury and Farquhar 1997; CLM5) */
#define PHOTO_FNPS          0.15      /* absorbed PAR not reaching photosystem II */
#define PHOTO_THETA_J       0.7       /* curvature of J */

/* Stomata */
#define PHOTO_VPD_MIN       0.05      /* kPa (CLM5, CABLE) */
#define PHOTO_BETA_MIN      1.0e-4    /* soil moisture factor floor (CABLE) */

/* Canopy radiation (de Pury and Farquhar 1997) and capacity profile (CLM5) */
#define PHOTO_G             0.5       /* leaf projection, spherical leaf angles */
#define PHOTO_SIGMA         0.15      /* leaf scattering of PAR */
#define PHOTO_KD            0.719     /* diffuse and scattered extinction */
#define PHOTO_RHO_CD        0.036     /* canopy reflection of diffuse PAR */
#define PHOTO_KN            0.30      /* extinction of Vcmax25 per unit LAI */

/* Temperature memories: cold acclimation (Kolari et al. 2007, pyAPES) and
   the 10-day mean air temperature (CLM5) */
#define PHOTO_TAU_ACCLIM    720000.0  /* s (200 h) */
#define PHOTO_ACCLIM_TBASE  -4.67     /* C */
#define PHOTO_ACCLIM_TMAX   15.0      /* C */
#define PHOTO_ACCLIM_FMIN   0.1
#define PHOTO_TAU_GROWTH    864000.0  /* s (10 d) */
#define PHOTO_UNSET         -273.15   /* C: a memory at or below takes the air temperature */

/* Kinetics of a step: the capacities per unit Vcmax25 at the leaf
   temperature, the Michaelis-Menten constant and the compensation point */
typedef struct {
  double Vcmax;      /* Vcmax / Vcmax25 */
  double Jmax;       /* Jmax / Vcmax25 */
  double Rd;         /* Rd / Vcmax25, in the dark */
  double Km;         /* Kc (1 + O/Ko), umol/mol */
  double GammaStar;  /* umol/mol */
} PHOTOKIN;

/* What the stomata of a canopy see */
typedef struct {
  double Ca;         /* atmospheric CO2, umol/mol */
  double Beta;       /* soil moisture factor */
  double G0;         /* intercept x Beta, mol H2O/m2/s */
  double Slope;      /* 1 + Beta g1 / sqrt(D) */
} PHOTOSTOMATA;

double PhotoArrhenius(double Tk, double Ha);
double PhotoHighTemperature(double Tk, double Hd, double Ds);
void PhotoKinetics(float Tair, float TGrowth, float Press, int Light, PHOTOKIN *K);
void PhotoStomata(double Ca, double VpdKpa, double Beta, double G1, double G0,
                  PHOTOSTOMATA *S);
double PhotoElectronTransport(double Jmax, double ParAbs);
double PhotoLeaf(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                 double ParAbs, int Light, double *Gs, double *Ci);
void PhotoSunShade(double Lai, double SinAlt, double Ib, double Id,
                   double *Kb, double *LSun, double *IcSun, double *IcSha);
double PhotoTwoLeafCanopy(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                          double Lai, double SinAlt, double VisBeam, double VisDiffuse,
                          double *An);
double PhotoBigLeafCanopy(const PHOTOKIN *K, const PHOTOSTOMATA *S, double Vcmax25,
                          double Lai, double VisAbs, double *An);
float PhotoMoistureStress(int NLayers, const float *RootFract, const float *Moist,
                          const float *WP, const float *FCap, float Threshold,
                          float *Share);
void PhotoTemperatureMemory(float Tair, float Dt, float *TAcclim, float *TGrowth);
float PhotoAcclimation(float TAcclim);
double PhotoMolarVolume(float Tair, float Press);

#endif
