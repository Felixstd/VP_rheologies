#!/usr/bin python

### SETTINGS AND DEFAULT PARAMETERS

### Parameters
e_d=2.0 ## Aspect ratio of the ellipse
SEAICE_strength_d=2.75e4
tnsFac_d=0.05 # adding tensile strength, proportinnal to P*
press0_d=SEAICE_strength_d # in this case, we assume A=1 and H=1
SEAICEpressReplFac=1.0 # using the replacement pressure

# TD parameters
q=0.5 #TD factor
SEAICEmcMuTD=1.0

## Morh-Coulomb parameters
SEAICEmcMu_d=0.7
tnsFac_mc_d=tnsFac_d
SEAICEmcC_d=SEAICEmcMu_d*tnsFac_mc_d*SEAICE_strength_d

# tnsFac_zeta=max(tnsFac,tnsFac_mc)

### Using smooth min and max to improve convergence
## For Delta
SEAICE_DELTA_SMOOTHREG=True
deltaMin=1e-10
deltaMinSq=deltaMin**2

## For Zeta
SEAICE_ZETA_SMOOTHREG=True
SEAICE_zetaMaxFac=2.5e8
ZMAX=SEAICE_zetaMaxFac*SEAICE_strength_d
ZMIN=1e-20

# for Eta
SEAICE_etaMaxFac=2.5e10
# EMAX=SEAICE_etaMaxFac*SEAICE_strength_d

## For the MC truncating
SEAICE_TMC_SMOOTHMIN=False