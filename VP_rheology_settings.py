#!/usr/bin python

### SETTINGS AND PARAMETERS

### Parameters
e=2.0 ## Aspect ratio of the ellipse
e_2=2.0 ## Aspect ratio of the plastic potential (flow rule)
recip_e2=1.0/(e**2.0) #e^-2
SEAICE_strength=2.75e4
tnsFac=0.05 # adding tensile strength, proportinnal to P*
press0=SEAICE_strength # in this case, we assume A=1 and H=1
SEAICEpressReplFac=1.0 # using the replacement pressure

# TD parameters
q=0.5 #TD factor
SEAICEmcMuTD=1.0

## Morh-Coulomb parameters
SEAICEmcMu=0.7
tnsFac_mc=tnsFac
SEAICEmcC=SEAICEmcMu*tnsFac_mc*SEAICE_strength

# tnsFac_zeta=max(tnsFac,tnsFac_mc)

### Using smooth min and max to improve convergence
## For Delta
SEAICE_DELTA_SMOOTHREG=True
deltaMin=1e-10
deltaMinSq=deltaMin**2

## For Zeta
SEAICE_ZETA_SMOOTHREG=True
SEAICE_zetaMaxFac=2.5e8
ZMAX=SEAICE_zetaMaxFac*SEAICE_strength
ZMIN=1e-20

# for Eta
SEAICE_etaMaxFac=2.5e10
# EMAX=SEAICE_etaMaxFac*SEAICE_strength

## For the MC truncating
SEAICE_TMC_SMOOTHMIN=False