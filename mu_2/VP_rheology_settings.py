#!/usr/bin python

### SETTINGS AND DEFAULT PARAMETERS

### Parameters
e_d=2.0 ## Aspect ratio of the ellipse
SEAICE_strength_d=2.75e4
tnsFac_d=0.05 # adding tensile strength, proportinnal to P*
press0_d=SEAICE_strength_d # in this case, we assume A=1 and H=1
SEAICEpressReplFac=0.0 # using the replacement pressure

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
SEAICE_etaMaxFac=2.5e10
# EMAX=SEAICE_etaMaxFac*SEAICE_strength_d

## For the MC truncating
SEAICE_TMC_SMOOTHMIN=False

## For the cap for MCE, MCTD and MCPL
mu_c_d = 4.0

## Default Values for mu(I) rheologies
## H22 refers to Hermam (2022) DOI:10.1098/rsta.2021.0260

# Maximal compaction
# phi_0_d = 1.1
phi_0_d = 1

# rate of drecrease of phi as function of I
# c_phi_d = 1e-2
c_phi_d = 1

# intrinsec density of floes
rho_d = 910

# average size of floes
# d_m_d = 300
d_m_d = 1e4

# mu(I) range for shear deformation
mu_0_d = 0.1 # H22 0.13
mu_i_d = 0.9# H33 0.4

# mu_0_d = 0.1 # H22 0.13
# mu_i_d = 0.9 # H33 0.4

# mu(I) range for dilatant deformatiom
mub_0_d = 5
mub_i_d = 1

# I threshold
# I_0_d = 6.8e-3 # H22 6.8e-3
I_0_d = 1e-3 # H22 6.8e-3

# cohesion (to add in equations)
coh_d = 0

# Constant for the computation of press0 as function of phi
Cstar_d = 7

# Constant for the computation of press0 as function of h
Pmax_d = 2.75e4

#--- Constants for the new friction pressure ---# 

K = 2
c_1 = 1e-1
c_2 = 1/4

# savefig = './VP_5x5_e2/'
savefig = '/aos/home/fstdenis/VP_rheologies/mu_2/'
# Pmax_d = e4