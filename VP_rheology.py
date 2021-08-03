#!/usr/bin python
import numpy as np
from smoothabs import smoothmin_good
from matplotlib.patches import Ellipse
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from matplotlib.colors import SymLogNorm
from math import copysign

#import functions defining the VP rheologies
from VP_rheology_functions import *

# import settings
from VP_rheology_settings import *

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

# create fake deformation data
e11, e22, e12, e21 = create_data(random=True,i=1e-3,j=0,plot=False,sym=False,s=21)

# compute simpler and additionnal variables
ep, em, e1, e2, eI, eII, e12s, e1s, e2s, eIs, eIIs = comp_sim_sr(e11,e22,e12,e21)


####################
# CHOICE OF RHEOLOGY
####################

ell= True
ell_rot= True

#######################
# COMPUTING VISCOSITIES
#######################

if ell :
    zeta_ell, eta_ell, press_ell = ellip(ep,em,e12s,e12s,e=e)

if ell_rot:
    zeta_ell_rot, eta_ell_rot, press_ell_rot = ellip(ep,em,e12,e21,e=e)

####################
# COMPUTING STRESSES
####################

if ell :

    s11_ell, s12_ell, s21_ell, s22_ell = compu_sigma(e11,e22,e12s,e12s,zeta_ell,eta_ell,press_ell,plot=False)
    sI_ell, sII_ell = comp_str_inva(s11_ell, s12_ell, s22_ell, s21_ell, plot=False) 

if ell_rot :

    s11_ell_rot, s12_ell_rot, s21_ell_rot, s22_ell_rot = compu_sigma(e11,e22,e12,e21,zeta_ell_rot,eta_ell_rot,press_ell_rot,plot=False)
    sI_ell_rot, sII_ell_rot = comp_str_inva(s11_ell_rot, s12_ell_rot, s22_ell_rot, s21_ell_rot, plot=False)


###################
# PLOTTING STRESSES
###################

fig1=plt.figure('stress states')
ax = fig1.gca()
plt.grid()
plt.axis('equal')

if ell:
    plot_inv(sI_ell,sII_ell,eIs,eIIs,opt='ellipse',ax=ax,carg='b.')

if ell_rot :
    plot_inv(sI_ell_rot,sII_ell_rot,eI,eII,opt='ellipse',ax=ax,carg='r.')

#################################
# TEST FOR THE EFFECT OF ROTATION
#################################

if ell and ell_rot :

    sI_diff=sI_ell_rot-sI_ell
    sII_diff=sII_ell_rot-sII_ell

    print('sI_diff mean', np.nanmean(sI_diff))
    print('sII_diff mean', np.nanmean(sII_diff))

    print('sI_diff median', np.nanmedian(sI_diff))
    print('sII_diff median', np.nanmedian(sII_diff))

    print('sI std diff ', np.nanstd(sI_ell_rot)-np.nanstd(sI_ell))
    print('sII std diff ', np.nanstd(sII_ell_rot)-np.nanstd(sII_ell))

    print('sI_diff mean sI>-0.5', np.nanmean(sI_diff[sI_ell>-0.5]))
    print('sII_diff mean sI<-0.5', np.nanmean(sII_diff[sI_ell<-0.5]))

    plt.figure('diff sI')
    plt.pcolormesh(sI_diff)
    plt.axis('equal')
    plt.colorbar()

    plt.figure('diff sII')
    plt.pcolormesh(sII_diff)
    plt.axis('equal')
    plt.colorbar()

    A1=[sI_ell.flatten(),sI_ell_rot.flatten()]
    B1=[sII_ell.flatten(),sII_ell_rot.flatten()]

    plt.figure('diff yield curve')
    plt.plot(A1,B1)
    plt.axis('equal')

    A2=[sI_ell.flatten(),sII_ell.flatten()]
    B2=[sI_ell_rot.flatten(),sII_ell_rot.flatten()]

    plt.figure('arrows diff yield curve')
    for c in range(len(A2[0])):
        # print(c)
        # print(A2[0][c])
        plt.arrow(A2[0][c], A2[1][c], B2[0][c] - A2[0][c], B2[1][c] - A2[1][c],head_width=0.01, length_includes_head=True,ec='k',fc='r')
    plt.axis('equal')

###########
# PLOT SHOW
###########

plt.show()