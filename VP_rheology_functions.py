#!/usr/bin python
import numpy as np
from smoothabs import smoothmin_good
from matplotlib.patches import Ellipse
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import copy
from matplotlib.colors import SymLogNorm
from math import copysign

# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')


from VP_rheology_settings import *

def TD(sI,p=1):
    '''
    Shape of the TD yield curve
    - sII/P as function as sI/P
    '''
    return -(sI/p-a)*(1+sI/p)**q

def create_data(random=True,i=1e-6,j=0,plot=False,sym=True,s=201):
    '''
    Creating fake random data
    note that this not the correct way to compute the strain rates on a C grid
    but it is irrelevant for the goal of this code
    '''
    if random :

        np.random.seed(42)
        u=np.random.random((s,s))*i+j
        np.random.seed(24)
        v=np.random.random((s,s))*i+j

        dxy=1000

        dudx=(u[1:,1:]-u[:-1,1:])/dxy
        dvdx=(v[1:,1:]-v[:-1,1:])/dxy
        dudy=(u[1:,1:]-u[1:,:-1])/dxy
        dvdy=(v[1:,1:]-v[1:,:-1])/dxy

        e11=dudx
        e22=dvdy
        if sym:
            e12=0.5*(dudy+dvdx)
            e21=e12
        else:
            e12=dudy
            e21=dvdx

    else:
        eg=np.mgrid[-i:i:200j,i:-i:200j]
        e11=eg[1,:,:]#*0.0
        e22=eg[0,:,:]#*-0.5
        e12=(1.0*e11+1.0*e22)*0.0
    

    if plot :
        plt.figure('e11')
        plt.pcolormesh(e11)
        plt.colorbar()
        plt.axis('equal')

        plt.figure('e22')
        plt.pcolormesh(e22)
        plt.colorbar()
        plt.axis('equal')

        plt.figure('e12')
        plt.pcolormesh(e12)
        plt.colorbar()
        plt.axis('equal')

        if not sym :
            plt.figure('e21')
            plt.pcolormesh(e21)
            plt.colorbar()
            plt.axis('equal')

        plt.show()

    return e11, e22, e12, e21

def comp_sim_sr(e11,e22,e12,e21):
    '''
    Computing simpler variables
    '''

    ep=e11+e22
    em=e11-e22
    e12s=0.5*(e12+e21)

    ### Computing deformations in principal stress
    # eTmp=np.sqrt(em**2+4*np.abs(e12*e21))
    eTmp=np.sqrt(np.abs(em**2+4*e12*e21))
    # eTmp=np.sqrt(em**2+4*e12*e21)

    e1=0.5*(ep+eTmp)
    e2=0.5*(ep-eTmp)

    ### Computing in stress invariant
    eI=0.5*(e1+e2)
    eII=0.5*(e1-e2)

    ### Computing deformations in principal stress
    eTmps=np.sqrt(em**2+4*e12s**2)
    e1s=0.5*(ep+eTmp)
    e2s=0.5*(ep-eTmp)

    ### Computing in stress invariant
    eIs=0.5*(e1+e2)
    eIIs=0.5*(e1-e2)

    return ep, em, e1, e2, eI, eII, e12s, e1s, e2s, eIs, eIIs


##########################
# RHEOLOGIES (VISCOSITIES)
##########################

def ellip(ep,em,e12,e21,e=2):
    '''
    ELLIPTICAL YIELD CURVE 
    '''

    recip_e2 = 1./(e*e) 

    ### Computing Delta
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*np.abs(e12*e21))
    deltaCsq=ep*ep+recip_e2*(np.abs(em*em+4.0*e12*e21))
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*e12*e21)

    deltaC=np.sqrt(deltaCsq)
    ## with a sqrt function
    if SEAICE_DELTA_SMOOTHREG :
      deltaCreg = np.sqrt(deltaCsq + deltaMinSq)
    ## with a sharp min function
    else:
      deltaCreg = np.sqrt(np.maximum(deltaCsq,deltaMinSq))

    ### Computing Zeta
    ## with a smooth tanh function
    if SEAICE_ZETA_SMOOTHREG:
      argTmp = np.exp(-1./(deltaCreg*SEAICE_zetaMaxFac))
      zeta=ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+tnsFac)
    ## with a sharp min function
    else:
      zeta=press0*(1.+tnsFac)/(2.*deltaCreg)
      zeta=np.minimum(zeta,ZMAX)
    zeta=np.maximum(zeta,ZMIN)

    ### Computing eta
    eta=zeta/e**2

    ### Computing pressure pressure
    press=(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+tnsFac))*(1.-tnsFac)

    return zeta, eta, press

##########
# STRESSES
##########

def compu_sigma(e11,e22,e12,e21,zeta,eta,press,plot=False):
    
    ep=e11+e22
    em=e11-e22

    sig11=zeta*ep+eta*em-press/2.0
    sig22=zeta*ep-eta*em-press/2.0
    sig12=2*e21*eta
    sig21=2*e12*eta

    if plot :
        plt.figure('sig11')
        plt.pcolormesh(sig11)
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig22')
        plt.pcolormesh(sig22)
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig12')
        plt.pcolormesh(sig12)
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig21')
        plt.pcolormesh(sig21)
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return sig11, sig12, sig21, sig22

def comp_princ_stress(sig11, sig12, sig22, sig21):

    sigp=sig11+sig22
    sigm=sig11-sig22
    # sigTmp=np.sqrt(sigm**2+4*np.abs(sig12*sig21))
    sigTmp=np.sqrt(np.abs(sigm**2+4*sig12*sig21))
    # sigTmp=np.sqrt(sigm**2+4*sig12*sig21)

    sig1=0.5*(sigp+sigTmp)*1/press0
    sig2=0.5*(sigp-sigTmp)*1/press0

    return sig1,sig2

def comp_str_inva(sig11, sig12, sig22, sig21, plot=False):

    sig1, sig2 = comp_princ_stress(sig11, sig12, sig22, sig21)

    sigI=0.5*(sig1+sig2)
    sigII=0.5*(sig1-sig2)

    if plot :
        plt.figure('sigI')
        plt.pcolormesh(sigI)
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sigII')
        plt.pcolormesh(sigII)
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return sigI, sigII

########## 
# PLOTTING
##########

def plot_inv(sigI,sigII,eI,eII,opt=None,arrows=False,ax=None,carg=None):
    '''
    Plotting the yield curve in invariant coordinate
    '''

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None : 
        ax.plot(sigI,sigII,carg)
    else:
        ax.plot(sigI,sigII)

    if arrows :
        qpfac=1
        eu=np.hypot(eI[::qpfac,::qpfac],eII[::qpfac,::qpfac])
        ax.quiver(sigI[::qpfac,::qpfac],sigII[::qpfac,::qpfac],eI[::qpfac,::qpfac]/eu,eII[::qpfac,::qpfac]/eu,scale=10)

    if opt=='ellipse':
        t=tnsFac
        f=SEAICE_strength
        f=1
        elli=Ellipse(xy=((-f+t)/2.,0), width=(f+t)/e, height=(f+t), angle=-90,edgecolor='b', fc='None', lw=0.5)
        ax.add_patch(elli)

    return None