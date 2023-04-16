#!/usr/bin python
import numpy as np
# from smoothabs import smoothmin_good
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

    return {'e11':e11, 'e22':e22, 'e12':e12, 'e21':e21}

def comp_sim_sr(data):
    '''
    Computing simpler variables
    '''

    data['ep'] = data['e11']+data['e22']
    data['em'] = data['e11']-data['e22']
    data['e12s'] = 0.5*(data['e12']+data['e21'])

    ### Computing deformations in principal stress
    # eTmp=np.sqrt(em**2+4*np.abs(e12*e21))
    data['eTmp'] = np.sqrt(np.abs(data['em']**2+4*data['e12']*data['e21']))
    # eTmp=np.sqrt(em**2+4*e12*e21)

    data['e1'] = 0.5*(data['ep']+data['eTmp'])
    data['e2'] = 0.5*(data['ep']-data['eTmp'])

    ### Computing in stress invariant
    data['eI'] = 0.5*(data['e1']+data['e2'])
    data['eII'] = 0.5*(data['e1']-data['e2'])

    ### Computing deformations in principal stress
    data['eTmps'] = np.sqrt(data['em']**2+4*data['e12s']**2)
    data['e1s'] = 0.5*(data['ep']+data['eTmp'])
    data['e2s'] = 0.5*(data['ep']-data['eTmp'])

    ### Computing in stress invariant
    data['eIs'] = 0.5*(data['e1']+data['e2'])
    data['eIIs'] = 0.5*(data['e1']-data['e2'])

    # return ep, em, e1, e2, eI, eII, e12s, e1s, e2s, eIs, eIIs
    return None


##########################
# RHEOLOGIES (VISCOSITIES)
##########################

def compute_visc(data={},rheos={}):

    for rheo_n in rheos:
        data['rheos'].append(rheo_n)
        data[rheo_n] = {}

        if rheos[rheo_n]['rheo_t'] == 'ell' :
            ellip(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'ell_rot' :
            ellip(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n, rot=True)
        elif rheos[rheo_n]['rheo_t'] == 'mce' :
            mce(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mcs' :
            mcs(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mcpl' :
            mctd(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)
        elif rheos[rheo_n]['rheo_t'] == 'mctd' :
            mcpl(data=data, rheo=rheos[rheo_n], rheo_n = rheo_n)

    return None

def ellip(data={}, rheo={}, rheo_n = '', rot=False):
    '''
    ELLIPTICAL YIELD CURVE
    '''

    # load data
    ep = data['ep']
    em = data['em']
    if rot :
        e12 = data['e12']
        e21 = data['e21']
    else:
        e12 = data['e12s']
        e21 = data['e12s']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'efr' in rheo:
        efr = rheo['efr']
    elif e != 2:
        efr = e
        rheo['efr'] = efr
    else:
        efr = e_d
        rheo['efr'] = efr

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    recip_e2 = 1./(e*e)
    recip_efr2 = 1./(efr*efr)
    efr2_recip_e4 = e**2/(efr**4)

    ### Computing Delta
    # deltaCsq=ep*ep+recip_e2*(em*em+4.0*np.abs(e12*e21))
    deltaCsq=ep*ep+efr2_recip_e4*(np.abs(em*em+4.0*e12*e21))
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
        zeta=ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+kt)
    ## with a sharp min function
    else:
        zeta=press0*(1.+kt)/(2.*deltaCreg)
        zeta=np.minimum(zeta,ZMAX)
    zeta=np.maximum(zeta,ZMIN)

    ### Computing eta
    eta=zeta*recip_efr2

    ### Computing pressure pressure
    press = 0.5*(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+kt))*(1.-kt)

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None

def mce(data={}, rheo={}, rheo_n = ''):
    '''
    MC-E rheology
    '''

    # load data
    ep = data['ep']
    em = data['em']
    e12 = data['e12s']
    eII  = data['eII']

    # load rheo parameters
    if 'e' in rheo:
        e = rheo['e']
    else:
        e = e_d
        rheo['e'] = e

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    recip_e2 = 1/e**2

    # compute delta
    deltaCsq = ep*ep+recip_e2*(np.abs(em*em+4.0*e12**2))

    deltaC = np.sqrt(deltaCsq)
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
        zeta = ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+kt)
    ## with a sharp min function
    else:
        zeta = press0*(1.+kt)/(2.*deltaCreg)
        zeta = np.minimum(zeta,ZMAX)
    zeta = np.maximum(zeta,ZMIN)

    ### Computing pressure pressure
    press = 0.5*(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+kt))*(1.-kt)

    ### Computing eta
    eta = mu*(press-zeta*(ep)+kt*press0)/(2*np.maximum(deltaMin,eII))

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def mcs(data={}, rheo={}, rheo_n = ''):
    '''
    MC-S rheology
    '''

    # load data
    ep = data['ep']
    em = data['em']
    e12 = data['e12s']
    eII = data['eII']

    if 'kt' in rheo:
        kt = rheo['kt']
    else:
        kt = tnsFac_d
        rheo['kt'] = kt

    if 'mu' in rheo:
        mu = rheo['mu']
    else:
        mu = SEAICEmcMu_d
        rheo['mu'] = mu

    if 'press0' in rheo:
        press0 = rheo['press0']
    else:
        press0 = press0_d
        rheo['press0'] = press0

    # compute the viscosities
    zeta = np.minimum(press0*(1+kt)/(2.*deltaMin),press0*(1+kt)/(2.*np.fabs(ep)))

    # press = (press0 * (1.-SEAICEpressReplFac) + SEAICEpressReplFac * 2 * zeta * np.fabs(ep)/(1.+kt))*(1.-kt)
    press = 0.5 * press0 * ( 1. - kt )

    eta = mu*(press-zeta*(ep)+kt*press0)/(2*np.maximum(deltaMin,eII))

    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def mctd(data={}, rheo={}, rheo_n = ''):
    '''
    MC-TD rheology
    '''

    # load data
    ep = data['ep']
    em = data['em']
    e12 = data['e12s']



    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


def mcpl(data={}, rheo={}, rheo_n = ''):
    '''
    MC-PL rheology
    '''

    # load data
    ep = data['ep']
    em = data['em']
    e12 = data['e12s']



    ### save in the dictionary
    data[rheo_n] = rheo
    data[rheo_n]['zeta'] = zeta
    data[rheo_n]['eta'] = eta
    data[rheo_n]['press'] = press

    return None


##########
# STRESSES
##########

def compute_stress(data={}):

    for rheo_n in data['rheos']:
        compu_sigma(data=data, rheo_n=rheo_n)
        comp_str_inva(data=data, rheo_n=rheo_n)

    return None

def compu_sigma(data={}, rheo_n='', plot=False):

    zeta = data[rheo_n]['zeta']
    eta = data[rheo_n]['eta']
    press = data[rheo_n]['press']
    ep = data['ep']
    em = data['em']
    e21 = data['e21']
    e12 = data['e12']

    data[rheo_n]['sig11'] = zeta*ep+eta*em-press
    data[rheo_n]['sig22'] = zeta*ep-eta*em-press
    data[rheo_n]['sig12'] = 2*e21*eta
    data[rheo_n]['sig21'] = 2*e12*eta

    if plot :
        plt.figure('sig11')
        plt.pcolormesh(data[rheo_n]['sig11'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig22')
        plt.pcolormesh(data[rheo_n]['sig22'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig12')
        plt.pcolormesh(data[rheo_n]['sig12'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sig21')
        plt.pcolormesh(data[rheo_n]['sig21'])
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return None

def comp_princ_stress(data={}, rheo_n=''):

    sig11 = data[rheo_n]['sig11']
    sig12 = data[rheo_n]['sig12']
    sig22 = data[rheo_n]['sig22']
    sig21 = data[rheo_n]['sig21']
    press0 = data[rheo_n]['press0']

    sigp=sig11+sig22
    sigm=sig11-sig22
    # sigTmp=np.sqrt(sigm**2+4*np.abs(sig12*sig21))
    sigTmp=np.sqrt(np.abs(sigm**2+4*sig12*sig21))
    # sigTmp=np.sqrt(sigm**2+4*sig12*sig21)

    data[rheo_n]['sig1']=0.5*(sigp+sigTmp)/press0
    data[rheo_n]['sig2']=0.5*(sigp-sigTmp)/press0

    return None

def comp_str_inva(data={}, rheo_n='', plot=False):

    comp_princ_stress(data=data, rheo_n=rheo_n)

    sig1 = data[rheo_n]['sig1']
    sig2 = data[rheo_n]['sig2']

    data[rheo_n]['sigI'] = 0.5*(sig1+sig2)
    data[rheo_n]['sigII'] = 0.5*(sig1-sig2)

    if plot :
        plt.figure('sigI')
        plt.pcolormesh(data[rheo_n]['sigI'])
        plt.axis('equal')
        plt.colorbar()

        plt.figure('sigII')
        plt.pcolormesh(data[rheo_n]['sigII'])
        plt.axis('equal')
        plt.colorbar()

        plt.show()

    return None

##########
# PLOTTING
##########

def plot_stress(data={}):

    fig1=plt.figure('stress states')
    ax = fig1.gca()
    plt.grid()
    plt.axis('equal')
    ax.set_ylabel('Sigma II')
    ax.set_xlabel('Sigma I')

    for rheo_n in data['rheos']:
        plot_inv(data=data, rheo_n=rheo_n, ax=ax, arrows=True)

    ax.legend()

    return None

# def plot_inv(sigI,sigII,eI,eII,opt=None,arrows=False,ax=None,carg=None):
def plot_inv(data={}, rheo_n='', opt=None, arrows=False, ax=None, carg=None):
    '''
    Plotting the yield curve in invariant coordinate
    '''

    sigI = data[rheo_n]['sigI']
    sigII = data[rheo_n]['sigII']
    eI = data['eI']
    eII = data['eII']

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None :
        ax.scatter(sigI, sigII, carg, label=rheo_n)
    else:
        p = ax.plot(sigI.ravel(), sigII.ravel(), '.', ms=1, label=rheo_n)
        carg = p[0].get_color()

    if arrows :
        qpfac=20
        eu=np.hypot(eI[::qpfac,::qpfac],eII[::qpfac,::qpfac])
        ax.quiver(sigI[::qpfac,::qpfac],sigII[::qpfac,::qpfac],eI[::qpfac,::qpfac]/eu,eII[::qpfac,::qpfac]/eu, scale=10, color=carg)

    if opt=='ellipse':
        t=tnsFac
        f=SEAICE_strength
        f=1
        elli=Ellipse(xy=((-f+t)/2.,0), width=(f+t)/e, height=(f+t), angle=-90,edgecolor='b', fc='None', lw=0.5)
        ax.add_patch(elli)

    return None


def plot_FR(data={}):

    fig1=plt.figure('flow rule')
    ax = fig1.gca()
    plt.grid()
    ax.set_xlabel('epsI/epsII')
    ax.set_ylabel('Sigma I')

    for rheo_n in data['rheos']:
        plot_sIFR(data=data, rheo_n=rheo_n, ax=ax)

    ax.legend()

    return None

def plot_sIFR(data={}, rheo_n='', ax=None, carg=None, opt=None):

    '''
    Plotting the the flow rule as function of sI
    '''

    sigI = data[rheo_n]['sigI']
    eI = data['eI']
    eII = data['eII']

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None :
        # ax.plot(sigI.ravel(),(eI/eII).ravel(),'.', color=carg, ms=4, label=rheo_n)
        ax.plot((eI/eII).ravel(),sigI.ravel(),'.', color=carg, ms=4, label=rheo_n, alpha=0.5)
    else:
        # p = ax.plot(sigI.ravel(),(eI/eII).ravel(),'.', ms=4, label=rheo_n)
        p = ax.plot((eI/eII).ravel(),sigI.ravel(),'.', ms=4, label=rheo_n, alpha=0.5)
        carg = p[0].get_color()

    return None


def plot_prAng(data={}):


    fig1=plt.figure('principal axes orientation')
    ax = fig1.gca()
    plt.grid()
    plt.axis('equal')
    ax.set_xlabel('principal stress orientation')
    ax.set_ylabel('principal strain rate orientation')

    for rheo_n in data['rheos']:
        plot_prAng_ori(data=data, rheo_n=rheo_n, ax=ax)

    ax.legend()

    return None

def plot_prAng_ori(data={}, rheo_n='', ax=None, carg=None, opt=None):

    '''
    Plotting the the flow rule as function of sI
    '''

    sig11 = data[rheo_n]['sig12']
    sig22 = data[rheo_n]['sig22']
    sig12 = data[rheo_n]['sig12']

    e11 = data['e11']
    e22 = data['e22']
    e12 = data['e12']

    psi_st = 0.5 * np.arctan2(2*sig12,(sig11-sig22)) * 180/np.pi
    psi_sr = 0.5 * np.arctan2(e12,(e11-e22)) * 180/np.pi

    if ax==None :
        fig1=plt.figure()
        ax = fig1.gca()
        plt.grid()
        plt.axis('equal')

    if carg != None :
        # ax.plot(sigI.ravel(),(eI/eII).ravel(),'.', color=carg, ms=4, label=rheo_n)
        ax.plot(psi_st.ravel(),psi_sr.ravel(),'.', color=carg, ms=4, label=rheo_n, alpha=0.5)
    else:
        # p = ax.plot(sigI.ravel(),(eI/eII).ravel(),'.', ms=4, label=rheo_n)
        p = ax.plot(psi_st.ravel(),psi_sr.ravel(),'.', ms=4, label=rheo_n, alpha=0.5)
        carg = p[0].get_color()

    return None




