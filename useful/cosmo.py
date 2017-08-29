from scipy import integrate
from scipy.optimize import brentq
from numpy.random import uniform
import numpy as np

def comovingDistance(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    dh = 3000. / h
    e = lambda zp: 1./np.sqrt(omegaM*(1+zp)**3 + omegaL)
    return dh*integrate.quad(e,0,z)[0]

def luminosityDistance(z):
     return (1+z)*comovingDistance(z)
    
def zAtLuminosityDistance(d):
    zero = lambda z: luminosityDistance(z) - d
    return brentq(zero,0,5)    
    
def lookbackTime(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    th = 9.78/h
    e = lambda zp: 1./(np.sqrt(omegaM*(1+zp)**3 + omegaL)*(1+zp))
    return th*integrate.quad(e,0,z)[0]

def zAtLookbackTime(t):
    zero = lambda z: lookbackTime(z) - t
    return brentq(zero,0,10)

def dVcdz(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    dh = 3000. / h
    e = np.sqrt(omegaM*(1+z)**3 + omegaL)
    return 4*np.pi*dh*comovingDistance(z)**2/e

def Vc(z):
    Dl = luminosityDistance(z)
    return 1.3333333334*np.pi*Dl**3 / (1+z)**3

def probZFlatInVc(z,Zmax):
    return dVcdz(z)/Vc(Zmax)

def zFlatInVc(zMax,size=1):
    VcMax = Vc(zMax)
    randomVol = uniform(0,VcMax,size=size)

    if size==1:
        zero = lambda z: Vc(z) - randomVol
        return brentq(zero,0,zMax)
    else:
        zs = []
        for r in randomVol:
            zero = lambda z: Vc(z) - r
            zs.append(brentq(zero,0,zMax))
        return array(zs)
    
