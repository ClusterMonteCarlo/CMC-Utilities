from scipy import integrate
from scipy.integrate import ode
from scipy.optimize import brentq
from numpy.random import uniform

def chirp_mass(m1,m2,z=0):
    return (1+z)*(m1*m2)**0.6 / (m1+m2)**0.2

def eta(m1,m2):
    return (m1*m2)/(m1+m2)**2

def deda_peters(a,e):
    num = 12*a*(1+(73./24)*e**2 + (37./96)*e**4)
    denom = 19*e*(1-e**2)*(1+(121./304)*e**2)
    return denom/num

def inspiral_time_peters(a0,e0,m1,m2,af=0):
    """
    Computes the inspiral time, in Gyr, for a binary
    a0 in Au, and masses in solar masses
    
    if different af is given, computes the time from a0,e0
    to that final semi-major axis
    
    for af=0, just returns inspiral time
    for af!=0, returns (t_insp,af,ef)
    """
    coef = 6.086768e-11 #G^3 / c^5 in au, gigayear, solar mass units
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        if not af == 0:
            print "ERROR: doesn't work for circular binaries"
            return 0
        return a0**4 / (4*beta)
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)
    
    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0)
        r.integrate(af)
        if not r.successful():
            print "ERROR, Integrator failed!"
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)
    
    if af==0:
        return integral * (12./19.) * c0**4. / beta
    else:
        return (integral * (12./19.) * c0**4. / beta,af,eFinal)

def a_at_fLow(m1,m2,fLow = 5):
    """
    Computes the semi-major axis at an orbital frequency of fLow
    Masses in solar masses, fLow in Hz
    """
    G = 3.9652611e-14 # G in au,solar mass, seconds
    quant = G*(m1+m2) / (4*pi**2 * fLow**2)
    return quant**(1./3)

def eccentricity_at_fLow(m1,m2,a_0,e_0,fLow=5):
    """
    Computes the eccentricity at a given fLow
    
    Masses are in solar masses, a_0 in AU
    
    NOTE!!! The frequency here is the ORBITAL frequency.  
    So if you want the eccentricity at a G-wave frequency of 10, set fLow=5
    (divide by 2)
    """
    a_low = a_at_fLow(m1,m2,fLow)
    
    r = ode(deda_peters)
    r.set_integrator('lsoda')
    r.set_initial_value(e_0,a_0)
    
    r.integrate(a_low)
    
    if not r.successful():
        print "ERROR, Integrator failed!"
    else:
        return r.y[0]

def period_to_au(p,m):
    """
    Returns the semi-major axis in AU for a binary
    with mass m (solar masses) and period p (days)
    """
    g = 2.96e-4 # G in days,AU,solar masses
    return (g*m*p**2 / (4*pi**2))**0.333333334

def au_to_period(a,m):
    """
    Returns the period (in days) for a binary
    with mass m (solar masses) and sem-major axis a (AU)
    """
    g = 2.96e-4 # G in days,AU,solar masses
    return sqrt(a**3 * 4*pi**2 / (g*m))

def au_to_vel(a,m):
    """
    Returns the orbital velocity (in ) for a binary
    with mass m (solar masses) and sem-major axis a (AU)
    """
    period = au_to_period(a,m)
    return 2*pi*a / period

def half_mass_relaxation_time(N,rh,m,gamma):
    """
    Returns the half-mass relaxation time of a system
    with N particles, a half-mass radius of rh (parsec),
    an average mass of m, and a Columb Log of gamma
    """
    g = 4499.72292 # G in gyr,pc,solar masses
    return 0.138*sqrt(N * rh**3 / (m*g)) / log(gamma*N)

def find_rtidal(mc,vg=220.,rg=8000.):
    """takes the cluster mass (mc) in solar mass, 
    the galactic circular velocity (vg) in km/s, 
    and the galactocentric distance (rg) in pc
    returns: the tidal radius of the cluster in pc"""
    return (6.67259e-8 * mc*1.99e33 / 2 / (vg*100000)**2.)**(1./3.) * (rg*3.086e18)**(2./3.) / 3.086e18


def comovingDistance(z):
    h = 0.679
    omegaM = 0.306
    omegaK = 0.
    omegaL = 1 - 0.306
    dh = 3000. / h
    e = lambda zp: 1./sqrt(omegaM*(1+zp)**3 + omegaL)
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
    e = lambda zp: 1./(sqrt(omegaM*(1+zp)**3 + omegaL)*(1+zp))
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
    e = sqrt(omegaM*(1+z)**3 + omegaL)
    return 4*pi*dh*comovingDistance(z)**2/e

def Vc(z):
    Dl = luminosityDistance(z)
    return 1.3333333334*pi*Dl**3 / (1+z)**3

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
    
