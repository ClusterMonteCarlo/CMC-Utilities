from scipy import integrate
from scipy.integrate import ode
import numpy as np

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
    quant = G*(m1+m2) / (4*np.pi**2 * fLow**2)
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
