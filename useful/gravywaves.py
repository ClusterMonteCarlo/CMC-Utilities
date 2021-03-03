from scipy import integrate
from scipy.integrate import ode
from astropy import units as U, constants as C
from . import rando as rn 
import numpy as np

def add_default_units(quant,unit):
    if type(quant) == U.quantity.Quantity:
        return quant
    else:
        return quant*unit

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

    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
    a0 = add_default_units(a0,U.AU)
    af = add_default_units(af,U.AU)

    coef = C.G**3 / C.c**5
    beta = (64./5.) * coef * m1 * m2 * (m1+m2)
    
    if e0 == 0:
        if not af == 0:
            print("ERROR: doesn't work for circular binaries")
            return 0
        return a0**4 / (4*beta)
    
    c0 = a0 * (1.-e0**2.) * e0**(-12./19.) * (1.+(121./304.)*e0**2.)**(-870./2299.)

    if af == 0:
        eFinal = 0.
    else:
        r = ode(deda_peters)
        r.set_integrator('lsoda')
        r.set_initial_value(e0,a0.value)
        r.integrate(af.value)
        if not r.successful():
            print("ERROR, Integrator failed!")
        else:
            eFinal = r.y[0]      
    
    time_integrand = lambda e: e**(29./19.)*(1.+(121./304.)*e**2.)**(1181./2299.) / (1.-e**2.)**1.5
    integral,abserr = integrate.quad(time_integrand,eFinal,e0)

    t_merger = (integral * (12./19.) * c0**4. / beta).to(U.Gyr)
    
    if af==0:
        return t_merger 
    else:
        return (t_merger,af,eFinal)

def a_at_fLow(m1,m2,fLow = 5):
    """
    Computes the semi-major axis at an orbital frequency of fLow
    Masses in solar masses, fLow in Hz
    """
    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
    fLow = add_default_units(fLow,1/U.s)

    quant = C.G*(m1+m2) / (4*np.pi**2 * fLow**2)
    return (quant**(1./3)).to(U.AU)

def eccentricity_at_a(m1,m2,a_0,e_0,a):
    """
    Computes the eccentricity at a given semi-major axis a 
    
    Masses are in solar masses, a_0 in AU
    
    """
    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
    a_0 = add_default_units(a_0,U.AU)
    a = add_default_units(a,U.AU)

    r = ode(deda_peters)
    r.set_integrator('lsoda')
    r.set_initial_value(e_0,a_0.value)
    
    r.integrate(a.value)
    
    if not r.successful():
        print("ERROR, Integrator failed!")
    else:
        return r.y[0]

def eccentricity_at_fLow(m1,m2,a_0,e_0,fLow=10):
    """
    Computes the eccentricity at a given fLow, assuming f_gw = 2*f_orb
    
    Masses are in solar masses, a_0 in AU
    
    NOTE!!! The frequency here is the gravitational-wave frequency.  
    """
    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
    a_0 = add_default_units(a_0,U.AU)
    fLow = add_default_units(fLow,1/U.s)

    a_low = a_at_fLow(m1,m2,fLow/2.)

    return eccentricity_at_a(m1,m2,a_0,e_0,a_low)
    

def eccentric_gwave_freq(a,m,e):
    """
    returns the gravitational-wave frequency for a binary at seperation a (AU), mass
    M (solar masses), and eccentricity e, using the formula from Wen 2003
    """

    m = add_default_units(m,U.solMass)
    a = add_default_units(a,U.AU)

    period = 1 / (rn.au_to_period(a,m))
    return (2*period*pow(1+e,1.1954)/pow(1-e*e,1.5)).to(1/U.s)

def eccentricity_at_eccentric_fLow(m1,m2,a_0,e_0,fLow=10,retHigh = False):
    """
    Computes the eccentricity at a given fLow using the peak frequency from Wen
    2003
    
    Masses are in solar masses, a_0 in AU.  The frequency here is the
    gravitational-wave frequency.

    Note that it is possible that, for binaries that merge in fewbody, there is
    no solution, since the binary was pushed into the LIGO band above 10Hz.  In
    this case, there is no a/e from that reference that will give you 10Hz, so
    this just returns e_0 by default, and 0.99 if retHigh is true
    """
    from scipy.optimize import brentq
    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
    a_0 = add_default_units(a_0,U.AU)
    fLow = add_default_units(fLow,1/U.s)

    ecc_at_a = lambda a: eccentricity_at_a(m1,m2,a_0,e_0,a)
    freq_at_a = lambda a: eccentric_gwave_freq(a,m1+m2,ecc_at_a(a))
    zero_eq = lambda a: freq_at_a(a).value - fLow.value #brentq doesn't handle 
                                                         # units well...

    lower_start = zero_eq(1e-10)
    upper_start = zero_eq(1)

    if (np.sign(lower_start) == np.sign(upper_start) or 
        np.isnan(lower_start) or np.isnan(upper_start)):
        if retHigh:
            return 0.99
        else:
            return e_0
    else:    
        a_low = brentq(zero_eq,1e-10,1)
        return ecc_at_a(a_low) 

def kick_this_yo(m1,m2,a1=0,a2=0,theta1=0,theta2=0,phi1=0,phi2=0):
    """
    Does exactly what it says on the jar

    returns (M_final, recoil kick, a_final)
    """

    from numpy.random import uniform

    m1 = add_default_units(m1,U.solMass)
    m2 = add_default_units(m2,U.solMass)
        
    csi = 145*np.pi / 180

    q_ratio = np.minimum(m1,m2)/np.maximum(m1,m2)
    eta = m1*m2 / (m1+m2)**2

    m1 = np.asarray(m1)
    if m1.ndim == 0:
        m1 = m1[None]  # Makes x 1D
        N = 1
    else:
        N = len(m1)
    m1 = add_default_units(m1,U.solMass)

    Theta = uniform(0,np.pi,N)

    theta1 = np.asarray(theta1)
    if theta1.ndim == 0:
        theta1 = theta1[None]
        theta1 = np.arccos(uniform(-1,1,N))
        theta2 = np.arccos(uniform(-1,1,N))

    phi1 = np.asarray(phi1)
    if phi1.ndim == 0:
        phi1 = phi1[None]
        phi1 = uniform(0,2*np.pi,N)
        phi2 = uniform(0,2*np.pi,N)
    
    a1 = np.asarray(a1)
    if a1.ndim == 0:
        a1 = a1[None]  # Makes x 1D
        chi1 = a1*np.ones(N)
        chi2 = a2*np.ones(N)
    else:
        N = len(a1)
        chi1 = a1
        chi2 = a2

    # This needs to be one line so it's done in place
    chi1, chi2 = np.choose(m1>m2,[chi2,chi1]), np.choose(m1>m2,[chi1,chi2])

    delta_x = (q_ratio*chi2*np.sin(theta2)*np.cos(phi2) - chi1*np.sin(theta1)*np.cos(phi1)) /(1+q_ratio)
    delta_y = (q_ratio*chi2*np.sin(theta2)*np.sin(phi2) - chi1*np.sin(theta1)*np.sin(phi1)) /(1+q_ratio)
    delta_z = (q_ratio*chi2*np.cos(theta2) - chi1*np.cos(theta1)) /(1+q_ratio)
               
    chi_x = (q_ratio**2*chi2*np.sin(theta2)*np.cos(phi2) + chi1*np.sin(theta1)*np.cos(phi1)) /(1+q_ratio)**2
    chi_y = (q_ratio**2*chi2*np.sin(theta2)*np.sin(phi2) + chi1*np.sin(theta1)*np.sin(phi1)) /(1+q_ratio)**2
    chi_z = (q_ratio**2*chi2*np.cos(theta2) + chi1*np.cos(theta1)) /(1+q_ratio)**2

    #Compute the approprite mass-weighted spin combinations and their
    #projections parallel and perpendicular to L
    delta_par = delta_z
    chi_par = chi_z
    delta_perp = np.sqrt(delta_y**2 + delta_x**2)
    chi_perp = np.sqrt(chi_y**2 + chi_x**2)
    
    z1 = 1 + (1-chi_par*chi_par)**0.3333333333*((1+chi_par)**0.3333333333+((1-chi_par)**0.3333333333))
    z2 = np.sqrt(3*chi_par*chi_par + z1*z1)
    rISCO = 3 + z2 - np.sign(chi_par)*np.sqrt((3-z1)*(3+z1+2*z2))
    eISCO = np.sqrt(1 - 2. / 3. / rISCO)

    mass_frac = 1 - eta*(1-4*eta)*(1-eISCO) - 16*eta*eta*(0.04827 + 4*0.01707*chi_par*(chi_par + 1))
    
    vm = 1.2e4*(U.km/U.s)* eta*eta * ((1-q_ratio) / (1+q_ratio)) * (1 - 0.93*eta)
    vs_perp = 6.9e3 *(U.km/U.s)* eta*eta * delta_par
    vs_par = 16.*(U.km/U.s)*eta*eta *(delta_perp*(3677.76 + 2*2481.21*chi_par + 4*1792.45*chi_par*chi_par + 8*1506.52*(chi_par**3.))
                             + 2.*chi_perp*delta_par*(1140. + 2*2481.*chi_par))*np.cos(Theta)
    vk = np.sqrt(vm*vm + 2*vm*vs_perp*np.cos(csi) + vs_perp*vs_perp + vs_par*vs_par)

    l = (2*np.sqrt(3) - 3.51712 * eta + 2.5763 * eta**2 - 0.1229 * (1 + q_ratio)**4 / 
            (1 + q_ratio**2)**2 * (chi_perp**2 + chi_par**2) + (0.4537 * eta - 2.8904 + 2.)
            *(1+q_ratio)**2 / (1+q_ratio**2) * chi_par) 

    afinal = np.minimum(1.,np.abs(((q_ratio**2*chi2*np.cos(theta2) + 
        chi1*np.cos(theta1))/((1+q_ratio)**2)) + q_ratio * l / ((1+q_ratio)**2.)));
    
    return mass_frac*(m1+m2),vk.to(U.km/U.s),afinal
