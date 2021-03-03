import numpy as np
from astropy import units as U, constants as C

def add_default_units(quant,unit):
	if type(quant) == U.quantity.Quantity:
		return quant
	else:
		return quant*unit

def period_to_au(p,m):
    """
    Returns the semi-major axis in AU for a binary
    with mass m (solar masses) and period p (days)
    """
    m = add_default_units(m,U.solMass)
    p = add_default_units(p,U.day)
    return ((C.G*m*p**2 / (4*np.pi**2))**(1/3)).to(U.AU)

def au_to_period(a,m):
    """
    Returns the period (in days) for a binary
    with mass m (solar masses) and sem-major axis a (AU)
    """
    m = add_default_units(m,U.solMass)
    a = add_default_units(a,U.AU)
    return np.sqrt(a**3 * 4*np.pi**2 / (C.G*m)).to(U.day)

def au_to_vel(a,m):
    """
    Returns the orbital velocity (in km/s) for a binary
    with mass m (solar masses) and sem-major axis a (AU)
    """
    a = add_default_units(a,U.AU)
    period = au_to_period(a,m)
    return (2*np.pi*a / period).to(U.km/U.s)

def half_mass_relaxation_time(N,rh,m,gamma):
    """
    Returns the half-mass relaxation time of a system
    with N particles, a half-mass radius of rh (parsec),
    an average mass of m, and a Columb Log of gamma
    """
    g = 4499.72292 # G in gyr,pc,solar masses
    return 0.138*np.sqrt(N * rh**3 / (m*g)) / np.log(gamma*N)

def find_rtidal(mc,vg=220.,rg=8000.):
    """takes the cluster mass (mc) in solar mass, 
    the galactic circular velocity (vg) in km/s, 
    and the galactocentric distance (rg) in pc
    returns: the tidal radius of the cluster in pc"""
    return (6.67259e-8 * mc*1.99e33 / 2 / (vg*100000)**2.)**(1./3.) * (rg*3.086e18)**(2./3.) / 3.086e18

def plummer_m_r(M,a,r):
    """
    Plummer mass internal to r
    """
    return M*(1+a**2 / r**2)**-1.5

def plummer_rho(M,a,r):
    """
    Plummer density at r
    """
    return 3*M*(1+r**2/a**2)**-2.5 / (4*np.pi*a**3)

def plummer_pot(M,a,r):
    """
    Plummer potential at r
    """
    return (-M/a) / np.sqrt(1+r**2/a**2)

def plummer_sigma(M,a,r):
    """
    Plummer 1D velocity dispersion at r 
    """
    return np.sqrt(-plummer_pot(M,a,r) / 6)
    
def plummer_vc(M,a,r):
    """
    Circular velocity for plummer sphere at r
    """
    return np.sqrt(plummer_m_r(M,a,r) / r)

def chirp_mass(m1,m2):
    """
        returns chirp mass from component masses
    """
    return (m1*m2)**0.6 / (m1+m2)**0.2

def eta(m1,m2):
    """
        returns symmetric mass ratio from component masses
    """
    return (m1*m2)/(m1+m2)**2

def m1_from_mchirp_eta(mchirp,eta):
    """
        returns larger component mass from chirp mass and eta 
    """
    prefac = mchirp*eta**(-0.6)/2.
    return prefac*(1+np.sqrt(1-4*eta))

def m2_from_mchirp_eta(mchirp,eta):
    """
        returns smaller component mass from chirp mass and eta 
    """
    prefac = mchirp*eta**(-0.6)/2.
    return prefac*(1-np.sqrt(1-4*eta))
