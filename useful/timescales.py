import astropy
import numpy as np
from astropy import units as u
from astropy import constants as c

def add_default_units(quant,unit):
    if type(quant) == astropy.units.quantity.Quantity:
        return quant
    else:
        return quant*unit

def perihelion_precession(a,e,M):
	"""
	Computes and returns the timescale for a single revolution of the periapse
	of a binary

	Time is in seconds, M in solar masses, a in AU
	"""
	G = 3.9652611e-14 # G in au,solar mass, seconds
	c = 0.0020039888 # C in au/sec

	return c**2 * a**2.5 * (1 - e**2) / (3*(G*M)**1.5)

def spinorbit_precession(a,e,M,chi_eff):
	"""
	Computes and returns the timescale for a single revolution of L about J for
	a spinning binary (the Lense-Thirring contribution)

	Time is in seconds, M in solar masses, a in AU. chi_eff is [-1,1]
	"""
	G = 3.9652611e-14 # G in au,solar mass, seconds
	c = 0.0020039888 # C in au/sec

	return c**3 * a**3 * (1 - e**2)**1.5 / (2*chi_eff*(G*M)**2)

def binary_single(a,M,sigma,n):
    """
    Computes the interaction timescale between binary and single objets

    default units are: 
        a is in AU, M is in solar masses, sigma in km/s, and n in pc^-3
    """

    a = add_default_units(a,u.AU)
    M = add_default_units(M,u.solMass)
    sigma = add_default_units(sigma,u.km/u.s)
    n = add_default_units(n,(u.pc)**-3)

    quantity = (n * a**2 * sigma * (1 + c.G * M/(2*a*sigma**2) ))**-1

    return quantity.to(u.yr)

def soft_binary_lifetime(a,sigma,rho):
    """
    Computes the expected lifetime for soft binaries in a cluster 
    Binney and Tremaine, 2nd Edition, Eqn 7.172

    default units are: 
        a ~ AU
        sigma ~ km/s
        rho ~ msun / pc^3
    """

    a = add_default_units(a,u.AU)
    sigma = add_default_units(sigma,u.km/u.s)
    rho = add_default_units(rho,u.solMass/(u.pc)**3)

    return (0.037 * sigma / (c.G * rho * a)).to(u.yr)

def kozai_timescale(m1,m2,m3,a1,a2,e1):
    """
    Standard Quadrupole-order Kozai timescale
    m1,m2 are binary, m3 is outer

    default units are:
        a1,a2 ~ AU
        m1,m2,m3 ~ msun
    """

    m1 = add_default_units(m1,u.solMass)
    m2 = add_default_units(m2,u.solMass)
    m3 = add_default_units(m3,u.solMass)
    a1 = add_default_units(a1,u.AU)
    a2 = add_default_units(a2,u.AU)

    t_kozai = ((m1+m2)/m3)*(a2/a1)**3*np.sqrt(a1**3*(1-e1**2)**3 / 
            c.G / (m1+m2))

    return t_kozai.to(u.yr)

def half_mass_relaxation_time(N,rh,m,gamma):
    """
    Returns the half-mass relaxation time of a system
    with N particles, a half-mass radius of rh (parsec),
    an average mass of m, and a Columb Log of gamma
    """
    rh = add_default_units(rh,u.pc)
    m = add_default_units(m,u.solMass)

    t_rh = 0.138*np.sqrt(N * rh**3 / (m*c.G)) / np.log(gamma*N) 

    return t_rh.to(u.Gyr)

#peri_valid = ((init[4]/init[3])**3 > (3 * c.c**2 * init[2] * init[3]*u.AU / (4*c.G*(init[0]+init[1])**2*u.solMass) *
