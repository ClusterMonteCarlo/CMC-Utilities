import numpy as np

def period_to_au(p,m):
	"""
	Returns the semi-major axis in AU for a binary
	with mass m (solar masses) and period p (days)
	"""
	g = 2.96e-4 # G in days,AU,solar masses
	return (g*m*p**2 / (4*np.pi**2))**0.333333334

def au_to_period(a,m):
	"""
	Returns the period (in days) for a binary
	with mass m (solar masses) and sem-major axis a (AU)
	"""
	g = 2.96e-4 # G in days,AU,solar masses
	return np.sqrt(a**3 * 4*np.pi**2 / (g*m))

def au_to_vel(a,m):
	"""
	Returns the orbital velocity (in ) for a binary
	with mass m (solar masses) and sem-major axis a (AU)
	"""
	period = au_to_period(a,m)
	return 2*np.pi*a / period

def half_mass_relaxation_time(N,rh,m,gamma):
	"""
	Returns the half-mass relaxation time of a system
	with N particles, a half-mass radius of rh (parsec),
	an average mass of m, and a Columb Log of gamma
	"""
	g = 4499.72292 # G in gyr,pc,solar masses
	return 0.138*np.sqrt(N * rh**3 / (m*g)) / log(gamma*N)

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
	return 3*M*(1+r**2/a**2)**-2.5 / (4*pi*a**3)

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
