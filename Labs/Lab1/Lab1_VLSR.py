
# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" folder by 5 PM on Jan 31st 2023

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of res (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants

def VLSR(R_0, mu=6.379, vsun=12.24*u.km/u.s):
    '''
    Compute VLSR from a given stellar radius with
        VLSR = 4.74*mu*R_0 - vsun
    Parameters:
        R_0: Solar Radius in kpc (Astropy Quantity)
        mu: Proper Motion of Sag A* in mas/yr (Astropy Quantity),
            default from Reid & Brunthaler 2004
        vpec: peculiar motion of the Sun in the v direction in km/s
              (Astropy Quantity), default from Schonrich 2010
    Returns:
        VLSR: An astropy quantity representing the velocity of the
              local standard of rest in km/s
    '''
    # Compute and return v_lsr with the given formula
    return 4.74*mu*(R_0/u.kpc)*(u.km/u.s) - vsun

print("part Aa:")
R_maser = 8.34*u.kpc # From Reid+2014
R_GravCollab = 8.178*u.kpc # From Abuter+2019 A&A 625
R_SG = 7.9*u.kpc# From Sparke & Gallagher
VLSR_Reid = np.around(VLSR(R_maser), 3)
VLSR_GravCollab = np.around(VLSR(R_GravCollab), 3)
VLSR_SG = np.around(VLSR(R_SG), 3)
print(f"VLSR Calcs:\nMaser: {VLSR_Reid}\nGRAVITY Collaboration: {VLSR_GravCollab}\nSparke & Gallagher: {VLSR_SG}")


# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s $\sim$ 1kpc/Gyr

def getPeriod(vsun, Rsun=8.178*u.kpc ):
    '''
    Gets an orbital period for the Sun at a given radius with
    a given velocity.
    Parameters:
        vsun: Astropy Quantity representing the velocity in km/s
        Rsun: Orbital distance of the sun (Astropy Quantity)
              Default from Abuter+2019
    Returns:
        Orbital Period of the Sun in Gyr as an Astropy Quantity
    '''
    vkpcGyr = vsun.to(u.kpc/u.Gyr) # convert to 
    return (2*np.pi*Rsun) / vkpcGyr # Orbital period


print("part Ab:")
vpec = 12.24*u.km/u.s
R_GravCollab = 8.178*u.kpc # From Abuter+2019 A&A 625
VLSR_GravCollab = np.around(VLSR(R_GravCollab), 3)
vsun = VLSR_GravCollab + vpec
period = getPeriod(vsun)
print("Orbital Period of the Sun:", np.around(period, 3))
print("Part Ac:")
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)
# numRev = Age / T_sun
numRevolutions = 13.8*u.Gyr / period
print("Number of Revolutions in Age of Universe:", np.around(numRevolutions, 2))






# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

# Mass = Integrate rho dV
#      = rho 4*pi*r^2 dr
#      = VLSR^2 / G/4*pi*r^2 * 4pir^2 dr
#      = VLSR**2 / G*r
def MassIso(r, VLSR):
    '''
    This function computes the mass enclosed at a given radius in the isothermal sphere
    model.
    Parameters:
        r: (Astropy Quantity) Distance from the Sun to the galactic center in kpc
        VLSR: (Astropy Quantity) Speed of local standard of rest in km/s
    Returns:
        M: (Astropy Quantity) Mass enclosed in within r in Solar Masses
    '''
    GravConst = const.G.to(u.kpc**3/(u.Gyr**2)/u.Msun) # G in desired units
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr)
    M = VLSRkpcGyr**2 / GravConst * r # Mass within r
    return M
print("part Ba: ")
MIsoSolar = MassIso(R_GravCollab, VLSR_GravCollab)
print(f"Solar Radius: {MIsoSolar:.2e}")
M260kpc = MassIso(260, VLSR_GravCollab)
print(f"260 kpc: {M260kpc:.2e}")






# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

#Potential from a Hernquist Profile: phi = -GM/(r+a)
# Then vesc^2 = 2GM/(r+a)
# So M = (r+a)(vesc^2)/2G
def MassFromVesc(vesc, r, a, ):
    '''
    This function determines the total mass needed for a given escape speed assuming
    a Hernquist Profile for the Dark Matter Halo,
    M = (r+a)(vesc^2)/2G
    Params:
        vesc: (Astropy Quantity) Escape Velocity at the given radius in km/s
        r: (Astropy Quantity) Distance from galactic center in kpc
        a: (Astropy Quantity) Hernquist scale length in kpc
    Returns:
        M: (Astropy Quantity) Total Mass within r in Msun
    '''
    GravConst = const.G.to(u.kpc**3/(u.Gyr**2)/u.Msun) # G in desired units
    vesckpcGyr = vesc.to(u.kpc/u.Gyr) # Converting vesc to proper units for G
    M = vesckpcGyr**2/2/GravConst*(r+a) # Required Mass to fit a Hernquist profile
    return M

Leo1v = 196*u.km/u.s # Leo 1 speed (Sohn+2013)
a = 30*u.kpc #Scale Radius for the Hernquist Halo
r = 260*u.kpc #Galactocentric distance of Leo 1
MLeo1 = MassFromVesc(Leo1v, r, a)
print(f"Hernquist Mass: {MLeo1:.2e}")
print(f"IsoSphere/Mleo1 = {M260kpc / MLeo1}")