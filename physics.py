units={}

g=9.81
units['g'] = 'm/s**2'

# Mass density of water
rho_w = 1031.
units['rho_w'] = 'kg/m**3'

cp_w =  4185.5
units['cp_w'] = 'J/(kg*K)'

# Expansion coefficient of water
alpha_w = 207e-6
units['alpha_w'] = '1/K'

# Dynamic viscosity of water
mu_w = 1.08e-3
units['mu_w'] = 'Pa*s'

# Kinematic viscocity of air
nu_a = 1.48e-5
units['nu_a'] = 'm**2/s'

def riseVelocity(diam, rho_d=859.870, rho_w=rho_w, mu=mu_w):
    """
    Calculates the droplet rise velocity for a quiescent fluid

    Parameters
    ----------
    diam: float
        diameter in micrometers
    rho_d: float
        mass density of dropplet
    rho_w: float
        mass density of fluid
    mu: float
        dynamic viscosity of fluid

    Returns
    -------
    rise_velocity: float
        in meters/second
    """
    delta_rho = rho_w - rho_d
    return 1e-12*delta_rho*g*(diam**2.)/(18.*mu)


def get_Ustokes(amp, omega, g=g):
    """
    Defines a Ustokes function of the depth z in meters

    Parameters
    ----------
    amp: float
        amplitude of the waves
    omega: float
        wavenumber

    Returns
    -------
    Ustokes: function
        drift velocity as function of depth
    """
    import numpy as np
    sigma = np.sqrt(g*omega)
    Us0 = sigma*omega*(amp**2.)
    def Ustokes(z, angle=None):
        """
        Returns stokes drift velocity at depth z (meters).
        """
        Us = Us0*np.exp(2.*omega*z)
        if angle:
            angle = np.radians(angle)
            return (Us0*np.cos(angle), Us0*np.sin(angle))
        else:
            return Us
    return Ustokes


def get_wT(Q, cp=cp_w, rho=1000.):
    """
    gets the heat flux in units of w*T (kinematic)
    """
    return Q/(rho*cp)

def get_Q(wt=None, simulation=None, cp=cp_w, rho=rho_w):
    return rho*cp*wt

def Hoennikker(simulation, alpha=alpha_w, g=g, rho=rho_w, cp=cp_w):
    """
    Calculates the Hoennikker (Ho) number
    """
    return None


def droppletTimeScale(diam, rho_d=859.870, mu=1.08e-3):
    return (rho_d + rho_w/2.)*diam**2./(18.*mu)


def get_zi(wT, xy_axis=(1,2), simulation=None):
    """ Calculates the inversion depth as a function of time """
    import numpy as np
    wT = np.asarray(wT)
    z_idx = np.argmin(wT.mean(axis=xy_axis), axis=1)
    if simulation==None:
        return z_idx
    else:
        return simulation.domain.z[z_idx]

def w_star(simulation=None, zi=None, wt_s=None, t_init=None):
    """Calculates the convective scale"""
    sim=simulation
    if type(zi)==type(None):
        zi=sim.inversion_depth
    if type(wt_s)==type(None):
        wt_s=sim.wt_s
    if type(t_init)==type(None):
        t_init=sim.t_init
    w_star = (g*wt_s*zi/t_init)**(1./3.)
    return w_star




