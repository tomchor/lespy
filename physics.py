units={}

g=9.81
units['g'] = 'm/s**2'

# Mass density of water
rho_w = 1031.
units['rho_w'] = 'kg/m**3'

# Mass density of oil
rho_oil=859.870
units['rho_oil'] = 'kg/m**3'

# Heat capacity
cp_w =  4185.5
units['cp_w'] = 'J/(kg*K)'

# Expansion coefficient of water
alpha_w = 2.0e-4
units['alpha_w'] = '1/K'

# Dynamic viscosity of water
mu_w = 1.08e-3
units['mu_w'] = 'Pa*s'

# Kinematic viscocity of air
nu_a = 1.48e-5
units['nu_a'] = 'm**2/s'


def get_R(d, rho=rho_w, delta_rho=(rho_w-rho_oil), g=9.81, mu=mu_w):
    ''' Gets the R variable, used to calculate rise velocity '''
    import numpy as np
    Nd = 4*rho_w*delta_rho*g*(d**3)/(3*mu_w**2)
    conds=[ Nd<=73, (73<Nd)*(Nd<=580), (580<Nd)*(Nd<=1.55e7) ]
    funcs=[]
    funcs.append( lambda Nd: Nd/24 - 1.7569e-4*(Nd**2) + 6.9252e-7*(Nd**3) - 2.3027e-10*(Nd**4) )
    funcs.append( lambda Nd: np.power(10, -1.7095 + 1.33438*np.log10(Nd) - 0.11591*np.log10(Nd)**2) )
    funcs.append( lambda Nd: np.power(10, -1.81391 + 1.34671*np.log10(Nd) - 0.12427*np.log10(Nd)**2 + 0.006344*np.log10(Nd)**3) )
    return np.piecewise(Nd, conds, funcs)


def termVelocity(d, mu=mu_w, rho=rho_w, delta_rho=(rho_w-rho_oil)):
    """
    Gets rise velocity of droplets according to Li Zheng and Poojitha D. Yapa, 2000.
    rho is the density o the fluid (in our case either water or air; default os water)
    """
    import numpy as np
    d=np.asarray(d)
    conds = [ d<=1e-3, (1e-3<d)*(d<15.e-3) ]
    funcs = []
    funcs.append( lambda d: get_R(d, mu=mu, rho=rho, delta_rho=delta_rho)*mu/(rho*d))
    funcs.append( None )
    return np.piecewise(d, conds, funcs)


def get_dropletSize(wr, mu=mu_w, rho=rho_w, delta_rho=(rho_w-rho_oil), nominal=False, nowarning=False, level=5):
    """
    Calculates the droplet size in micrometers for a given set of terminal velocities
    """
    import numpy as np
    D=np.linspace(1e-8, 1e-3, 500)
    Wr = termVelocity(D, mu=mu, rho=rho, delta_rho=delta_rho)
    diam = np.interp(wr, Wr, D)
    if nominal:
        if not nowarning: print('In micrometers!')
        return np.around(10**(level)*diam, decimals=0)*10**(6-level)
    else:
        return diam



def get_Ustokes(amp, omega, g=g, sigma=None):
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
    if type(sigma)==type(None):
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
        t_init=1/alpha_w
    w_star = (g*wt_s*zi/t_init)**(1./3.)
    return w_star




