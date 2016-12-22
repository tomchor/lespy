g=9.81


def riseVelocity(diam, rho_d=859.870, rho_w=1031., mu=1.08e-3):
    """
    Calculates the droplet rise velocity for a quiescent fluid
    """
    delta_rho = rho_w - rho_d
    return 1e-12*delta_rho*g*(diam**2.)/(18.*mu)


def droppletTimeScale(diam, rho_d=859.870, mu=1.08e-3):
    return (rho_d + rho_w/2.)*diam**2./(18.*mu)
