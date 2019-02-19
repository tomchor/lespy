def get_W(sim, kappa=0.41):
    """ Get big W based on 2018 GRL paper """
    A_L = 0.816
    A_c = 1.170
    if sim.wt_s>=0:
        W3 = (sim.u_star**3)*((kappa**3) + (A_L**3)/sim.La_t**2) + (A_c**3)*sim.w_star**3
    else:
        phi = 1 + A_L**3/kappa**3/sim.La_t**2 + (A_c**3/kappa**3) * abs(sim.w_star**3)/sim.u_star**3
        W3 = kappa**3 * sim.u_star**3 / phi**3
    return W3**(1/3)

