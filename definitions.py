def get_dicts(sim, nz=None):
    if nz is None:
        nz=sim.nz_tot

    vdict = dict()
    vdict["u"] = dict(name="u", 
                      dims=["x", "y", "z"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_u[:nz],)),

    vdict["v"] = dict(name="v", 
                      dims=["x", "y", "z"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_u[:nz],)),

    vdict["w"] = dict(name="w", 
                      dims=["x", "y", "z"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_w[:nz],)),

    vdict["θ"] = dict(name="θ", 
                      dims=["x", "y", "z"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_u[:nz],)),

    vdict["b"] = dict(name="b", 
                      dims=["x", "y", "z"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_u[:nz],)),

    vdict["c"] = dict(name="c", 
                      dims=["x", "y", "z", "index"],
                      coords=dict(x=sim.domain.x,
                                  y=sim.domain.y,
                                  z=sim.domain.z_u[:nz]
                                  index=np.arange(sim.n_con),)),


    return vdict
