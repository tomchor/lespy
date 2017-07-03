#!/usr/bin/python
import lespy as lp
import matplotlib.pyplot as plt
import numpy as np
sim_s=lp.Simulation('/data/1/tomaschor/LES/homo4/sml_param.nml')
sim_b=lp.Simulation('/data/1/tomaschor/LES/homo4/big_param.nml')

us,vs,ws = lp.routines.readBinary2('/data/1/tomaschor/LES/homo4/restart_s/vel_tt01000000.out', simulation=sim_s)
ub,vb,wb = lp.routines.readBinary2('vel_tt01000000.out', simulation=sim_b)

Ts = lp.routines.readBinary2('/data/1/tomaschor/LES/homo4/restart_s/temp_tt01000000.out', simulation=sim_s)
Tb = lp.routines.readBinary2('temp_tt01000000.out', simulation=sim_b)

Zs = np.linspace(0,-120, 61)
Zb = np.linspace(0,-120, 121)

if 1:
    plt.plot((us*ws).mean(axis=(0,1)), Zs, label='small')
    plt.plot((ub*wb).mean(axis=(0,1)), Zb, label='large')
    plt.plot((vs*ws).mean(axis=(0,1)), Zs, label='small')
    plt.plot((vb*wb).mean(axis=(0,1)), Zb, label='large')
    plt.legend()
    plt.show()
else:
    plt.figure(0)
    plt.pcolormesh(ws[:-2,0,:])
    plt.colorbar()
    plt.figure(1)
    plt.pcolormesh(wb[:-2,0,:])
    plt.colorbar()
    plt.show()

