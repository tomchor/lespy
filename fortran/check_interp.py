#!/usr/bin/python
import lespy as lp
import matplotlib.pyplot as plt
sim_s=lp.Simulation('/data/1/tomaschor/LES/homo4/sml_param.nml')
sim_b=lp.Simulation('/data/1/tomaschor/LES/homo4/big_param.nml')

us,vs,ws = lp.routines.readBinary2('/data/1/tomaschor/LES/homo4/restart/vel_tt00002000.out', simulation=sim_s)
ub,vb,wb = lp.routines.readBinary2('vel_tt00002000.out', simulation=sim_b)

Ts = lp.routines.readBinary2('/data/1/tomaschor/LES/homo4/restart/temp_tt00002000.out', simulation=sim_s)
Tb = lp.routines.readBinary2('temp_tt00002000.out', simulation=sim_b)

plt.figure(0)
plt.pcolormesh(ws[:-2,0,:])
plt.figure(1)
plt.pcolormesh(wb[:-2,0,:])
plt.show()

