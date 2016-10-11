import lespy as lp

sim_a = lp.simulation('/data/1/tomaschor/LES/LV3_test/param.nml')

avgs ='/data/1/tomaschor/LES/LV3_test/output/postprocessing/avgs.csv'
out = lp.readMeanPP(avgs)

lp.postProcessAvgs('/data/1/tomaschor/LES/LV3_test/output', simulation=sim_a, t_ini=49198)

