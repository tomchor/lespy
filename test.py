import lespy as lp
import numpy as np
import pandas as pd
pd.options.display.max_rows=10

sim_a = lp.simulation('/data/1/tomaschor/LES/LV3_test/output/codebkp/param.nml')
print(sim_a.timelength)
exit()
sim_a=timelength=720000
print(sim_a.timelength)

avgs ='/data/1/tomaschor/LES/LV3_test/output/postprocessing/avgs.csv'
out = lp.readMeanPP(avgs)

arr = lp.postProcess2D('/data/1/tomaschor/LES/LV3_test/output', simulation=sim_a, t_ini=50000, t_end=100000)
a = pd.DataFrame(arr.T)

