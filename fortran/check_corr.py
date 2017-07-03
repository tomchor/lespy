import matplotlib.pyplot as plt
import numpy as np
import numfor as nf
from time import time

if 0:
    a = np.arange(3**2*3*2).reshape(2,3,3,-1)
    print(a.mean(axis=(2,3)))
    corr = nf.correlate4d(a, 1, 1, 1, 1)


[X, Y] = np.meshgrid(2 * np.pi * np.arange(100) / 12,
                             2 * np.pi * np.arange(100) / 34)
S = np.sin(X) + np.cos(Y) + np.random.uniform(0, 1, X.shape)
S = np.stack([S]*500)
S = S.reshape(5,100,100,100)
print(S.shape)

#nf.correlate4d(
corr = nf.correlate4d(S, 50, 50, 5, 100)
corr = corr[0,0]
print(corr.shape)

plt.pcolormesh(corr)
plt.show()
