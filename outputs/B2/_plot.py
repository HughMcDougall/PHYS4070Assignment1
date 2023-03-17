import numpy as np
import matplotlib.pylab as plt

A = np.loadtxt("./waves_s.txt")
B = np.loadtxt("./waves_l.txt")
X = A[:,0]

fig,ax = plt.subplots(2,1,sharex=True)
for i in range(1,4):
    ax[0].plot(X, A[:,i]**2,label="s-wave %i" %i)
    ax[1].plot(X, B[:,i]**2,label="l-wave %i" %i)

ax[0].set_xlim([0,20])

plt.legend(loc='best')
plt.show()
