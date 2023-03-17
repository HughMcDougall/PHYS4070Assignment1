import numpy as np
import matplotlib.pylab as plt

A = np.loadtxt("./waves.txt")
X = A[:,0]

for i in range(1,4):
    plt.plot(X, A[:,i],label="wave %i" %i)
plt.xlim([0,20])
plt.axhline(0)
plt.legend(loc='best')
plt.show()
