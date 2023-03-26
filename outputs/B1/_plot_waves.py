import numpy as np
import matplotlib.pylab as plt

WAVES_S = np.loadtxt("./waves_s.txt")
WAVES_L = np.loadtxt("./waves_l.txt")

POTENTIALS = np.loadtxt("./potential.txt")

R = WAVES_S[:,0]

#Plot first 3 orbitals
fig,ax = plt.subplots(2,1,sharex=True)
fig.tight_layout()
for i in range(1,4):
    ax[0].plot(R, WAVES_S[:,i]**2)
    ax[1].plot(R, WAVES_L[:,i]**2)

ax[0].set_xlim([0,50])
ax[0].legend(["1s", "2s", "3s"], loc='best')
ax[1].legend(["2p", "3p", "4p"], loc='best')
ax[0].set_ylabel("$|P_{ns}(r)|^2$")
ax[1].set_ylabel("$|P_{nl}(r)|^2$")
ax[1].set_xlabel("Radius (au)")

fig.tight_layout()
plt.savefig("Waves-subplot.png",fmt="png")

#All on same plot
fig = plt.figure(figsize=(4,2))

plt.plot(R, WAVES_S[:,1]**2)
plt.plot(R, WAVES_S[:,2]**2)
plt.plot(R, WAVES_L[:,1]**2)

plt.xlim([0,10])
plt.legend(["1s", "2s", "2p"], loc='best')
plt.ylabel("$|P_{ns}(r)|^2$")
plt.xlabel("Radius (au)")

fig.tight_layout()
plt.savefig("Waves.png",fmt="png")
