import numpy as np
import matplotlib.pylab as plt

POTENTIALS = np.loadtxt("./potential.txt")
R = POTENTIALS[:,0]

#Plot Potentials
fig,ax = plt.subplots(2,1,sharex=True)

ax[0].plot(R, POTENTIALS[:,1],label="l=0 Raw Z,L Potential")
ax[0].plot(R, POTENTIALS[:,3],label="l=0 Potential /w Greens Function")
ax[0].plot(R, POTENTIALS[:,5],label="l=0 Potential /w Full Hartree")

ax[1].plot(R, POTENTIALS[:,2],label="l=1 Raw Z,L Potential")
ax[1].plot(R, POTENTIALS[:,4],label="l=1 Potential /w Greens Function")
ax[1].plot(R, POTENTIALS[:,6],label="l=1 Potential /w Full Hartree")

ax[0].set_xlim([0,10])
ax[0].set_ylim([-2,2 ])
ax[1].set_ylim([-2,2 ])

ax[0].legend(loc='best')
ax[1].legend(loc='best')
ax[0].set_ylabel("$- \frac{Z}/{r} + \frac{l(l+1)}{2/r^2}$")
ax[1].set_ylabel("$- \frac{Z}/{r} + \frac{l(l+1)}{2/r^2}$")
ax[1].set_xlabel("Radius (au)")
ax[0].axhline(0,c='k',ls='--')
ax[1].axhline(0,c='k',ls='--')

fig.tight_layout()
plt.savefig("Potentials.png",fmt="png")
