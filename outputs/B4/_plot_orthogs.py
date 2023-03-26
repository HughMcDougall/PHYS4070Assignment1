import numpy as np
import matplotlib.pylab as plt

WAVES_S = np.loadtxt("./waves_s.txt")
WAVES_L = np.loadtxt("./waves_l.txt")

R = WAVES_S[:,0]
dr = (R[-1]-R[0]) / (len(R)-1)
S2 = WAVES_S[:,2]
P2 = WAVES_L[:,1]

Rss = dr * np.sum(R*S2*S2)
Rpp = dr * np.sum(R*P2*P2)
Rsp = dr * np.sum(R*P2*P2)

print(Rss,Rpp,Rsp)
print(Rsp / (Rss * Rpp)**0.5 )
