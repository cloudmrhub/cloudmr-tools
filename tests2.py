import pynico_eros_montin.pynico as pn

import twixtools as tt


# a=pn.Pathable('cmtools/testdata.pkl')

# A=a.readPkl()

# S=A[0]["signal"]
# N=A[0]["noise"]

import numpy as np


A=tt.map_twix('/data/MYDATA/siemensrawdataexamples/26092019_brain/26092019_meas_MID01481_FID318481_gre_snr.dat')
S=A[0]["image"]
S=S[:][0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:]
S=np.transpose(S,[0,2,1])


B=tt.map_twix('/data/MYDATA/siemensrawdataexamples/26092019_brain/26092019_meas_MID01493_FID318493_gre_noise.dat')
N=B[0]["image"]
N=N[:][0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:]
N=np.transpose(N,[0,2,1])


from cmtools import cm,cm2D
import numpy as np
import matplotlib.pyplot as plt
nx,ny,nz=S.shape


S=S[nx//2,:,1]
W=np.zeros_like(S)
O=np.zeros_like(S)
ACL=48
window=np.hamming(ACL)
s=ny//2-ACL//2


W[s:s+ACL]=window
O[s:s+ACL]=1

WS=(S*(W/np.sum(W)))
OS=(S*(O/np.sum(O)))

plt.plot(W/np.max(W),'g')
plt.plot(np.abs(np.abs(S)/np.max(np.abs(S))),'b')
plt.plot(O,'k')
plt.plot(np.abs(WS)/np.max(np.abs(WS)),'r')
plt.xlim([ny//4,ny//4*3])
plt.grid(True)
plt.legend(['Hamming Window','Signal','Window','Hamming Windowed Signal'])
plt.savefig('/g/signal_windowed.png')
plt.close()

plt.plot(np.abs(np.abs(S)/np.max(np.abs(S))),'b')
plt.xlim([ny//4,ny//4*3])
plt.legend(['Signal'])
plt.savefig('/g/signal.png')
plt.close()

plt.plot(np.abs(np.abs(S)/np.max(np.abs(S))),'b')
plt.plot(O,'k')
plt.legend(['Signal','Window'])
plt.xlim([ny//4,ny//4*3])
plt.grid(True)
plt.savefig('/g/signal_trunk.png')
plt.close()

plt.plot(np.abs(np.abs(S)/np.max(np.abs(S))),'b')
plt.plot(O,'k')
plt.grid(True)
plt.plot(W/np.max(W),'g')
plt.xlim([ny//4,ny//4*3])
plt.legend(['Signal','Window','hamming Window'])
plt.savefig('/g/signal_trunk_window_.png')

plt.close()




plt.plot(O,'k')
plt.plot(np.abs(OS)/np.max(np.abs(OS)),'g')
plt.xlim([ny//4,ny//4*3])
plt.grid(True)
plt.legend(['Window','Windowed Signal'])
plt.savefig('/g/original_signal.png')
plt.close()


plt.plot(W/np.max(W),'k')
plt.plot(np.abs(WS)/np.max(np.abs(WS)),'g')
plt.xlim([ny//4,ny//4*3])
plt.legend(['Hamming Window','Windowed Signal'])
plt.grid(True)
plt.savefig('/g/windowed_signals.png')

plt.close()




