import pynico_eros_montin.pynico as pn

import twixtools as tt


# a=pn.Pathable('cmtools/testdata.pkl')

# A=a.readPkl()

# S=A[0]["signal"]
# N=A[0]["noise"]

import numpy as np

S=np.load('/data/PROJECTS/fsqu/DATA/15042019_MR/KS.npz')['signal']  # (r,Nx, Ny, nc)
N=np.load('/data/PROJECTS/fsqu/DATA/15042019_MR/KN.npz')['noise']  # ( Nc, Nx, Ny)

from cmtools import cm,cm2D
import numpy as np
import matplotlib.pyplot as plt


x,y=4,4

NU=1000

L=[]
M=[]
M2=[]
mask=None
for a in range(1,(x*y)+1):
    K=np.mean(S[:a],axis=0)
    K=np.fft.ifft2(K, axes=(0,1))
    SC=np.sqrt(np.prod(np.array(K.shape[0:2])))
    K=K*SC
    K*=1000
    O=np.sqrt(np.sum(np.abs(K), axis=-1))
    O=np.fft.ifftshift(O)# (Nx, Ny)
    if mask is None:
        mask=np.abs(O)>np.percentile(np.abs(O),40)
        mask=mask.astype(np.float64)
        mask[mask==0]=np.nan
    O=O*mask
    o=np.nanmean(np.abs(O))
    plt.subplot(x,y,a)
    plt.imshow(np.abs(O),cmap='gray')
    plt.colorbar()
    plt.tight_layout()
    plt.title(f'NAvg({a}) {o:.3f}')
    
plt.suptitle(f'Reconstructed Images vs Number of Averages max multipied by{NU}', fontsize=20)
S1=[]
INI=None
for a in range(1,S.shape[0]):
    K=np.mean(S[:a],axis=0)
    K=np.fft.ifft2(K, axes=(0,1))
    SC=np.sqrt(np.prod(np.array(K.shape[0:2])))
    K=K*SC
    O=np.sqrt(np.sum(np.abs(K)**2, axis=-1))
    O=np.fft.ifftshift(O)# (Nx, Ny)
    NON=np.abs(O[~np.isnan(mask)])
    ON=np.abs(O[np.isnan(mask)])
    
    o=NON.max()
    print(o)
    o2=np.mean(ON)
    o3=np.std(ON)
    L.append(o)
    M.append(o2)
    M2.append(o3)
    if INI is None:
        INI=o2
        print(o2)
    print(a,INI,np.sqrt(a),INI/np.sqrt(a))
    S1.append(INI/np.sqrt(a))
    


plt.figure(figsize=(15,15))
plt.plot(L)
plt.plot(S1,'*')
plt.title('Max of Reconstructed Image vs Number of Averages', fontsize=20)


plt.figure(figsize=(15,15))

plt.plot(M)
plt.plot(S1,'*')
plt.title('Mean of Reconstructed Image vs Number of Averages', fontsize=20)

plt.show()
plt.clf()