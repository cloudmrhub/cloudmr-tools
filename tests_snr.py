import pynico_eros_montin.pynico as pn

import twixtools as tt


# a=pn.Pathable('cmtools/testdata.pkl')

# A=a.readPkl()

# S=A[0]["signal"]
# N=A[0]["noise"]

import numpy as np

S=np.load('/data/PROJECTS/fsqu/DATA/15042019_MR/KS.npz')['signal']  # (r,Nc, Nx, Ny)
N=np.load('/data/PROJECTS/fsqu/DATA/15042019_MR/KN.npz')['noise']  # ( Nc, Nx, Ny)


def estimate_noise_cov(no_rf_complex: np.ndarray, ncoils: int) -> np.ndarray:
    """no_rf_complex: array with coil as last axis or reshape-able to (-1, ncoils)."""
    X = no_rf_complex.reshape(-1, ncoils)          # (Nsamples_noise, Ncoils)
    A = X.T                                        # (Ncoils, Nsamples_noise)
    Sigma = (A @ A.conj().T) / A.shape[1]
    Sigma = 0.5 * (Sigma + Sigma.conj().T)
    # regularize if needed
    lam = 1e-8 * np.trace(Sigma) / Sigma.shape[0]
    Sigma = Sigma + lam * np.eye(Sigma.shape[0], dtype=Sigma.dtype)
    return Sigma


from cmtools import cm,cm2D
import numpy as np
import matplotlib.pyplot as plt


S=S[:,:, :, 0:1]  # (R, Nx, Ny, Nc)
N=N[ :, :, 0:1]  # ( Nc, Nx, Ny)
R=cm2D.cm2DReconRSS()
R.setNoiseKSpace(N)
MR=cm2D.cm2DSignalToNoiseRatioMultipleReplicas()
MR.setReconstructor(R)
for a in range(S.shape[0]):
    MR.add2DKspace(S[a])
O=MR.getOutput()


plt.subplot(3,3,1)
plt.imshow(np.abs(O),cmap='gray')
plt.title(f'MREP SNR Map over Replicas {np.max(np.abs(O)):.1f}')
plt.colorbar()
A=cm2D.cm2DKellmanRSS()
A.setNoiseKSpace(N)
A.setSignalKSpace(S[0])
O=A.getOutput()

plt.subplot(3,3,2)
plt.imshow(np.abs(O),cmap='gray')
plt.colorbar()
plt.title(f'AC SNR Map over 1 Replicas {np.max(np.abs(O)):.1f}')



plt.subplot(3,3,3)
A=cm2D.cm2DKellmanRSS()
A.setNoiseKSpace(N)
A.setSignalKSpace(np.mean(S,axis=0))
O=A.getOutput()

plt.imshow(np.abs(O),cmap='gray')
plt.colorbar()
plt.title(f'AC SNR Map over avg Replicas {np.max(np.abs(O)):.1f}')




# plt.subplot(3,3,4)
# A=cm2D.cm2DKellmanRSS()
# A.setNoiseCovariance(estimate_noise_cov(np.mean(S,0), S.shape[1]))
# A.setSignalKSpace(np.mean(S,axis=0))
# O=A.getOutput()

# plt.imshow(np.abs(O),cmap='gray')
# plt.colorbar()
# plt.title(f'AC SNR Map over avg Replicas {np.max(np.abs(O)):.1f}')



for a in range(5,35,5):
    A=cm2D.cm2DKellmanRSS()
    A.setNoiseKSpace(N)
    A.setSignalKSpace(np.mean(S[:a],axis=0))
    O=A.getOutput()
    plt.subplot(3,3,3+a//5)
    plt.imshow(np.abs(O),cmap='gray')
    plt.colorbar()
    plt.title(f'AC SNR Map over avg({a}) Replicas {np.mean(np.abs(O)):.1f}')






def concat_replicas_along_freq(S, replica_axis=0, freq_axis=1, phase_axis=2, coil_axis=3):
    """
    Return array with shape (n_phase, nreplica * n_freq, n_coils)

    S: 4-D k-space with axes (replica, freq, phase, coil) in any order;
       default assumed (R, Nx, Ny, Nc).
    replica_axis, freq_axis, phase_axis, coil_axis: indices of those axes in S.
    """
    # move axes to (replica, freq, phase, coil)
    S4 = np.moveaxis(S, (replica_axis, freq_axis, phase_axis, coil_axis), (0, 1, 2, 3))
    R, Nx, Ny, Nc = S4.shape
    # (phase, replica, freq, coil)
    S_t = S4.transpose(2, 0, 1, 3)
    # reshape to (phase, replica*freq, coil)
    out = S_t.reshape(Ny, R * Nx, Nc)
    return out

A=cm2D.cm2DKellmanRSS()

meanS = np.mean(S, axis=0, keepdims=True)   # shape (1, Nx, Ny, Nc)
residuals = S - meanS  


A.setSignalKSpace(np.mean(S,axis=0))
A.setNoiseKSpace(concat_replicas_along_freq(residuals))
O=A.getOutput()
plt.figure()
plt.imshow(np.abs(O),cmap='gray')
plt.colorbar()
plt.title(f'AC SNR Map over avg Replicas with estimated noise covariance using S-mean(S) {np.max(np.abs(O)):.1f}')



kmean = np.mean(S, axis=0)     # (Nx, Ny, Nc) after fix above
img = np.fft.ifftshift(np.fft.ifft2(kmean, axes=(0,1))  )
img_rss = np.sqrt(np.sum(np.abs(img)**2, axis=-1))   # (Nx, Ny)
maskPlot = img_rss > (0.1 * np.max(img_rss))   # boolean mask for low-signal areas
outside_coords = np.nonzero(~maskPlot)          # tuple (xs, ys)
n_loc = len(outside_coords[0])


OK_img = np.fft.ifft2(S, axes=(1, 2))
K_out_img = OK_img[:, outside_coords[0], outside_coords[1], :]
# mean across replicas for each location -> shape (n_loc, Nc)
Kmean_loc = K_out_img.mean(axis=0)

# residuals: noise-only approximation in image space (R, n_loc, Nc)
residuals = K_out_img - Kmean_loc[None, ...]
# flatten samples -> shape (n_samples, Nc) where n_samples = R * n_loc
res_flat = residuals.reshape(-1, residuals.shape[-1])
n_samples = res_flat.shape[0]

# covariance (Nc x Nc) unbiased estimator from image-space residuals
noise_cov_est = (res_flat.conj().T @ res_flat) / (n_samples - 1)
# enforce Hermitian symmetry and small regularization for stability
noise_cov_est = 0.5 * (noise_cov_est + noise_cov_est.conj().T)
eps = 1e-6 * np.trace(noise_cov_est) / noise_cov_est.shape[0]
noise_cov_est += np.eye(noise_cov_est.shape[0]) * eps

        
A.setSignalKSpace(np.mean(S,axis=0))
A.setNoiseKSpace(noise_cov_est)
O=A.getOutput()
plt.figure()
plt.imshow(np.abs(O),cmap='gray')
plt.colorbar()
plt.title(f'AM noise_cov_set estimated {np.max(np.abs(O)):.1f}')



def SNRIM(K,mask):
    outside_coords = np.nonzero(~mask)          # tuple (xs, ys)
    K=img_rss[outside_coords[0], outside_coords[1]]  # (n_loc, Nc)
    img_rss/=np.abs(np.std(K))

plt.figure()


plt.imshow(np.abs(img_rss),cmap='gray')
plt.colorbar()
plt.title(f'AC im/std(N) {np.max(np.abs(img_rss)):.1f}')



plt.show()







class SNRIMsingle():
    def __init__(self,kspace,mask=None):
        self.kspace = kspace
        self.NR=kspace.shape[0]
        self.mask = mask
    def setMask(self,mask):
        self.mask = mask
    def getOutputs(self):
        K=np.mean(self.kspace,axis=0)  # (Nx, Ny, Nc)
        K_img = np.fft.ifft2(K, axes=(0, 1))
        outside_coords = np.nonzero(self.mask)          # tuple (xs, ys)
 
        