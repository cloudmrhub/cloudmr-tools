import twixtools
import json
from pynico_eros_montin import pynico as pn
import numpy as np
import scipy

from cm import *
#----------------------------------------------------------------------------------------------------------
def adapt_array_2d(yn, rn=None, normalize=False):
    """Adaptive reconstruction of MRI array data.

    Args:
        yn: Array data to be combined.
        rn: Data covariance matrix.
        normalize: Normalize image intensity (False)

    Returns:
        recon: Reconstructed image.
        cmap: Estimated coil sensitivity maps.
        wfull: Coefficients to combine coil signals.
    """

    yn = np.transpose(yn, [2,0,1])
    nc, ny, nx = yn.shape
    if rn is None:
        rn = np.eye(nc)

    maxcoil = np.argmax(np.sum(np.sum(np.abs(yn), axis=1),axis=1))
    bs1 = 2
    bs2 = 2
    st = 1

    wsmall = np.zeros((nc, round(ny / st), round(nx / st)),dtype=np.complex64)
    cmapsmall = np.zeros((nc, round(ny / st), round(nx / st)),dtype=np.complex64)

    for x in np.arange(0, nx, st):
        for y in np.arange(0, ny, st):
            ymin1 = max([y - bs1 // 2, 0])
            xmin1 = max([x - bs2 // 2, 0])
            ymax1 = min([(y + bs1 // 2)+1, ny])
            xmax1 = min([(x + bs2 // 2)+1, nx])

            ly1 = len(np.arange(ymin1,ymax1))
            lx1 = len(np.arange(xmin1,xmax1))

            m1 = yn[:, ymin1:ymax1, xmin1:xmax1].reshape(nc, lx1 * ly1)

            m = m1 @ m1.conj().T

            v, e = np.linalg.eigh(np.linalg.inv(rn) @ m, UPLO='U')

            
            ind = np.argmax(v)

            mf = e[:, ind]

            mf = np.divide(mf, mf.conj().T @ np.linalg.inv(rn) @ mf)
            normmf = e[:, ind]

            mf = np.dot(mf,np.exp(-1j * np.angle(mf[maxcoil])))
            normmf = np.dot(normmf, np.exp(-1j * np.angle(normmf[maxcoil])))

            wsmall[:, y, x] = mf
            cmapsmall[:, y, x] = normmf

    recon = np.zeros((ny, nx),dtype=complex)
    wfull = np.zeros((nc,ny,nx),dtype=complex)
    
    # for i in range(nc):
    #     wfull[i, :, :] = 
    #         np.array(
                
    #             scipy.misc.imresize(
    #                 np.abs(wsmall[i, :, :]), (ny, nx), order=1, mode="nearest"
    #             )
    #         )
    #         * np.exp(1j * np.array(
    #             scipy.ndimage.interpolation.zoom(
    #                 np.angle(wsmall[i, :, :]), (ny, nx), order=0, mode="nearest"
    #             )
    #         )
    #     )
    #     recon += wfull[i, :, :] * yn[:, :, i]

    # if normalize:
    #     recon *= np.sum(np.abs(cmap)) ** 2

    # cmap = np.transpose(cmap,[2,1,0])
    # return cmap


# def calculate_simple_sense_sensitivitymaps_acl(K,autocalibrationF,autocalibrationP,mask=None):
#     """Calculates the coil sensitivity maps using the simple SENSE method on undersampled Kspace.

#     Args:
#         K: freq,phase,coil numpy array of the kspace 


#     Returns:
#         The coil sensitivity matrix.
#     """
#     return calculate_simple_sense_sensitivitymaps(getAutocalibrationsLines2DKSpaceZeroPadded(K,autocalibrationF,autocalibrationP),mask)

