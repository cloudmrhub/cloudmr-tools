import pynico_eros_montin.pynico as pn

import twixtools as tt


# a=pn.Pathable('cmtools/testdata.pkl')

# A=a.readPkl()

# S=A[0]["signal"]
# N=A[0]["noise"]

import numpy as np


A=tt.map_twix('/data/MYDATA/siemensrawdataexamples/26092019_brain/26092019_meas_MID01481_FID318481_gre_snr.dat')
S=A[0]["image"]
S.flags['remove_os'] = True  # activate automatic os removal
S.flags['regrid'] = True  # activate ramp sampling regridding
S.flags['average']['Rep'] = True  # average all repetitions

S=S[:][0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:]
S=np.transpose(S,[0,2,1])


B=tt.map_twix('/data/MYDATA/siemensrawdataexamples/26092019_brain/26092019_meas_MID01493_FID318493_gre_noise.dat')
N=B[0]["image"]
N=N[:][0,0,0,0,0,0,0,0,0,0,0,0,0,:,:,:]
N=np.transpose(N,[0,2,1])


from cmtools import cm,cm2D
import numpy as np
import matplotlib.pyplot as plt





for PA in [1,2,3,4]:
    FA=1
    ACL=24

    # Filters=[None,'hanning','hamming','blackman','bartlett','kaiser']

    Filters=['blackman']

    for f in Filters:
        US,REF=cm.mimicAcceleration2D(S,[FA,PA],ACL=[np.nan,ACL],filter=f)
        L2=cm2D.cm2DKellmanSENSE()
        L2.setAcceleration([FA,PA])
        L2.setSignalKSpace(US)
        L2.setNoiseKSpace(N)
        L2.setReferenceKSpace(REF)
        L2.setAutocalibrationLines(ACL)
        L2.setMaskCoilSensitivityMatrixBasedOnEspirit(6)
        L2.setMaskCoilSensitivityMatrixDefault()
        L2.setMaskCoilSensitivityMatrixBasedOnEspirit(6)
        # plt.figure()
        # plt.imshow(np.abs(L2.getOutput()))
        # plt.colorbar()
        # if f is None:
        #     plt.title(f'SENSE Reconstruction {FA}x{PA} No Filter')
        # else:
        #     plt.title(f'SENSE Reconstruction {FA}x{PA} {f}')
        # plt.savefig(f'/g/SENSE_{FA}x{PA}_{f}.png')
        # plt.close()
        
        #make a whole window figure
        plt.figure(figsize=(20,10))
        L2.plotCoilSensitivityMatrix(f'/g/coil_{FA}x{PA}_{f}.png',title_addition=f'\n{FA}x{PA} {f}')
        plt.close()
        
        
