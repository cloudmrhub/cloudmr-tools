
from cmtools.cm2D import cm2DKellmanRSS, cm2DReconRSS, cm2DKellmanB1,cm2DReconB1,cm2DReconSENSE,cm2DKellmanSENSE,cm2DGFactorSENSE,cm2DReconGRAPPA
import cmtools.cm as cm

import numpy as np
import matplotlib.pyplot as plt




S,N,nc=cm.getMRoptimumTestData(figure='cloudmr-tools/eros.jpg')

L=cm2DReconB1()
L.setSignalKSpace(S)
L.setNoiseKSpace(N)
L.setReferenceKSpace(S)
L.prepareCoilSensitivityMatrixPlot(title='B1 Coil Sensitivity Maps')
plt.figure()
plt.imshow(np.abs(L.getOutput()))
plt.title('B1 Recon')
cm.saveMatlab('/g/B1.mat',[{"name":'B1',"data":L.getOutput()}])
plt.show()
