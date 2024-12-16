# Cloudmr-tools

```python
from cmrtools.cm2D import cm2DReconB1,cm2DReconRSS,cm2DReconSENSE,cm2DGFactorSENSE

#S= your multi coil K-Space 2D signal
#N=your multi coil K-Space 2D noise or Noise covariance

L0=cm2DReconRSS()
L0.setSignalKSpace(S)
L0.setNoiseKSpace(N)
plt.figure()
plt.imshow(np.abs(L0.getOutput()))
plt.colorbar()
plt.title('RSS Reconstruction')

```
# Installation
```
#create an environment 
python3 -m venv CMT
source CMT/bin/activate
pip install git+https://github.com/cloudmrhub/cloudmr-tools.git
```
# Live Example

https://colab.research.google.com/drive/1WIEwRrNy9rpo_2X_zVdwRHRtH4ImPaan?usp=sharing

# Cite Us
- Montin E, Lattanzi R. Seeking a Widely Adoptable Practical Standard to Estimate Signal-to-Noise Ratio in Magnetic Resonance Imaging for Multiple-Coil Reconstructions. J Magn Reson Imaging. 2021 Dec;54(6):1952-1964. doi: 10.1002/jmri.27816. Epub 2021 Jul 4. PMID: 34219312; PMCID: PMC8633048.



[*Dr. Eros Montin, PhD*]
(http://me.biodimensional.com)\
