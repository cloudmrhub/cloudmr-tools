# Cloudmr-tools
![License](https://img.shields.io/github/license/cloudmrhub/cloudmr-tools)
![GitHub last commit](https://img.shields.io/github/last-commit/cloudmrhub/cloudmr-tools)
![GitHub issues](https://img.shields.io/github/issues/cloudmrhub/cloudmr-tools)

![GitHub forks](https://img.shields.io/github/forks/cloudmrhub/cloudmr-tools)
![GitHub stars](https://img.shields.io/github/stars/cloudmrhub/cloudmr-tools)

**Cloudmr-tools** provides tools for advanced multi-coil reconstruction methods for Magnetic Resonance Imaging (MRI). Designed for researchers and developers in the field of MRI, this package supports streamlined implementation of reconstruction techniques like RSS, SENSE, and G-Factor calculation.


## Quickstart
```python
from cmtools.cm2D import cm2DReconB1,cm2DReconRSS,cm2DReconSENSE,cm2DGFactorSENSE

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
If you use **Cloudmr-tools** in your research, please cite:

Montin E, Lattanzi R. Seeking a Widely Adoptable Practical Standard to Estimate Signal-to-Noise Ratio in Magnetic Resonance Imaging for Multiple-Coil Reconstructions. J Magn Reson Imaging. 2021 Dec;54(6):1952-1964. doi: 10.1002/jmri.27816. Epub 2021 Jul 4. PMID: 34219312; PMCID: PMC8633048.


# **Versioning**

The **Cloudmr-tools** package has two versions:

### **V1 (Deprecated)**
- **Name:** `cloudmrhub`
- **Status:** Deprecated, but still functional for backward compatibility. (v1 branch)
- **Details:** This version is no longer actively maintained and will not receive updates or bug fixes.


## **Version 2 (Current)**
- **Name:** `cloudmr-tools`
- **Status:** Actively maintained (main branch).
- **Details:** This is the recommended version for new projects. It includes updated functionality and better support for advanced features.

---

## **Key Differences**
| Feature                 | Version 1 (`cloudmrhub`)      | Version 2 (`cloudmr-tools`)  |
|-------------------------|------------------------------|-----------------------------|
| Maintenance             | Deprecated                  | Actively maintained         |
| Compatibility           | Legacy projects             | New and legacy projects     |
| Features                | Limited                     | Updated and expanded        |
---

## **Migration**

If you're currently using **Version 1** of the library, consider migrating to **Version 2** to take advantage of the latest features and updates.

If you need to continue using the Version 1 code, simply change the import path from `cloudmrhub` to `cmtools`. For example:

**Original (Version 1):**
```python
import numpy as np
import cloudmrhub.cm2D as cm2D
```
**Modified version (Version 2)**
```python
import numpy as np
import cmtools.cm2D as cm2D
```

# Contributors
[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)\
[![GitHub](https://img.shields.io/badge/GitHub-erosmontin-blue)](https://github.com/erosmontin)\
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--1773--0064-green)](https://orcid.org/0000-0002-1773-0064)\
[![Scopus](https://img.shields.io/badge/Scopus-35604121500-orange)](https://www.scopus.com/authid/detail.uri?authorId=35604121500)


[*Prof. Riccardo Lattanzi, PhD*](https://med.nyu.edu/faculty/riccardo-lattanzi)\
[![GitHub](https://img.shields.io/badge/GitHub-rlattanzi-blue)](https://github.com/rlattanzi)\
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--8240--5903-green)](https://orcid.org/0000-0002-8240-5903)\
[![Scopus](https://img.shields.io/badge/Scopus-6701330033-orange)](https://www.scopus.com/authid/detail.uri?authorId=6701330033)
