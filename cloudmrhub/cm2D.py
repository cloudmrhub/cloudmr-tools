import numpy as np
try:
    import cm
except:
    import cloudmrhub.cm as cm
    
import matplotlib.pyplot as plt
import scipy


class cm2DRecon(cm.cmOutput):
    """
    Python implementation of the cm2DRecon MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.


    :Attributes:
        SignalKSpace (cm.k2d): The signal k-space data.
        
        NoiseKSpace (cm.k2d): The noise k-space data.
                
        NoiseCovariance (np.ndarray): The noise covariance matrix.
        
        InverseNoiseCovariance (np.ndarray): The Inverse Noise Covariance matrix.
        
        SignalPrewhitened (np.ndarray): The prewhitened signal k-space data.
        
        HasSensitivity (bool): Whether the reconstruction contains sensitivity information.
        
        HasAcceleration (bool): Whether the reconstruction has accelerations.
   
    """
    def __init__(self):
        """
        Initializes the cm2DRecon object.
        """
        self.SignalKSpace = cm.k2d()        
        self.NoiseKSpace = cm.k2d()
        self.NoiseCovariance = np.array([])
        self.InverseNoiseCovariance = np.array([])
        self.SignalPrewhitened = cm.k2d()
        self.HasSensitivity = False
        self.HasAcceleration = False
        self.NoiseBandWidth = None

    def setSignalKSpace(self, signalKSpace):
        """
        Sets the signal k-space data.
        
        :param signalkspace: The signal k-space data
        :type: np.ndarray
        """
        
        self.SignalKSpace.set(signalKSpace)

    def getSignalKSpace(self):
        """
        Gets the signal k-space data.
        
        Returns:
            _type_: nd.array(f,p,c)
        """
        return self.SignalKSpace.get()
    def getSignalKSpaceSize(self):
        """
        Gets the signal k-space size.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.SignalKSpace.getSize()
    def getsignalNCoils(self):
        """
        Gets the number of coils in the signal k-space data.
        
        Returns:
            int: ncoils
        """
        return self.SignalKSpace.getNCoils()
    def setNoiseKSpace(self, noiseKSpace):
        """ Sets the noise k-space data.

        :param noiseKspace: The noise k-space data
        :type: np.ndarray(f,p,c)
        """
        self.NoiseKSpace.set(noiseKSpace)
    def getNoiseKSpace(self):
        """
        Gets the noise k-space data.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.NoiseKSpace.get()
    def getNoiseKSpaceSize(self):
        """
        Gets the noise k-space size.
        
        Returns:
            NoiseKspace: nd.array(f,p,c)
        """
        return self.NoiseKSpace.getSize()
    def getNoiseNCoils(self):
        """
        Gets the number of coils in the noise k-space data.
        
        Returns:
            int: ncoils
        """
        return self.NoiseKSpace.getNCoils()
    
    def getNoiseBandWidth(self):
            if(self.NoiseBandWidth is None):
                noise_bandwidth = cm.mrir_noise_bandwidth(self.getNoiseKSpace());
                self.NoiseBandWidth=noise_bandwidth
            try:
                return  self.NoiseBandWidth
            except:
                print("---problem in the noisebadwidth----")
                self.appendLog("---problem in the noisebadwidth----","warning")
                return  1.0
            
    def getInverseNoiseCovariancePrewhitened(self):
        """when prewhitened the correlation is an eye

        """
        return np.eye(self.getsignalNCoils())
    def setNoiseCovariance(self, noiseCovariance):
        """
        Sets the noise covariance matrix.

        :param noiseCovariance: The noise covariance matrix.
        :type: mp.ndarray(c,c)

        """
        self.NoiseCovariance = noiseCovariance
        self.InverseNoiseCovariance = np.linalg.inv(noiseCovariance)
    def getNoiseCovariance(self):
        """
        Return the covariance matrix
        
        Returns:
            np.ndarray(c,c): Covariance Matrix
        """
        if not self.NoiseCovariance.any():
            self.NoiseCovariance = cm.calculate_covariance_matrix(self.getNoiseKSpace(),self.getNoiseBandWidth())

        return self.NoiseCovariance
    def getInverseNoiseCovariance(self):
        """
        Gets the inverse noise covariance matrix.

        Returns:
            np.ndarray(c,c): inverse of the covariance matrix
        """
        return self.InverseNoiseCovariance

    def getPrewhitenedSignal(self):
        """
        Gets the prewhitened signal.

        Returns:
            np.ndarray(f,p,c): prewhitened signal
        """
        if self.SignalPrewhitened.isEmpty():
            self.SignalPrewhitened.set(cm.prewhiteningSignal(self.getSignalKSpace(), self.getNoiseCovariance()))
        return self.SignalPrewhitened.get()
    def setPrewhitenedSignal(self, prewhitenedSignal):
        self.SignalPrewhitened.set(prewhitenedSignal)
    def plotImageAfterTest(self,IM,tit):
        # Create the figure and subplots
        fig, axarr = plt.subplots(2, 1)
        axarr[0].imshow(IM)
        axarr[0].set_title(tit)
        axarr[1].imshow(abs(IM))
        axarr[1].set_title(tit + ' abs')

        ha = plt.axes([0, 0, 1, 1], frameon=False, visible=False,
              xlim=[0, 1], ylim=[0, 1], aspect='equal')
        # Add text to the axes object
        text = ha.text(0.5, 0.98, f'{IM.shape[0]}x{IM.shape[1]}',
                    transform=ha.transAxes, horizontalalignment='center',
                    verticalalignment='top')
        plt.show()
    def test(self):
        TEST = self.testrecon()
        self.plotImageAfterTest(TEST.getOutput(), "recon")
        return TEST
    def get2DKSIFFT(self):
        SC=np.sqrt(np.prod(np.array(self.getSignalKSpaceSize())))
        return cm.MRifft(self.getPrewhitenedSignal(),[0,1])*SC

class cm2DReconRSS(cm2DRecon):
    """
    Python implementation of the cm2DReconRSS MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=False
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        im = np.sqrt(np.sum(np.abs(img_matrix)**2,axis=-1))
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconRSS()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST

class cm2DReconSS(cm2DReconRSS):
    """
    Python implementation of the cm2DReconSS MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        im = np.sum(img_matrix**2,axis=-1)
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconSS()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST
class cm2DReconRSSunAbs(cm2DReconSS):
    """
    Python implementation of the cm2DReconRSSunAbs MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2023
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    def getOutput(self):
        im = np.sqrt(super().getOutput())
        return im

    @staticmethod
    def testrecon():
        TEST = cm2DReconRSSunAbs()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST


class cm2DKellmanRSS(cm2DReconRSS):
    def __init__(self):
        super().__init__()
    def getOutput(self):
        nf,nph =self.getSignalKSpaceSize()
        img_matrix = self.get2DKSIFFT()
        snr = np.zeros((nf,nph))
        for irow in range(nf):
            for icol in range(nph):
                        B=np.expand_dims(img_matrix[irow,icol],axis=-1)
                        A=B.conj().T
                        snr[irow,icol] = np.abs(np.sqrt(2*(A @ B)))
        return snr
    


class cm2DReconWithSensitivity(cm2DRecon):
    """
    Python implementation of the cm2DReconWithSensitivity MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=True
        self.CoilSensitivityMatrixCalculationMethod='inner'
        self.CoilSensitivityMatrix=cm.k2d()
        self.CoilSensitivityMatrixSource=cm.k2d()
        self.CoilSensitivityMatrixSourcePrewhitened=cm.k2d()
        self.MaskCoilSensitivityMatrix='ref'

    
    def setCoilSensitivityMatrix(self, S):
        self.CoilSensitivityMatrix.set(S)

    def resetCoilSensitivityMatrix(self):
        self.setSensitivityMatrix.reset()

    def calculateCoilSensitivityMatrix(self):
        if self.CoilSensitivityMatrix.isEmpty():
            if ((self.getCoilSensitivityMatrixCalculationMethod() == 'simplesense') or (self.getCoilSensitivityMatrixCalculationMethod() =='inner')):
                s=self.getCoilSensitivityMatrixSimpleSense()
        else:
            s=self.getCoilSensitivityMatrix()
        return s

    def getCoilSensitivityMatrix(self):
        if self.CoilSensitivityMatrix.isEmpty():
            coilsens_set = self.calculateCoilSensitivityMatrix()
            self.setCoilSensitivityMatrix(coilsens_set)
        return self.CoilSensitivityMatrix.get()

    def setCoilSensitivityMatrixSourcePrewhitened(self, x):
        self.CoilSensitivityMatrixSourcePrewhitened.set(x)

    def getCoilSensitivityMatrixSource(self):
        return self.CoilSensitivityMatrixSource.get()
    
    def getCoilSensitivityMatrixSourcePrewhitened(self):
        if self.CoilSensitivityMatrixSourcePrewhitened.isEmpty():
            S = self.getCoilSensitivityMatrixSource()
            Rn = self.getNoiseCovariance()
            pw_S = cm.prewhiteningSignal(S, Rn)
            self.setCoilSensitivityMatrixSourcePrewhitened(pw_S)
        return pw_S

    def setCoilSensitivityMatrixSource(self, IM):
        self.CoilSensitivityMatrixSource.set(IM)

    def setMaskCoilSensitivityMatrix(self, x):
        self.MaskCoilSensitivityMatrix = x

    def getMaskCoilSensitivityMatrix(self):
        return self.MaskCoilSensitivityMatrix

    def setCoilSensitivityMatrixCalculationMethod(self, x):
        self.CoilSensitivityMatrixCalculationMethod = x

    def getCoilSensitivityMatrixCalculationMethod(self):
        return self.CoilSensitivityMatrixCalculationMethod

    def getCoilSensitivityMatrixSimpleSense(self):
        # MASK 
        if self.CoilSensitivityMatrix.isEmpty():       
            if self.CoilSensitivityMatrixSource.isEmpty():
                s=self.getSignalKSpace()
            else:
                s=self.getCoilSensitivityMatrixSource()

            self.setCoilSensitivityMatrix(cm.prewhiteningSignal(cm.calculate_simple_sense_sensitivitymaps(s,self.MaskCoilSensitivityMatrix), self.getNoiseCovariance() ))
        return self.CoilSensitivityMatrix.get()


class cm2DReconB1(cm2DReconWithSensitivity):
    """
    Python implementation of the cm2DReconB1 MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=False
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        pw_sensmap=self.getCoilSensitivityMatrix()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        nf,nph =self.getSignalKSpaceSize()
        im = np.zeros((nf,nph))
        for irow in range(nf):
            for icol in range(nph):
                        s_matrix=pw_sensmap[irow,icol,:]
                        if s_matrix.sum() !=0:
                            im[irow,icol] = s_matrix.conj().T @ invRn @ img_matrix[irow,icol,:]
        return im
    

    @staticmethod
    def testrecon():
        TEST = cm2DReconB1()
        [K, N, nc] = cm.getMRoptimumTestData()
        TEST.setSignalKSpace(K)
        TEST.setNoiseCovariance(nc)
        return TEST

class cm2DKellmanB1(cm2DReconB1):
    """
    Python implementation of the cm2DReconB1 MATLAB class
   
    :author:
        Dr. Eros Montin, Ph.D. <eros.montin@gmail.com>
    :date:
        16/06/2003
    :note:
        This work was supported in part by the National Institute of Biomedical Imaging and Bioengineering (NIBIB) of the National Institutes of Health under Award Number R01 EB024536 and P41 EB017183. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.
    """
    
    def __init__(self):
        """
        Initializes the RSS reconstruction.
        
        """
        super().__init__()
        self.HasAcceleration= False
        self.HasSensitivity=False
    def getOutput(self):
        img_matrix = self.get2DKSIFFT()
        pw_sensmap=self.getCoilSensitivityMatrix()
        invRn=self.getInverseNoiseCovariancePrewhitened()
        nf,nph =self.getSignalKSpaceSize()
        im = np.zeros((nf,nph))
        SR=np.sqrt(2.0)
        for irow in range(nf):
            for icol in range(nph):
                        s_matrix=pw_sensmap[irow,icol,:]
                        if s_matrix.sum() !=0:
                            S=s_matrix
                            ST=s_matrix.conj().T
                            I=img_matrix[irow,icol,:]
                            num=np.dot(SR,np.abs(ST @ invRn @ I))
                            den=np.sqrt(np.abs(ST @ invRn @ S))
                            im[irow,icol] = np.divide(num,den)
        return im
    