
import numpy as np
from cm2D import cm2DReconWithSensitivityAutocalibrated, cm2DReconRSS,cm2DReconWithSensitivity,cm2DReconB1

# import devcm as cm
import matplotlib.pyplot as plt
import devcm as cm
from  pyable_eros_montin import imaginable as ima
from cm import k2d,i3d

class f3dmc(k2d):
    """Field 3D multicoil

    Args:
        k2d (_type_): _description_
    """
    def __init__(self, k=None):
        super().__init__(k)
    def getSize(self):
        return self.k.shape[:-2]
    def getNumberOfElementPerVoxel(self):
        return self.k.shape[-2]

class t3d(f3dmc):
    """tensor 3D
    [x,y,z,v]
    """
    def __init__(self, k=None):
        super().__init__(k)
    def getSize(self):
        return self.k.shape[:-1]
    def getNumberOfElementPerVoxel(self):
        return self.k.shape[:-1]
    

class Field(f3dmc):
    """Class Field
        the 5d image of a field [x,y,z,v,C]
    """
    def __init__(self, k=None):
        super().__init__(k)
      
    

# class HField(Field):
#     def __init__(self,F=None):
#         super().__init__(F)
#         self.B1plus =t3d()
#         self.B1minus =t3d()

#     def reset(self): #override
#         #we want to recalculate the minus and plus in case H was changed 
#         self.B1plus =t3d()
#         self.B1minus =t3d()

#     def setB1plus(self,F):
#         self.B1plus.set(F)


#     def getB1plus(self):
#         if self.B1plus.isEmpty():
#             b=self.__calcB1plusminus__(minus=False)
#             self.setB1plus(b)
#             return self.B1plus
#         else:
#             return self.B1plus.get()
            
#     def __calcB1plusminus__(self,minus=True):
#         """Calculate the B1 minus or plus 

#         Args:
#             minus (bool, optional): 
#             - B1minus -> minus=True
#             - B1plus -> minus=Flase

#         Returns:
#             _type_: im.Imaginable with b1 plus or minus
#         """        
#         H=self.k.get()
#         X,Y,Z,V,C=self.k.geteSize()
#         R=np.array(H.getImageSpacing())
#         O=np.array(H.getImageOrigin())
#         d=H.getImageDimension()
#         D=np.reshape(np.array(H.getImageDirections()),[d, d])
#         DO=D[0:4,0:4]
#         #numpy is zyx while itk is xyz!!
#         # IM=np.zeros((X,Y,Z,C),dtype=np.complex128)
#         IM=np.zeros((C,Z,Y,X),dtype=np.complex128)
#         A=H.getImageArray() #[C,V,Z,Y,X]
#         if (minus):
#             IM=cnt.mu_0*A[:,0,:,:,:]-A[:,1,:,:,:]*1j
#         else:
#             IM=cnt.mu_0*A[:,0,:,:,:]+A[:,1,:,:,:]*1j
#         B1_SITK=im.createSITKImagefromArray(IM,R[np.r_[0,1,2,4]],O[np.r_[0,1,2,4]],DO.flatten())
#         O=im.Imaginable()
#         O.setImage(B1_SITK)
#         return O

# def getSNRScaling(self):
#     return np.random.random(1)
# def getNoiseCovarianceMatrix(self):
#     if self.NoiseCovarianceMatrix is not None:
#         return self.NoiseCovarianceMatrix
#     else:  
#         self.setNoiseCovarianceMatrix(self.__calculateNCM__(self.operator,self.solver))
#         return self.NoiseCovarianceMatrix
        
# def getSNR(self):
#     if self.SNR.isImageSet():
#         return self.SNR
#     else:
#         NC=self.E.getNoiseCovarianceMatrix()
#         invPsi = np.linalg.inv(NC)
#         B=self.H.getB1minus()
#         X,Y,Z,C=B.getImageSize()
#         scaling_snr=self.getSNRScaling()
#         SNR=np.zeros((Z,Y,X),dtype=np.complex128)
#         A=B.getImageArray() #[C,Z,Y,X]
#         # I know i hate this kind of loops, i miss itk and its iterator classes :(
#         for x in range(X):
#             for y in range(Y):
#                 for z in range(Z):
#                     Svox=np.matrix(A[:,z,y,x])
#                     Svox=Svox.T
#                     Spsiss=Svox.T*invPsi*Svox
#                     SNR[z,y,x]=scaling_snr*np.sqrt(Spsiss)
#         R=np.array(B.getImageSpacing())
#         O=np.array(B.getImageOrigin())
#         d=B.getImageDimension()
#         D=np.reshape(np.array(B.getImageDirections()),[d, d])
#         DO=D[0:3,0:3]
#         SNR_SITK=im.createSITKImagefromArray(SNR,R[np.r_[0,1,2]],O[np.r_[0,1,2]],DO.flatten())
#         O=im.Imaginable()
#         O.setImage(SNR_SITK)
#         self.EHSNR=O
#         # self.Log.append("Calculating Fields SNR","ok")
#         return O



