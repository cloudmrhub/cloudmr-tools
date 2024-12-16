% test gfactor

SIGNALFILE='../data/meas_MID00024_FID188178_Multislice.dat';
NOISEFILE='../data/meas_MID00027_FID188181_Multislice_no_RF.dat';



N=readRawDataFilter('noise'); %we need to inform the class that it's going to read a noise data so it won't remove the fewquency oversample
N.setFilename(NOISEFILE); %set the filename
nReader=N.getOutput(); %the ouptut is a rawdata reader class (cm2DRawDataReader)


%instantiate the readingFilter
S=readRawDataFilter('signal'); %we need to inform the class that it's going to read a signal data so it will reomove the frequency overesampling
S.setFilename(SIGNALFILE); %set the filename
sReader=S.getOutput(); %the ouptut is a rawdata reader class (cm2DRawDataReader)


%select a slice and reconstruct
average=1;
contrast=1;
repetition=1;
slice=1

KS=sReader.getRawDataImageKSpaceSlice(average,contrast,repetition,slice);
KN=nReader.getRawDataImageKSpaceSlice(average,contrast,repetition,slice);


R1=1;
R2=4;
ACL=24;
RE=cm2DGFactorSENSE();
%mimik a grappa acquisition from a fully sample
TITLE='SENSE ' ;
ACCK=RE.mimicmSenseDataFromFullysampledZeroPadded(KS,R1,R2,ACL);
RE.setSignalKSpace(ACCK)
RE.setNoiseKSpace(KN)
RE.setCoilSensitivityMatrixCalculationMethod('innerACL');
%don't need to set the sensitivity source since it will be the signal for simple
RE.setAccelerationFrequency(R1)
RE.setAccelerationPhase(R2);
RE.setAutocalibration(ACL);
OUT=RE.getOutput();

figure()
imshow(OUT,[])
colorbar
