
import uuid
import boto3
from pynico_eros_montin import pynico as pn
import shutil
try:
    from .cm import saveMatlab
except:
    from cmtools.cm import saveMatlab


#FILESDICT example
# Initialize an S3 client
# s3
# filedict = {
#     "type": "file",
#     "id": -1,
#     "options": {
#         "type": "s3",
#         "filename": "file.nii.gz",
#         "bucket": "my-s3-bucket",
#         "key": "aaaaaxxxcdcdcdcdcdcd.nii.gz",
#         "options": {}
#     }
# }


# local
# filedict = {
#     "type": "file",
#     "id": -1,
#     "options": {
#         "type": "local",
#         "filename": "/data/MYDATA/TESStestData/Density.nii.gz",
#         "options": {}
#     }
# }


def getAwsCredentials(credentialsfn='~/.aws/credentials'):
    with open(credentialsfn) as f:
        lines = f.readlines()
    AWS_ACCESS_KEY = lines[1].strip()[len('aws_access_key_id')+1:]
    AWS_SECRET_KEY = lines[2].strip()[len('aws_secret_access_key')+1:]
    AWS_SESSION_TOKEN = lines[3].strip()[len('aws_session_token')+1:]   
    # aws_session_token = lines[3].strip()
    return AWS_ACCESS_KEY, AWS_SECRET_KEY,AWS_SESSION_TOKEN

def getS3Resource(aws_access_key_id, aws_secret_access_key,aws_session_token):
    return boto3.resource('s3',
                             aws_access_key_id=aws_access_key_id,
                             aws_secret_access_key=aws_secret_access_key,
                             aws_session_token=aws_session_token)

def getS3ResourceFromCredentials(credentialsfn='~/.aws/credentials'):
    AWS_ACCESS_KEY, AWS_SECRET_KEY,AWS_SESSION_TOKEN = getAwsCredentials(credentialsfn)
    return getS3Resource(AWS_ACCESS_KEY, AWS_SECRET_KEY,AWS_SESSION_TOKEN)



def s3FileTolocal(J, s3=None, pt="/tmp"):
    key = J["key"]
    bucket = J["bucket"]
    filename = J["filename"]
    O = pn.Pathable(pt)
    O.addBaseName(filename)
    O.changeFileNameRandom()
    f = O.getPosition()
    downloadFileFromS3(bucket,key,f, s3)
    J["filename"] = f
    J["type"] = "local"
    return J

def getBucketAndKeyIdFromUplaodEvent(event):
    """
    Args:
        event (_type_): _description_
    Returns:
        bucket_name (str): bucket name
        file_key (str): file key identifier
    """
    bucket_name = event["Records"][0]["s3"]["bucket"]["name"]
    file_key = event["Records"][0]["s3"]["object"]["key"]
    return bucket_name, file_key

def downloadFileFromS3(bucket_name,file_key,outfile=None, s3=None):
    """
    Args:
        bucket_name (str): bucket name
        file_key (str): file key identifier
        outfile (str): output file name
        s3 (_type_): the s3 resource
    Returns:
        filename (str): output file name
    """
    if s3 == None:
        s3 = boto3.resource("s3")
    if outfile == None:
        outfile = pn.createRandomTemporaryPathableFromFileName(file_key).getPosition()
    s3.Bucket(bucket_name).download_file(file_key, outfile)
    return outfile

def uploadFiletoS3(filename,bucket_name,file_key=None, s3=None):
    if s3 == None:
        s3 = boto3.resource("s3")
    if file_key == None:
        file_key = pn.Pathable(filename).addSuffix(f'-{str(uuid.uuid1())}').getBaseName()
    
    s3.Bucket(bucket_name).upload_file(filename, file_key)
    return {"bucket": bucket_name, "key": file_key}



def getCMRFile(filedict,s3=None):
    """    
    Args:
        s (dict): {
                "type": "file",
                "id": -1,
                "options": {
                    "type": "local",
                    "filename": "/data/MYDATA/TESStestData/Density.nii.gz",
                    "options": {}
                }
            }
    Returns:
      fn (str): position of the file in the local filesystem
    """
    s=filedict["options"]
    if (s["type"].lower()=='local'):
        return s["filename"]
    elif (s["type"].lower()=='s3'):
        T=pn.createRandomTemporaryPathableFromFileName(s["filename"]).getPosition()
        T=downloadFileFromS3(s["bucket"],s["key"],T,s3=s3)
        return T
    else:
        raise Exception("I can't get this file modality")

import pyable_eros_montin.imaginable as ima      
import numpy as np


    
import os

class cmrOutput:
    def __init__(self,app=None,outputfilename=None,tmppath=None,s3=None):
            self.out={"headers":{"options":{}},"data":[]}
            
            if outputfilename!=None:
                self.outputfilenamepathable=pn.Pathable(outputfilename)
            else:
                self.outputfilenamepathable=pn.createRandomTemporaryPathableFromFileName("a.zip")
            self.outputfilenamepathable.ensureDirectoryExistence()
            
            if tmppath==None:
                self.tmppathable=pn.createRandomTemporaryPathableFromFileName(self.outputfilenamepathable.getBaseName()).appendPathRandom()
            else:
                if not tmppath.endswith("/"):
                    tmppath=tmppath+"/"
                self.tmppathable=pn.Pathable(tmppath)
            self.tmppathable.ensureDirectoryExistence()
            self.outputfilenamepathable.ensureDirectoryExistence()
            self.setApp(app)
            #self.tmppathable wil, be the place where i stored the directory before making the zip
            self.forkable=self.tmppathable.fork()
            self.forkable.addBaseName("data")
            self.savematlab=True          

    
    def addAbleFromFilename(self,filename,id,name,type="output"):
        L=ima.Imaginable(filename=filename)
        N=pn.Pathable(filename)
        o=self.addAble(L,id,name,type,N.getBaseName())
        o["filename"]=filename
        return o
        
    def addAble(self,L,id,name,type="output",basename=None):
        pixeltype='real'
        im=L.getImageAsNumpy()

        if np.iscomplexobj(im):
            pixeltype='complex'
            L.setImageFromNumpy(im.astype(np.singlecomplex))
        if basename==None:
            basename=pn.createRandomTemporaryPathableFromFileName("a.nii.gz").getBaseName()
        o={'filename':None,
           'basename':basename,
           'able':L,
                'id':id,
                'dim':L.getImageDimension(),
                'name':name,
                'type':type,
                'numpyPixelType':im.dtype.name,
                'pixelType':pixeltype}
        
        self.out["data"].append(o)
        return o
    
    def setHeader(self,a={}):
        for k in a.keys():
            self.out["headers"][k]=a[k]
    def setPipeline(self,a):
        if "options" not in self.out["headers"].keys():
            self.out["headers"]["options"]={}
        self.out["headers"]["options"]["pipeline"]=a
        self.out["headers"]["options"]["pipelineid"]=a

    def setToken(self,a):
        if "options" not in self.out["headers"].keys():
            self.out["headers"]["options"]={}
        self.out["headers"]["options"]["token"]=a
    def setLog(self,a):
        l=self.forkable.fork()
        l.changeBaseName("log.json")
        a.writeLogAs(l.getPosition())
        self.out["log"]=a.getLog()

    def setApp(self,a):
        self.out["app"]=a
    def setTask(self,a):
        l=self.forkable.fork()
        l.changeBaseName("task.json")
        l.writeJson(a)

    
    def setOptions(self,a):
        for k in a.keys():
            self.out["headers"]["options"][k]=a[k]
        l=self.forkable.fork()
        l.changeBaseName("options.json")
        l.writeJson(a)

    def setEvent(self,a):
        l=self.forkable.fork()
        l.changeBaseName("event.json")
        l.writeJson(a)

    def exportResults(self):
        tmpdirectory=self.tmppathable.getPath()
        #check if the data are in the expected directory
        J=[]
        INFO=self.out.copy()
        INFO["data"]=[]
        for d in self.out["data"]:
            relativename="data/"+d["basename"]
            tmpfile= tmpdirectory+"/" + relativename
            pn.Pathable(tmpfile).ensureDirectoryExistence()
            if d["filename"]==None:
                d["able"].writeImageAs(tmpfile)
            elif d["filename"]!=tmpfile:
                shutil.copy(d["filename"],tmpfile)
            
            if self.savematlab:
                J.append({"name":d["name"],"data":d["able"].getImageAsNumpy()})
            # if "able" in d.keys():
            #     del d["able"]
            # if "basename" in d.keys():
            #     del d["basename"]
            info=d.copy()
            info["filename"]=relativename
            if "able" in info.keys():
                del info["able"]
            if "basename" in info.keys():
                del info["basename"]
            INFO["data"].append(info)
        
        #write the json file
        OUT=self.forkable.fork()
        OUT.changeBaseName("info.json")
        OUT.writeJson(INFO)
        if self.savematlab:
            OUT.changeBaseName("matlab.mat")
            saveMatlab(OUT.getPosition(),J)
        return tmpdirectory

    def changeOutputPath(self,path):
        #copy the data to the new path
        shutil.copytree(self.outputpath.getPosition(),path)
        shutil.rmtree(self.outputpath.getPosition())          
        self.outputpath=pn.Pathable(path)
        self.outputpath.ensureDirectoryExistence()
        self.forkable=self.outputpath.fork()
        self.forkable.addBaseName("data.nii.gz")
        
        return self.outputpath.getPosition()


    def exportAndZipResults(self,outzipfile=None,deletetemporarydirectory=False):
        tmpdirectory=self.exportResults()
        # in case the output file is not specified, i will use the outputfilenamepathable
        if outzipfile==None:
            outzipfile=self.outputfilenamepathable.getPosition()
        ext=pn.Pathable(outzipfile).getExtension()
        fi=outzipfile.replace(f'.{ext}',"")
        shutil.make_archive(fi,ext , tmpdirectory)
        if deletetemporarydirectory:
            shutil.rmtree(tmpdirectory)
        return outzipfile
    
    def exportAndZipResultsToS3(self,bucket,key=None,outzipfile=None,deletetemporarydirectory=False,deleteoutputzip=False,s3=None):
        p=self.exportAndZipResults(outzipfile=outzipfile,deletetemporarydirectory=deletetemporarydirectory)
        O= uploadFiletoS3(p,bucket,key,s3=s3)
        if deleteoutputzip:
            os.remove(p)
        return O
    
    
    
if __name__=="__main__":
    
#     filedict = {
#     "options": {
#     "type": "local",
#     "filename": "/data/PROJECTS/mroptimum/_data/signal.dat",
#     "options": {},
#     "multiraid": False,
#     "vendor": "Siemens"
# }
#     }
#     getCMRFile(filedict["options"])

    
    A=cmrOutput()
    A.addAbleFromFilename("/data/garbage/dataMYDATASRIKARPCFT1173original.nii.gz",1,"signal")
    
    P=A.exportResults()
    print(P)
    P=A.exportAndZipResults()
    print(P)
    s3=getS3ResourceFromCredentials("/home/eros/.aws/credentials")
    O=A.exportAndZipResultsToS3("mytestcmr",s3=s3,deletetemporarydirectory=True,deleteoutputzip=True)
    print(O)
    print("done")