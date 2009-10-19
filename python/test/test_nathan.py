from nose.tools import assert_equal, raises
import numpy as np
import vigranumpycmodule as vm
import arraytypes as at

#in the hope, that functions are tested in C++, we basicly test return types

#image=vm.readImage("/export/home/nhuesken/sas/experiments/testdata/bmpmultilabel.bmp")
image=at.RGBImage(np.random.rand(100,100,3)*255,dtype=np.float32)
scalar_image=at.ScalarImage(np.random.rand(100,100)*255,dtype=np.float32)
volume256=at.Volume(np.random.rand(100,100,100)*255,dtype=np.uint8)
volumeBin=at.Volume(np.random.rand(100,100,100))>0.5

def checkImages(i1,i2):
    assert(i1.shape==i2.shape)
    assert(np.sum(i1==i2)!=0)

def checkAboutSame(i1,i2):
    assert(i1.shape==i2.shape)
    difference=np.sum(np.abs(i1-i2))/float(np.size(i1))
    assert(difference<5)

def test_simpleRotate():
    i2=vm.rotateImageSimple(image,orientation=vm.RotationDirection.CLOCKWISE);
    #simple basic test
    assert(i2.shape[0] == image.shape[1])
    assert(i2.shape[1] == image.shape[0])
    #test, that they are compatible
    i3=vm.rotateImageSimple(i2,orientation=vm.RotationDirection.COUNTER_CLOCKWISE);
    checkImages(image,i3)
    i2=vm.rotateImageSimple(image,orientation=vm.RotationDirection.UPSIDE_DOWN);
    i3=vm.rotateImageSimple(i2,orientation=vm.RotationDirection.COUNTER_CLOCKWISE);
    i2=vm.rotateImageSimple(i3,orientation=vm.RotationDirection.COUNTER_CLOCKWISE);
    checkImages(image,i2)

def test_rotate():
    i2=vm.rotateImageSimple(image)
    i3=i2
    vm.rotateImageDegree(image,degree=90,out=i3)
    checkImages(i2,i3)

def test_resample():
    #just testing the size
    i2=vm.resampleImage(image,factor=0.5)  
    assert(i2.shape[0]==image.shape[0]*0.5)
    
    
def test_resize():
    i2=vm.resizeImageNoInterpolation(image,destSize=(image.shape[0]+10,image.shape[1]+10))
    i2=vm.resizeImageNoInterpolation(image,destSize=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=vm.resizeImageCatmullRomInterpolation(image,destSize=(image.shape[0]+10,image.shape[1]+10))
    i2=vm.resizeImageCatmullRomInterpolation(image,destSize=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=vm.resizeImageCoscotInterpolation( image,destSize=(image.shape[0]+10,image.shape[1]+10))
    i2=vm.resizeImageCoscotInterpolation( image,destSize=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=vm.resizeImageLinearInterpolation( image,destSize=(image.shape[0]+10,image.shape[1]+10))
    i2=vm.resizeImageLinearInterpolation( image,destSize=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=vm.resizeImageSplineInterpolation(
        image,destSize=(image.shape[0]+10,image.shape[1]+10),splineOrder=4)
    i2=vm.resizeImageSplineInterpolation(
        image,destSize=(image.shape[0],image.shape[1]), splineOrder=4)
    checkAboutSame(i2,image)

def test_2DMorphology():
    i2=vm.discErosion(image.astype(np.uint8),radius=2)
    i3=(255-vm.discDilation((256-image).astype(np.uint8),radius=2))
    checkImages(i2,i3)
    i2=vm.discOpening(image.astype(np.uint8),radius=2)
    i3=(255-vm.discDilation((256-image).astype(np.uint8),radius=2))
    checkImages(i2,i3)


def test_3DMorphologyBinary():
    i3=0
    i2=vm.multiBinaryClosing(volumeBin,radius=2)
    i3=vm.multiBinaryOpening(volumeBin==False,radius=2)
    i3=i3==False
    checkImages(i2,i3)

def test_3DMorphologyGrayscale():
    i2=vm.multiGrayscaleErosion(volume256,radius=2)
    i3=(256-vm.multiGrayscaleDilation(256-volume256,radius=2))
    checkImages(i2,i3)
    i2=vm.multiGrayscaleOpening(volume256,radius=2)
    i3=(256-vm.multiGrayscaleClosing(256-volume256,radius=2))
    checkImages(i2,i3)

def test_Noise():
    #I just test, if the things run or not
    vm.noiseVarianceEstimation(scalar_image)
    vm.noiseVarianceClustering(scalar_image)
    vm.nonparametricNoiseNormalization(image)
    vm.quadraticNoiseNormalizationEstimated(image)
    vm.linearNoiseNormalizationEstimated(image)
    vm.quadraticNoiseNormalization(image,1.0,1.0,1.0)
    vm.linearNoiseNormalization(image,1.0,1.0)
    
def test_Kernel1D():
    #we just test functions, that where not directly imported
    k1=vm.Kernel1D()
    k1.initSetExplicitly(-3,3,np.array([-1,2,3,5,3,2,-1],dtype=np.float64))
    assert(k1[-3]==-1)
    assert(k1[3]==-1)
    k1[-2]=5
    assert(k1[-2]==5)

def test_Kernel2D():
    k2=vm.Kernel2D()
    k2.initSetExplicitly((-1,-1),(1,1),np.array([[0,1,2],[1,2,3],[2,3,4]],dtype=np.float64))
    assert(k2[(-1,0)]==1)
    k2[(0,-1)]=-5
    assert(k2[(0,-1)]==-5)

