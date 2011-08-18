#######################################################################
#
#         Copyright 2009-2010 by Ullrich Koethe
#
#    This file is part of the VIGRA computer vision library.
#    The VIGRA Website is
#        http://hci.iwr.uni-heidelberg.de/vigra/
#    Please direct questions, bug reports, and contributions to
#        ullrich.koethe@iwr.uni-heidelberg.de    or
#        vigra@informatik.uni-hamburg.de
#
#    Permission is hereby granted, free of charge, to any person
#    obtaining a copy of this software and associated documentation
#    files (the "Software"), to deal in the Software without
#    restriction, including without limitation the rights to use,
#    copy, modify, merge, publish, distribute, sublicense, and/or
#    sell copies of the Software, and to permit persons to whom the
#    Software is furnished to do so, subject to the following
#    conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the
#    Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#    OTHER DEALINGS IN THE SOFTWARE.
#
#######################################################################

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

from nose.tools import assert_equal, raises
import numpy as np
from vigra import *
from vigra.filters import *
from vigra.sampling import *
from vigra.noise import *
import vigra.arraytypes as at

#in the hope, that functions are tested in C++, we basicly test return types

#image=readImage("/export/home/nhuesken/sas/experiments/testdata/bmpmultilabel.bmp")
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
    i2=rotateImageSimple(image,orientation=RotationDirection.CLOCKWISE);
    #simple basic test
    assert(i2.shape[0] == image.shape[1])
    assert(i2.shape[1] == image.shape[0])
    #test, that they are compatible
    i3=rotateImageSimple(i2,orientation=RotationDirection.COUNTER_CLOCKWISE);
    checkImages(image,i3)
    i2=rotateImageSimple(image,orientation=RotationDirection.UPSIDE_DOWN);
    i3=rotateImageSimple(i2,orientation=RotationDirection.COUNTER_CLOCKWISE);
    i2=rotateImageSimple(i3,orientation=RotationDirection.COUNTER_CLOCKWISE);
    checkImages(image,i2)

def test_rotate():
    i2=rotateImageSimple(image)
    i3=i2
    rotateImageDegree(image,degree=90,out=i3)
    checkImages(i2,i3)

def test_resample():
    #just testing the size
    i2=resampleImage(image,factor=0.5)
    assert(i2.shape[0]==image.shape[0]*0.5)
    
    
def test_resize():
    i2=resizeImageNoInterpolation(image,shape=(image.shape[0]+10,image.shape[1]+10))
    i2=resizeImageNoInterpolation(image,shape=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=resizeImageCatmullRomInterpolation(image,shape=(image.shape[0]+10,image.shape[1]+10))
    i2=resizeImageCatmullRomInterpolation(image,shape=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=resizeImageCoscotInterpolation( image,shape=(image.shape[0]+10,image.shape[1]+10))
    i2=resizeImageCoscotInterpolation( image,shape=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=resizeImageLinearInterpolation( image,shape=(image.shape[0]+10,image.shape[1]+10))
    i2=resizeImageLinearInterpolation( image,shape=(image.shape[0],image.shape[1]))
    checkAboutSame(i2,image)

    i2=resizeImageSplineInterpolation(
        image,shape=(image.shape[0]+10,image.shape[1]+10),order=4)
    i2=resizeImageSplineInterpolation(
        image,shape=(image.shape[0],image.shape[1]), order=4)
    checkAboutSame(i2,image)

def test_2DMorphology():
    i2=discErosion(image.astype(np.uint8),radius=2)
    i3=(255-discDilation((256-image).astype(np.uint8),radius=2))
    checkImages(i2,i3)
    i2=discOpening(image.astype(np.uint8),radius=2)
    i3=(255-discDilation((256-image).astype(np.uint8),radius=2))
    checkImages(i2,i3)


def test_3DMorphologyBinary():
    i3=0
    i2=multiBinaryClosing(volumeBin,radius=2)
    i3=multiBinaryOpening(volumeBin==False,radius=2)
    i3=i3==False
    checkImages(i2,i3)

def test_3DMorphologyGrayscale():
    i2=multiGrayscaleErosion(volume256,sigma=2)
    i3=(256-multiGrayscaleDilation(256-volume256,sigma=2))
    checkImages(i2,i3)
    i2=multiGrayscaleOpening(volume256,sigma=2)
    i3=(256-multiGrayscaleClosing(256-volume256,sigma=2))
    checkImages(i2,i3)

def test_Noise():
    # ATM, we only test that these things run
    noiseVarianceEstimation(scalar_image)
    # noiseVarianceClustering(scalar_image) # This test is not expected to work whith the kind of image used here.
    nonparametricNoiseNormalization(image)
    quadraticNoiseNormalizationEstimated(image)
    linearNoiseNormalizationEstimated(image)
    quadraticNoiseNormalization(image,1.0,1.0,1.0)
    linearNoiseNormalization(image,1.0,1.0)
    
def test_Kernel1D():
    # we just test functions that were not directly imported
    contents = np.array([-1,2,3,5,3,2,-1], dtype=np.float64)
    k1 = Kernel1D()
    k1.initExplicitly(-3,3, contents)
    for k in xrange(-3, 4):
        assert(k1[k]==contents[k+3])
    k1[-2]=5
    assert(k1[-2]==5)

def test_Kernel2D():
    contents = np.array([[0,1,2],[3,4,5],[6,7,8]],dtype=np.float64)
    k2=Kernel2D()
    k2.initExplicitly((-1,-1),(1,1), contents)
    for i in xrange(-1, 2):
        for j in xrange(-1, 2):
            assert(k2[i,j]==contents[i+1, j+1])
    k2[0,-1]=-5
    assert(k2[0,-1]==-5)

