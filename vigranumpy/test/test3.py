import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

from nose.tools import assert_equal, raises
from vigra import numpy as np
from vigra import *
from vigra import arraytypes as at
from vigra.filters import *
from vigra.analysis import *

img_rgb_f = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.float32)
img_scalar_f = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.float32)
img_multi_f = at.Vector4Image(np.random.rand(100,200,4)*255,dtype=np.float32)
img_3_f = at.Vector3Image(np.random.rand(100,200,3)*255,dtype=np.float32)
 
img_rgb_i = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.int32)
img_scalar_i = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.int32)
img_multi_i = at.Vector4Image(np.random.rand(100,200,4)*255,dtype=np.int32)

vol_rgb_f = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.float32)
vol_scalar_f = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.float32)
vol_multi_f = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.float32)
 
vol_rgb_i = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.int32)
vol_scalar_i = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.int32)
vol_multi_i = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.int32)

def checkImages(i1,i2):
    assert(i1.shape==i2.shape)
    assert(np.sum(i1==i2)!=0)

def checkAboutSame(i1,i2):
    assert(i1.shape==i2.shape)
    difference=np.sum(np.abs(i1-i2))/float(np.size(i1))
    assert(difference<5)
	
def test_watersheds():
	res = watersheds(img_scalar_f)
	assert(res[0].shape==img_scalar_f.shape)
    
def test_structureTensor():
    res = structureTensor(img_scalar_f,1.0,2.0, out=img_3_f)
    res = structureTensor(img_scalar_f,1.0,2.0)
    res = structureTensor(img_rgb_f,1.0,2.0)
    assert(res.shape==img_scalar_f.shape + (3,))
	
def test_simpleSharpening():
    res = simpleSharpening(img_scalar_f)
    
def test_gaussianSharpening():
    res = gaussianSharpening(img_scalar_f)
    
def test_convolution():
    krnl = gaussianKernel(0.5)

    
    k2_ = Kernel2D()
    k2_.initDisk(10)
    
    #k3 = gaussianDerivativeKernel(sigma, order)
    #guassianDerivative(img, sx,ox,sy,oy,sz,oz)
    #guassianDerivative(img, (sx,sy, sz), (ox, oy,oz))
    
    k2=Kernel2D()
    k2.initExplicitly((-1,-1),(1,1),np.array([[0,1,2],[1,2,3],[2,3,4]],dtype=np.float64))
    res = convolve(img_scalar_f, k2)


def test_multiconvolution():
   vol = vol_scalar_f
   k=Kernel1D()
   k.initGaussian(0.5)
   a=convolveOneDimension(vol, 0, k)
   a=convolveOneDimension(a, 1, k)
   a=convolveOneDimension(a, 2, k)
         
   b=convolve(vol, k)
     
   c=convolve(vol, (k, k, k))
   checkAboutSame(a, b)
   checkAboutSame(a, c)
    
