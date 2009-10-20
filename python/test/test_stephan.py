execfile('set_paths.py')

import sys

from nose.tools import assert_equal, raises
import numpy as np
import vigranumpycmodule as vm
import arraytypes as at

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
	
def tttest_watersheds():
	print >> sys.stderr, " asdfa sasf "
	res = vm.watersheds(img_scalar_f)
	assert(res.shape==img_scalar_f.shape)
    
def test_structureTensor():
    print >> sys.stderr, "shape:"
    print >> sys.stderr, img_rgb_f.shape
    res = vm.structureTensor(img_scalar_f,1.0,2.0, out=img_3_f)
    res = vm.structureTensor(img_scalar_f,1.0,2.0)
    #assert(res.shape==img_scalar_f.shape)
	
def test_simpleSharpening():
    res = vm.simpleSharpening(img_scalar_f)
    
def test_gaussianSharpening():
    res = vm.gaussianSharpening(img_scalar_f)
    
def test_convolution():
    k2_ = vm.Kernel2D()
    k2_.initWithFactoryKernel(vm.Kernel2D.kernelDisk(10))
    
    #k3 = vm.gaussianDerivativeKernel(sigma, order)
    #guassianDerivative(img, sx,ox,sy,oy,sz,oz)
    #guassianDerivative(img, (sx,sy, sz), (ox, oy,oz))
    
    k2=vm.Kernel2D()
    k2.initSetExplicitly((-1,-1),(1,1),np.array([[0,1,2],[1,2,3],[2,3,4]],dtype=np.float64))
    res = vm.convolveImage(img_scalar_f, k2)
