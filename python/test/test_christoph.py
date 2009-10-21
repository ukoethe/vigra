import sys
print >> sys.stderr, "executing test file", __file__
execfile('set_paths.py')

from nose.tools import assert_equal, raises, assert_raises
import numpy as np
import vigranumpycmodule as vm
import arraytypes as at

img_rgb_f = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.float32)
img_scalar_f = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.float32)
img_multi_f = at.Vector3Image(np.random.rand(100,200,3)*255,dtype=np.float32)
 
img_rgb_i = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.int32)
img_scalar_i = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.int32)
img_scalar_i64 = at.ScalarImage(np.random.rand(199,199)*4095,dtype=np.int64)
img_scalar_ui8 = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.uint8)
img_multi_i = at.Vector3Image(np.random.rand(100,200,3)*255,dtype=np.int32)

vol_rgb_f = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.float32)
vol_scalar_f = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.float32)
vol_multi_f = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.float32)
 
vol_rgb_i = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.int32)
vol_scalar_i = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.int32)
vol_multi_i = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.int32)

def checkShape(i1,i2):
    assert(i1.shape==i2.shape)
def checkType(i,type):
	assert(i.dtype == type)

def checkAboutSame(i1,i2):
    assert(i1.shape==i2.shape)
    difference=np.sum(np.abs(i1-i2))/float(np.size(i1))
    assert(difference<5)
	
def test_watersheds():

	res = vm.watersheds(img_scalar_f)
	checkShape(img_scalar_f,res)

	res = vm.watersheds(img_scalar_f,4)
	checkShape(img_scalar_f,res)

	vm.watersheds(img_scalar_f,8,img_scalar_i)
	checkShape(img_scalar_f,img_scalar_i)
	checkType(img_scalar_i,np.int32)

	assert_raises(RuntimeError, vm.watersheds, img_scalar_f, 5, img_scalar_i)
	
	# 3D
	res = vm.watersheds(vol_scalar_f,6)
	checkShape(vol_scalar_f,res)

	vm.watersheds(vol_scalar_f,26,vol_scalar_i)
	checkShape(vol_scalar_f,res)
	
	assert_raises(RuntimeError, vm.watersheds, img_scalar_f, 17, img_scalar_i)

def test_region_growing():
	res = vm.seededRegionGrowingSeeded2D(img_scalar_f,img_scalar_i)
	checkShape(img_scalar_f,res)

def test_seededRegionGrowing2D():
	res = vm.seededRegionGrowing2D(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,np.int32)
def test_MinimaMaxima():
	res = vm.localMinima2D(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)
	
	res = vm.extendedLocalMinima2D(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)
	
	res = vm.localMaxima2D(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)
	
	res = vm.extendedLocalMaxima2D(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)
	
	res = vm.labelImage(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)
	
	res = vm.labelImageWithBackground(img_scalar_f)
	checkShape(img_scalar_f,res)
	checkType(res,img_scalar_f.dtype)

def test_Region2Creack():
	res = vm.regionImageToCrackEdgeImage(img_scalar_i)
	assert(img_scalar_f.shape[0]*2-1 == res.shape[0])
	assert(img_scalar_f.shape[1]*2-1 == res.shape[1])
	checkType(res,res.dtype)

	vm.regionImageToCrackEdgeImage(img_scalar_i64[0:100,0:100], 1, img_scalar_i64)
	checkType(img_scalar_i64, np.int64)
	
	res = vm.regionImageToEdgeImage(img_scalar_i)
	checkShape(res,img_scalar_i)
	
def test_transforms():
	res = vm.distanceTransform2D(img_scalar_i)
	checkShape(img_scalar_i, res)
	#print >> sys.stderr, res.dtype,
	checkType(res, np.float32)
	
	res = vm.distanceTransform2D(img_scalar_ui8)
	checkShape(img_scalar_ui8, res)
	checkType(res, np.float32)
	
	res = vm.radialSymmetryTransform2D(img_scalar_f,1)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)

def test_cornerss():
	res = vm.cornerResponseFunction2D(img_scalar_f,1)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.foerstnerCornerDetector2D(img_scalar_f,2)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.rohrCornerDetector2D(img_scalar_f,0.5)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.beaudetCornerDetector2D(img_scalar_f,1)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)

def test_edges():	
	res = vm.cannyEdgeImage(img_scalar_f, 1,128,255)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.cannyEdgeImageWithThinning(img_scalar_f, 1,128,255)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.shenCastanEdgeImage(img_scalar_f, 1,128,255)
	checkShape(img_scalar_f, res)
	checkType(res, np.float32)
	
	res = vm.shenCastanCrackEdgeImage(img_scalar_f, 1,128,255)
	assert(img_scalar_f.shape[0]*2-1 == res.shape[0])
	assert(img_scalar_f.shape[1]*2-1 == res.shape[1])
	
	res = vm.removeShortEdges(img_scalar_i, 10, 0)
	checkShape(img_scalar_f, res)
	checkType(res, np.int32)
	
	res = vm.beautifyCrackEdgeImage(img_scalar_i,  1, 0)
	checkShape(img_scalar_f, res)
	checkType(res, np.int32)
	
	res = vm.closeGapsInCrackEdgeImage(img_scalar_i, 4)
	checkShape(img_scalar_f, res)
	checkType(res, np.int32)
	
	res = vm.boundaryTensor2D(img_scalar_f, 1)
	assert(img_scalar_f.shape[0]== res.shape[0])
	assert(img_scalar_f.shape[1] == res.shape[1])
	assert(res.shape[2] == 3)	
	checkType(res, np.float32)
	
	res = vm.hourGlassFilter2D(img_multi_f, 1, 2)
	assert(img_multi_f.shape[0]== res.shape[0])
	assert(img_multi_f.shape[1] == res.shape[1])
	checkType(res, np.float32)
	
	res = vm.tensorEigenRepresentation2D(img_multi_f)
	assert(img_multi_f.shape[0]== res.shape[0])
	assert(img_multi_f.shape[1] == res.shape[1])
	checkType(res, np.float32)
	
	res = vm.tensorTrace2D(img_multi_f)
	assert(img_multi_f.shape[0]== res.shape[0])
	assert(img_multi_f.shape[1] == res.shape[1])
	checkType(res, np.float32)
	
	res = vm.rieszTransformOfLOG2D(img_scalar_f, 1, 1, 1)
	assert(img_multi_f.shape[0]== res.shape[0])
	assert(img_multi_f.shape[1] == res.shape[1])
	checkType(res, np.float32)
	
	
def ok_():
	print >> sys.stderr, ".",
		
