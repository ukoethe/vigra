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

from nose.tools import assert_equal, raises, assert_raises
import vigra
import numpy as np
from vigra.analysis import *
from vigra.filters import *
import vigra.arraytypes as at

img_rgb_f = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.float32)
img_scalar_f = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.float32)
img_multi_f = at.Vector3Image(np.random.rand(100,200,3)*255,dtype=np.float32)
 
img_rgb_i = at.RGBImage(np.random.rand(100,200,3)*255,dtype=np.int32)
img_scalar_i = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.uint32)
img_scalar_i64 = at.ScalarImage(np.random.rand(199,199)*4095,dtype=np.uint64)
img_scalar_ui8 = at.ScalarImage(np.random.rand(100,200)*255,dtype=np.uint8)
img_multi_i = at.Vector3Image(np.random.rand(100,200,3)*255,dtype=np.int32)

vol_rgb_f = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.float32)
vol_scalar_f = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.float32)
vol_multi_f = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.float32)
 
vol_rgb_i = at.RGBVolume(np.random.rand(100,200,60,3)*255,dtype=np.int32)
vol_scalar_i = at.ScalarVolume(np.random.rand(100,200,50)*255,dtype=np.int32)
vol_multi_i = at.Vector6Volume(np.random.rand(100,200,50,6)*255,dtype=np.int32)

def checkShape(shape1, shape2):
    if isinstance(shape1, np.ndarray):
        shape1 = shape1.shape
    if isinstance(shape2, np.ndarray):
        shape2 = shape2.shape
    assert_equal(shape1, shape2)

def checkType(i,type):
    assert_equal(i.dtype, type)

def checkAboutSame(i1,i2):
    compare_shapes(i1.shape, i2.shape)
    difference=np.sum(np.abs(i1-i2))/float(np.size(i1))
    assert(difference<5)
    
def test_watersheds():

    res = watershedsUnionFind(img_scalar_f)
    checkShape(img_scalar_f,res[0])
    checkType(res[0], np.uint32)

    res = watershedsUnionFind(img_scalar_f, 4)
    checkShape(img_scalar_f, res[0])
    checkType(res[0], np.uint32)

    img_scalar_i = at.ScalarImage(img_scalar_f.shape, dtype=np.uint32)
    res = watershedsUnionFind(img_scalar_f, 8, img_scalar_i)
    assert(res[0] is img_scalar_i)

    assert_raises(RuntimeError, watershedsUnionFind, img_scalar_f, 5, img_scalar_i)
    
    res = watersheds(img_scalar_f, seeds=img_scalar_i, method="RegionGrowing")
    checkShape(img_scalar_f, res[0])
    checkType(res[0], np.uint32)

    # 3D
    res = watershedsUnionFind(vol_scalar_f, 6)
    checkShape(vol_scalar_f, res[0])

    vol_scalar_i = at.ScalarVolume(vol_scalar_f.shape, dtype=np.uint32)
    res = watershedsUnionFind(vol_scalar_f, 26, vol_scalar_i)
    assert(res[0] is vol_scalar_i)
    
    assert_raises(RuntimeError, watershedsUnionFind, img_scalar_f, 17, img_scalar_i)

def test_MinimaMaxima():
    #2D Case
    res = localMinima(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,img_scalar_f.dtype)
    
    res = extendedLocalMinima(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,img_scalar_f.dtype)
    
    res = localMaxima(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,img_scalar_f.dtype)
    
    res = extendedLocalMaxima(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,img_scalar_f.dtype)
    
    res = labelImage(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,np.uint32)
    
    res = labelImageWithBackground(img_scalar_f)
    checkShape(img_scalar_f,res)
    checkType(res,np.uint32)
    
    #3D Case
    res = localMinima3D(vol_scalar_f)
    checkShape(vol_scalar_f,res)
    checkType(res,vol_scalar_f.dtype)
    
    res = extendedLocalMinima3D(vol_scalar_f)
    checkShape(vol_scalar_f,res)
    checkType(res,vol_scalar_f.dtype)
    
    res = localMaxima3D(vol_scalar_f)
    checkShape(vol_scalar_f,res)
    checkType(res,vol_scalar_f.dtype)
    
    res = extendedLocalMaxima3D(vol_scalar_f)
    checkShape(vol_scalar_f,res)
    checkType(res,vol_scalar_f.dtype)
    
    data = vigra.taggedView(np.zeros((100,200,50),dtype=np.float32), 
                            vigra.defaultAxistags(3, noChannels=True))

    data[10,5,10]=1
    data[98,10,10]=1
    data[90,50,30]=1
    data[10,10,2]=1
    data[10,80,8]=1
    data[10,10,8]=1
    data[10,150,2]=1
    data[80,10,10]=1
    data[70,120,40]=1
    
    res = localMaxima3D(data,neighborhood=26)
    np.testing.assert_array_equal(res, data, "Error in the calculation of the Maxima")

    res = extendedLocalMaxima3D(data,neighborhood=26)
    np.testing.assert_array_equal(res, data, "Error in the calculation of the Extended Maxima")
    
    res = localMinima3D(1-data,neighborhood=26)
    np.testing.assert_array_equal(res, data, "Error in the calculation of the Minima")
    
    res = extendedLocalMinima3D(1-data,neighborhood=26)
    np.testing.assert_array_equal(res, data, "Error in the calculation of the Extended Minima")

def test_Region2Crack():
    res = regionImageToCrackEdgeImage(img_scalar_i)
    assert_equal(img_scalar_f.shape[0]*2-1, res.shape[0])
    assert_equal(img_scalar_f.shape[1]*2-1, res.shape[1])
    checkType(res,res.dtype)

    regionImageToCrackEdgeImage(img_scalar_i64[0:100,0:100], 1, img_scalar_i64)
    checkType(img_scalar_i64, np.uint64)
    
    res = regionImageToEdgeImage(img_scalar_i)
    checkShape(res,img_scalar_i)
    
def test_transforms():
    res = distanceTransform2D(img_scalar_f)
    checkShape(img_scalar_i, res)
    #print >> sys.stderr, res.dtype,
    checkType(res, np.float32)
    
    res = distanceTransform2D(img_scalar_ui8)
    checkShape(img_scalar_ui8, res)
    checkType(res, np.float32)
    
    res = radialSymmetryTransform2D(img_scalar_f,1)
    checkShape(img_scalar_f, res)
    checkType(res, np.float32)

def test_cornerss():
    res = cornernessHarris(img_scalar_f,1)
    checkShape(img_scalar_f, res)
    checkType(res, np.float32)
    
    res = cornernessFoerstner(img_scalar_f,2)
    checkShape(img_scalar_f, res)
    checkType(res, np.float32)
    
    res = cornernessRohr(img_scalar_f,0.5)
    checkShape(img_scalar_f, res)
    checkType(res, np.float32)
    
    res = cornernessBeaudet(img_scalar_f,1)
    checkShape(img_scalar_f, res)
    checkType(res, np.float32)

def test_edges():    
    res = cannyEdgeImage(img_scalar_f, 1,128,255)
    checkShape(img_scalar_f, res)
    checkType(res, np.uint8)
    
    res = cannyEdgeImageWithThinning(img_scalar_f, 1,128,255)
    checkShape(img_scalar_f, res)
    checkType(res, np.uint8)
    
    res = shenCastanEdgeImage(img_scalar_f, 1,128,255)
    checkShape(img_scalar_f, res)
    checkType(res, np.uint8)
    
    res = shenCastanCrackEdgeImage(img_scalar_f, 1,128,255)
    assert_equal(img_scalar_f.shape[0]*2-1, res.shape[0])
    assert_equal(img_scalar_f.shape[1]*2-1, res.shape[1])
    
    res1 = beautifyCrackEdgeImage(res,  1, 0)
    checkShape(res1, res)
    checkType(res1, np.uint8)
    
    res1 = closeGapsInCrackEdgeImage(res, 4)
    checkShape(res1, res)
    checkType(res1, np.uint8)
    
    res = removeShortEdges(img_scalar_ui8, 10, 0)
    checkShape(img_scalar_ui8, res)
    checkType(res, np.uint8)
    
    res = boundaryTensor2D(img_scalar_f, 1)
    assert_equal(img_scalar_f.shape[0], res.shape[0])
    assert_equal(img_scalar_f.shape[1], res.shape[1])
    assert_equal(res.shape[2], 3)    
    checkType(res, np.float32)
    
    res = hourGlassFilter2D(img_multi_f, 1, 2)
    assert_equal(img_multi_f.shape[0], res.shape[0])
    assert_equal(img_multi_f.shape[1], res.shape[1])
    checkType(res, np.float32)
    
    res = tensorEigenRepresentation2D(img_multi_f)
    assert_equal(img_multi_f.shape[0], res.shape[0])
    assert_equal(img_multi_f.shape[1], res.shape[1])
    checkType(res, np.float32)
    
    res = tensorTrace(img_multi_f)
    assert_equal(img_multi_f.shape[0], res.shape[0])
    assert_equal(img_multi_f.shape[1], res.shape[1])
    checkType(res, np.float32)
    
    res = rieszTransformOfLOG2D(img_scalar_f, 1, 1, 1)
    assert_equal(img_multi_f.shape[0], res.shape[0])
    assert_equal(img_multi_f.shape[1], res.shape[1])
    checkType(res, np.float32)
    
    
def ok_():
    print >> sys.stderr, ".",
        
