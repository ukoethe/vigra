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

# run with a simple 'nosetests' in this directory
# (and nose installed, i.e. 'easy_install nose')

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

import vigra.arraytypes as arraytypes
import vigra.ufunc as ufunc
import numpy, copy
import vigranumpytest as vt
from nose.tools import assert_equal, raises

from vigra.arraytypes import AxisTags, AxisInfo

numpyHasComplexNegateBug = numpy.version.version.startswith('1.0')

try:
    vt.testAny()
except Exception, e:
    ArgumentError = type(e)
    
allTests = set()
for n, f in vt.__dict__.items():
    if n.startswith('test'):
        allTests.add(n)
 
def checkShape(shape1, shape2):
    assert_equal(shape1, shape2)

def checkStride(stride1, stride2):
    assert_equal(stride1, stride2)

def checkArray(cls, channels, dim):
        def testCopy(img):
            b = cls(img, order='A')
            assert_equal(sys.getrefcount(b), 2)
            assert b.__class__ is img.__class__
            assert_equal(b.shape, img.shape)
            # print b.shape, img.shape, b.strides, img.strides
            # assert False
            assert_equal(b.strides, img.strides)
            assert_equal(b.order, img.order)
            assert_equal(b.flags.c_contiguous, img.flags.c_contiguous)
            assert_equal(b.flags.f_contiguous, img.flags.f_contiguous)
            assert (b == img).all()
            assert not numpy.may_share_memory(b, img)
            b = img.copy(order='A')
            assert_equal(sys.getrefcount(b), 2)
            assert_equal(b.shape, img.shape)
            assert_equal(b.strides, img.strides)
            assert_equal(b.order, img.order)
            assert_equal(b.flags.c_contiguous, img.flags.c_contiguous)
            assert_equal(b.flags.f_contiguous, img.flags.f_contiguous)
            assert (b == img).all()
            assert not numpy.may_share_memory(b, img)
        
        shape = (20, 10, 5, channels)[-(dim+1):]
        
        # figure out expected strides for various memory layouts
        d = dim+1
        rshape = shape[:d]
        fcstrides = (200*channels, 20*channels, 4*channels, 4)[-d:]
        bcstrides = (50*channels, 5*channels, channels, 1)[-d:]
        ffstrides = (4, shape[0]*4, shape[0]*shape[1]*4, shape[0]*shape[1]*shape[2]*4)[:d]
        bfstrides = (1, shape[0], shape[0]*shape[1], shape[0]*shape[1]*shape[2])[:d]
        fvstrides = (channels*4, channels*shape[0]*4, channels*shape[0]*shape[1]*4)[:d-1]+(4,)
        bvstrides = (channels, channels*shape[0], channels*shape[0]*shape[1])[:d-1]+(1,)
        caxistags = AxisTags([AxisInfo.z, AxisInfo.y, AxisInfo.x][3-dim:] + [AxisInfo.c])
#        caxistags = AxisTags([AxisInfo.x, AxisInfo.y, AxisInfo.z][:dim] + [AxisInfo.c])
        vaxistags = AxisTags([AxisInfo.x, AxisInfo.y, AxisInfo.z][:dim] + [AxisInfo.c])
        faxistags = AxisTags([AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z][:dim+1])
        
        value = 1 if channels == 1 else range(1,channels+1)

        # test type
        img = cls(shape, order="V")
#        assert type(img) is cls
        assert isinstance(img, numpy.ndarray)
        assert_equal(img.dtype, numpy.float32)
        assert_equal(sys.getrefcount(img), 2)
        
        # test shape
        checkShape(img.shape, rshape)
        assert_equal(img.width, rshape[0])
        assert_equal(img.height, rshape[1])
        if dim == 3:
            assert_equal(img.depth, rshape[2])
        assert_equal(img.channels, channels)
        assert_equal(img.spatialDimensions, dim)
       
        # test strides and order
        checkStride(img.strides, fvstrides)
#        assert_equal(img.order, "F" if channels == 1 else "V")
        assert_equal(img.order, "V")
        assert not img.flags.c_contiguous
#        assert_equal(img.flags.f_contiguous, True if channels == 1 else False)
        assert_equal(img.flags.f_contiguous, False)
        
        # test axistags
        assert_equal(img.axistags, vaxistags)
        # FIXME: add more tests

        # test initialization and assignment
        assert_equal(img.min(), 0.0)
        assert_equal(img.max(), 0.0)
        img.flat[:] = range(img.size)
        assert_equal(img.flatten().tolist(), range(img.size))
        img[1,2] = value
        assert (img[1,2]==value).all()

        # test that copy and ufuncs preserve memory layout
        testCopy(img)
        assert_equal(img.shape, (-img).shape)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.shape, (img+img).shape)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.shape, (img*2).shape)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (img*2).axistags)
        
        if channels == 1 or getattr(cls, 'channels', 0) > 0:
            img = cls(shape)
            checkShape(img.shape, rshape)
            checkStride(img.strides, fvstrides)

        # test shape, strides, and copy for 'F' order
        img = cls(shape, order='F')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, rshape)
        checkStride(img.strides, ffstrides)
        assert_equal(img.axistags, faxistags)
        assert_equal(img.order, "F")
        assert not img.flags.c_contiguous
        assert img.flags.f_contiguous
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)

        # test shape, strides, and copy for 'A' order (should be equal to 'V' order)
        img = cls(shape, order='A')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, rshape)
        checkStride(img.strides, fvstrides)
#        assert_equal(img.order, "F" if channels == 1 else "V")
        assert_equal(img.order, "V")
        assert_equal(img.axistags, vaxistags)
        assert not img.flags.c_contiguous
#        assert_equal(img.flags.f_contiguous, True if channels == 1 else False)
        assert_equal(img.flags.f_contiguous, False)
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)

        # test shape, strides, and copy for 'C' order
        img = cls(shape, order='C')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, rshape)
        checkStride(img.strides, fcstrides)
        assert_equal(img.axistags, caxistags)
        assert_equal(img.order, "C")
        assert img.flags.c_contiguous
        assert not img.flags.f_contiguous
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)

        value = 10 if channels == 1 else range(10,channels+10)
        zero = 0 if channels == 1 else (0,)*channels
        
        # test shape, strides, and copy for dtype uint8
        b = cls(img, dtype=numpy.uint8, order='V')
        assert_equal(sys.getrefcount(b), 2)
        assert_equal(b.dtype, numpy.uint8)
        checkShape(b.shape, rshape)
        checkStride(b.strides, bvstrides)
#        assert_equal(b.order, "F" if channels == 1 else "V")
        assert_equal(b.order, "V")
        assert_equal(b.axistags, img.axistags)
        assert not b.flags.c_contiguous
#        assert_equal(b.flags.f_contiguous, True if channels == 1 else False)
        assert_equal(b.flags.f_contiguous, False)
        assert (b==img).all()
        b[2,1] = value
        assert (b[2,1]==value).all()
        assert (img[2,1]==zero).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)
        assert_equal(b.axistags, (-b).axistags)
        assert_equal(b.axistags, (b+b).axistags)
        assert_equal(b.axistags, (b*2).axistags)

        b = cls(img, dtype=numpy.uint8, order='C')
        assert_equal(sys.getrefcount(b), 2)
        checkShape(b.shape, rshape)
        checkStride(b.strides, bcstrides)
        assert_equal(b.axistags, img.axistags)
        assert_equal(b.order, "C")
        assert b.flags.c_contiguous
        assert not b.flags.f_contiguous
        assert (b==img).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)
        assert_equal(b.axistags, (-b).axistags)
        assert_equal(b.axistags, (b+b).axistags)
        assert_equal(b.axistags, (b*2).axistags)
        
        b = cls(img, dtype=numpy.uint8, order='F')
        assert_equal(sys.getrefcount(b), 2)
        checkShape(b.shape, rshape)
        checkStride(b.strides, bfstrides)
        assert_equal(b.axistags, img.axistags)
        assert_equal(b.order, "F")
        assert not b.flags.c_contiguous
        assert b.flags.f_contiguous
        assert (b==img).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)
        assert_equal(b.axistags, (-b).axistags)
        assert_equal(b.axistags, (b+b).axistags)
        assert_equal(b.axistags, (b*2).axistags)
        
        value = 100 if channels == 1 else range(100,channels+100)

        # test ndarray view
        v = img.view(numpy.ndarray)
        assert type(v) is numpy.ndarray
        assert (v==img).all()
        assert numpy.may_share_memory(v, img)
        v[3,4] = value
        assert (v[3,4]==value).all()
        assert (v==img).all()

def checkFailure(obj, n):
    f = getattr(vt, n)
    try:
        f(obj)
    except ArgumentError:
        return
    raise AssertionError, "%r did not throw ArgumentError as expected when passed a %r with shape %s, stride %s, axistags '%s'" % (n, type(obj), str(obj.shape), str(obj.strides), str(obj.axistags))

def checkCompatibility(obj, compatible):
    for n in compatible:
        try:
            f = getattr(vt, n)
            res = f(obj)
            assert_equal(obj.shape, res)
        except Exception:
            print "exception in %s with shape %s strides %s tags (%s)" % (n, obj.shape, obj.strides, 
                                            getattr(obj, "axistags", "none"))
            raise
        
        
    incompatible = allTests.difference(compatible)
    
    for n in incompatible:
        try:
            checkFailure(obj, n)
        except Exception:
            print "exception in",n
            raise

def testImage1():
    checkArray(arraytypes.Image, 1, 2)
    
    shape = (20, 10)
    rshape = (10, 20)

    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]

    checkCompatibility(arraytypes.Image(rshape, order='C'), c)
    
    # FIXME: this should behave almost like a C-order array
#    checkCompatibility(numpy.ndarray(rshape), c)
    
    checkCompatibility(arraytypes.Image(shape, order='V'), c)

    checkCompatibility(arraytypes.Image(shape, order='F'), c)
    
    img = arraytypes.Image(rshape, order='C')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.Image(shape, order='V')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.Image(shape, order='F')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[:,0,0], 1)
 
def testImage2():
    checkArray(arraytypes.Image, 2, 2)
    
    shape = (10, 20, 2)
    cshape = (20, 10, 2)
    fshape = (2, 10, 20)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
         
    checkCompatibility(arraytypes.Image(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.Image(shape, order='V'), c)
    
    checkCompatibility(arraytypes.Image(fshape, order='F'), c)
    
    img = arraytypes.Image(cshape, order='C')
    assert_equal(vt.viewArray3Unstrided(img), fshape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Image(shape, order='V')
    assert_equal(vt.viewArray3Unstrided(img), fshape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Image(fshape, order='F')
    assert_equal(vt.viewArray3Strided(img), fshape)
    assert (img[:,0,0]==(1,0)).all()
    
    img = arraytypes.Image(cshape, order='C')
    assert_equal(vt.viewImageVector2Strided(img), shape[:-1])
    assert (img[0,0]==(1,1)).all()
    
    img = arraytypes.Image(shape, order='V')
    assert_equal(vt.viewImageVector2Unstrided(img), shape[:-1])
    assert (img[0,0]==(1,1)).all()
 
def testScalarImage():
    checkArray(arraytypes.ScalarImage, 1, 2)
    
    cshape = (10, 20)
    shape = (20, 10)

    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
         
    checkCompatibility(arraytypes.ScalarImage(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.ScalarImage(shape, order='V'), c)
    
    checkCompatibility(arraytypes.ScalarImage(shape, order='F'), c)
    
    img = arraytypes.ScalarImage(cshape, order='C')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.ScalarImage(shape, order='V')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.ScalarImage(shape, order='F')
    checkShape(vt.viewArray2Strided(img), shape)
    # FIXME: remove third index (strip it in the ScalarImage factory function, and change
    #        checkArray() accordingly)
    assert_equal(img[0,0,0], 1)
 
def testRGBImage():
    checkArray(arraytypes.RGBImage, 3, 2)
    
    cshape = (10, 20)
    shape = (20, 10)
    rshape = (3, 20, 10)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageRGBUnstrided", "testImageRGBStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
         
    checkCompatibility(arraytypes.RGBImage(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.RGBImage(shape, order='V'), c)

    checkCompatibility(arraytypes.RGBImage(shape, order='F'), c)
    
    img = arraytypes.RGBImage(cshape, order='C')
    assert_equal(vt.viewArray3Unstrided(img), rshape)
    assert (img[0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage(shape, order='V')
    assert_equal(vt.viewArray3Unstrided(img), rshape)
    assert (img[0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage(shape, order='F')
    assert_equal(vt.viewArray3Strided(img), rshape)
    assert (img[:,0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage(cshape, order='C')
    assert_equal(vt.viewImageRGBUnstrided(img), shape)
    assert (img[0,0]==(1,1,1)).all()
    
    img = arraytypes.RGBImage(shape, order='V')
    assert_equal(vt.viewImageRGBStrided(img), shape)
    assert (img[0,0]==(1,1,1)).all()
 
def testVector2Image():
    checkArray(arraytypes.Vector2Image, 2, 2)
    
    cshape = (10, 20)
    shape = (20, 10)
    rshape = (2, 20, 10)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
         
    checkCompatibility(arraytypes.Vector2Image(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.Vector2Image(shape, order='V'), c)
    
    checkCompatibility(arraytypes.Vector2Image(shape, order='F'), c)
    
    img = arraytypes.Vector2Image(cshape, order='C')
    assert_equal(vt.viewArray3Unstrided(img), rshape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image(shape, order='V')
    assert_equal(vt.viewArray3Unstrided(img), rshape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image(shape, order='F')
    assert_equal(vt.viewArray3Strided(img), rshape)
    assert (img[:,0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image(cshape, order='C')
    assert_equal(vt.viewImageVector2Unstrided(img), shape)
    assert (img[0,0]==(1,1)).all()
    
    img = arraytypes.Vector2Image(shape, order='V')
    assert_equal(vt.viewImageVector2Unstrided(img), shape)
    assert (img[0,0]==(1,1)).all()
 
def testVector3Image():
    checkArray(arraytypes.Vector3Image, 3, 2)
 
def testVector4Image():
    checkArray(arraytypes.Vector4Image, 4, 2)
 
def testVolume1():
    checkArray(arraytypes.Volume, 1, 3)
    
    shape = (20, 10, 5)
    rshape = (5, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
         
    checkCompatibility(arraytypes.Volume(rshape, order='C'), c)
    
    checkCompatibility(arraytypes.Volume(shape, order='V'), c)
    
    checkCompatibility(arraytypes.Volume(shape, order='F'), c)
    
    vol = arraytypes.Volume(rshape, order='C')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)
    
    vol = arraytypes.Volume(shape, order='V')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)
    
    vol = arraytypes.Volume(shape, order='F')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[:,0,0,0], 1)
 
def testVolume2():
    checkArray(arraytypes.Volume, 2, 3)
    
    shape = (5, 10, 20, 2)
    cshape = (20, 10, 5, 2)
    fshape = (2, 5, 10, 20)
        
    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeVector2Unstrided", "testVolumeVector2Strided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
         
    checkCompatibility(arraytypes.Volume(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.Volume(shape, order='V'), c)
    
    checkCompatibility(arraytypes.Volume(fshape, order='F'), c)
    
    vol = arraytypes.Volume(cshape, order='C')
    assert_equal(vt.viewArray4Unstrided(vol), fshape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume(shape, order='V')
    assert_equal(vt.viewArray4Unstrided(vol), fshape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume(fshape, order='F')
    assert_equal(vt.viewArray4Strided(vol), fshape)
    assert (vol[:,0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume(cshape, order='C')
    assert_equal(vt.viewVolumeVector2Unstrided(vol), shape[:-1])
    assert (vol[0,0,0]==(1,1)).all()
    
    vol = arraytypes.Volume(shape, order='V')
    assert_equal(vt.viewVolumeVector2Unstrided(vol), shape[:-1])
    assert (vol[0,0,0]==(1,1)).all()
 
def testScalarVolume():
    checkArray(arraytypes.ScalarVolume, 1, 3)
    
    cshape = (20, 10, 5)
    shape = (5, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
         
    checkCompatibility(arraytypes.ScalarVolume(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.ScalarVolume(shape, order='V'), c)
    
    checkCompatibility(arraytypes.ScalarVolume(shape, order='F'), c)
    
    vol = arraytypes.ScalarVolume(cshape, order='C')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)
    
    vol = arraytypes.ScalarVolume(shape, order='V')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)
 
def testRGBVolume():
    checkArray(arraytypes.RGBVolume, 3, 3)
    
    cshape = (20, 10, 5)
    shape = (5, 10, 20)
    rshape = (3, 5, 10, 20)

    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeRGBUnstrided", "testVolumeRGBStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.RGBVolume(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.RGBVolume(shape, order='V'), c)
    
    checkCompatibility(arraytypes.RGBVolume(shape, order='F'), c)
    
    vol = arraytypes.RGBVolume(cshape, order='C')
    checkShape(vt.viewArray4Unstrided(vol), rshape)
    assert (vol[0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume(shape, order='V')
    checkShape(vt.viewArray4Unstrided(vol), rshape)
    assert (vol[0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume(shape, order='F')
    checkShape(vt.viewArray4Strided(vol), rshape)
    assert (vol[:,0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume(cshape, order='C')
    checkShape(vt.viewVolumeRGBUnstrided(vol), shape)
    assert (vol[0,0,0]==(1,1,1)).all()
    
    vol = arraytypes.RGBVolume(shape, order='V')
    checkShape(vt.viewVolumeRGBUnstrided(vol), shape)
    assert (vol[0,0,0]==(1,1,1)).all()
 
def testVector2Volume():
    checkArray(arraytypes.Vector2Volume, 2, 3)
    
    cshape = (20, 10, 5)
    shape = (5, 10, 20)
    rshape = (2, 5, 10, 20)

    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeVector2Unstrided", "testVolumeVector2Strided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
         
    checkCompatibility(arraytypes.Vector2Volume(cshape, order='C'), c)
    
    checkCompatibility(arraytypes.Vector2Volume(shape, order='V'), c)
    
    checkCompatibility(arraytypes.Vector2Volume(shape, order='F'), c)
    
    vol = arraytypes.Vector2Volume(cshape, order='C')
    checkShape(vt.viewArray4Unstrided(vol), rshape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume(shape, order='V')
    checkShape(vt.viewArray4Unstrided(vol), rshape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume(shape, order='F')
    checkShape(vt.viewArray4Strided(vol), rshape)
    assert (vol[:,0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume(cshape, order='C')
    checkShape(vt.viewVolumeVector2Unstrided(vol), shape)
    assert (vol[0,0,0]==(1,1)).all()
    
    vol = arraytypes.Vector2Volume(shape, order='V')
    checkShape(vt.viewVolumeVector2Unstrided(vol), shape)
    assert (vol[0,0,0]==(1,1)).all()
 
def testVector3Volume():
    checkArray(arraytypes.Vector3Volume, 3, 3)
 
def testVector4Volume():
    checkArray(arraytypes.Vector4Volume, 4, 3)
 
def testDeepcopy():
        a = arraytypes.Image(numpy.random.random((10, 4)), order='C')
        b = copy.deepcopy(a)
        assert numpy.all(a == b)
        assert_equal(b.flags.c_contiguous, a.flags.c_contiguous)
        assert_equal(b.flags.f_contiguous, a.flags.f_contiguous)
        a[0,0] = 42
        assert b[0,0] != 42

def testDeepcopyWithAttributes():
        a = arraytypes.Image((320, 200), order='C')
        a.myCustomAttribute = 42
        b = copy.deepcopy(a)
        assert hasattr(b, "myCustomAttribute")
        assert_equal(b.myCustomAttribute, 42)

def testDeepcopyWithCyclicReference():
        a = arraytypes.Image((320, 200), order='C')
        b = arraytypes.Image((320, 200), order='C')
        a.myCustomAttribute = b
        b.backLink = a
        c = copy.deepcopy(a)
        assert hasattr(c, "myCustomAttribute")
        assert c.myCustomAttribute.backLink is c

def testUfuncs():
    from numpy import bool, int8, uint8, int16, uint16, int32, uint32, int64, uint64
    from numpy import float32, float64, longdouble, complex64, complex128, clongdouble
    integers = [ int8, uint8, int16, uint16, int32, uint32, int64, uint64]
    floats = [float32, float64, longdouble]
    compl = [complex64, complex128, clongdouble]
    types = integers + floats + compl
    
    arrays, ones  = {}, {}
    for t in types:
        arrays[t] = arraytypes.ScalarImage((2,2), t, value=2)
        ones[t] = arraytypes.ScalarImage((1,1), t, value=1)
    for t, a in arrays.iteritems():
        b = -a
        assert_equal(t, b.dtype)
        if not numpyHasComplexNegateBug or t is not clongdouble:
            assert (b == -t(2)).all()
        b = a + a
        assert_equal(t, b.dtype)
        assert (b == t(4)).all()
        b = a == a
        assert_equal(bool, b.dtype)
        assert (b == True).all()
        b = ones[t] + a
        assert_equal(a.shape, b.shape)
        assert_equal(t, b.dtype)
        assert (b == t(3)).all()
        b = a + ones[t]
        assert_equal(a.shape, b.shape)
        assert_equal(t, b.dtype)
        assert (b == t(3)).all()
        b = 3 + a
        assert_equal(t, b.dtype)
        assert (b == t(5)).all()
        b = a + 3
        assert_equal(t, b.dtype)
        assert (b == t(5)).all()
    for i in integers:
        a1 = arrays[i]
        for j in integers:
            if i == j:
                continue
            a2 = arrays[j]
            b = a1 * a2
            if a1.dtype.itemsize < 8 and a2.dtype.itemsize < 8:
                assert_equal(int32, b.dtype)
            else:
                assert_equal(int64, b.dtype)
            assert (b == 4).all()
            b = a1 <= a2
            assert_equal(bool, b.dtype)
            assert (b == True).all()
        for j in floats + compl:
            a2 = arrays[j]
            b = a1 * a2
            assert_equal(j, b.dtype)
            assert (b == 4).all()
            b = a2 * a1
            assert_equal(j, b.dtype)
            assert (b == 4).all()
            b = a1 >= a2
            assert_equal(bool, b.dtype)
            assert (b == True).all()
            b = a2 > a1
            assert_equal(bool, b.dtype)
            assert (b == False).all()
        b = a1 + 1
        assert (b == 3).all()
        assert_equal(a1.dtype, b.dtype)
        b = a1 + 1.0
        assert (b == 3.0).all()
        if a1.dtype.itemsize < 8:
            assert_equal(float32, b.dtype)
        else:
            assert_equal(float64, b.dtype)
    for i in floats:
        a1 = arrays[i]
        for j in compl:
            if i == j:
                continue
            a2 = arrays[j]
            b = a1 * a2
            assert_equal(j, b.dtype)
            assert (b == 4).all()
            b = a2 * a1
            assert_equal(j, b.dtype)
            assert (b == 4).all()
            b = a1 >= a2
            assert_equal(bool, b.dtype)
            assert (b == True).all()
            b = a2 > a1
            assert_equal(bool, b.dtype)
            assert (b == False).all()
        b = a1 + 1.0
        assert (b == 3.0).all()
        assert_equal(a1.dtype, b.dtype)
    assert_equal(float64, (arrays[float32]+arrays[float64]).dtype)
    assert_equal(float64, (arrays[float64]+arrays[float32]).dtype)
    assert_equal(longdouble, (arrays[float32]+arrays[longdouble]).dtype)
    assert_equal(longdouble, (arrays[longdouble]+arrays[float32]).dtype)
    assert_equal(longdouble, (arrays[float64]+arrays[longdouble]).dtype)
    assert_equal(longdouble, (arrays[longdouble]+arrays[float64]).dtype)

    assert_equal(complex128, (arrays[complex64]+arrays[complex128]).dtype)
    assert_equal(complex128, (arrays[complex128]+arrays[complex64]).dtype)
    assert_equal(clongdouble, (arrays[complex64]+arrays[clongdouble]).dtype)
    assert_equal(clongdouble, (arrays[clongdouble]+arrays[complex64]).dtype)
    assert_equal(clongdouble, (arrays[complex128]+arrays[clongdouble]).dtype)
    assert_equal(clongdouble, (arrays[clongdouble]+arrays[complex128]).dtype)

    b = abs(arrays[complex64])
    assert (b == 2.0).all()
    assert_equal(float32, b.dtype)
    b = abs(arrays[complex128])
    assert (b == 2.0).all()
    assert_equal(float64, b.dtype)
    b = abs(arrays[clongdouble])
    assert (b == 2.0).all()
    assert_equal(longdouble, b.dtype)

    a = arraytypes.ScalarImage((2,2), uint8, value=255)
    b = a + a
    assert (b == 254).all()
    b = a + 255
    assert (b == 254).all()
    b = 255 + a
    assert (b == 254).all()
    b = arraytypes.ScalarImage((2,2), int32, value=0)
    ufunc.add(a, a, b)
    assert (b == 510).all()
    b = arraytypes.ScalarImage((2,2), int32, value=0)
    ufunc.add(a, 255, b)
    assert (b == 510).all()
    b = arraytypes.ScalarImage((2,2), int32, value=0)
    ufunc.add(255, a, b)
    assert (b == 510).all()
