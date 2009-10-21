# run with a simple 'nosetests' in this directory
# (and nose installed, i.e. 'easy_install nose')

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

import arraytypes, ufunc, numpy, copy
import vigranumpytest as vt
from nose.tools import assert_equal, raises

numpyHasComplexNegateBug = numpy.version.version.startswith('1.0')

try:
    vt.testAny()
except Exception, e:
    ArgumentError = type(e)
    
allTests = set()
for n, f in vt.__dict__.items():
    if n.startswith('test'):
        allTests.add(n)
 
def checkArray(cls, channels, dim):
        def testCopy(img):
            b = cls(img, order='A')
            assert_equal(sys.getrefcount(b), 2)
            assert_equal(b.shape, img.shape)
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
        if channels > 1:
            d = dim+1
            rshape = shape[:d]
            fcstrides = (200*channels, 20*channels, 4*channels, 4)[-d:]
            bcstrides = (50*channels, 5*channels, channels, 1)[-d:]
            ffstrides = (4, shape[0]*4, shape[0]*shape[1]*4, shape[0]*shape[1]*shape[2]*4)[:d]
            bfstrides = (1, shape[0], shape[0]*shape[1], shape[0]*shape[1]*shape[2])[:d]
            fvstrides = (channels*4, channels*shape[0]*4, channels*shape[0]*shape[1]*4)[:d-1]+(4,)
            bvstrides = (channels, channels*shape[0], channels*shape[0]*shape[1])[:d-1]+(1,)
        else:
            d = dim
            rshape = shape[:d]
            fcstrides = (200, 20, 4)[-d:]
            bcstrides = (50, 5, 1)[-d:]
            ffstrides = (4, shape[0]*4, shape[0]*shape[1]*4)[:d]
            bfstrides = (1, shape[0], shape[0]*shape[1])[:d]
            fvstrides = ffstrides
            bvstrides = bfstrides
            
        if cls.channels > 0 or channels == 1:
            shape = shape[:-1]
        
        value = 1 if channels == 1 else range(1,channels+1)

        # test type
        img = cls(shape, order="V")
        assert type(img) is cls
        assert isinstance(img, numpy.ndarray)
        assert_equal(img.dtype, numpy.float32)
        assert_equal(sys.getrefcount(img), 2)
        
        # test shape
        assert_equal(img.shape, rshape)
        assert_equal(img.width, rshape[0])
        assert_equal(img.height, rshape[1])
        if dim == 3:
            assert_equal(img.depth, rshape[2])
        assert_equal(img.channels, channels)
        assert_equal(img.spatialDimensions, dim)
       
        # test strides and order
        assert_equal(img.strides, fvstrides)
        assert_equal(img.order, "F" if channels == 1 else "V")
        assert not img.flags.c_contiguous
        assert_equal(img.flags.f_contiguous, True if channels == 1 else False)

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
        assert_equal(img.shape, (img+img).shape)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.shape, (img*2).shape)
        assert_equal(img.strides, (img*2).strides)
        
        if cls.channels > 0 or channels == 1:
            img = cls(shape + (channels,))
            assert_equal(img.shape, rshape)
            assert_equal(img.strides, fvstrides)

        # test shape, strides, and copy for 'F' order
        img = cls(shape, order='F')
        assert_equal(sys.getrefcount(img), 2)
        assert_equal(img.shape, rshape)
        assert_equal(img.strides, ffstrides)
        assert_equal(img.order, "F")
        assert not img.flags.c_contiguous
        assert img.flags.f_contiguous
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)

        # test shape, strides, and copy for 'A' order (should be equal to 'V' order)
        img = cls(shape, order='A')
        assert_equal(sys.getrefcount(img), 2)
        assert_equal(img.shape, rshape)
        assert_equal(img.strides, fvstrides)
        assert_equal(img.order, "F" if channels == 1 else "V")
        assert not img.flags.c_contiguous
        assert_equal(img.flags.f_contiguous, True if channels == 1 else False)
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)

        # test shape, strides, and copy for 'C' order
        img = cls(shape, order='C')
        assert_equal(sys.getrefcount(img), 2)
        assert_equal(img.shape, rshape)
        assert_equal(img.strides, fcstrides)
        assert_equal(img.order, "C")
        assert img.flags.c_contiguous
        assert not img.flags.f_contiguous
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)

        value = 10 if channels == 1 else range(10,channels+10)
        zero = 0 if channels == 1 else (0,)*channels
        
        # test shape, strides, and copy for dtype uint8
        b = cls(img, dtype=numpy.uint8, order='V')
        assert_equal(sys.getrefcount(b), 2)
        assert_equal(b.dtype, numpy.uint8)
        assert_equal(b.shape, rshape)
        assert_equal(b.strides, bvstrides)
        assert_equal(b.order, "F" if channels == 1 else "V")
        assert not b.flags.c_contiguous
        assert_equal(b.flags.f_contiguous, True if channels == 1 else False)
        assert (b==img).all()
        b[2,1] = value
        assert (b[2,1]==value).all()
        assert (img[2,1]==zero).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)

        b = cls(img, dtype=numpy.uint8, order='C')
        assert_equal(sys.getrefcount(b), 2)
        assert_equal(b.shape, rshape)
        assert_equal(b.strides, bcstrides)
        assert_equal(b.order, "C")
        assert b.flags.c_contiguous
        assert not b.flags.f_contiguous
        assert (b==img).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)
        
        b = cls(img, dtype=numpy.uint8, order='F')
        assert_equal(sys.getrefcount(b), 2)
        assert_equal(b.shape, rshape)
        assert_equal(b.strides, bfstrides)
        assert_equal(b.order, "F")
        assert not b.flags.c_contiguous
        assert b.flags.f_contiguous
        assert (b==img).all()
        assert_equal(b.strides, (-b).strides)
        assert_equal(b.strides, (b+b).strides)
        assert_equal(b.strides, (b*2).strides)
        
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
    raise AssertionError, "%r did not throw ArgumentError as expected when passed a %r with shape %s and stride %s" % (n, type(obj), str(obj.shape), str(obj.strides))

def checkCompatibility(obj, compatible):
    for n in compatible:
        f = getattr(vt, n)
        assert_equal(obj.shape, f(obj))
        
    incompatible = allTests.difference(compatible)
    
    for n in incompatible:
        checkFailure(obj, n)

def testImage1():
    checkArray(arraytypes.Image, 1, 2)
    
    c = ["testAny",
         "testArray2Strided",
         "testArray3Strided",
         "testImageSinglebandStrided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10), order='C'), c)
    
    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10), order='F'), c)
    
    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10), order='V'), c)
    
    img = arraytypes.Image((20, 10), order='F')
    assert_equal(vt.viewArray2Unstrided(img), img.shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.Image((20, 10), order='V')
    assert_equal(vt.viewArray2Unstrided(img), img.shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.Image((20, 10), order='C')
    assert_equal(vt.viewArray2Strided(img), img.shape)
    assert_equal(img[0,0], 1)
 
def testImage2():
    checkArray(arraytypes.Image, 2, 2)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageVector2Strided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10, 2), order='C'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10, 2), order='F'), c)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Image((20, 10, 2), order='V'), c)
    
    img = arraytypes.Image((20, 10, 2), order='F')
    assert_equal(vt.viewArray3Unstrided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Image((20, 10, 2), order='V')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Image((20, 10, 2), order='C')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Image((20, 10, 2), order='V')
    assert_equal(vt.viewImageVector2Unstrided(img), (20, 10))
    assert (img[0,0]==(1,1)).all()
    
    img = arraytypes.Image((20, 10, 2), order='C')
    assert_equal(vt.viewImageVector2Strided(img), (20,10))
    assert (img[0,0]==(1,1)).all()
 
def testScalarImage():
    checkArray(arraytypes.ScalarImage, 1, 2)
    
    c = ["testAny",
         "testArray2Strided",
         "testArray3Strided",
         "testImageSinglebandStrided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.ScalarImage((20, 10), order='C'), c)
    
    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.ScalarImage((20, 10), order='F'), c)
    
    checkCompatibility(arraytypes.ScalarImage((20, 10), order='V'), c)
    
    img = arraytypes.ScalarImage((20, 10), order='F')
    assert_equal(vt.viewArray2Unstrided(img), img.shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.ScalarImage((20, 10), order='V')
    assert_equal(vt.viewArray2Unstrided(img), img.shape)
    assert_equal(img[0,0], 1)
    
    img = arraytypes.ScalarImage((20, 10), order='C')
    assert_equal(vt.viewArray2Strided(img), img.shape)
    assert_equal(img[0,0], 1)
 
def testRGBImage():
    checkArray(arraytypes.RGBImage, 3, 2)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageRGBStrided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.RGBImage((20, 10), order='C'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.RGBImage((20, 10), order='F'), c)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageRGBUnstrided", "testImageRGBStrided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.RGBImage((20, 10), order='V'), c)
    
    img = arraytypes.RGBImage((20, 10), order='F')
    assert_equal(vt.viewArray3Unstrided(img), img.shape)
    assert (img[0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage((20, 10), order='V')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage((20, 10), order='C')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0,0)).all()
    
    img = arraytypes.RGBImage((20, 10), order='V')
    assert_equal(vt.viewImageRGBUnstrided(img), (20, 10))
    assert (img[0,0]==(1,1,1)).all()
    
    img = arraytypes.RGBImage((20, 10), order='C')
    assert_equal(vt.viewImageRGBStrided(img), (20,10))
    assert (img[0,0]==(1,1,1)).all()
 
def testVector2Image():
    checkArray(arraytypes.Vector2Image, 2, 2)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageVector2Strided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Image((20, 10), order='C'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Image((20, 10), order='F'), c)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Image((20, 10), order='V'), c)
    
    img = arraytypes.Vector2Image((20, 10), order='F')
    assert_equal(vt.viewArray3Unstrided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image((20, 10), order='V')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image((20, 10), order='C')
    assert_equal(vt.viewArray3Strided(img), img.shape)
    assert (img[0,0]==(1,0)).all()
    
    img = arraytypes.Vector2Image((20, 10), order='V')
    assert_equal(vt.viewImageVector2Unstrided(img), (20, 10))
    assert (img[0,0]==(1,1)).all()
    
    img = arraytypes.Vector2Image((20, 10), order='C')
    assert_equal(vt.viewImageVector2Strided(img), (20,10))
    assert (img[0,0]==(1,1)).all()
 
def testVector3Image():
    checkArray(arraytypes.Vector3Image, 3, 2)
 
def testVector4Image():
    checkArray(arraytypes.Vector4Image, 4, 2)
 
def testVolume1():
    checkArray(arraytypes.Volume, 1, 3)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5), order='C'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5), order='F'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5), order='V'), c)
    
    vol = arraytypes.Volume((20, 10, 5), order='V')
    assert_equal(vt.viewArray3Unstrided(vol), vol.shape)
    assert_equal(vol[0,0,0], 1)
    
    vol = arraytypes.Volume((20, 10, 5), order='C')
    assert_equal(vt.viewArray3Strided(vol), vol.shape)
    assert_equal(vol[0,0,0], 1)
 
def testVolume2():
    checkArray(arraytypes.Volume, 2, 3)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5, 2), order='C'), c)
    
    c = ["testAny",
         "testArray4Unstrided", 
		 "testArray4Strided",
         "testVolumeMultibandUnstrided", 
		 "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5, 2), order='F'), c)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Unstrided", 
		 "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Volume((20, 10, 5, 2), order='V'), c)
    
    vol = arraytypes.Volume((20, 10, 5, 2), order='F')
    assert_equal(vt.viewArray4Unstrided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume((20, 10, 5, 2), order='V')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume((20, 10, 5, 2), order='C')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Volume((20, 10, 5, 2), order='V')
    assert_equal(vt.viewVolumeVector2Unstrided(vol), (20, 10, 5))
    assert (vol[0,0,0]==(1,1)).all()
    
    vol = arraytypes.Volume((20, 10, 5, 2), order='C')
    assert_equal(vt.viewVolumeVector2Strided(vol), (20, 10, 5))
    assert (vol[0,0,0]==(1,1)).all()
 
def testScalarVolume():
    checkArray(arraytypes.ScalarVolume, 1, 3)
    
    c = ["testAny",
         "testArray3Strided",
         "testArray4Strided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.ScalarVolume((20, 10, 5), order='C'), c)
    
    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.ScalarVolume((20, 10, 5), order='F'), c)
    
    checkCompatibility(arraytypes.ScalarVolume((20, 10, 5), order='V'), c)
    
    vol = arraytypes.ScalarVolume((20, 10, 5), order='V')
    assert_equal(vt.viewArray3Unstrided(vol), vol.shape)
    assert_equal(vol[0,0,0], 1)
    
    vol = arraytypes.ScalarVolume((20, 10, 5), order='C')
    assert_equal(vt.viewArray3Strided(vol), vol.shape)
    assert_equal(vol[0,0,0], 1)
 
def testRGBVolume():
    checkArray(arraytypes.RGBVolume, 3, 3)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeRGBStrided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.RGBVolume((20, 10, 5), order='C'), c)
    
    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.RGBVolume((20, 10, 5), order='F'), c)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeRGBUnstrided", "testVolumeRGBStrided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.RGBVolume((20, 10, 5), order='V'), c)
    
    vol = arraytypes.RGBVolume((20, 10, 5), order='F')
    assert_equal(vt.viewArray4Unstrided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume((20, 10, 5), order='V')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume((20, 10, 5), order='C')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0,0)).all()
    
    vol = arraytypes.RGBVolume((20, 10, 5), order='V')
    assert_equal(vt.viewVolumeRGBUnstrided(vol), (20, 10, 5))
    assert (vol[0,0,0]==(1,1,1)).all()
    
    vol = arraytypes.RGBVolume((20, 10, 5), order='C')
    assert_equal(vt.viewVolumeRGBStrided(vol), (20, 10, 5))
    assert (vol[0,0,0]==(1,1,1)).all()
 
def testVector2Volume():
    checkArray(arraytypes.Vector2Volume, 2, 3)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Volume((20, 10, 5), order='C'), c)
    
    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Volume((20, 10, 5), order='F'), c)
    
    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Unstrided", "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]
    checkCompatibility(arraytypes.Vector2Volume((20, 10, 5), order='V'), c)
    
    vol = arraytypes.Vector2Volume((20, 10, 5), order='F')
    assert_equal(vt.viewArray4Unstrided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume((20, 10, 5), order='V')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume((20, 10, 5), order='C')
    assert_equal(vt.viewArray4Strided(vol), vol.shape)
    assert (vol[0,0,0]==(1,0)).all()
    
    vol = arraytypes.Vector2Volume((20, 10, 5), order='V')
    assert_equal(vt.viewVolumeVector2Unstrided(vol), (20, 10, 5))
    assert (vol[0,0,0]==(1,1)).all()
    
    vol = arraytypes.Vector2Volume((20, 10, 5), order='C')
    assert_equal(vt.viewVolumeVector2Strided(vol), (20, 10, 5))
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
