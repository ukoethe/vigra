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

# import vigra  # FIXME: without this line, C++ constructors don't find VigraArray
import vigra.arraytypes as arraytypes
import vigra.ufunc as ufunc
import numpy, copy
import vigranumpytest as vt
from nose.tools import assert_equal, raises, assert_true

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

def computeFStrides(shape):
    return tuple(numpy.cumprod((1,)+shape[:-1]))

def computeCStrides(shape):
    return tuple(reversed(computeFStrides(tuple(reversed(shape)))))

def computeVStrides(shape, hasChannelAxis):
    if not hasChannelAxis:
        return computeFStrides(shape)
    stride = computeFStrides(shape[-1:] + shape[:-1])
    return stride[1:] + stride[0:1]

def checkArray(cls, channels, dim, hasChannelAxis=True):
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

        shape = (channels, 5, 10, 20)
        axistags = [AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z]
        axistags3 = AxisTags(AxisInfo.y, AxisInfo.z, AxisInfo.x)
        axistags4 = AxisTags(AxisInfo.y, AxisInfo.z, AxisInfo.x, AxisInfo.c)
        axistags5 = AxisTags(AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z, AxisInfo.t)

        # figure out expected strides and axistags
        s = 0 if hasChannelAxis else 1
        d = dim + 1
        fshape = shape[s:d]
        fstrides = computeFStrides(fshape)
        faxistags = AxisTags(axistags[s:d])

        cshape = tuple(reversed(fshape))
        cstrides = computeCStrides(cshape)
        caxistags = AxisTags(list(reversed(faxistags)))

        vshape = fshape[1-s:d] + fshape[:1-s]
        vstrides = computeVStrides(vshape, hasChannelAxis)
        vaxistags = AxisTags(axistags[1:d] + axistags[:1-s])

        fcstrides = tuple([k*4 for k in cstrides])
        fvstrides = tuple([k*4 for k in vstrides])
        ffstrides = tuple([k*4 for k in fstrides])

        value = 1 if channels == 1 else range(1,channels+1)

        # test type
        img = cls(vshape, order="V")
#        assert type(img) is cls
        assert isinstance(img, numpy.ndarray)
        assert_equal(img.dtype, numpy.float32)
        assert_equal(sys.getrefcount(img), 2)

        # test shape
        checkShape(img.shape, vshape)
        assert_equal(img.width, vshape[0])
        assert_equal(img.height, vshape[1])
        if dim == 3:
            assert_equal(img.depth, vshape[2])
        assert_equal(img.channels, channels)
        assert_equal(img.spatialDimensions, dim)

        # test strides and order
        checkStride(img.strides, fvstrides)
        if channels > 1:
            assert_equal(img.order, "V" if hasChannelAxis else "F")
        else:
            assert_true(img.order in ['V', 'F'])
        assert_equal(img.flags.c_contiguous, False)

        # test axistags
        assert_equal(img.axistags, vaxistags)
        assert_equal(img.view5D('F').axistags, axistags5)
        assert_equal(img.withAxes('y', 'z', 'x', 'c').axistags, axistags4)
        assert_equal(img.withAxes('yzxc').axistags, axistags4)
        assert_equal(img.withAxes(axistags4).axistags, axistags4)
        assert_true(img.withAxes(img.axistags) is img)
        array = img.noTags()
        assert_equal(type(array), numpy.ndarray)
        assert_equal(arraytypes.taggedView(array, vaxistags).axistags, vaxistags)
        assert_equal(arraytypes.taggedView(array, vaxistags.keys()).axistags, vaxistags)
        assert_equal(arraytypes.taggedView(array, ''.join(vaxistags.keys())).axistags, vaxistags)
        if img.channels == 1:
            assert_equal(img.withAxes('y', 'z', 'x').axistags, axistags3)
        else:
            try:
                img.withAxes('y', 'z', 'x')
                raise AssertionError, "img.withAxes() failed to throw on non-singleton channel."
            except RuntimeError:
                pass
        # FIXME: add more tests

        # test initialization and assignment
        assert_equal(img.min(), 0.0)
        assert_equal(img.max(), 0.0)
        img = cls(vshape, order="V", value=99.0)
        assert_equal(img.min(), 99.0)
        assert_equal(img.max(), 99.0)
        img.flat[:] = range(img.size)
        assert_equal(img.flatten().tolist(), range(img.size))
        img[1,2] = value
        assert_equal((img[1,2]==value).all(), True)

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

        # test shape, strides, and copy for 'F' order
        img = cls(fshape, order='F')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, fshape)
        checkStride(img.strides, ffstrides)
        assert_equal(img.axistags, faxistags)
        assert_equal(img.order, "F")
        assert_equal(img.flags.c_contiguous, False)
        assert_equal(img.flags.f_contiguous, True)

        assert_equal(img.min(), 0.0)
        assert_equal(img.max(), 0.0)
        img = cls(fshape, order="F", value=99.0)
        assert_equal(img.min(), 99.0)
        assert_equal(img.max(), 99.0)

        if dim == 2:
            img[...,1,2] = value
        else:
            img[...,0,1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)
        assert_equal(img.view5D('F').axistags, axistags5)
        assert_equal(img.withAxes('y', 'z', 'x', 'c').axistags, axistags4)

        # test shape, strides, and copy for 'A' order (should be equal to 'V' order)
        img = cls(vshape, order='A')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, vshape)
        checkStride(img.strides, fvstrides)
        if channels > 1:
            assert_equal(img.order, "V" if hasChannelAxis else "F")
        else:
            assert_true(img.order in ['V', 'F'])
        assert_equal(img.flags.c_contiguous, False)
        assert_equal(img.axistags, vaxistags)
        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)

        # test shape, strides, and copy for 'C' order
        img = cls(cshape, order='C')
        assert_equal(sys.getrefcount(img), 2)
        checkShape(img.shape, cshape)
        checkStride(img.strides, fcstrides)
        assert_equal(img.axistags, caxistags)
        assert_equal(img.order, "C")
        assert_equal(img.flags.c_contiguous, True)
        assert_equal(img.flags.f_contiguous, False)

        assert_equal(img.min(), 0.0)
        assert_equal(img.max(), 0.0)
        img = cls(cshape, order="C", value=99.0)
        assert_equal(img.min(), 99.0)
        assert_equal(img.max(), 99.0)

        img[1,2] = value
        testCopy(img)
        assert_equal(img.strides, (-img).strides)
        assert_equal(img.strides, (img+img).strides)
        assert_equal(img.strides, (img*2).strides)
        assert_equal(img.axistags, (-img).axistags)
        assert_equal(img.axistags, (img+img).axistags)
        assert_equal(img.axistags, (img*2).axistags)
        assert_equal(img.view5D('F').axistags, axistags5)
        assert_equal(img.withAxes('y', 'z', 'x', 'c').axistags, axistags4)

        value = 10 if channels == 1 else range(10,channels+10)
        zero = 0 if channels == 1 else (0,)*channels

        # test shape, strides, and copy for dtype uint8
        img = cls(vshape, order="V")
        b = cls(img, dtype=numpy.uint8, order='V')
        assert_equal(sys.getrefcount(b), 2)
        assert_equal(b.dtype, numpy.uint8)
        checkShape(b.shape, img.shape)
        checkStride(b.strides, computeVStrides(b.shape, hasChannelAxis))
        if channels > 1:
            assert_equal(img.order, "V" if hasChannelAxis else "F")
        else:
            assert_true(img.order in ['V', 'F'])
        assert_equal(b.axistags, img.axistags)
        assert_equal(b.flags.c_contiguous, False)
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

        img = cls(cshape, order="C")
        b = cls(img, dtype=numpy.uint8, order='C')
        assert_equal(sys.getrefcount(b), 2)
        checkShape(b.shape, img.shape)
        checkStride(b.strides, computeCStrides(b.shape))
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

        img = cls(fshape, order="F")
        b = cls(img, dtype=numpy.uint8, order='F')
        assert_equal(sys.getrefcount(b), 2)
        checkShape(b.shape, img.shape)
        checkStride(b.strides, computeFStrides(b.shape))
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
        img = cls(vshape, order="V")
        v1 = img.view(numpy.ndarray)
        v2 = img.view(numpy.ndarray)
        assert type(v1) is numpy.ndarray
        assert (v1==v2).all()
        assert numpy.may_share_memory(v1, img)
        v1[3,4] = value
        assert (v2[3,4]==value).all()
        assert (v1==v2).all()

def checkFailure(obj, n):
    f = getattr(vt, n)
    try:
        f(obj)
    except ArgumentError:
        return
    raise AssertionError, "%r did not throw ArgumentError as expected when passed a %r with shape %s, stride %s, axistags '%s'" % (n, type(obj), str(obj.shape), str(obj.strides), repr(getattr(obj, "axistags", "none")))

def checkCompatibility(obj, compatible):
    for n in compatible:
        try:
            f = getattr(vt, n)
            shape, acopy, default_ordering, same_ordering = f(obj)


            assert_equal(obj.shape, shape)

            assert_equal(obj.__class__, acopy.__class__)
            assert_equal(obj.shape, acopy.shape)
            if hasattr(obj, 'axistags'):
                assert_equal(obj.axistags, acopy.axistags)
            else:
                assert(not hasattr(acopy, 'axistags'))

            if n != "testAny":
                assert_equal(obj.shape, same_ordering.shape)
                assert(obj.view(numpy.ndarray) == same_ordering.view(numpy.ndarray)).all()

                if not hasattr(obj, 'axistags'):
                    assert_equal(numpy.ndarray, same_ordering.__class__)
                    assert(not hasattr(same_ordering, 'axistags'))

                    if n.startswith("testArray"):
                        assert_equal(numpy.ndarray, default_ordering.__class__)
                        assert_equal(obj.shape, default_ordering.shape)
                        assert(obj.view(numpy.ndarray) == default_ordering.view(numpy.ndarray)).all()
                        assert(not hasattr(default_ordering, 'axistags'))
                    else:
                        assert_equal(arraytypes.VigraArray, default_ordering.__class__)
                        assert_equal(default_ordering.axistags,
                                     arraytypes.VigraArray.defaultAxistags(default_ordering.ndim))

                        if obj.ndim == default_ordering.ndim:
                            assert_equal(obj.shape, default_ordering.shape)
                            assert(obj.view(numpy.ndarray) == default_ordering.view(numpy.ndarray)).all()
                        else:
                            assert_equal(obj.shape + (1,), default_ordering.shape)
                            assert(obj.view(numpy.ndarray) == default_ordering[...,0].view(numpy.ndarray)).all()

                else:
                    assert_equal(arraytypes.VigraArray, same_ordering.__class__)
                    assert_equal(obj.axistags, same_ordering.axistags)

                    if n.startswith("testArray"):
                        assert_equal(numpy.ndarray, default_ordering.__class__)
                        fobj = obj.transposeToNormalOrder()
                        fshape = fobj.shape
                        assert_equal(fshape, default_ordering.shape)
                        assert(fobj.view(numpy.ndarray) == default_ordering.view(numpy.ndarray)).all()
                        assert(not hasattr(default_ordering, 'axistags'))
                    else:
                        assert_equal(arraytypes.VigraArray, default_ordering.__class__)
                        dobj = obj.transposeToOrder(arraytypes.VigraArray.defaultOrder)
                        dshape = dobj.shape
                        assert_equal(default_ordering.axistags,
                                     arraytypes.VigraArray.defaultAxistags(default_ordering.ndim))
                        if obj.ndim == default_ordering.ndim:
                            assert_equal(dshape, default_ordering.shape)
                            assert(dobj.view(numpy.ndarray) == default_ordering.view(numpy.ndarray)).all()
                        else:
                            assert_equal(dshape + (1,), default_ordering.shape)
                            assert(fobj.view(numpy.ndarray) == default_ordering[...,0].view(numpy.ndarray)).all()
        except Exception:
            print "exception in %s with shape %s strides %s tags (%s)" % (n, obj.shape, obj.strides,
                                            repr(getattr(obj, "axistags", "none")))
            raise

    incompatible = allTests.difference(compatible)

    for n in incompatible:
        try:
            checkFailure(obj, n)
        except Exception:
            print "exception in %s with shape %s strides %s tags (%s)" % (n, obj.shape, obj.strides,
                                            repr(getattr(obj, "axistags", "none")))
            raise

def testAxisTags():
    axistags = AxisTags(AxisInfo.c(description="RGB"),
                        AxisInfo.ft(3.0, "time frequency"),
                        AxisInfo.y(0.5),
                        AxisInfo.z(4, "confocal depth"))

    json = '''{
  "axes": [
    {
      "key": "c",
      "typeFlags": 1,
      "resolution": 0,
      "description": "RGB"
    },
    {
      "key": "t",
      "typeFlags": 24,
      "resolution": 3,
      "description": "time frequency"
    },
    {
      "key": "y",
      "typeFlags": 2,
      "resolution": 0.5,
      "description": ""
    },
    {
      "key": "z",
      "typeFlags": 2,
      "resolution": 4,
      "description": "confocal depth"
    }
  ]
}'''
    assert_equal(axistags.toJSON(), json)

    readBack = AxisTags.fromJSON(json)
    assert_equal(axistags, readBack)
    assert_equal(readBack[0].description, "RGB")
    assert_equal(readBack[1].description, "time frequency")
    assert_equal(readBack[2].description, "")
    assert_equal(readBack[3].description, "confocal depth")
    assert_equal(readBack[0].resolution, 0)
    assert_equal(readBack[1].resolution, 3)
    assert_equal(readBack[2].resolution, 0.5)
    assert_equal(readBack[3].resolution, 4)

    import pickle
    s = pickle.dumps(axistags)
    unpickled = pickle.loads(s)
    assert_equal(axistags, unpickled)
    assert_equal(unpickled[0].description, "RGB")
    assert_equal(unpickled[1].description, "time frequency")
    assert_equal(unpickled[2].description, "")
    assert_equal(unpickled[3].description, "confocal depth")
    assert_equal(unpickled[0].resolution, 0)
    assert_equal(unpickled[1].resolution, 3)
    assert_equal(unpickled[2].resolution, 0.5)
    assert_equal(unpickled[3].resolution, 4)

    # FIXME: add more tests here
    defaultTags = arraytypes.VigraArray.defaultAxistags('cxyt')
    assert_equal(defaultTags.permutationToOrder('A'), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationToOrder('F'), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationToOrder('C'), (3, 2, 1, 0))
    assert_equal(defaultTags.permutationToOrder('V'), (1, 2, 3, 0))
    assert_equal(defaultTags.permutationToNormalOrder(), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationToNumpyOrder(), (3, 2, 1, 0))
    assert_equal(defaultTags.permutationToVigraOrder(), (1, 2, 3, 0))
    assert_equal(defaultTags.permutationFromNormalOrder(), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationFromNumpyOrder(), (3, 2, 1, 0))
    assert_equal(defaultTags.permutationFromVigraOrder(), (3, 0, 1, 2))

    defaultTags = arraytypes.AxisTags(4)
    assert_equal(defaultTags.permutationToOrder('A'), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationToOrder('F'), (0, 1, 2, 3))
    assert_equal(defaultTags.permutationToOrder('C'), (3, 2, 1, 0))
    assert_equal(defaultTags.permutationToOrder('V'), (0, 1, 2, 3))

    assert_equal(arraytypes.VigraArray.defaultAxistags('cxyz'),
                 arraytypes.VigraArray.defaultAxistags(4, order='F'))
    assert_equal(arraytypes.VigraArray.defaultAxistags('zyxc'),
                 arraytypes.VigraArray.defaultAxistags(4, order='C'))
    assert_equal(arraytypes.VigraArray.defaultAxistags('xyzc'),
                 arraytypes.VigraArray.defaultAxistags(4, order='V'))
    assert_equal(arraytypes.VigraArray.defaultAxistags('xyzc'),
                 arraytypes.VigraArray.defaultAxistags(4, order='A'))

    for order in 'VCF':
        defaultTags = arraytypes.VigraArray.defaultAxistags(3, order=order)
        assert_equal(defaultTags.permutationToOrder(order), (0, 1, 2))
        assert (defaultTags.channelIndex == 0 if order == 'F' else 2)

        defaultTags.transpose(defaultTags.permutationToOrder('V'))
        assert_equal(defaultTags.permutationToVigraOrder(), (0, 1, 2))
        assert (defaultTags.channelIndex == 2)

        defaultTags = arraytypes.VigraArray.defaultAxistags(3, order=order, noChannels=True)
        assert_equal(defaultTags.permutationToOrder(order), (0, 1, 2))
        assert (defaultTags.channelIndex == 3)

def testImage1():
    checkArray(arraytypes.Image, 1, 2)

    shape = (10, 20)
    rshape = (20, 10)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]

    checkCompatibility(arraytypes.Image(rshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Image(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Image(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray2Strided",
         "testImageSinglebandStrided",
         "testImageMultibandStrided"]

    a = numpy.ndarray(rshape, dtype=numpy.float32)
    a[...] = 2
    checkCompatibility(a, c)

    img = arraytypes.Image(rshape, order='C').dropChannelAxis()
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)

    img = arraytypes.Image(shape, order='V').dropChannelAxis()
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)

    img = arraytypes.Image(shape, order='F').dropChannelAxis()
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)

def testImage2():
    checkArray(arraytypes.Image, 2, 2)

    shape = (10, 20, 2)
    cshape = (20, 10, 2)
    fshape = (2, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandStrided"]

    checkCompatibility(arraytypes.Image(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Image(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Image(fshape, order='F', value=2), c)

    c = ["testAny",
         "testArray3Strided",
         "testImageMultibandStrided",
         "testImageVector2Strided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]

    a = numpy.ndarray(cshape, dtype=numpy.float32)
    a[...] = 2
    checkCompatibility(a, c)

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
    checkArray(arraytypes.ScalarImage, 1, 2, False)

    shape = (10, 20)
    cshape = (20, 10)

    c = ["testAny",
         "testArray2Unstrided", "testArray2Strided",
         "testImageSinglebandUnstrided", "testImageSinglebandStrided",
         "testImageMultibandUnstrided", "testImageMultibandStrided"]

    checkCompatibility(arraytypes.ScalarImage(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.ScalarImage(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.ScalarImage(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray2Strided",
         "testImageSinglebandStrided",
         "testImageMultibandStrided"]

    checkCompatibility(arraytypes.ScalarImage(cshape, order='C', value=2).view(numpy.ndarray), c)

    img = arraytypes.ScalarImage(cshape, order='C')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)

    img = arraytypes.ScalarImage(shape, order='V')
    checkShape(vt.viewArray2Unstrided(img), shape)
    assert_equal(img[0,0], 1)

    img = arraytypes.ScalarImage(shape, order='F')
    checkShape(vt.viewArray2Strided(img), shape)
    assert_equal(img[0,0], 1)

def testRGBImage():
    checkArray(arraytypes.RGBImage, 3, 2)

    cshape = (20, 10)
    shape = (10, 20)
    rshape = (3, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageRGBUnstrided", "testImageRGBStrided",
         "testImageMultibandStrided"]

    checkCompatibility(arraytypes.RGBImage(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.RGBImage(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.RGBImage(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray3Strided",
         "testImageMultibandStrided",
         "testImageRGBStrided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.RGBImage(cshape, order='C', value=2).view(numpy.ndarray), c)

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

    cshape = (20, 10)
    shape = (10, 20)
    rshape = (2, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testImageVector2Unstrided", "testImageVector2Strided",
         "testImageMultibandStrided"]

    checkCompatibility(arraytypes.Vector2Image(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Vector2Image(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Vector2Image(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray3Strided",
         "testImageMultibandStrided",
         "testImageVector2Strided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.Vector2Image(cshape, order='C', value=2).view(numpy.ndarray), c)

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

    shape = (5, 10, 20)
    rshape = (20, 10, 5)

    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.Volume(rshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Volume(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Volume(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray3Strided",
         "testImageMultibandStrided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]

    a = numpy.ndarray(rshape, dtype=numpy.float32)
    a[...] = 2
    checkCompatibility(a, c)

    vol = arraytypes.Volume(rshape, order='C').dropChannelAxis()
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)

    vol = arraytypes.Volume(shape, order='V').dropChannelAxis()
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)

    vol = arraytypes.Volume(shape, order='F').dropChannelAxis()
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)

def testVolume2():
    checkArray(arraytypes.Volume, 2, 3)

    shape = (5, 10, 20, 2)
    cshape = (20, 10, 5, 2)
    fshape = (2, 5, 10, 20)

    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeVector2Unstrided", "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.Volume(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Volume(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Volume(fshape, order='F', value=2), c)

    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]

    a = numpy.ndarray(cshape, dtype=numpy.float32)
    a[...] = 2
    checkCompatibility(a, c)

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
    checkArray(arraytypes.ScalarVolume, 1, 3, False)

    cshape = (20, 10, 5)
    shape = (5, 10, 20)

    c = ["testAny",
         "testArray3Unstrided", "testArray3Strided",
         "testVolumeSinglebandUnstrided", "testVolumeSinglebandStrided",
         "testVolumeMultibandUnstrided", "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.ScalarVolume(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.ScalarVolume(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.ScalarVolume(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray3Strided",
         "testImageMultibandStrided",
         "testVolumeSinglebandStrided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.ScalarVolume(cshape, order='C', value=2).view(numpy.ndarray), c)

    vol = arraytypes.ScalarVolume(cshape, order='C')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)

    vol = arraytypes.ScalarVolume(shape, order='V')
    checkShape(vt.viewArray3Unstrided(vol), shape)
    assert_equal(vol[0,0,0], 1)

    vol = arraytypes.ScalarVolume(shape, order='F')
    checkShape(vt.viewArray3Strided(vol), shape)
    assert_equal(vol[0,0,0], 1)

def testRGBVolume():
    checkArray(arraytypes.RGBVolume, 3, 3)

    cshape = (20, 10, 5)
    shape = (5, 10, 20)
    rshape = (3, 5, 10, 20)

    c = ["testAny",
         "testArray4Unstrided", "testArray4Strided",
         "testVolumeRGBUnstrided", "testVolumeRGBStrided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.RGBVolume(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.RGBVolume(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.RGBVolume(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray4Strided",
         "testVolumeRGBStrided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.RGBVolume(cshape, order='C', value=2).view(numpy.ndarray), c)

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
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.Vector2Volume(cshape, order='C', value=2), c)

    checkCompatibility(arraytypes.Vector2Volume(shape, order='V', value=2), c)

    checkCompatibility(arraytypes.Vector2Volume(shape, order='F', value=2), c)

    c = ["testAny",
         "testArray4Strided",
         "testVolumeVector2Strided",
         "testVolumeMultibandStrided"]

    checkCompatibility(arraytypes.Vector2Volume(cshape, order='C', value=2).view(numpy.ndarray), c)

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

def testTaggedShape():
    a = arraytypes.makeAxistags(4)
    assert_equal(repr(a), 'x y z c')
    a = arraytypes.makeAxistags(4, order='C')
    assert_equal(repr(a), 'z y x c')
    a = arraytypes.makeAxistags(4, order='F')
    assert_equal(repr(a), 'c x y z')
    a = arraytypes.makeAxistags(4, order='C', noChannels=True)
    assert_equal(repr(a), 't z y x')
    a = arraytypes.makeAxistags(4, noChannels=True)
    assert_equal(repr(a), 'x y z t')

    aa = arraytypes.makeAxistags('xyc')
    a = arraytypes.makeAxistags(aa)
    assert_equal(repr(a), 'x y c')
    assert(a is not aa)
    a = arraytypes.makeAxistags(aa, order='C')
    assert_equal(repr(a), 'y x c')
    assert(a is not aa)
    a = arraytypes.makeAxistags(aa, order='F')
    assert_equal(repr(a), 'c x y')
    assert(a is not aa)

    a = arraytypes.makeAxistags(aa, noChannels=True)
    assert_equal(repr(a), 'x y')
    assert(a is not aa)
    a = arraytypes.makeAxistags(aa, order='V', noChannels=True)
    assert_equal(repr(a), 'x y')
    assert(a is not aa)
    a = arraytypes.makeAxistags(aa, order='C', noChannels=True)
    assert_equal(repr(a), 'y x')
    assert(a is not aa)
    a = arraytypes.makeAxistags(aa, order='F', noChannels=True)
    assert_equal(repr(a), 'x y')
    assert(a is not aa)

    a = arraytypes.makeAxistags('xyc', order='V')
    assert_equal(repr(a), 'x y c')
    a = arraytypes.makeAxistags('xyc', order='C')
    assert_equal(repr(a), 'y x c')
    a = arraytypes.makeAxistags('xyc', order='F')
    assert_equal(repr(a), 'c x y')
    a = arraytypes.makeAxistags('xyc', order='F', noChannels=True)
    assert_equal(repr(a), 'x y')
    a = arraytypes.makeAxistags('xyc', order='V', noChannels=True)
    assert_equal(repr(a), 'x y')
    a = arraytypes.makeAxistags('xyc', order='C', noChannels=True)
    assert_equal(repr(a), 'y x')

    a = arraytypes.Image((20,10))
    a.axistags.setChannelDescription("in")
    res = vt.checkTaggedShapeMultiband(a)

    assert_equal(res[0].shape, (20,10,1))
    assert_equal(res[0].axistags, a.axistags)
    assert_equal(res[0].axistags[2].description, "in")

    assert_equal(res[1].shape, (20,10,1))
    assert_equal(res[1].axistags, a.axistags)
    assert_equal(res[1].axistags[2].description, "res2")

    assert_equal(res[2].shape, (20,10,1))
    assert_equal(res[2].axistags, a.axistags)
    assert_equal(res[2].axistags[2].description, "res3")

    assert_equal(res[3].shape, (20,10,3))
    assert_equal(res[3].axistags, a.axistags)
    assert_equal(res[3].axistags[2].description, "res4")

    assert_equal(res[4].shape, (20,10,1))
    assert_equal(res[4].axistags, a.axistags)
    assert_equal(res[4].axistags[2].description, "res5")

    assert_equal(res[5].shape, (20,10,3))
    assert_equal(res[5].axistags, a.axistags)
    assert_equal(res[5].axistags[2].description, "res6")

    res = vt.checkTaggedShapeSingleband(a)

    assert_equal(res[0].shape, (20,10,1))
    assert_equal(res[0].axistags, a.axistags)
    assert_equal(res[0].axistags[2].description, "in")

    assert_equal(res[1].shape, (20,10,1))
    assert_equal(res[1].axistags, a.axistags)
    assert_equal(res[1].axistags[2].description, "res2")

    assert_equal(res[2].shape, (20,10,1))
    assert_equal(res[2].axistags, a.axistags)
    assert_equal(res[2].axistags[2].description, "res3")

    assert_equal(res[3].shape, (20,10,3))
    assert_equal(res[3].axistags, a.axistags)
    assert_equal(res[3].axistags[2].description, "res4")

    assert_equal(res[4].shape, (20,10,1))
    assert_equal(res[4].axistags, a.axistags)
    assert_equal(res[4].axistags[2].description, "res5")

    assert_equal(res[5].shape, (20,10,3))
    assert_equal(res[5].axistags, a.axistags)
    assert_equal(res[5].axistags[2].description, "res6")

    a = arraytypes.Image((20,10,2))
    a.axistags.setChannelDescription("in")
    res = vt.checkTaggedShapeMultiband(a)

    assert_equal(res[0].shape, (20,10,2))
    assert_equal(res[0].axistags, a.axistags)
    assert_equal(res[0].axistags[2].description, "in")

    assert_equal(res[1].shape, (20,10,2))
    assert_equal(res[1].axistags, a.axistags)
    assert_equal(res[1].axistags[2].description, "res2")

    assert_equal(res[2].shape, (20,10,1))
    assert_equal(res[2].axistags, a.axistags)
    assert_equal(res[2].axistags[2].description, "res3")

    assert_equal(res[3].shape, (20,10,3))
    assert_equal(res[3].axistags, a.axistags)
    assert_equal(res[3].axistags[2].description, "res4")

    assert_equal(res[4].shape, (20,10,1))
    assert_equal(res[4].axistags, a.axistags)
    assert_equal(res[4].axistags[2].description, "res5")

    assert_equal(res[5].shape, (20,10,3))
    assert_equal(res[5].axistags, a.axistags)
    assert_equal(res[5].axistags[2].description, "res6")

    a = arraytypes.ScalarImage((20,10))
    a.axistags.setChannelDescription("in")
    resaxistags = copy.copy(a.axistags)
    resaxistags.insertChannelAxis()

    res = vt.checkTaggedShapeMultiband(a)

    assert_equal(res[0].shape, (20,10))
    assert_equal(res[0].axistags, a.axistags)
    assert_equal(len(res[0].axistags), 2)

    assert_equal(res[1].shape, (20,10))
    assert_equal(res[1].axistags, a.axistags)
    assert_equal(len(res[1].axistags), 2)

    assert_equal(res[2].shape, (20,10))
    assert_equal(res[2].axistags, a.axistags)
    assert_equal(len(res[1].axistags), 2)

    assert_equal(res[3].shape, (20,10,3))
    assert_equal(res[3].axistags, resaxistags)
    assert_equal(res[3].axistags[2].description, "res4")

    assert_equal(res[4].shape, (20,10))
    assert_equal(res[4].axistags, a.axistags)
    assert_equal(len(res[4].axistags), 2)

    assert_equal(res[5].shape, (20,10,3))
    assert_equal(res[5].axistags, resaxistags)
    assert_equal(res[5].axistags[2].description, "res6")

    res = vt.checkTaggedShapeSingleband(a)

    assert_equal(res[0].shape, (20,10))
    assert_equal(res[0].axistags, a.axistags)
    assert_equal(len(res[0].axistags), 2)

    assert_equal(res[1].shape, (20,10))
    assert_equal(res[1].axistags, a.axistags)
    assert_equal(len(res[1].axistags), 2)

    assert_equal(res[2].shape, (20,10))
    assert_equal(res[2].axistags, a.axistags)
    assert_equal(len(res[1].axistags), 2)

    assert_equal(res[3].shape, (20,10,3))
    assert_equal(res[3].axistags, resaxistags)
    assert_equal(res[3].axistags[2].description, "res4")

    assert_equal(res[4].shape, (20,10))
    assert_equal(res[4].axistags, a.axistags)
    assert_equal(len(res[4].axistags), 2)

    assert_equal(res[5].shape, (20,10,3))
    assert_equal(res[5].axistags, resaxistags)
    assert_equal(res[5].axistags[2].description, "res6")

    a = numpy.zeros((3,4,5))
    at = AxisTags(AxisInfo.x, AxisInfo.y, AxisInfo.z)

    r = arraytypes.taggedView(a, at)
    assert_equal(r.shape, (3,4,5))
    assert_equal(r.axistags, at)
    assert(r.axistags is not at)

    r = arraytypes.taggedView(a, 'xyz')
    assert_equal(r.shape, (3,4,5))
    assert_equal(r.axistags, at)

    r = arraytypes.taggedView(a, 'cyx')
    assert_equal(r.shape, (3,4,5))
    assert_equal(repr(r.axistags), 'c y x')

    try:
        r = arraytypes.taggedView(a, 'cxy', order='C')
        raise AssertionError, "arraytypes.taggedView() failed to throw."
    except RuntimeError:
        pass

    a = arraytypes.taggedView(a, 'zyx')
    r = arraytypes.taggedView(a, order='C')
    assert_equal(r.shape, (3,4,5))
    assert_equal(repr(r.axistags), 'z y x')
    r = arraytypes.taggedView(a, order='V')
    assert_equal(r.shape, (5,4,3))
    assert_equal(repr(r.axistags), 'x y z')
    r = arraytypes.taggedView(a, order='F')
    assert_equal(r.shape, (5,4,3))
    assert_equal(repr(r.axistags), 'x y z')

    a = a[0,...]
    r = arraytypes.taggedView(a, 'xcy')
    assert_equal(r.shape, (5,1,4))
    assert_equal(repr(r.axistags), 'x c y')
    r = arraytypes.taggedView(a, 'yxc')
    assert_equal(r.shape, (4,5,1))
    assert_equal(repr(r.axistags), 'y x c')

    try:
        r = arraytypes.taggedView(a, 'xcz')
        raise AssertionError, "arraytypes.taggedView() failed to throw."
    except RuntimeError:
        pass

    try:
        r = arraytypes.taggedView(a, 'xcz', force=True)
        raise AssertionError, "arraytypes.taggedView() failed to throw."
    except RuntimeError:
        pass

    r = arraytypes.taggedView(a, 'xz', force=True)
    assert_equal(r.shape, (4,5))
    assert_equal(repr(r.axistags), 'x z')

    a = a[..., arraytypes.newaxis('c')]
    r = arraytypes.taggedView(a, order='V')
    assert_equal(r.shape, (5, 4, 1))
    assert_equal(repr(r.axistags), 'x y c')

    r = arraytypes.taggedView(a, order='V', noChannels=True)
    assert_equal(r.shape, (5, 4))
    assert_equal(repr(r.axistags), 'x y')


def testDeepcopy():
    a = arraytypes.RGBImage(numpy.random.random((10, 4, 3)), order='C')
    b = copy.deepcopy(a)
    assert numpy.all(a == b)
    assert_equal(b.flags.c_contiguous, a.flags.c_contiguous)
    assert_equal(b.flags.f_contiguous, a.flags.f_contiguous)
    assert_equal(b.axistags, a.axistags)
    a[0,0,0] += 42
    assert b[0,0,0] != a[0,0,0]

    a = arraytypes.RGBImage(numpy.random.random((4, 10, 3)), order='V')
    b = copy.deepcopy(a)
    assert numpy.all(a == b)
    assert_equal(b.flags.c_contiguous, a.flags.c_contiguous)
    assert_equal(b.flags.f_contiguous, a.flags.f_contiguous)
    assert_equal(b.axistags, a.axistags)
    a[0,0,0] += 42
    assert b[0,0,0] != a[0,0,0]

    a = arraytypes.RGBImage(numpy.random.random((3, 4, 10)), order='F')
    b = copy.deepcopy(a)
    assert numpy.all(a == b)
    assert_equal(b.flags.c_contiguous, a.flags.c_contiguous)
    assert_equal(b.flags.f_contiguous, a.flags.f_contiguous)
    assert_equal(b.axistags, a.axistags)
    a[0,0,0] += 42
    assert b[0,0,0] != a[0,0,0]

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

def testPickle():
    import pickle
    a = arraytypes.RGBImage(numpy.random.random((10, 4, 3)), order='C')
    s = pickle.dumps(a)
    b = pickle.loads(s)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

    a = arraytypes.RGBImage(numpy.random.random((4, 10, 3)), order='V')
    s = pickle.dumps(a)
    b = pickle.loads(s)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

    a = arraytypes.RGBImage(numpy.random.random((3, 4, 10)), order='F')
    s = pickle.dumps(a)
    b = pickle.loads(s)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

def testZMQ():
    try:
        import zmq
        ctx = zmq.Context.instance()
        sender = zmq.Socket(ctx, zmq.PUSH)
        receiver = zmq.Socket(ctx, zmq.PULL)
        sender.bind('inproc://a')
        receiver.connect('inproc://a')
    except:
        return

    a = arraytypes.RGBImage(numpy.random.random((10, 4, 3)), order='C')
    a.sendSocket(sender, copy=False)
    b = arraytypes.VigraArray.receiveSocket(receiver, copy=False)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

    a = arraytypes.RGBImage(numpy.random.random((4, 10, 3)), order='V')
    a.sendSocket(sender, copy=False)
    b = arraytypes.VigraArray.receiveSocket(receiver, copy=False)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

    a = arraytypes.RGBImage(numpy.random.random((3, 4, 10)), order='F')
    a.sendSocket(sender, copy=False)
    b = arraytypes.VigraArray.receiveSocket(receiver, copy=False)
    assert_equal(b.shape, a.shape)
    assert_equal(b.strides, a.strides)
    assert_equal(b.axistags, a.axistags)
    assert numpy.all(a == b)

def testSlicing():
    a = arraytypes.Vector2Volume((5,4,3))
    a.flat[...] = xrange(a.size)

    tags = arraytypes.VigraArray.defaultAxistags('xyzc')
    assert_equal(tags, a.axistags)

    b = a[...]
    assert_true((a==b).all())
    assert_equal(tags, b.axistags)

    b = a[...,0]
    assert_equal(b.shape, a.shape[:-1])
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('xyz'))
    assert_equal(b[3,2,1], a[3,2,1,0])

    b = a[1,...]
    assert_equal(b.shape, a.shape[1:])
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('yzc'))
    assert_equal(b[3,2,1], a[1,3,2,1])

    b = a[:,2,...]
    assert_equal(b.shape, (5,3,2))
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('xzc'))
    assert_equal(b[3,2,1], a[3,2,2,1])

    b = a[:,1,2,...]
    assert_equal(b.shape, (5,2))
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('xc'))
    assert_equal(b[2,1], a[2,1,2,1])

    b = a[2:4, :, 2, ...]
    assert_equal(b.shape, (2, 4, 2))
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('xyc'))
    assert_equal(b[0,2,1], a[2,2,2,1])

    b = a[1:4, :, arraytypes.newaxis(arraytypes.AxisInfo.t), 1, ...]
    assert_equal(b.shape, (3, 4, 1, 2))
    assert_equal(b.axistags, arraytypes.VigraArray.defaultAxistags('xytc'))
    assert_equal(b[0,2,0,0], a[1,2,1,0])

    b = a[..., None, :,1]
    assert_equal(b.shape, (5, 4, 1, 3))
    rtags = arraytypes.AxisTags(arraytypes.AxisInfo.x, arraytypes.AxisInfo.y, arraytypes.AxisInfo(), arraytypes.AxisInfo.z)
    assert_equal(b.axistags, rtags)
    assert_equal(b[0,3,0,1], a[0,3,1,1])

    b = a.subarray((4,3,2))
    assert_equal(b.shape, (4,3,2,2))
    assert_true((a[:4,:3,:2,:]==b).all())
    assert_equal(tags, b.axistags)

    b = a.subarray((1,1,1),(4,3,2))
    assert_equal(b.shape, (3,2,1,2))
    assert_true((a[1:4,1:3,1:2]==b).all())
    assert_equal(tags, b.axistags)

    b = a.subarray((1,1,1,1),(4,3,2,2))
    assert_equal(b.shape, (3,2,1,1))
    assert_true((a[1:4,1:3,1:2,1:]==b).all())
    assert_equal(tags, b.axistags)

def testMethods():
    a = arraytypes.ScalarImage((20, 30))
    ones = arraytypes.ScalarImage((20, 30), value=1)

    a.ravel()[...] = range(a.size)

    for k, i in zip(a.flat, xrange(a.size)):
        assert_equal(k, i)

    assert (a.flatten() == range(a.size)).all()

    assert (a >= 0).all()
    assert not (a == 0).all()

    assert (a == 0).any()
    assert not (a == -1).any()

    assert_equal(a.argmax(), a.size-1)
    assert (a.argmax(axis='y') == a.shape[1]-1).all()

    assert_equal(a.argmin(), 0)
    assert (a.argmin(axis='y') == 0).all()

    assert (ones.cumsum()-1 == a.ravel()).all()
    oc = ones.cumsum(axis='x')-1
    for s in oc.sliceIter('y'):
        assert (s == range(a.shape[0])).all()

    assert (ones.cumprod() == 1).all()
    assert (ones.cumprod(axis='x') == 1).all()

    assert_equal(a.max(), a.size-1)
    assert (a.max(axis='y') == range(a.size-a.shape[0], a.size)).all()

    assert_equal(a.mean(dtype=numpy.longdouble), (a.size - 1.0) / 2.0)
    assert (a.mean(axis='y', dtype=numpy.longdouble) ==
            range((a.size-a.shape[0])/2, (a.size+a.shape[0])/2)).all()

    assert_equal(a.min(), 0)
    assert (a.min(axis='y') == range(a.shape[0])).all()

    n = arraytypes.ScalarImage(numpy.array([[1, 0, 0],[0, 1, 1],[1, 0, 1]]))
    nz = n.nonzero()
    assert (nz[0] == [0, 1, 1, 2, 2]).all()
    assert_equal(nz[0].axistags, n.defaultAxistags('x'))
    assert (nz[1] == [0, 1, 2, 0, 2]).all()
    assert_equal(nz[1].axistags, n.defaultAxistags('y'))

    assert_equal(ones.prod(), 1.0)
    assert (ones.prod(axis='y') == [1]*ones.shape[0]).all()

    assert_equal(a.ptp(), a.size-1)
    assert (a.ptp(axis='x') == [a.shape[0]-1]*a.shape[1]).all()

    r = arraytypes.ScalarImage((2,2))
    r.ravel()[...] = range(4)

    assert (r.repeat(1) == r.ravel()).all()
    assert (r.repeat(2) == reduce(lambda x,y: x+[y,y], range(4), [])).all()
    assert (r.repeat([0,1,2,3]) == [1,2,2,3,3,3]).all()
    assert (r.repeat(2, axis='y').ravel() == [0,1,0,1,2,3,2,3]).all()
    assert (r.repeat([1,2], axis='y').ravel() == [0,1,2,3,2,3]).all()

    s = a[arraytypes.AxisInfo.c,:,arraytypes.AxisInfo.z,:,arraytypes.AxisInfo.t]
    assert_equal(s.shape, (1, a.shape[0], 1, a.shape[1], 1))
    assert_equal(s.axistags,a.defaultAxistags('cxzyt'))
    ss = s.squeeze()
    assert_equal(ss.shape, a.shape)
    assert_equal(ss.axistags,a.axistags)

    assert_equal(ones.std(dtype=numpy.longdouble), 0.0)
    assert (ones.std(axis='x', dtype=numpy.longdouble) == [0.0]*a.shape[1]).all()

    assert_equal(ones.sum(dtype=numpy.longdouble), ones.size)
    assert (ones.sum(axis='x', dtype=numpy.longdouble) == [a.shape[0]]*a.shape[1]).all()

    b = a.swapaxes(0, 1)
    assert_equal(b.shape, (a.shape[1], a.shape[0]))
    assert_equal(len(b.axistags), 2)
    assert_equal(b.axistags[0], a.axistags[1])
    assert_equal(b.axistags[1], a.axistags[0])

    b = a.swapaxes(0, 1, keepTags=True)
    assert_equal(b.shape, (a.shape[1], a.shape[0]))
    assert_equal(len(b.axistags), 2)
    assert_equal(b.axistags, a.axistags)

    rt = r.take([1,2])
    assert (rt == [1,2]).all()
    assert_equal(rt.axistags, arraytypes.AxisTags(1))
    rt = r.take([0,1], axis='y')
    assert (rt == r).all()
    assert_equal(rt.axistags, rt.axistags)

    assert_equal(ones.var(dtype=numpy.longdouble), 0.0)
    assert (ones.var(axis='x', dtype=numpy.longdouble) == [0.0]*a.shape[1]).all()

    a = arraytypes.Image((5,4,3))
    b = a.transpose()
    assert_equal(b.shape, (3,4,5))
    assert_equal(len(b.axistags), len(a.axistags))
    assert_equal(b.axistags[0], a.axistags[2])
    assert_equal(b.axistags[1], a.axistags[1])
    assert_equal(b.axistags[2], a.axistags[0])
    b = a.transpose((1,2,0))
    assert_equal(b.shape, (4,3,5))
    assert_equal(len(b.axistags), len(a.axistags))
    assert_equal(b.axistags[0], a.axistags[1])
    assert_equal(b.axistags[1], a.axistags[2])
    assert_equal(b.axistags[2], a.axistags[0])
    b = a.transpose(keepTags=True)
    assert_equal(b.shape, (3,4,5))
    assert_equal(len(b.axistags), len(a.axistags))
    assert_equal(b.axistags, a.axistags)
    b = a.transpose((1,2,0), keepTags=True)
    assert_equal(b.shape, (4,3,5))
    assert_equal(len(b.axistags), len(a.axistags))
    assert_equal(b.axistags, a.axistags)

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
        assert_equal(a.axistags, b.axistags)
        if not numpyHasComplexNegateBug or t is not clongdouble:
            assert (b == -t(2)).all()
        b = a + a
        assert_equal(t, b.dtype)
        assert_equal(a.axistags, b.axistags)
        assert (b == t(4)).all()
        b = a == a
        assert_equal(bool, b.dtype)
        assert_equal(a.axistags, b.axistags)
        assert (b == True).all()
        b = ones[t] + a
        assert_equal(a.shape, b.shape)
        assert_equal(a.axistags, b.axistags)
        assert_equal(t, b.dtype)
        assert (b == t(3)).all()
        b = a + ones[t]
        assert_equal(a.shape, b.shape)
        assert_equal(a.axistags, b.axistags)
        assert_equal(t, b.dtype)
        assert (b == t(3)).all()
        b = 3 + a
        assert_equal(t, b.dtype)
        assert_equal(a.axistags, b.axistags)
        assert (b == t(5)).all()
        b = a + 3
        assert_equal(t, b.dtype)
        assert_equal(a.axistags, b.axistags)
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
            assert_equal(a1.axistags, b.axistags)
            assert (b == 4).all()
            b = a1 <= a2
            assert_equal(bool, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == True).all()
        for j in floats + compl:
            a2 = arrays[j]
            b = a1 * a2
            assert_equal(j, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == 4).all()
            b = a2 * a1
            assert_equal(j, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == 4).all()
            b = a1 >= a2
            assert_equal(bool, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == True).all()
            b = a2 > a1
            assert_equal(bool, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == False).all()
        b = a1 + 1
        assert (b == 3).all()
        assert_equal(a1.dtype, b.dtype)
        assert_equal(a1.axistags, b.axistags)
        b = a1 + 1.0
        assert_equal(a1.axistags, b.axistags)
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
            assert_equal(a1.axistags, b.axistags)
            assert (b == 4).all()
            b = a2 * a1
            assert_equal(j, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == 4).all()
            b = a1 >= a2
            assert_equal(bool, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == True).all()
            b = a2 > a1
            assert_equal(bool, b.dtype)
            assert_equal(a1.axistags, b.axistags)
            assert (b == False).all()
        b = a1 + 0.5
        assert (b == 2.5).all()
        assert_equal(a1.dtype, b.dtype)
        assert_equal(a1.axistags, b.axistags)

        fractional, integral = ufunc.modf(b)
        assert (fractional == 0.5).all()
        assert (integral == 2.0).all()
        assert_equal(a1.dtype, fractional.dtype)
        assert_equal(a1.axistags, fractional.axistags)
        assert_equal(a1.dtype, integral.dtype)
        assert_equal(a1.axistags, integral.axistags)

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
    assert_equal(arrays[complex64].axistags, b.axistags)
    b = abs(arrays[complex128])
    assert (b == 2.0).all()
    assert_equal(float64, b.dtype)
    assert_equal(arrays[complex128].axistags, b.axistags)
    b = abs(arrays[clongdouble])
    assert (b == 2.0).all()
    assert_equal(longdouble, b.dtype)
    assert_equal(arrays[clongdouble].axistags, b.axistags)

    a = arraytypes.ScalarImage((2,2), uint8, value=255)
    b = a + a
    assert (b == 254).all()
    b = a + 255
    assert (b == 254).all()
    b = 255 + a
    assert (b == 254).all()

    b = arraytypes.ScalarImage((2,2), int32, value=0)
    bb = ufunc.add(a, a, b)
    assert bb is b
    assert (b == 510).all()
    b = arraytypes.ScalarImage((2,2), int32, value=0)
    bb = ufunc.add(a, 255, b)
    assert bb is b
    assert (b == 510).all()
    b = arraytypes.ScalarImage((2,2), int32, value=0)
    bb = ufunc.add(255, a, b)
    assert bb is b
    assert (b == 510).all()
