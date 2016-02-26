from __future__ import division, print_function
import sys
print("\nexecuting test file", __file__, file=sys.stderr)
exec(compile(open('set_paths.py', "rb").read(), 'set_paths.py', 'exec'))

import numpy
import vigra

def _impl_test_labelMultiArray(dtype):
    a = numpy.zeros( (10,10,10,10,2), dtype=dtype )
    a[3:4, 5:8, 5:8, 5:8] = 100
    a[4:5, 8:10, 5:8, 5:8] = 100 # touches above object on corners only.
    
    a[3:4, 1:4, 5:8, 5:8] = 200
    
    labeled = vigra.analysis.labelMultiArray(a)
    assert labeled.dtype == numpy.uint32
    assert labeled.max() == 4

    labeled = vigra.analysis.labelMultiArrayWithBackground(a)
    assert labeled.dtype == numpy.uint32
    assert labeled.max() == 3

    labeled = vigra.analysis.labelMultiArrayWithBackground(a, neighborhood='indirect')
    assert labeled.dtype == numpy.uint32
    assert labeled.max() == 2
    
    labeled = vigra.analysis.labelMultiArrayWithBackground(a, background_value=100)
    assert labeled.dtype == numpy.uint32
    assert labeled.max() == 2

def test_labelMultiArray():
    _impl_test_labelMultiArray(numpy.uint8)
    _impl_test_labelMultiArray(numpy.uint32)
    _impl_test_labelMultiArray(numpy.float32)


def _impl_test_applyMapping(dtype):
    original = numpy.arange(100, dtype=dtype ).reshape(10,10)
    mapping = dict( zip( original.flat[:], original.flat[:] + 100 ) )

    # Not in-place
    remapped = vigra.analysis.applyMapping(original, mapping)
    assert (remapped == original+100).all()

    # in-place
    original_copy = original.copy()
    vigra.analysis.applyMapping(original_copy, mapping, out=original_copy)
    assert (original_copy == original+100).all()

    # Different dtypes
    mapping = dict( zip( original.flat[:], (original.flat[:] + 100).astype(numpy.uint64) ) )

    result = numpy.zeros_like( original, dtype=numpy.uint64 )
    vigra.analysis.applyMapping(original, mapping, out=result)
    assert (result == original+100).all()

    mapping = dict( zip( original.flat[:], (original.flat[:] + 100).astype(numpy.uint8) ) )

    result = numpy.zeros_like( original, dtype=numpy.uint8 )
    vigra.analysis.applyMapping(original, mapping, out=result)
    assert (result == original+100).all()

    # Incomplete mapping
    for i in range(10):
        del mapping[i]

    remapped = vigra.analysis.applyMapping(original, mapping, allow_incomplete_mapping=True)
    assert (remapped[0] == original[0]).all()
    assert (remapped[1:] == original[1:]+100).all()
    
    try:
        remapped = vigra.analysis.applyMapping(original, mapping, allow_incomplete_mapping=False)
    except IndexError:
        pass
    else:
        assert False, "Expected to get an exception due to the incomplete mapping!"

def test_applyMapping():
    _impl_test_applyMapping(numpy.uint8)
    _impl_test_applyMapping(numpy.uint32)
    _impl_test_applyMapping(numpy.uint64)


def _impl_test_unique(dtype):
    a = numpy.array([2,3,5,7,11,13,17,19,23,29] + [2,3,5,7,11,13,17,19,23,29], dtype=dtype)
    u = vigra.analysis.unique(a)
    assert set(u) == set([2,3,5,7,11,13,17,19,23,29])

def test_unique():
    _impl_test_unique(numpy.uint8)
    _impl_test_unique(numpy.uint32)
    _impl_test_unique(numpy.uint64)

def ok_():
    print(".", file=sys.stderr)
