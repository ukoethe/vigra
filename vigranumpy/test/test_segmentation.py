import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

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
    
def ok_():
    print >> sys.stderr, ".",
