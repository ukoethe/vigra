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

import copy
import numpy
import ufunc
import sys

from numpy import int8, uint8, int16, uint16, int32, uint32, int64, uint64
from numpy import float32, float64, longdouble, complex64, complex128, clongdouble

try:
    import qimage2ndarray
except:
    import vigra
    vigra._fallbackModule('qimage2ndarray',
    '''    It can be obtained at
    http://pypi.python.org/pypi/qimage2ndarray/.''')
    import qimage2ndarray

def qimage2array(q):
    '''Create a view to the given array with the appropriate type.

       q.format() must be QImage.Format_Indexed8, QImage.Format_RGB32, or
       QImage.Format_ARGB32, and you will get ScalarImage, RGBImage, or
       Vector4Image respectively, all with dtype=uint8. The channels in
       a Vector4Image will be ordered as [alpha, red, green, blue].
    '''
    import PyQt4.QtGui as qt
    import qimage2ndarray
    if q.format() == qt.QImage.Format_Indexed8:
        width, height = q.width(), q.height()
        return qimage2ndarray.byte_view(q).swapaxes(0,1).reshape((width, height)).view(ScalarImage)
    if q.format() == qt.QImage.Format_RGB32:
        return qimage2ndarray.rgb_view(q).swapaxes(0,1).view(RGBImage)
    if q.format() == qt.QImage.Format_ARGB32:
        return qimage2ndarray.byte_view(q, 'big').swapaxes(0,1).view(Vector4Image)
    raise RuntimeError("qimage2array(): q.format() must be Format_Indexed8, Format_RGB32, or Format_ARGB32")

class classproperty(object):
    def __get__(self, instance, cls):
            if self.__instance_method is not None and instance is not None:
                return self.__instance_method(instance)
            else:
                return self.__class_method(cls)
    def __init__(self, class_method, instance_method = None):
            self.__class_method = class_method
            self.__instance_method = instance_method

def _array_docstring_(name, shape, compat):
    return '''
    Constructor:
    
    .. method:: %(name)s(obj, dtype=numpy.float32, order='V', init = True, value = None)

        :param obj: a data or shape object (see below)
        :param dtype: desired element type
        :param order: desired memory layout (see below)
        :param init: True: initialize the image with zeros; False: do not initialize the image
        :type init: boolean
        :param value: initialize the image with this value (overrides init)
        :type value: convertible to dtype
 
        **obj** may be one of the following

        * If obj is a vigra.%(name)s or a subclass, a copy of obj with the
          given dtype and order is created, and obj's class is transferred.
        * If obj is another subtype of numpy.ndarray with compatible shape,
          a *transposed* copy of obj with the given dtype, order and class
          vigra.%(name)s is created. Transposition changes the order of the
          spatial dimensions (and therefore the index order for element access)
          from [x,y] or [z,y,x] to the VIGRA convention [x,y] and [x,y,z]
          respectively. The channel dimension is assumed to be the last
          dimension and remains in that position. (Use numpy.rollaxis()
          to adjust the input if necessary.)
        * If obj is a sequence, it is interpreted as a shape. When
          the shape is compatible, a new vigra.%(name)s with the given
          dtype and order is created.
        * Otherwise, or if the shape is not compatible, an exception
          is raised.
      
        %(shape)s
      
        **order** can be 'C' (C order), 'F' (Fortran order), 'V' (vector-valued
        order), and 'A'.

          'C' and 'F' order:
            have the usual numpy meaning

          'V' order:
            is an interleaved memory layout that simulates vector-
            valued pixels or voxels: while the spatial dimensions are arranged
            as in Fortran order, the major memory-aligned dimension is the
            channel (i.e. last) dimension. Arrays in 'V'-order are compatible
            with vector-valued NumpyArrays. For example, an RGBImage((4,3), uint8)
            has strides (3, 12, 1) and is compatible with
            NumpyArray<2, RGBValue<UInt8>, UnstridedArrayTag>.

          'A' order:
            defaults to 'V' when a new array is created, and means
            'preserve order' when an existing array is copied.
    
        In particular, the following compatibility rules apply (Note that
        compatibility with 'UnstridedArrayTag' implies compatibility with
        'StridedArrayTag'. Due to their loop order, VIGRA algorithms are
        generally more efficient when the memory layout is compatible with
        'UnstridedArrayTag'. T is the array's dtype.):
        
%(compat)s
    ''' % {'name': name, 'shape': shape, 'compat': compat}

##################################################################

class _VigraArray(numpy.ndarray):
    """
This base class ensures that arrays created in Python are
compatible with the memory layout requirements of
VIGRA's NumpyArray family of C++ views. Do always use
this class via its subclasses!
    """
    def __new__(cls, obj, dtype=numpy.float32, order='V', init = True, value = None):
        from vigranumpycore import constructNumpyArray
        if isinstance(obj, numpy.ndarray) and not isinstance(obj, _VigraArray):
            obj = obj.swapaxes(0, cls.spatialDimensions-1)
        if value is not None:
            res = constructNumpyArray(cls, obj, cls.spatialDimensions, cls.channels, dtype, order, False)
            res[...] = value
            return res
        else:
            return constructNumpyArray(cls, obj, cls.spatialDimensions, cls.channels, dtype, order, init)
    
    @property
    def order(self):
        if self.flags.c_contiguous:
            return 'C'
        elif self.flags.f_contiguous:
            return 'F'
        elif self.channels > 1 and self.itemsize == self.strides[-1] and \
             reduce(lambda x, y: y if y >= x and x >= 0 else -1, self.strides[:-1], 0) >= 0:
            return 'V'
        return 'A'
    
    def astype(self, dtype):
        return self.__class__(self, dtype=dtype)
    
    def copy(self, order = 'A'):
        return self.__class__(self, dtype=self.dtype, order=order)
    
    def __copy__(self, order = 'A'):
        return self.copy(order)
    
    def __deepcopy__(self, memo):
        result = self.__class__(self, dtype=self.dtype, order="A")
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        return result
    
    def flatten(self, order = 'C'):
        return self.view(numpy.ndarray).swapaxes(0, self.spatialDimensions-1).flatten(order)
    
    def __str__(self, separator = ' ', precision=2, suppress_small=True):
        return numpy.array2string(self.T, separator = separator, precision=precision, suppress_small=suppress_small)
    
    def __repr__(self):
        return "%s(dtype=%s, shape=%s, data=\n%s)" % \
          (self.__class__.__name__, str(self.dtype), str(self.shape), self.__str__(', '))
          
    def bands(self):
        if len(self.shape) == self.spatialDimensions:
            return 1
        else:
            return self.shape[-1]

    channels = classproperty(lambda cls: 0, bands)
    
    @property
    def flat(self):
        return self.view(numpy.ndarray).swapaxes(0, self.spatialDimensions-1).flat
    
    __array_priority__ = 10.0
    
    # we reimplement the numerical operators in order to make sure that array order is preserved
    def __abs__(self):
        return ufunc.absolute(self)
    
    def __add__(self, other):
        return ufunc.add(self, other)
        
    def __and__(self, other):
        return ufunc.bitwise_and(self, other)
        
    def __div__(self, other):
        return ufunc.divide(self, other)
    
    def __divmod__(self, other):
        return ufunc.floor_divide(self, other), ufunc.remainder(self, other)
    
    def __eq__(self, other):
        return ufunc.equal(self, other)
    
    def __floordiv__(self, other):
        return ufunc.floor_divide(self, other)
    
    def __ge__(self, other):
        return ufunc.greater_equal(self, other)
    
    def __gt__(self, other):
        return ufunc.greater(self, other)
        
    def __invert__(self):
        return ufunc.invert(self)
    
    def __le__(self, other):
        return ufunc.less_equal(self, other)
    
    def __lshift__(self, other):
        return ufunc.left_shift(self, other)
    
    def __lt__(self, other):
        return ufunc.less(self, other)
    
    def __mod__(self, other):
        return ufunc.remainder(self, other)
    
    def __mul__(self, other):
        return ufunc.multiply(self, other)
    
    def __ne__(self, other):
        return ufunc.not_equal(self, other)
    
    def __neg__(self):
        return ufunc.negative(self)
    
    def __or__(self, other):
        return ufunc.bitwise_or(self, other)
    
    def __pos__(self):
        return self
    
    def __pow__(self, other):
        return ufunc.power(self, other)
    
    def __radd__(self, other):
        return ufunc.add(other, self)
    
    def __radd__(self, other):
        return ufunc.add(other, self)
    
    def __rand__(self, other):
        return ufunc.bitwise_and(other, self)
    
    def __rdiv__(self, other):
        return ufunc.divide(other, self)
    
    def __rdivmod__(self, other):
        return ufunc.floor_divide(other, self), ufunc.remainder(other, self)
    
    def __rfloordiv__(self, other):
        return ufunc.floor_divide(other, self)
    
    def __rlshift__(self, other):
        return ufunc.left_shoft(other, self)
    
    def __rmod__(self, other):
        return ufunc.remainder(other, self)
    
    def __rmul__(self, other):
        return ufunc.multiply(other, self)
    
    def __ror__(self, other):
        return ufunc.bitwise_or(other, self)
    
    def __rpow__(self, other):
        return ufunc.power(other, self)
    
    def __rrshift__(self, other):
        return ufunc.right_shift(other, self)
    
    def __rshift__(self, other):
        return ufunc.right_shift(self, other)

    def __rsub__(self, other):
        return ufunc.subtract(other, self)
    
    def __rtruediv__(self, other):
        return ufunc.true_divide(other, self)
    
    def __rxor__(self, other):
        return ufunc.bitwise_xor(other, self)
    
    def __sub__(self, other):
        return ufunc.subtract(self, other)
    
    def __truediv__(self, other):
        return ufunc.true_divide(self, other)
    
    def __xor__(self, other):
        return ufunc.bitwise_xor(self, other)

##################################################################

class Image(_VigraArray):
    __doc__ = _array_docstring_('Image', '''A shape is compatible when it has two dimensions (width, height) or three dimensions (width, height, channels).''', """
          'C':
             | NumpyArray<2, T, StridedArrayTag> (if channels=1),
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag> (if channels>1),
             | NumpyArray<2, TinyVector<T, M>, StridedArrayTag> (if channels=M),
             | NumpyArray<2, RGBValue<T>, StridedArrayTag> (if channels=3),
             | NumpyArray<2, Singleband<T>, StridedArrayTag> (if channels=1),
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<2, T, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag> (if channels>1),
             | NumpyArray<2, Singleband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<2, T, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, T, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, T, StridedArrayTag> (if channels>1),
             | NumpyArray<4, T, StridedArrayTag> (if channels>1),
             | NumpyArray<2, Singleband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<2, TinyVector<T, M>, UnstridedArrayTag> (if channels=M),
             | NumpyArray<2, RGBValue<T>, UnstridedArrayTag> (if channels=3),
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, Multiband<T>, StridedArrayTag> (if channels>1)
           
""")
            
    def write(self, filename, dtype = '', compression = ''):
        "Write an image to a file. Consult :func:`vigra.impex.writeImage` for detailed documentation"
        import vigra.impex
        vigra.impex.writeImage(self, filename, dtype, compression)
            
    def writeHDF5(self, filename, pathInFile, dtype = ''):
        "Write an image to a HDF5 file. Consult :func:`vigra.impex.writeImageToHDF5` for detailed documentation"
        import vigra.impex
        vigra.impex.writeImageToHDF5(self, filename, pathInFile, dtype)

    def show(self, normalize = True):
        '''
        Display this image in a vigra.pyqt.ImageWindow.
        
        The parameter `normalize` can be used to normalize an image's
        value range to 0..255:

        `normalize` = (nmin, nmax):
          scale & clip image values from nmin..nmax to 0..255

        `normalize` = nmax:
          lets nmin default to zero, i.e. scale & clip the range 0..nmax
          to 0..255

        `normalize` = True: (default)
          scale the image's actual range min()..max() to 0..255

        `normalize` = False:
          don't scale the image's values
           
        '''
        from pyqt.imagewindow import showImage
        return showImage(self, normalize)

    def qimage(self, normalize = True):
        """
        Convert this image to a Qt QImage (mainly for display purposes).
        The present image must have 1, 2, 3, or 4 channels, and the resulting
        QImage will have QImage.Format_Indexed8 iff there was only one
        channel and QImage.Format_[A]RGB32 otherwise (with the last of
        2/4 channels being used as alpha channel).
        
        The parameter `normalize` can be used to normalize an image's
        value range to 0..255:

        `normalize` = (nmin, nmax):
          scale & clip image values from nmin..nmax to 0..255

        `normalize` = nmax:
          lets nmin default to zero, i.e. scale & clip the range 0..nmax
          to 0..255

        `normalize` = True: (default)
          scale the image's actual range min()..max() to 0..255

        `normalize` = False:
          don't scale the image's values
           
        """
        import qimage2ndarray

        yxImage = self.swapaxes(0, 1)

        if self.channels == 1:
            q = qimage2ndarray.gray2qimage(yxImage, normalize)
        else:
            q = qimage2ndarray.array2qimage(yxImage, normalize)

        return q
        
    @property
    def width(self):
        """the image's width"""
        return self.shape[0]
    
    @property
    def height(self):
        "the image's height"
        return self.shape[1]

    @classproperty
    def spatialDimensions(cls):
        "number of spatial dimensions (useful for distinguishing RGBImage and ScalarVolume)"
        return 2

class ScalarImage(Image):
    __doc__ = _array_docstring_('ScalarImage', '''A shape is compatible when it has two dimensions (width, height) or three dimensions (width, height, 1).''', """
          'C':
             | NumpyArray<2, T, StridedArrayTag>,
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<2, Singleband<T>, StridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<2, T, UnstridedArrayTag>,
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<2, Singleband<T>, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | like 'F'""")

    channels = classproperty(lambda cls: 1, Image.bands)
        
class Vector2Image(Image):
    __doc__ = _array_docstring_('Vector2Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 2).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 2>, StridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 2>, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 2, Image.bands)

class Vector3Image(Image):
    __doc__ = _array_docstring_('Vector3Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 3).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 3, Image.bands)

class Vector4Image(Image):
    __doc__ = _array_docstring_('Vector4Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 4).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 4>, StridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 4>, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 4, Image.bands)

class RGBImage(Vector3Image):
    __doc__ = _array_docstring_('RGBImage', '''A shape is compatible when it has two dimensions (width, height) or three dimensions (width, height, 3).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, RGBValue<T>, StridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<2, RGBValue<T>, UnstridedArrayTag>,
             | NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
             | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

#################################################################

class Volume(_VigraArray):
    __doc__ = _array_docstring_('Volume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, channels).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag> (if channels=1),
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, M>, StridedArrayTag> (if channels=M),
             | NumpyArray<3, RGBValue<T>, StridedArrayTag> (if channels=3),
             | NumpyArray<3, Singleband<T>, StridedArrayTag> (if channels=1),
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag> (if channels=1),
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Singleband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<3, T, UnstridedArrayTag> (if channels=1),
             | NumpyArray<4, T, UnstridedArrayTag> (if channels=1),
             || NumpyArray<4, T, StridedArrayTag> (if channels>1),
             | NumpyArray<3, Singleband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<3, TinyVector<T, M>, UnstridedArrayTag> (if channels=M),
             | NumpyArray<3, RGBValue<T>, UnstridedArrayTag> (if channels=3),
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag> (if channels=1),
             | NumpyArray<4, Multiband<T>, StridedArrayTag> (if channels>1)""")
            
    def write(self, filename_base, filename_ext, dtype = '', compression = ''):
        "Write a volume to a sequence of files. Consult :func:`vigra.impex.writeVolume` for detailed documentation.\n"
        import vigra.impex
        vigra.impex.writeVolume(self, filename_base, filename_ext, dtype, compression)
            
    def writeHDF5(self, filename, pathInFile, dtype = ''):
        "Write a volume to a HDF5 file. Consult :func:`vigra.impex.writeVolumeToHDF5` for detailed documentation.\n"
        import vigra.impex
        vigra.impex.writeVolumeToHDF5(self, filename, pathInFile, dtype)
    
    @classproperty
    def spatialDimensions(cls): return 3
        
    @property
    def width(self):
        return self.shape[0]
    
    @property
    def height(self):
        return self.shape[1]
    
    @property
    def depth(self):
        return self.shape[2]

class ScalarVolume(Volume):
    __doc__ = _array_docstring_('ScalarVolume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 1).''', """
          'C':
             | NumpyArray<3, T, StridedArrayTag>,
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, Singleband<T>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<3, T, UnstridedArrayTag>,
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<3, Singleband<T>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | like 'F'""")

    channels = classproperty(lambda cls: 1, Volume.bands)

class Vector2Volume(Volume):
    __doc__ = _array_docstring_('Vector2Volume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 2).''', """
          'C':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 2>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 2>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 2, Volume.bands)

class Vector3Volume(Volume):
    __doc__ = _array_docstring_('Vector3Volume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 3).''', """
          'C':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 3, Volume.bands)

class Vector4Volume(Volume):
    __doc__ = _array_docstring_('Vector4Volume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 4).''', """
          'C':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 4>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 4>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 4, Volume.bands)
    
class Vector6Volume(Volume):
    __doc__ = _array_docstring_('Vector4Volume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 6).''', """
          'C':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 6>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 6>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 6, Volume.bands)
    
class RGBVolume(Vector3Volume):
    __doc__ = _array_docstring_('RGBVolume', '''
    A shape is compatible when it has three dimensions (width, height,
    depth) or four dimensions (width, height, depth, 3).''', """
          'C':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, RGBValue<T>, StridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>
          'F':
             | NumpyArray<4, T, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          'V':
             | NumpyArray<4, T, StridedArrayTag>,
             | NumpyArray<3, RGBValue<T>, UnstridedArrayTag>,
             | NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
             | NumpyArray<4, Multiband<T>, StridedArrayTag>""")


#################################################################

class ImagePyramid(list):
    def __init__(self, image, copyImagedestLevel = 0, lowestLevel = 0, highestLevel = 0):
        ''' Create a new pyramid.
            The new pyramid levels range from 'lowestLevel' to 'highestLevel' (inclusive),
            and the given 'image' is copied to 'copyImagedestLevel'. The images at other
            levels are filled with zeros and sized so that the shape is reduced by half
            when going up (to higher levels), and doubled when going down.
        '''
        if lowestLevel > copyImagedestLevel or highestLevel < copyImagedestLevel:
            raise ValueError('ImagePyramid(): copyImagedestLevel must be between lowestLevel and highestLevel (inclusive)')
        
        list.__init__(self, [image.__class__(image, dtype=image.dtype)])
        self._lowestLevel = copyImagedestLevel
        self._highestLevel = copyImagedestLevel
        self.createLevel(lowestLevel)
        self.createLevel(highestLevel)

    @property
    def lowestLevel(self):
        '''The pyramids lowest level.
        '''
        return self._lowestLevel
    
    @property
    def highestLevel(self):
        '''The pyramids highest level (inclusive).
        '''
        return self._highestLevel
    
    def __getitem__(self, level):
        '''Get the image at 'level'.
           Raises IndexError when the level does not exist.
        '''
        if level < self.lowestLevel or level > self.highestLevel:
            raise IndexError("ImagePyramid[level]: level out of range.")
        return list.__getitem__(self, level - self.lowestLevel)
    
    def __setitem__(self, level, image):
        '''Copy the data of the given 'image' to the image at 'level'.
           Raises IndexError when the level does not exist.
        '''
        self[level][...] = image[...]
        
    def expandImpl(self, src, dest, centerValue):
        import filters
        
        ss, ds = src.shape, dest.shape
        s = [ss[k] if 2*ss[k] == ds[k] else -1 for k in range(len(ss))]
    
        smooth1 = filters.explicitlyKernel(-1, 1, numpy.array([0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue]))
        smooth2 = filters.explicitlyKernel(-1, 0, numpy.array([0.5, 0.5]));

        filters.convolve(src, (smooth1, smooth1), out=dest[::2,::2])
        filters.convolve(src[:,:s[1]], (smooth1, smooth2), out=dest[::2,1::2])
        filters.convolve(src[:s[0],:], (smooth2, smooth1), out=dest[1::2,::2])
        filters.convolve(src[:s[0],:s[1]], (smooth2, smooth2), out=dest[1::2,1::2])
    
    def reduce(self, srcLevel, destLevel, centerValue = 0.42):
        '''Reduce the image at 'srcLevel' to 'destLevel', using the Burt smoothing filter
           with the given 'centerValue'. srcLevel must be smaller than destLevel.
           
           For more details, see pyramidReduceBurtFilter_ in the C++ documentation.
        '''
        # FIXME: This should be implemented in C++
        # FIXME: This should be implemented for arbitrary dimensions
        import filters
        
        if srcLevel > destLevel:
            raise RuntimeError("ImagePyramid::reduce(): srcLevel <= destLevel required.")
        if srcLevel < self.lowestLevel or srcLevel > self.highestLevel:
            raise RuntimeError("ImagePyramid::reduce(): srcLevel does not exist.")
        self.createLevel(destLevel)
        
        smooth = filters.burtFilterKernel(0.25 - 0.5*centerValue)
        for k in range(srcLevel, destLevel):
            i = filters.convolve(self[k], smooth)
            self[k+1] = i[::2,::2]

    def expand(self, srcLevel, destLevel, centerValue = 0.42):
        '''Expand the image at 'srcLevel' to 'destLevel', using the Burt smoothing filter
           with the given 'centerValue'. srcLevel must be larger than destLevel.
           
           For more details, see pyramidExpandBurtFilter_ in the C++ documentation.
        '''
        # FIXME: This should be implemented in C++
        # FIXME: This should be implemented for arbitrary dimensions
        if srcLevel < destLevel:
            raise RuntimeError("ImagePyramid::expand(): srcLevel >= destLevel required.")
        if srcLevel < self.lowestLevel or srcLevel > self.highestLevel:
            raise RuntimeError("ImagePyramid::expand(): srcLevel does not exist.")
        self.createLevel(destLevel)

        for k in range(srcLevel, destLevel, -1):
            self.expandImpl(self[k], self[k-1], centerValue)

    def reduceLaplacian(self, srcLevel, destLevel, centerValue = 0.42):
        '''Reduce the image at 'srcLevel' to 'destLevel', using the Burt smoothing filter
           with the given 'centerValue', and compute Laplacian images for the levels
           srcLevel ... destLevel-1. srcLevel must be smaller than destLevel.
           
           For more details, see pyramidReduceBurtLaplacian_ in the C++ documentation.
        '''
        # FIXME: This should be implemented in C++
        # FIXME: This should be implemented for arbitrary dimensions
        import filters
        
        if srcLevel > destLevel:
            raise RuntimeError("ImagePyramid::reduceLaplacian(): srcLevel <= destLevel required.")
        if srcLevel < self.lowestLevel or srcLevel > self.highestLevel:
            raise RuntimeError("ImagePyramid::reduceLaplacian(): srcLevel does not exist.")
        self.createLevel(destLevel)

        smooth = filters.burtFilterKernel(0.25 - 0.5*centerValue)
        for k in range(srcLevel, destLevel):
            i = filters.convolve(self[k], smooth)
            self[k+1] = i[::2,::2]
            self.expandImpl(self[k+1], i, centerValue)
            self[k] = i - self[k]

    def expandLaplacian(self, srcLevel, destLevel, centerValue = 0.42):
        '''Expand the image at 'srcLevel' to 'destLevel', using the Burt smoothing filter
           with the given 'centerValue', and reconstruct the images for the levels
           srcLevel-1 ... destLevel from their Laplacian images. srcLevel must be larger than destLevel.
           
           For more details, see pyramidExpandBurtLaplacian_ in the C++ documentation.
        '''
        # FIXME: This should be implemented in C++
        # FIXME: This should be implemented for arbitrary dimensions
        import filters
        
        if srcLevel < destLevel:
            raise RuntimeError("ImagePyramid::expandLaplacian(): srcLevel >= destLevel required.")
        if srcLevel < self.lowestLevel or srcLevel > self.highestLevel:
            raise RuntimeError("ImagePyramid::expandLaplacian(): srcLevel does not exist.")
        self.createLevel(destLevel)

        smooth1 = filters.explicitlyKernel(-1, 1, numpy.array([0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue]))
        smooth2 = filters.explicitlyKernel(-1, 0, numpy.array([0.5, 0.5]));
        for k in range(srcLevel, destLevel, -1):
            i = self[k-1].__class__(self[k-1].shape, dtype = self[k-1].dtype)
            self.expandImpl(self[k], i, centerValue)
            self[k-1] = i - self[k-1]

    def createLevel(self, level):
        ''' Make sure that 'level' exists. If 'level' is outside the current range of levels,
            empty images of the appropriate shape are inserted into the pyramid.
        '''
        if level > self.highestLevel:
            for i in range(self.highestLevel, level):
                image = list.__getitem__(self, -1)
                newShape = [int((k + 1) / 2) for k in image.shape]
                self.append(image.__class__(newShape, dtype=image.dtype))
            self._highestLevel = level
        elif level < self.lowestLevel:
            image = list.__getitem__(self, 0)
            for i in range(self.lowestLevel, level, -1):
                newShape = [2*k-1 for k in image.shape]
                self.insert(0, image.__class__(newShape, dtype=image.dtype))
            self._lowestLevel = level
             
#################################################################

def _registerArrayTypes():
    from vigranumpycore import registerPythonArrayType
    
    def checkImage(obj):
        return (type(obj) is numpy.ndarray) or (obj.spatialDimensions == 2)
    def checkVolume(obj):
        return (type(obj) is numpy.ndarray) or (obj.spatialDimensions == 3)

    registerPythonArrayType("NumpyArray<2, Singleband<*> >", ScalarImage, checkImage)
    registerPythonArrayType("NumpyArray<2, RGBValue<*> >", RGBImage, checkImage)
    registerPythonArrayType("NumpyArray<2, TinyVector<*, 2> >", Vector2Image, checkImage)
    registerPythonArrayType("NumpyArray<2, TinyVector<*, 3> >", Vector3Image, checkImage)
    registerPythonArrayType("NumpyArray<2, TinyVector<*, 4> >", Vector4Image, checkImage)
    registerPythonArrayType("NumpyArray<3, Multiband<*> >", Image, checkImage)
    registerPythonArrayType("NumpyArray<3, Singleband<*> >", ScalarVolume, checkVolume)
    registerPythonArrayType("NumpyArray<3, RGBValue<*> >", RGBVolume, checkVolume)
    registerPythonArrayType("NumpyArray<3, TinyVector<*, 2> >", Vector2Volume, checkVolume)
    registerPythonArrayType("NumpyArray<3, TinyVector<*, 3> >", Vector3Volume, checkVolume)
    registerPythonArrayType("NumpyArray<3, TinyVector<*, 4> >", Vector4Volume, checkVolume)
    registerPythonArrayType("NumpyArray<3, TinyVector<*, 6> >", Vector6Volume, checkVolume)
    registerPythonArrayType("NumpyArray<4, Multiband<*> >", Volume, checkVolume)

_registerArrayTypes()
del _registerArrayTypes
