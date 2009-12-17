import copy
import numpy
import ufunc
import sys

from numpy import int8, uint8, int16, uint16, int32, uint32, int64, uint64
from numpy import float32, float64, longdouble, complex64, complex128, clongdouble

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
        * If obj is a numpy.ndarray with compatible shape, a copy
          of obj with the given dtype, order and class vigra.%(name)s is 
          created.
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
    '''
This base class ensures that arrays created in Python are
compatible with the memory layout requirements of
VIGRA's NumpyArray family of C++ views. Do always use
this class via its subclasses!
    '''
    def __new__(cls, obj, dtype=numpy.float32, order='V', init = True, value = None):
        from vigranumpycmodule import constructNumpyArray
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
    
    def flatten(self):
        return self.view(numpy.ndarray).T.flatten()
    
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
        return self.view(numpy.ndarray).T.flat
    
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
    __doc__ = _array_docstring_('Image', '''
        A shape is compatible when it has two dimensions (width, height)
        or three dimensions (width, height, channels).''', """

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
            
    def write(self, filename, export_type = '', compression = ''):
        "consult :func:`vigra.writeImage` for detailed documentation"
        import vigranumpycmodule as vn
        vn.writeImage(self, filename, export_type, compression)

    def qimage(self, normalize = True):
        '''
        Convert this image to a Qt QImage (mainly for display purposes).
        The present image must have 1 or 3 channels, and the resulting
        QImage will have QImage.Format_Indexed8 or QImage.Format_RGB32
        respectively.
        
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
        import PyQt4.QtGui as qt
        import qimage2ndarray

        if self.channels not in [1,3]:
            raise RuntimeError("Image.qimage(): channels == 1 or channels == 3 required.")
        
        if normalize is None or normalize is False:
            nmin, nmax = 0.0, 255.0
        elif normalize is True:
            nmin, nmax = float(self.min()), float(self.max())
        else:
            try:
                nmax = float(normalize)
                nmin = 0.0
            except:
                nmin, nmax = map(float, normalize)
        if nmax < nmin:
            raise RuntimeError("Image.qimage(): invalid normalization (nmax < nmin).")
        
        if self.channels == 1:
            q = qt.QImage(self.width, self.height, qt.QImage.Format_Indexed8)
            for i in range(256):
                q.setColor(i, qt.QColor(i,i,i).rgba())
            if nmax == nmin:
                q.fill(0)
            else:
                # FIXME: use proper rounding
                ufunc.multiply(self - nmin, 255.0 / (nmax - nmin),
                               qimage2ndarray.byte_view(q).swapaxes(0,1).reshape(self.shape))
        else:
            q = qt.QImage(self.width, self.height, qt.QImage.Format_RGB32)
            if nmax == nmin:
                q.fill(0)
            else:
                ufunc.multiply(self - nmin, 255.0 / (nmax - nmin),
                               qimage2ndarray.rgb_view(q).swapaxes(0,1))
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
    __doc__ = _array_docstring_('ScalarImage', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 1).''', """
    * 'C': NumpyArray<2, T, StridedArrayTag>,
           NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<2, Singleband<T>, StridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<2, T, UnstridedArrayTag>,
           NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<2, Singleband<T>, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, UnstridedArrayTag>
    * 'V': like 'F'""")

    channels = classproperty(lambda cls: 1, Image.bands)
        
class Vector2Image(Image):
    __doc__ = _array_docstring_('Vector2Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 2).''', """
    * 'C': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 2>, StridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 2>, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 2, Image.bands)

class Vector3Image(Image):
    __doc__ = _array_docstring_('Vector3Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 3).''', """
    * 'C': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 3, Image.bands)

class Vector4Image(Image):
    __doc__ = _array_docstring_('Vector4Image', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 4).''', """
    * 'C': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 4>, StridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 4>, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 4, Image.bands)

class RGBImage(Vector3Image):
    __doc__ = _array_docstring_('RGBImage', '''
    A shape is compatible when it has two dimensions (width, height)
    or three dimensions (width, height, 3).''', """
    * 'C': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, RGBValue<T>, StridedArrayTag>,
           NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<2, RGBValue<T>, UnstridedArrayTag>,
           NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
           NumpyArray<3, Multiband<T>, StridedArrayTag>""")

#################################################################

class Volume(_VigraArray):
    __doc__ = _array_docstring_('Volume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, channels).''', """
    * 'C': NumpyArray<3, T, StridedArrayTag> (if channels=1),
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, M>, StridedArrayTag> (if channels=M),
           NumpyArray<3, RGBValue<T>, StridedArrayTag> (if channels=3),
           NumpyArray<3, Singleband<T>, StridedArrayTag> (if channels=1),
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag> (if channels=1),
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Singleband<T>, UnstridedArrayTag> (if channels=1),
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<3, T, UnstridedArrayTag> (if channels=1),
           NumpyArray<4, T, UnstridedArrayTag> (if channels=1),
           NumpyArray<4, T, StridedArrayTag> (if channels>1),
           NumpyArray<3, Singleband<T>, UnstridedArrayTag> (if channels=1),
           NumpyArray<3, TinyVector<T, M>, UnstridedArrayTag> (if channels=M),
           NumpyArray<3, RGBValue<T>, UnstridedArrayTag> (if channels=3),
           NumpyArray<4, Multiband<T>, UnstridedArrayTag> (if channels=1),
           NumpyArray<4, Multiband<T>, StridedArrayTag> (if channels>1)""")
    
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
    * 'C': NumpyArray<3, T, StridedArrayTag>,
           NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, Singleband<T>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<3, T, UnstridedArrayTag>,
           NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<3, Singleband<T>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': like 'F'""")

    channels = classproperty(lambda cls: 1, Volume.bands)

class Vector2Volume(Volume):
    __doc__ = _array_docstring_('Vector2Volume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, 2).''', """
    * 'C': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 2>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 2>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 2, Volume.bands)

class Vector3Volume(Volume):
    __doc__ = _array_docstring_('Vector3Volume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, 3).''', """
    * 'C': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 3, Volume.bands)

class Vector4Volume(Volume):
    __doc__ = _array_docstring_('Vector4Volume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, 4).''', """
    * 'C': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 4>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 4>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 4, Volume.bands)
    
class Vector6Volume(Volume):
    __doc__ = _array_docstring_('Vector4Volume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, 6).''', """
    * 'C': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 6>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 6>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    channels = classproperty(lambda cls: 6, Volume.bands)
    
class RGBVolume(Vector3Volume):
    __doc__ = _array_docstring_('RGBVolume', '''
    A shape is compatible when it has three dimensions (width, height, 
    depth) or four dimensions (width, height, depth, 3).''', """
    * 'C': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, RGBValue<T>, StridedArrayTag>,
           NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>
    * 'F': NumpyArray<4, T, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, UnstridedArrayTag>
    * 'V': NumpyArray<4, T, StridedArrayTag>,
           NumpyArray<3, RGBValue<T>, UnstridedArrayTag>,
           NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
           NumpyArray<4, Multiband<T>, StridedArrayTag>""")

def _registerArrayTypes():
    from vigranumpycmodule import registerPythonArrayType
    
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
    registerPythonArrayType("NumpyArray<4, Multiband<*> >", Volume, checkVolume)

_registerArrayTypes()
