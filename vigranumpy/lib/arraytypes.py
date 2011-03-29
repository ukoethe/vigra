#######################################################################
#
#         Copyright 2009-2011 by Ullrich Koethe
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

"""
Vigranumpy can work directly on numpy.ndarrays. However, plain ndarrays do not carry
any information about the semantics of the different coordinate axes. For example,
one cannot distinguish a 2-dimensional RGB image from a scalar volume data set that
happens to contain only three slices. In order to distinguish between arrays that
have the same structure but different interpretation, vigra.arraytypes provides the
following array classes:

    numpy.ndarray
        Image
            ScalarImage
            Vector2Image
            Vector3Image
                RGBImage
            Vector4Image
        Volume
            ScalarVolume
            Vector2Volume
            Vector3Volume
                RGBVolume
            Vector4Volume
            Vector6Volume

    list
        ImagePyramid
"""

import sys
import copy
import numpy
import ufunc
import vigranumpycore

from tagged_array import TaggedArray, preserve_doc

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

from vigranumpycore import AxisType, AxisInfo, AxisTags

def dropChannelAxis(array):
    try:
        return array.dropChannelAxis()
    except:
        return array
                
# How to construct a VigraArray
#
# case 1: from shape and order or axistags
# conventions: - shape has explicit channel axis
#              - 'A' order defaults to 'V' order
#              - order implies axistags and vice versa, you cannot provide both
# * look up the array type. If it is a plain ndarray, skip axistags
# * construct array according to order, optionally init with a constant
# * create and assign normalized axistags, if not explicitly given
# * optionally remove a singleton channel dimension (while we know where it is)
# * permute the array by the inverse normalization
# * assign axistags, if explicitly given (check compatibility)
#
# case 2: from another array
# * if taget order is 'A' or source and target order are equal, copy as is (including axistags)
# * otherwise, normalize the shape according to target order and
#   remember the normalizing permutation
# * construct array in normalized order
# * permute the array by the inverse normalization
# * copy original data and axistags

_constructArrayFromAxistags = vigranumpycore.constructArrayFromAxistags 
    
def _constructArrayFromOrder(cls, shape, dtype, order, init):
    axistags = VigraArray.defaultAxistags(len(shape), order)
    return _constructArrayFromAxistags(cls, shape, dtype, axistags, init)
    
def _constructArrayFromArray(cls, obj, dtype, order, init, axistags):
    if order is None:
        order = 'A'
    if order == 'A':
        # we cannot use ndarray.copy('A') here, because this only preserves 'C' and 'F'
        # order, whereas any other order is silently transformed into 'C'
        
        # we must also make sure that a singleton channel index has the smallest stride
        # (otherwise, strides may not exactly match in the copy)
        strides = list(obj.strides)
        try:
            channelIndex = obj.axistags.index('c')
            if channelIndex < obj.ndim and obj.shape[channelIndex] == 1:
                strides[channelIndex] = 0
        except:
            pass
        permutation = list(numpy.array(strides).argsort())
        norm_shape = tuple(numpy.array(obj.shape)[permutation])
        inverse_permutation = list(numpy.array(permutation).argsort())
        array = numpy.ndarray.__new__(cls, norm_shape, dtype, order='F')
        array = array.transpose(inverse_permutation)
    else:
        array = _constructArrayFromOrder(cls, obj.shape, dtype, order, False)
        
    if init:
        array[...] = obj
    if cls is not numpy.ndarray:
        if axistags is not None:
            array.axistags = axistags
        elif hasattr(array, 'axistags'):
            del array.axistags
    return array
    
def _array_docstring_(name, shape, compat):
    return '''
    Constructor:
    
    .. method:: %(name)s(obj, dtype=numpy.float32, order=None, init = True, value = None)

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

class VigraArray(TaggedArray):
    """
This base class ensures that arrays created in Python are
compatible with the memory layout requirements of
VIGRA's NumpyArray family of C++ views. Do always use
this class via its subclasses!
    """
    
    # IMPORTANT: do not remove this function, it is called from C++
    @classproperty
    def defaultOrder(cls):
        return 'V'

    # IMPORTANT: do not remove this function, it is called from C++
    @staticmethod
    def defaultAxistags(ndim, order=None):
        if order is None:
            order = VigraArray.defaultOrder
        if order == 'F':
            tags = [AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z, AxisInfo()][:ndim]
        elif order == 'C':
            tags = [AxisInfo(), AxisInfo.z, AxisInfo.y, AxisInfo.x, AxisInfo.c][-ndim:]
        else: # order in ['A', 'V']:
            tags = [AxisInfo.x, AxisInfo.y, AxisInfo.z, AxisInfo()][:ndim-1] + [AxisInfo.c]
        return AxisTags(tags)

    # IMPORTANT: do not remove this function, it is called from C++
    @staticmethod
    def copyValuesImpl(target, source):
        try:
            target = target.transposeToNumpyOrder()
            if target.ndim > source.ndim:
                target = target.dropChannelAxis()
        except:
            pass

        try:
            source = source.transposeToNumpyOrder()
            if target.ndim < source.ndim:
                source = source.dropChannelAxis()
        except:
            pass
        
        try:
            compatible = source.axistags.compatible(target.axistags)
        except:
            compatible = True

        if not compatible:
            raise RuntimeError("VigraArray.copyValuesImpl(): incompatible axistags")

        target[...] = source
    
    def __new__(cls, obj, dtype=numpy.float32, order=None, init=True, value=None, axistags=None):
        if value is not None:
            init = False
        if isinstance(obj, numpy.ndarray):
            if axistags is None:
                if hasattr(obj, 'axistags'):
                    axistags = copy.copy(obj.axistags)
                else:
                    raise RuntimeError("VigraArray(): axistags must be given when constructing from plain array.")
            elif obj.ndim != len(axistags):
                raise RuntimeError("VigraArray(): axistags have wrong length.")
            if order is None:
                res = _constructArrayFromAxistags(cls, obj.shape, dtype, axistags, init)
                if init:
                    res[...] = obj
            else:
                res = _constructArrayFromArray(cls, obj, dtype, order, init, axistags)
        else:
            if axistags is None:
                if order is None:
                    order = VigraArray.defaultOrder
            elif len(axistags) != len(obj):
                raise RuntimeError("VigraArray(): axistags have wrong length.")
            if order is None:
                res = _constructArrayFromAxistags(cls, obj, dtype, axistags, init)
            else:
                res = _constructArrayFromOrder(cls, obj, dtype, order, init)
                if cls is not numpy.ndarray and axistags is not None:
                    res.axistags = axistags
        if value is not None:
            res.fill(value)
        return res

    __array_priority__ = 15.0

    def __copy__(self, order='A'):
        return self.copy(order)
        
    @preserve_doc
    def __deepcopy__(self, memo):
        # numpy.ndarray.__deepcopy__ always creates C-order arrays =>
        #   transpose self accordingly, and transpose back after the copy
        result = numpy.ndarray.__deepcopy__(self.transposeToNumpyOrder(), memo)
        result = result.transpose(self.permutationFromNumpyOrder())
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        return result
    
    def copy(self, order='A'):
        return self.__class__(self, dtype=self.dtype, order=order)
    
    def copyValues(self, other):
        self.copyValuesImpl(self, other)
    
    @property
    def channelIndex(self):
        return self.axistags.channelIndex
    
    @property
    def majorNonchannelIndex(self):
        return self.axistags.majorNonchannelIndex
    
    @property
    def channels(self):
        i = self.channelIndex
        if i < self.ndim:
            return self.shape[i]
        else:
            return 1
            
    @property
    def width(self):
        i = self.axistags.index('x')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.width(): axistag 'x' does not exist.")
    
    @property
    def height(self):
        i = self.axistags.index('y')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.height(): axistag 'y' does not exist.")
    
    @property
    def depth(self):
        i = self.axistags.index('z')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.depth(): axistag 'z' does not exist.")
    
    @property
    def duration(self):
        i = self.axistags.index('t')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.duration(): axistag 't' does not exist.")
            
    @property
    def spatialDimensions(self):
        return self.axistags.axisTypeCount(AxisType.Space)
    
    @staticmethod
    def empty_axistags(ndim):
        '''Create an axistags object with non-informative entries.
        '''
        return AxisTags(ndim)
    
    def transpose_axistags(self, permutation = None):
        if hasattr(self, 'axistags'):
            res = copy.copy(self.axistags)
            res.transpose(permutation)
            return res
        else:
            return self.empty_axistags(self.ndim)

    def transform_axistags(self, index):
        if hasattr(self, 'axistags'):
            return self.axistags.transform(index, self.ndim)
        else:
            return self.empty_axistags(self.ndim)

    def permutationToNormalOrder(self, types=AxisType.AllAxes):
        return list(self.axistags.permutationToNormalOrder(types))
    
    def permutationToNumpyOrder(self):
        return list(reversed(self.axistags.permutationToNormalOrder()))
    
    def permutationFromNormalOrder(self):
        return list(self.axistags.permutationFromNormalOrder())
    
    def permutationFromNumpyOrder(self):
        return map(lambda x: self.ndim-1-x, self.axistags.permutationFromNormalOrder())
    
    @property
    def order(self):
        if self.flags.c_contiguous:
            return 'C'
        elif self.flags.f_contiguous:
            return 'F'
        elif self.channelIndex == self.ndim-1 and self.itemsize == self.strides[-1] and \
             reduce(lambda x, y: y if y >= x and x >= 0 else -1, self.strides[:-1], 0) >= 0:
            return 'V'
        return 'A'
    
    @property
    def flat(self):
        return self.transposeToNumpyOrder().view(numpy.ndarray).flat
    
    def flatten(self, order='C'):
        return self.transposeToNumpyOrder().view(numpy.ndarray).flatten(order)        

    def ravel(self, order='C'):
        return self.transposeToNumpyOrder().view(numpy.ndarray).ravel(order)        

    def __str__(self):
        return str(self.transposeToVigraOrder().transpose().view(numpy.ndarray))
    
    def bindAxis(self, which, value=0):
        if type(which) == str:
            which = self.axistags.index(which)
        return self[(slice(None),)*which + (value,) + (slice(None),)*(self.ndim-which-1)]
    
    def dropChannelAxis(self):
        ci = self.channelIndex
        if ci == self.ndim:
            return self
        
        if self.shape[ci] != 1:
            raise RuntimeError("dropChannelAxis(): only allowed when there is a single channel.")
        return self.bindAxis(ci, 0)
    
    def insertChannelAxis(self, order=None):
        ci = self.channelIndex
        if ci != self.ndim:
            return self
        
        if order == 'F':
            res = self[numpy.newaxis,...]
            res.axistags[0] = AxisInfo.c
        else:
            res = self[..., numpy.newaxis]
            res.axistags[-1] = AxisInfo.c
        return res
    
    def transposeToOrder(self, order):
        if order == 'A':
            return self
        permutation = self.permutationToNormalOrder()
        if order == 'C':
            permutation.reverse()
        elif order == 'V':
            if self.channelIndex < self.ndim:
                permutation = permutation[1:] + [permutation[0]]
        return self.transpose(permutation)
    
    def transposeToDefaultOrder(self):
        return self.transposeToOrder(self.defaultOrder)

    def transposeToVigraOrder(self):
        return self.transposeToOrder('V')

    def transposeToNumpyOrder(self):
        return self.transposeToOrder('C')
    
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

    def writeImage(self, filename, dtype = '', compression = ''):
        "Write an image to a file. Consult :func:`vigra.impex.writeImage` for detailed documentation"
        import vigra.impex

        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim != 2:
            raise RuntimeError("VigraArray.writeImage(): array must have 2 non-channel axes.")

        vigra.impex.writeImage(self, writeImage, dtype, compression)
            
    def writeSlices(self, filename_base, filename_ext, dtype = '', compression = ''):
        "Write a volume to a sequence of files. Consult :func:`vigra.impex.writeVolume` for detailed documentation.\n"
        import vigra.impex

        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim != 3:
            raise RuntimeError("VigraArray.writeSlices(): array must have 3 non-channel axes.")

        vigra.impex.writeVolume(self, filename_base, filename_ext, dtype, compression)
            
    def writeHDF5(self, filename, pathInFile, dtype = ''):
        "Write an image or volume to a HDF5 file. Consult :func:`vigra.impex.writeImageToHDF5` and :func:`vigra.impex.writeVolumeToHDF5` for detailed documentation"
        import vigra.impex

        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim == 2:
            vigra.impex.writeImageToHDF5(self, filename, pathInFile, dtype)
        elif ndim == 3:
            vigra.impex.writeVolumeToHDF5(self, filename, pathInFile, dtype)
        else:
            raise RuntimeError("VigraArray.writeHDF5(): array must have 2 or 3 non-channel axes.")

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

        ndim = self.ndim
        channelIndex = self.channelIndex
        if channelIndex < ndim:
            if self.shape[channelIndex] not in [1, 3]:
                raise RuntimeError("VigraArray.show(): array must have 1 or 3 channels.")
            ndim -= 1
        if ndim != 2:
            raise RuntimeError("VigraArray.show(): array must have 2 non-channel axes.")

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
        
        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim != 2:
            raise RuntimeError("VigraArray.qimage(): array must have 2 non-channel axes.")

        yxImage = self.transposeToNumpyOrder()

        if self.channels == 1:
            q = qimage2ndarray.gray2qimage(yxImage.dropChannelAxis(), normalize)
        else:
            q = qimage2ndarray.array2qimage(yxImage, normalize)

        return q

##################################################################

# channelCount == None: array must not have channels
# channelCount == 0:    array can have arbitrary number of channels (including None)
def _adjustShape(shape, order, spatialDimensions, channelCount, axistags, name):
    if order is None:
        order = VigraArray.defaultOrder
    if len(shape) == spatialDimensions:
        if channelCount is not None and channelCount == 0:
            channelCount = 1
        if channelCount:
            if order == 'F':
                shape = (channelCount,) + shape
            else:
                shape = shape + (channelCount,)
    else:        
        if channelCount is None or len(shape) != spatialDimensions + 1:
            raise RuntimeError("%s: input shape has wrong length." % name)
        if channelCount > 0:
            if order == 'F':
                if shape[0] != channelCount:
                    raise RuntimeError("%s: input shape has wrong number of channels." % name)
            else:
                if shape[-1] != channelCount:
                    raise RuntimeError("%s: input shape has wrong number of channels." % name)
    if axistags is None:
        axistags = VigraArray.defaultAxistags(spatialDimensions+1, order)
    if len(shape) == spatialDimensions:
        axistags.dropChannelAxis()
    if len(shape) != len(axistags):
        raise RuntimeError("%s: size mismatch between shape and axistags." % name)
    return shape, axistags

def _adjustArray(array, order, spatialDimensions, channelCount, axistags, name):
    if order is None:
        order = VigraArray.defaultOrder
    if array.ndim == spatialDimensions:
        if channelCount is not None and channelCount > 1:
            raise RuntimeError("%s: input array has too few dimensions." % name)
        if hasattr(array, 'axistags'):
            if array.channelIndex != array.ndim:
                raise RuntimeError("%s: input array has too few non-channel axes." % name)
        if channelCount:
            if hasattr(array, 'axistags'):
                array = array.insertChannelAxis(order)
            elif order == 'F':
                array = array[numpy.newaxis,...]
            else:
                array = array[...,numpy.newaxis]
    else:
        if channelCount is None or array.ndim != spatialDimensions+1:
            raise RuntimeError("%s: input array has wrong number of dimensions." % name)
        if hasattr(array, 'axistags'):
            channelIndex = array.channelIndex
            if channelIndex == array.ndim:
                raise RuntimeError("%s: input array has no channel axis." % name)
            if channelCount > 0 and array.shape[channelIndex] != channelCount:
                raise RuntimeError("%s: input array has wrong number of channels." % name)
    if axistags is None:
        if hasattr(array, 'axistags'):
            axistags = copy.copy(array.axistags)
        else:
            axistags = VigraArray.defaultAxistags(spatialDimensions+1, order)
    if array.ndim == spatialDimensions:
        axistags.dropChannelAxis()
    if array.ndim != len(axistags):
        raise RuntimeError("%s: axistags have wrong number of axes." % name)
    return array, axistags

def _adjustInput(obj, order, spatialDimensions, channelCount, axistags, name):
    if isinstance(obj, numpy.ndarray):
        return _adjustArray(obj, order, spatialDimensions, channelCount, axistags, name)
    else:
        return _adjustShape(obj, order, spatialDimensions, channelCount, axistags, name)

#################################################################

def Image(obj, dtype=numpy.float32, order=None, 
          init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, 0, axistags, "vigra.Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)
        
def ScalarImage(obj, dtype=numpy.float32, order=None, 
                init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, None, axistags, "vigra.ScalarImage()")
    return VigraArray(obj, dtype, None, init, value, axistags)
        
def Vector2Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, 2, axistags, "vigra.Vector2Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector3Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, 3, axistags, "vigra.Vector3Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector4Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, 4, axistags, "vigra.Vector4Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def RGBImage(obj, dtype=numpy.float32, order=None, 
             init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 2, 3, axistags, "vigra.RGBImage()")
    res = VigraArray(obj, dtype, None, init, value, axistags)
    res.axistags.setChannelDescription('RGB')
    return res

#################################################################

def Volume(obj, dtype=numpy.float32, order=None, 
           init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 0, axistags, "vigra.Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def ScalarVolume(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, None, axistags, "vigra.ScalarVolume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector2Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 2, axistags, "vigra.Vector2Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector3Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 3, axistags, "vigra.Vector3Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector4Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 4, axistags, "vigra.Vector4Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector6Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 6, axistags, "vigra.Vector6Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def RGBVolume(obj, dtype=numpy.float32, order=None, 
              init=True, value=None, axistags=None):
    obj, axistags = _adjustInput(obj, order, 3, 3, axistags, "vigra.RGBVolume()")
    res = VigraArray(obj, dtype, None, init, value, axistags)
    res.axistags.setChannelDescription('RGB')
    return res

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
             