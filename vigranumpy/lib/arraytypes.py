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

import copy
import numpy
import ufunc
import vigranumpycore

from tagged_array import TaggedArray

from vigranumpycore import AxisType, AxisInfo, AxisTags

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

def AxisInfo_copy(axisinfo):
    return AxisInfo(axisinfo)

def AxisInfo_deepcopy(axisinfo, memo):
    result = AxisInfo(axisinfo)
    memo[id(axisinfo)] = result
    return result

AxisInfo.__copy__ = AxisInfo_copy
AxisInfo.__deepcopy__ = AxisInfo_deepcopy

del AxisInfo_copy
del AxisInfo_deepcopy

def AxisTags_copy(axistags):
    return AxisTags(axistags)

def AxisTags_deepcopy(axistags, memo):
    taglist = [k for k in axistags]
    taglist = copy.deepcopy(taglist, memo)
    result = AxisTags(taglist)
    memo[id(axistags)] = result
    return result

AxisTags.__copy__ = AxisTags_copy
AxisTags.__deepcopy__ = AxisTags_deepcopy

del AxisTags_copy
del AxisTags_deepcopy

defaultOrder = 'V'

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

def normalizeShape(shape, order):
    if order == 'F':
        norm_shape = shape
        inverse_permutation = range(len(shape))
    elif order == 'C':
        norm_shape = tuple(reversed(shape))
        inverse_permutation = range(len(shape)-1, -1, -1)
    else: # order in ['A', 'V']:
        norm_shape = (shape[-1],) + shape[:-1]
        inverse_permutation = range(1, len(shape)) + [0]
    return norm_shape, inverse_permutation
    
def constructArrayFromOrder(cls, shape, dtype, order, init):
    norm_shape, inverse_permutation = normalizeShape(shape, order)
    
    axistags = AxisTags([AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z][:len(shape)])    
    array = TaggedArray.__new__(cls, norm_shape, dtype, order='F', axistags=axistags)
    if init:
        array.fill(0)
    return array.transpose(inverse_permutation)
    
def constructArrayFromAxistags(cls, shape, dtype, axistags, init):
    permutation = list(numpy.array(map(lambda x: ord(x.key[-1]), axistags)).argsort())
    norm_shape = tuple(numpy.array(shape)[permutation])
    inverse_permutation = list(numpy.array(permutation).argsort())

    array = numpy.ndarray.__new__(cls, norm_shape, dtype, order='F')
    array = array.transpose(inverse_permutation)
    if init:
        array.fill(0)
    if cls is not numpy.ndarray:
        array.axistags = axistags
    return array
    
def constructArrayFromArray(cls, obj, dtype, order, init):
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
        array = constructArrayFromOrder(cls, obj.shape, dtype, order, False)
        
    if init:
        array[...] = obj
    if hasattr(obj, 'axistags') and cls is not numpy.ndarray:
        array.axistags = copy.copy(obj.axistags)
    return array
    
def dropChannelDimension(array):
    try:
        channelIndex = array.axistags.index('c')
    except:
        return array
    
    if channelIndex < array.ndim:
        if array.shape[channelIndex] != 1:
            raise RuntimeError("dropChannelDimension(): only allowed when there is a single channel.")
        return array[(slice(None),)*channelIndex + (0,) + (slice(None),)*(array.ndim-channelIndex-1)]
    else:
        return array
        
def addChannelDimension(shape, channels, order):
    if order == 'F':
        return (channels,) + shape
    else:
        return shape + (channels,)
        
def checkChannelCount(shape, channels, order):
    if order == 'F':
        return shape[0] == channels
    else:
        return shape[-1] == channels
        
def _array_docstring_(name, shape, compat):
    return '''
    Constructor:
    
    .. method:: %(name)s(obj, dtype=numpy.float32, order=defaultOrder, init = True, value = None)

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

def constructNumpyArrayOld(cls, obj, spatialDimensions, channels, dtype, order, init):
    if isinstance(obj, numpy.ndarray):
        shape = list(obj.shape)
        strideOrdering = list(numpy.array(obj.strides).argsort().argsort())
    else:
        shape = list(obj)
        strideOrdering = None
    if channels == 0: # if the requested number of channels is not given ...
        # ... deduce it
        if len(shape) == spatialDimensions:
            channels = 1
        else:
            channels = shape[-1]

    # if we have only one channel, no explicit channel dimension should be in the shape
    shapeSize = spatialDimensions if channels == 1 else spatialDimensions + 1
    shape.append(0)

    # create the shape object with optional channel dimension
    pshape = shape[:shapeSize]
    if shapeSize > spatialDimensions:
        pshape[-1] = channels

    # order "A" means "preserve order" when an array is copied, and
    # defaults to "V" when a new array is created without explicit strideOrdering
    if order == "A":
        if strideOrdering is None:
            order = "V"
        elif len(strideOrdering) > shapeSize:
            # make sure that strideOrdering length matches shape length
            pstride = strideOrdering[:shapeSize]

            # adjust the ordering when the channel dimension has been dropped because channel == 1
            if strideOrdering[shapeSize] == 0:
                pstride = [k-1 for k in pstride]
            strideOrdering = pstride
        elif len(strideOrdering) < shapeSize:
            # make sure that strideOrdering length matches shape length
            # adjust the ordering when the channel dimension has been dropped because channel == 1
            strideOrdering = [k+1 for k in strideOrdering]
            strideOrdering.append(0)

    # create the appropriate strideOrdering objects for the other memory orders
    # (when strideOrdering already contained data, it is ignored because order != "A")
    if order == "C":
        strideOrdering = range(len(pshape)-1, -1, -1)
    elif order == "F" or (order == "V" and channels == 1):
        strideOrdering = range(len(pshape))
    elif order == "V":
        strideOrdering = range(1, len(pshape)+1)
        strideOrdering[-1] = 0
        
    ppshape = [0]*len(pshape)
    for k in xrange(len(pshape)):
        ppshape[strideOrdering[k]] = pshape[k]
    
    res = TaggedArray.__new__(cls, ppshape, dtype, order='F')
    res = res.transpose(strideOrdering)

    if init:
        if isinstance(obj, numpy.ndarray):
            res[...] = obj
        else:
            res[...] = 0
    return res
        
def constructNumpyArray(cls, obj, spatialDimensions, channels, dtype, order, init):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim == spatialDimensions:
            if order == 'F':
                obj = obj[numpy.newaxis, ...]
            else:
                obj = obj[..., numpy.newaxis]
        shape = list(obj.shape)
    else:
        shape = list(obj)
        if order == 'A':
            order = 'V'
        if len(shape) == spatialDimensions:
            if order == 'F':
                shape.insert(0, 1 if channels == 0 else channels)
            else:
                shape.append(1 if channels == 0 else channels)

    # create the appropriate strideOrdering objects
    if order == "C":
        strideOrdering = range(len(shape)-1, -1, -1)
    elif order == "F":
        strideOrdering = range(len(shape))
    elif order == "V":
        strideOrdering = range(1, len(shape)+1)
        strideOrdering[-1] = 0
    elif order == "A":
        # this is only reached when obj is an ndarray
        strideOrdering = numpy.array(obj.strides).argsort()
        # make sure that the channel dimension is listed first if it has minimum stride
        # (this is not always automatically the case if there is only 1 channel)
        if obj.strides[strideOrdering[0]] == obj.strides[strideOrdering[1]]:
            strideOrdering[0], strideOrdering[1] = (max(strideOrdering[0], strideOrdering[1]), 
                                                     min(strideOrdering[0], strideOrdering[1]))
        strideOrdering = list(strideOrdering.argsort())
        
    pshape = [0]*len(shape)
    for k in xrange(len(shape)):
        pshape[strideOrdering[k]] = shape[k]
    
    # we construct the array in 'F' order because strideOrdering is ascending
    res = TaggedArray.__new__(cls, pshape, dtype, order='F')
    res = res.transpose(strideOrdering)

    if init:
        if isinstance(obj, numpy.ndarray):
            res[...] = obj
        else:
            res[...] = 0
    return res
        

class VigraArray(TaggedArray):
    """
This base class ensures that arrays created in Python are
compatible with the memory layout requirements of
VIGRA's NumpyArray family of C++ views. Do always use
this class via its subclasses!
    """
    # def __new__(cls, obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None, axistags=None):
        # if axistags is None and hasattr(obj, 'axistags'):
            # axistags = obj.axistags
        # if axistags is None:
            # try:
                # spatialDimensions = obj.ndim
            # except:
                # spatialDimensions = len(obj)
        # else:
            # spatialDimensions = axistags.axisTypeCount(AxisType.Space)
        # if isinstance(obj, numpy.ndarray) and not isinstance(obj, VigraArray):
            # obj = obj.swapaxes(0, spatialDimensions-1)
        # channels = 0
        # if value is not None:
            # res = constructNumpyArray(cls, obj, spatialDimensions, channels, dtype, order, False)
            # res[...] = value
        # else:
            # res = constructNumpyArray(cls, obj, spatialDimensions, channels, dtype, order, init)
        # # FIXME: vigranumpycore.constructNumpyArray() should lead to the same results and should
        # #        be preferred for the sake of consistence (currently, we get different strides for 
        # #        singleton axes)
        # # if value is not None:
            # # res = vigranumpycore.constructNumpyArray(cls, obj, spatialDimensions, channels, dtype, order, False)
            # # res[...] = value
        # # else:
            # # res = vigranumpycore.constructNumpyArray(cls, obj, spatialDimensions, channels, dtype, order, init)
        # if len(axistags) != res.ndim:
            # raise RuntimeError("VigraArray(): len(axistags) must match ndim.")
        # if axistags is not None:
            # res.axistags = copy.copy(axistags)
        # return res

    def __new__(cls, obj, dtype=numpy.float32, order=None, init=True, value=None, axistags=None):
        if order is not None and axistags is not None:
            raise RuntimeError("VigraArray(): You cannot provide both order and axistags.")
        if value is not None:
            init = False
        if isinstance(obj, numpy.ndarray):
            if order is None:
                order = 'A'
            res = constructArrayFromArray(cls, obj, dtype, order, init)
            if axistags is not None and cls is not numpy.ndarray:
                if len(axistags) != obj.ndim:
                    raise RuntimeError("VigraArray(): axistags have wrong length.")
                res.axistags = axistags
        elif axistags is not None:
            if len(axistags) != len(obj):
                raise RuntimeError("VigraArray(): axistags have wrong length.")
            res = constructArrayFromAxistags(cls, obj, dtype, axistags, init)
        else:
            if order is None:
                order = defaultOrder
            res = constructArrayFromOrder(cls, obj, dtype, order, init)
        if value is not None:
            res.fill(value)
        return res

    __array_priority__ = 15.0

    def __copy__(self, order='A'):
        return self.copy(order)
    
    def copy(self, order='A'):
        return self.__class__(self, dtype=self.dtype, order=order)
    
    @property
    def channelIndex(self):
        return self.axistags.index('c')
    
    @property
    def majorNonchannelIndex(self):
        # FIXME: this must be generalized to the case when 'x' is not present.
        return self.axistags.index('x')
    
    @property
    def canonicalOrdering(self):
        # FIXME: this must be generalized to arbitrary axistags.
        #        (implement as a method in AxisTags)
        return list(numpy.array(map(lambda x: ord(x.key[-1]), self.axistags)).argsort())
    
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
    def timeSteps(self):
        i = self.axistags.index('t')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.timeSteps(): axistag 't' does not exist.")
            
    @property
    def spatialDimensions(self):
        return self.axistags.axisTypeCount(AxisType.Space)

    def default_axistags(self):
        '''Create an axistags object with non-informative entries.
        '''
        return AxisTags(self.ndim)
    
    def transform_axistags(self, index):
        if hasattr(self, 'axistags'):
            return self.axistags.transform(index, self.ndim)
        else:
            return self.default_axistags()

    @property
    def order(self):
        if self.flags.c_contiguous:
            return 'C'
        elif self.flags.f_contiguous:
            return 'F'
        # FIXME: this should use axistags
        elif self.itemsize == self.strides[-1] and \
             reduce(lambda x, y: y if y >= x and x >= 0 else -1, self.strides[:-1], 0) >= 0:
            return 'V'
        return 'A'
    
    @property
    # FIXME: this should depend on axistags
    def flat(self):
        return self.view(numpy.ndarray).swapaxes(0, self.spatialDimensions-1).flat
    
    # FIXME: this should depend on axistags
    def flatten(self, order='C'):
        return self.view(TaggedArray).swapaxes(0, self.spatialDimensions-1).flatten(order)        

    # FIXME: this should depend on axistags
    def ravel(self, order='C'):
        return self.view(TaggedArray).swapaxes(0, self.spatialDimensions-1).ravel(order)        

    # FIXME: to be implemented
    # def __str__(self):
    
    # def __repr__(self):
    
    def transposeToOrder(self, order = 'C'):
        if order == 'A':
            return self
        permutation = [int(k) for k in numpy.array(self.strides).argsort()]
        if order == 'C':
            permutation.reverse()
        elif order == 'V':
            if hasattr(self, 'axistags'):
                permutation = self.axistags.canonicalOrdering()
            else:
                permutation.reverse()
                d = self.spatialDimensions - 1
                permutation[0], permutation[d] = permutation[d], permutation[0]
        return self.transpose(permutation)
    
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

##################################################################

class Image(VigraArray):
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
            
    def __new__(cls, obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None, axistags=None):
        if isinstance(obj, numpy.ndarray):
            if obj.ndim != 2 and obj.ndim != 3:
                raise RuntimeError("Image(): shape mismatch")
        elif len(obj) == 2:
            if order == 'F':
                obj = (1,) + obj
            else:
                obj += (1,)
        return VigraArray.__new__(cls, obj, dtype, order, init, value, axistags)
    
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
        
    # @property
    # def width(self):
        # """the image's width"""
        # return self.shape[0]
    
    # @property
    # def height(self):
        # "the image's height"
        # return self.shape[1]

    # @classproperty
    # def spatialDimensions(cls):
        # "number of spatial dimensions (useful for distinguishing RGBImage and ScalarVolume)"
        # return 2

# class ScalarImage(Image):
    # __doc__ = _array_docstring_('ScalarImage', '''A shape is compatible when it has two dimensions (width, height) or three dimensions (width, height, 1).''', """
          # 'C':
             # | NumpyArray<2, T, StridedArrayTag>,
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<2, Singleband<T>, StridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<2, T, UnstridedArrayTag>,
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<2, Singleband<T>, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | like 'F'""")

    # channels = classproperty(lambda cls: 1, Image.bands)
        
def ScalarImage(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 2 and obj.ndim != 3:
            raise RuntimeError("ScalarImage(): shape mismatch")
    elif len(obj) == 2:
        if order == 'F':
            obj = (1,) + obj
        else:
            obj += (1,)
    res = Image(obj, dtype, order, init, value)
    # FIXME: activate this after test refactoring
    # res = dropChannelDimension(res)
    return res
        
# class Vector2Image(Image):
    # __doc__ = _array_docstring_('Vector2Image', '''
    # A shape is compatible when it has two dimensions (width, height)
    # or three dimensions (width, height, 2).''', """
          # 'C':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 2>, StridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 2>, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 2, Image.bands)

def Vector2Image(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 3:
            raise RuntimeError("Vector2Image(): shape mismatch")
    elif len(obj) == 2:
        if order == 'F':
            obj = (2,) + obj
        else:
            obj += (2,)
    return Image(obj, dtype, order, init, value)

# class Vector3Image(Image):
    # __doc__ = _array_docstring_('Vector3Image', '''
    # A shape is compatible when it has two dimensions (width, height)
    # or three dimensions (width, height, 3).''', """
          # 'C':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 3, Image.bands)

def Vector3Image(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 3:
            raise RuntimeError("Vector3Image(): shape mismatch")
    elif len(obj) == 2:
        if order == 'F':
            obj = (3,) + obj
        else:
            obj += (3,)
    return Image(obj, dtype, order, init, value)

# class Vector4Image(Image):
    # __doc__ = _array_docstring_('Vector4Image', '''
    # A shape is compatible when it has two dimensions (width, height)
    # or three dimensions (width, height, 4).''', """
          # 'C':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 4>, StridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 4>, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 4, Image.bands)

def Vector4Image(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 3:
            raise RuntimeError("Vector4Image(): shape mismatch")
    elif len(obj) == 2:
        if order == 'F':
            obj = (4,) + obj
        else:
            obj += (4,)
    return Image(obj, dtype, order, init, value)

# class RGBImage(Vector3Image):
    # __doc__ = _array_docstring_('RGBImage', '''A shape is compatible when it has two dimensions (width, height) or three dimensions (width, height, 3).''', """
          # 'C':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, RGBValue<T>, StridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 3>, StridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<2, RGBValue<T>, UnstridedArrayTag>,
             # | NumpyArray<2, TinyVector<T, 3>, UnstridedArrayTag>,
             # | NumpyArray<3, Multiband<T>, StridedArrayTag>""")

def RGBImage(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 3:
            raise RuntimeError("RGBImage(): shape mismatch")
    elif len(obj) == 2:
        if order == 'F':
            obj = (3,) + obj
        else:
            obj += (3,)
    return Image(obj, dtype, order, init, value)

#################################################################

class Volume(VigraArray):
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
            
    def __new__(cls, obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None, axistags=None):
        if isinstance(obj, numpy.ndarray):
            if obj.ndim != 3 and obj.ndim != 4:
                raise RuntimeError("Volume(): shape mismatch")
        elif len(obj) == 3:
            if order == 'F':
                obj = (1,) + obj
            else:
                obj += (1,)
        return VigraArray.__new__(cls, obj, dtype, order, init, value, axistags)
    
    def write(self, filename_base, filename_ext, dtype = '', compression = ''):
        "Write a volume to a sequence of files. Consult :func:`vigra.impex.writeVolume` for detailed documentation.\n"
        import vigra.impex
        vigra.impex.writeVolume(self, filename_base, filename_ext, dtype, compression)
            
    def writeHDF5(self, filename, pathInFile, dtype = ''):
        "Write a volume to a HDF5 file. Consult :func:`vigra.impex.writeVolumeToHDF5` for detailed documentation.\n"
        import vigra.impex
        vigra.impex.writeVolumeToHDF5(self, filename, pathInFile, dtype)
    
    # @classproperty
    # def spatialDimensions(cls): return 3
        
    # @property
    # def width(self):
        # return self.shape[0]
    
    # @property
    # def height(self):
        # return self.shape[1]
    
    # @property
    # def depth(self):
        # return self.shape[2]

# class ScalarVolume(Volume):
    # __doc__ = _array_docstring_('ScalarVolume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 1).''', """
          # 'C':
             # | NumpyArray<3, T, StridedArrayTag>,
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, Singleband<T>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<3, T, UnstridedArrayTag>,
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<3, Singleband<T>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | like 'F'""")

    # channels = classproperty(lambda cls: 1, Volume.bands)

def ScalarVolume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 3 and obj.ndim != 4:
            raise RuntimeError("ScalarVolume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (1,) + obj
        else:
            obj += (1,)
    res = Volume(obj, dtype, order, init, value)
    # FIXME: activate this after test refactoring
    # res = dropChannelDimension(res)
    return res

# class Vector2Volume(Volume):
    # __doc__ = _array_docstring_('Vector2Volume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 2).''', """
          # 'C':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 2>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 2>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 2, Volume.bands)

def Vector2Volume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 4:
            raise RuntimeError("Vector2Volume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (2,) + obj
        else:
            obj += (2,)
    return Volume(obj, dtype, order, init, value)

# class Vector3Volume(Volume):
    # __doc__ = _array_docstring_('Vector3Volume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 3).''', """
          # 'C':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 3, Volume.bands)

def Vector3Volume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 4:
            raise RuntimeError("Vector3Volume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (3,) + obj
        else:
            obj += (3,)
    return Volume(obj, dtype, order, init, value)

# class Vector4Volume(Volume):
    # __doc__ = _array_docstring_('Vector4Volume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 4).''', """
          # 'C':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 4>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 4>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 4, Volume.bands)
    
def Vector4Volume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 4:
            raise RuntimeError("Vector4Volume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (4,) + obj
        else:
            obj += (4,)
    return Volume(obj, dtype, order, init, value)


# class Vector6Volume(Volume):
    # __doc__ = _array_docstring_('Vector4Volume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 6).''', """
          # 'C':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 6>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 6>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>""")

    # channels = classproperty(lambda cls: 6, Volume.bands)
    
def Vector6Volume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 4:
            raise RuntimeError("Vector6Volume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (6,) + obj
        else:
            obj += (6,)
    return Volume(obj, dtype, order, init, value)

# class RGBVolume(Vector3Volume):
    # __doc__ = _array_docstring_('RGBVolume', '''
    # A shape is compatible when it has three dimensions (width, height,
    # depth) or four dimensions (width, height, depth, 3).''', """
          # 'C':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, RGBValue<T>, StridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 3>, StridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>
          # 'F':
             # | NumpyArray<4, T, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, UnstridedArrayTag>
          # 'V':
             # | NumpyArray<4, T, StridedArrayTag>,
             # | NumpyArray<3, RGBValue<T>, UnstridedArrayTag>,
             # | NumpyArray<3, TinyVector<T, 3>, UnstridedArrayTag>,
             # | NumpyArray<4, Multiband<T>, StridedArrayTag>""")


def RGBVolume(obj, dtype=numpy.float32, order=defaultOrder, init=True, value=None):
    if isinstance(obj, numpy.ndarray):
        if obj.ndim != 4:
            raise RuntimeError("RGBVolume(): shape mismatch")
    elif len(obj) == 3:
        if order == 'F':
            obj = (3,) + obj
        else:
            obj += (3,)
    return Volume(obj, dtype, order, init, value)

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

# def _registerArrayTypes():
    # from vigranumpycore import registerPythonArrayType
    
    # def checkImage(obj):
        # return (type(obj) is numpy.ndarray) or (obj.spatialDimensions == 2)
    # def checkVolume(obj):
        # return (type(obj) is numpy.ndarray) or (obj.spatialDimensions == 3)

    # registerPythonArrayType("NumpyArray<2, Singleband<*> >", ScalarImage, checkImage)
    # registerPythonArrayType("NumpyArray<2, RGBValue<*> >", RGBImage, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 2> >", Vector2Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 3> >", Vector3Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 4> >", Vector4Image, checkImage)
    # registerPythonArrayType("NumpyArray<3, Multiband<*> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<3, Singleband<*> >", ScalarVolume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, RGBValue<*> >", RGBVolume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 2> >", Vector2Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 3> >", Vector3Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 4> >", Vector4Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 6> >", Vector6Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<4, Multiband<*> >", Volume, checkVolume)

# _registerArrayTypes()
# del _registerArrayTypes

# def _registerArrayTypes():
    # from vigranumpycore import registerPythonArrayType
    
    # def checkImage(obj):
        # return True
    # def checkVolume(obj):
        # return True

    # registerPythonArrayType("NumpyArray<2, Singleband<*> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, RGBValue<*> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 2> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 3> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<2, TinyVector<*, 4> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<3, Multiband<*> >", Image, checkImage)
    # registerPythonArrayType("NumpyArray<3, Singleband<*> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, RGBValue<*> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 2> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 3> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 4> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<3, TinyVector<*, 6> >", Volume, checkVolume)
    # registerPythonArrayType("NumpyArray<4, Multiband<*> >", Volume, checkVolume)

# _registerArrayTypes()
# del _registerArrayTypes
