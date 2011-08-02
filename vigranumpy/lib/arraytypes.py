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

import sys
import copy
import numpy
import ufunc
import vigranumpycore

from vigranumpycore import AxisType, AxisInfo, AxisTags

def _preserve_doc(f):
    f.__doc__ = eval('numpy.ndarray.%s.__doc__' % f.__name__)
    return f

# a decorator to finalize the return value of a 
# dimension-reducing function (e.g. array.max())
def _finalize_reduce_result(f):
    def new_f(self, axis=None, *args, **kw):
        if type(axis) == str:
            axis = self.axistags.index(axis)
        res = f(self, axis, *args, **kw)
        if axis is not None:
            res.axistags = res._copy_axistags()
            del res.axistags[axis]
        else:
            # this 'else' is necessary because numpy 1.6.0 gives 
            #     type(res) == type(self)
            # instead of the desired
            #     type(res) == self.dtype
            # when res is scalar and self is a subclass of ndarray 
            # (this is probably a bug in numpy, since it works correctly 
            #  when self is a plain ndarray)
            res = res.dtype.type(res)
        return res
    return new_f

class classproperty(object):
    def __get__(self, instance, cls):
            if self.__instance_method is not None and instance is not None:
                return self.__instance_method(instance)
            else:
                return self.__class_method(cls)
    def __init__(self, class_method, instance_method=None):
            self.__class_method = class_method
            self.__instance_method = instance_method

def newaxis(axisinfo=AxisInfo()):
    '''
    Create a new singleton axis via the indexing operator. This works similar to 
    `numpy.newaxis`, but allows to provide an AxisInfo object for the new axis. 
    For example::
    
        >>> s = vigra.ScalarImage((width, height))
        >>> s.axistags  # no channel axis
        x y
        >>> t = s[..., numpy.newaxis]
        >>> t.axistags  # with unknown axis type
        x y ?
        >>> t = s[..., vigra.newaxis(vigra.AxisInfo.c)]
        >>> t.axistags  # with channel axis 
        x y c
    '''
    return axisinfo

def taggedView(array, axistags):
    '''
    Create a view to the given array with type :class:`~vigra.VigraArray` and the 
    given axistags. This is a shorthand for::
    
        >>> view = array.view(vigra.VigraArray)
        >>> view.axistags = copy.copy(axistags)
    '''
    if array.ndim != len(axistags):
        raise RuntimeError('vigra.taggedView(): array.ndim must match len(axistags).')
    res = array.view(VigraArray)
    res.axistags = copy.copy(axistags)
    return res

def dropChannelAxis(array):
    '''
    Return the view created by ``array.``:meth:`~vigra.VigraArray.dropChannelAxis` if 
    the given array supports that function, or return ``array`` unchanged otherwise.
    '''
    try:
        return array.dropChannelAxis()
    except:
        return array

# FIXME: This is a workaround for the disabled C++ function for the same purpose.
#        Enable the C++ version when boost 1.41 is available on all relevant platforms.
def _AxisTags_fromJSON(json_rep):
    '''
        Construct a new AxisTags object from the given JSON representation.
        This is mainly used to reconstruct arrays from HDF5 datasets with
        a suitable axistags attribute (see :func:`~vigra.impex.readHDF5`).
    '''
    tag_dict = eval(json_rep)
    tag_list = []
    for tags in tag_dict['axes']:
        tags['typeFlags'] = AxisType(tags['typeFlags'])
        tag_list.append(AxisInfo(**tags))
    return AxisTags(tag_list)

AxisTags.fromJSON = staticmethod(_AxisTags_fromJSON)
AxisTags.fromJSON.__doc__ = _AxisTags_fromJSON.__doc__
    
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
        # (otherwise, strides in the copy may not exactly match those in obj)
        strides = list(obj.strides)
        try:
            channelIndex = obj.channelIndex
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
    
##################################################################

class VigraArray(numpy.ndarray):
    '''
    This class extends numpy.ndarray with the concept of **axistags** 
    which encode the semantics of the array's axes. VigraArray overrides all
    numpy.ndarray methods in order to handle axistags in a sensible way.
    In particular, operations acting on two arrays simultaneously (e.g.
    addition) will first transpose the arguments such that their axis 
    ordering matches.

    Constructor:
    
    .. method:: VigraArray(obj, dtype=numpy.float32, order=None, init=True, value=None, axistags=None)

        :param obj: an array or shape object (see below)
        :param dtype: desired element type
        :param order: desired memory layout (see below)
        :param init: True: initialize the image with zeros; False: do not 
                     initialize the image
        :type init: boolean
        :param value: initialize the image with this value (overrides init)
        :type value: convertible to dtype
        :param axistags: the AxisTags object of the new array. The length of 
                         axistags must match the array's shape. It axistags=None, 
                         obj.axistags is used if it exists. Otherwise, a new 
                         axistags object is created by a call to 
                         :meth:`~vigra.VigraArray.defaultAxistags`.
 
        **obj** may be one of the following

        * If obj is a numpy.ndarray or a subtype, a copy of obj with the given 
          dtype, order and resulting class VigraArray is created. If obj.axistags
          exists, the new array will have these axistags as well, unless new
          axistags are explicitly passed to the constructor.
        * If obj is a sequence, it is interpreted as a shape. 
        * Otherwise, or if shape and axistags are incompatible, an exception
          is raised.
          
        **order** can be 'C' (C order), 'F' (Fortran order), 'V' (VIGRA
        order), 'A' (any), or None. This parameter controls the order of strides 
        and axistags (unless axistags are explicit passed into the constructor). 
        See the :ref:`order definitions <array-order-parameter>` for details. If 
        'order=None', the order is determined by :attr:`VigraArray.defaultOrder`.
    '''
    
    ###############################################################
    #                                                             #
    #       a number of helper functions related to axistags      #
    #                                                             #
    ###############################################################    
    
    # a number of helper functions related to axistags
    
    # IMPORTANT: do not remove or rename this function, it is called from C++
    @classproperty
    def defaultOrder(cls):
        '''
        Get the default axis ordering, currently 'V' (VIGRA order)
        '''
        return 'V'

    # IMPORTANT: do not remove or rename this function, it is called from C++
    @staticmethod
    def defaultAxistags(tagSpec, order=None, noChannels=False):
        '''
        Get default axistags for the given specification 'tagSpec'. TagSpec can be the 
        number of dimensions of the array (``array.ndim``, must be <= 5) or a string 
        containing a sequence of axis keys (only the default keys 'x', 'y', 'z', 't', 
        and 'c' are currently supported). The 'order' parameter determines the axis 
        ordering, see the :ref:`order definitions <array-order-parameter>` for details.
        If 'noChannels' is True, there will be no channel axis. Examples::
        
            >>> vigra.VigraArray.defaultAxistags(3)
            x y c
            >>> vigra.VigraArray.defaultAxistags(4)
            x y z c
            >>> vigra.VigraArray.defaultAxistags(5)
            x y z t c
            >>> vigra.VigraArray.defaultAxistags(3, order='C')
            y x c
            >>> vigra.VigraArray.defaultAxistags(2, noChannels=True)
            x y z
            >>> vigra.VigraArray.defaultAxistags(3, noChannels=True)
            x y z
            >>> vigra.VigraArray.defaultAxistags(4, noChannels=True)
            x y z t
            >>> vigra.VigraArray.defaultAxistags('xty')
            x t y
            >>> vigra.VigraArray.defaultAxistags('xty', order='V')
            x y t
        '''
        if type(tagSpec) == str:
            taglist = [eval('AxisInfo.' + k) for k in tagSpec]
        else:
            start = 1 if noChannels else 0
            end = start + tagSpec
            taglist = [AxisInfo.c, AxisInfo.x, AxisInfo.y, AxisInfo.z, AxisInfo.t][start:end]
            if order is None or order == 'A':
                order = VigraArray.defaultOrder
        tags = AxisTags(taglist)
        if order is not None:
            tags.transpose(tags.permutationToOrder(order))
        return tags

    # IMPORTANT: do not remove or rename this function, it is called from C++
    @staticmethod
    def _copyValuesImpl(target, source):
        try:
            target = target.squeeze()
            target = target.transposeToNumpyOrder()
        except:
            pass

        try:
            source = source.squeeze()
            source = source.transposeToNumpyOrder()
        except:
            pass
        
        try:
            compatible = source.axistags.compatible(target.axistags)
        except:
            compatible = True

        if not compatible:
            raise RuntimeError("VigraArray._copyValuesImpl(): incompatible axistags")

        target[...] = source
    
    # IMPORTANT: do not remove or rename this function, it is called from C++
    @staticmethod
    def _empty_axistags(ndim):
        '''Create an axistags object with non-informative entries. 
           That is, all axisinfo objects are '?'.
        '''
        return AxisTags(ndim)
    
    def _copy_axistags(self):
        '''Create a copy of 'self.axistags'. If the array doesn't have axistags, _empty_axistags() 
           will be returned.
        '''
        return copy.copy(getattr(self, 'axistags', self._empty_axistags(self.ndim)))
        
    def _transform_axistags(self, index):
        if hasattr(self, 'axistags'):
            return self.axistags.transform(index, self.ndim)
        else:
            return self._empty_axistags(self.ndim)

    def _transpose_axistags(self, permutation = None):
        '''Create a copy of self.axistags with transposed entries.
        '''
        if hasattr(self, 'axistags'):
            res = copy.copy(self.axistags)
            res.transpose(permutation)
            return res
        else:
            return self._empty_axistags(self.ndim)

    ###############################################################
    #                                                             #
    #                   standard array functions                  #
    #                                                             #
    ###############################################################    
    
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

    def __array_finalize__(self, obj):
        if hasattr(obj, 'axistags'):
            self.axistags = obj.axistags

    def __copy__(self, order='A'):
        result = numpy.ndarray.__copy__(self, order)
        result.axistags = result._copy_axistags()
        return result
        
    @_preserve_doc
    def __deepcopy__(self, memo):
        # numpy.ndarray.__deepcopy__ always creates C-order arrays =>
        #   transpose self accordingly, and transpose back after the copy
        result = numpy.ndarray.__deepcopy__(self.transposeToNumpyOrder(), memo)
        result = result.transpose(self.permutationFromNumpyOrder())
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        return result
    
    def __repr__(self):
        return "%s(shape=%s, axistags=%s, dtype=%s, data=\n%s)" % \
          (self.__class__.__name__, str(self.shape), repr(self.axistags), str(self.dtype), str(self))
          
    def __str__(self):
        try:
            self = self.transposeToVigraOrder().transpose()
        except:
            pass
        return str(self.view(numpy.ndarray))
    
    ###############################################################
    #                                                             #
    #                     array I/O and display                   #
    #                                                             #
    ###############################################################    
    
    def writeImage(self, filename, dtype = '', compression = ''):
        '''Write an image to a file. 
        Consult :func:`vigra.impex.writeImage` for detailed documentation'''
        import vigra.impex

        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim != 2:
            raise RuntimeError("VigraArray.writeImage(): array must have 2 non-channel axes.")

        vigra.impex.writeImage(self, writeImage, dtype, compression)
            
    def writeSlices(self, filename_base, filename_ext, dtype = '', compression = ''):
        '''Write a volume to a sequence of files.         
        Consult :func:`vigra.impex.writeVolume` for detailed documentation.
        '''
        import vigra.impex

        ndim = self.ndim
        if self.channelIndex < ndim:
            ndim -= 1
        if ndim != 3:
            raise RuntimeError("VigraArray.writeSlices(): array must have 3 non-channel axes.")

        vigra.impex.writeVolume(self, filename_base, filename_ext, dtype, compression)
            
    def writeHDF5(self, filenameOurGroup, pathInFile):
        '''Write the array to a HDF5 file. 
           This is just a shortcut for :func:`vigra.impex.writeHDF5`
        '''
        import vigra.impex

        vigra.impex.writeHDF5(self, filenameOurGroup, pathInFile)

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
        '''
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
           
        '''
        try:
            import qimage2ndarray
        except:
            from vigra import _fallbackModule
            _fallbackModule('qimage2ndarray',
            '''    It can be obtained at
            http://pypi.python.org/pypi/qimage2ndarray/.''')
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

    ###############################################################
    #                                                             #
    #           new functionality enabled by axistags             #
    #                                                             #
    ###############################################################    
    
    def copyValues(self, other):
        '''
        Copy the values of an array to another one. This is similar to::
        
            self[...] = other
            
        but will first transpose both arrays so that axistags are aligned. If
        there is no valid alignment, RuntimeError will be raised.
        '''
        self._copyValuesImpl(self, other)
    
    # IMPORTANT: do not remove or rename this property, it is called from C++
    @property
    def channelIndex(self):
        '''
        The index of the channel axis according to the axistags.        
        For example, when axistags are 'x y c', the channel index is 2. 
        If the axistags contain no channel axis, self.ndim is returned.
        '''
        return self.axistags.channelIndex
    
    # IMPORTANT: do not remove or rename this property, it is called from C++
    @property
    def innerNonchannelIndex(self):
        '''
        The index of the innermost non-channel axis according to the axistags.
        The innermost axis is determined by the AxisInfo sorting rules (see 
        the :ref:`order definitions <array-order-parameter>` for details).
        For example, when axistags are 'x y c', the innerNonchannelIndex is 0.
        '''
        return self.axistags.innerNonchannelIndex
    
    @property
    def channels(self):
        '''
        The number of channels in this array (shape of the 'c' axis).
        If the axistags contain no channel axis, the number of channels is implicitly 1.
        '''
        i = self.channelIndex
        if i < self.ndim:
            return self.shape[i]
        else:
            return 1
            
    @property
    def width(self):
        '''
        The width of the array (shape of the 'x' axis).
        If the axistags contain no 'x' axis, RuntimeError will be raised.
        '''
        i = self.axistags.index('x')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.width(): axistag 'x' does not exist.")
    
    @property
    def height(self):
        '''
        The height of the array (shape of the 'y' axis).
        If the axistags contain no 'y' axis, RuntimeError will be raised.
        '''
        i = self.axistags.index('y')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.height(): axistag 'y' does not exist.")
    
    @property
    def depth(self):
        '''
        The depth of the array (shape of the 'z' axis).
        If the axistags contain no 'z' axis, RuntimeError will be raised.
        '''
        i = self.axistags.index('z')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.depth(): axistag 'z' does not exist.")
    
    @property
    def duration(self):
        '''
        The number of time steps in the array (shape of the 't' axis).
        If the axistags contain no 't' axis, RuntimeError will be raised.
        '''
        i = self.axistags.index('t')
        if i < self.ndim:
            return self.shape[i]
        else:
            raise RuntimeError("VigraArray.duration(): axistag 't' does not exist.")
            
    @property
    def spatialDimensions(self):
        '''
        The number of spatial axes in the array.
        That is, the number of entries in the axistags where the flag 'AxisType.Space'
        is set.
        '''
        return self.axistags.axisTypeCount(AxisType.Space)
            
    def iterImpl(self, type):
        axes = [k for k in xrange(self.ndim) if self.axistags[k].isType(type)]
        if axes:
            axes.sort(key=lambda x: self.axistags[x], reverse=True)
            slices = [slice(None)]*self.ndim
            for point in numpy.ndindex(*(self.shape[k] for k in axes)):
                for j in xrange(len(point)):
                    slices[axes[j]] = point[j]
                yield self[tuple(slices)]
        else:
            yield self
            
    def channelIter(self):
        '''
        Create an iterator over the channels of the array.
        In each iteration, you get the array corresponding to a single channel.
        If the axistags contain no channel axis, there is only one iteration
        which yields the entire array. Example::
        
            >>> rgb = vigra.RGBImage((200, 100))
            >>> rgb.axistags
            x y c
            >>> red, green, blue = rgb.channelIter()
            >>> red.axistags
            x y
            >>> red.shape
            (200, 100)
        '''
        return self.iterImpl(AxisType.Channels)
            
    def spaceIter(self):
        '''
        Create an iterator over all the spatial coordinates in the array.
        In each iteration, you get the value corresponding to a single 
        coordinate location. If the axistags contain no spatial axes, 
        there is only one iteration which yields the entire array. Example::
        
            >>> s = vigra.ScalarImage((2,2))
            >>> s.ravel()[...] = range(4)
            >>> for p in s.spaceIter():
            ....    print p
            0.0
            1.0
            2.0
            3.0
        '''
        return self.iterImpl(AxisType.Space)
        
    def timeIter(self):
        '''
        Create an iterator over the time points of the array.
        In each iteration, you get the array corresponding to a single time point.
        If the axistags contain no time axis, there is only one iteration
        which yields the entire array. Example::
        
            >>> from vigra import *
            >>> axistags = AxisTags(AxisInfo.t, AxisInfo.x, AxisInfo.y)
            >>> timesteps, width, height = 2, 200, 100
            >>> image_sequence = Image((timesteps, width, height), axistags=axistags)
            >>> step1, step2 = image_sequence.timeIter()
        '''
        return self.iterImpl(AxisType.Time)
        
    def sliceIter(self, key='z'):
        '''
        Create an iterator over a single spatial axis of the array.
        In each iteration, you get the array corresponding to one coordinate
        along the axis given by 'key'. For example, to iterate along the z-axis
        to get all x-y-slices in turn, you write::
        
            >>> volume = vigra.Volume((width, height, depth))
            >>> for slice in volume.sliceIter('z'):
            ...     processSlice(slice)
        '''
        i = self.axistags.index(key)
        if i < self.ndim:
            if not self.axistags[i].isSpatial():
                raise RuntimeError("VigraArray.sliceIter(): %s is not a spatial axis." % key)
            for k in xrange(self.shape[i]):
                yield self.bindAxis(i, k)
        else:
            yield self
    
    def bindAxis(self, which, index=0):
        '''
        Bind the axis identified by 'which' to the given 'index'.
        This is similar to::
        
            array[:, index, ...]
            
        but you don not need to know the position of the axis when you use the 
        axis key (according to axistags). For example, to get the green channel
        of an RGBImage, you write::
        
            >>> rgb = vigra.RGBImage((200, 100))
            >>> green = rgb.bindAxis('c', 1)
            
        This gives the correct result irrespective of the axis ordering.
        '''
        if type(which) == str:
            which = self.axistags.index(which)
        return self[(slice(None),)*which + (index,) + (slice(None),)*(self.ndim-which-1)]
    
    def dropChannelAxis(self, ignoreMultiChannel=False):
        '''
        Drop the channel axis when it is a singleton.
        This function is for easy transformation of an array shaped 
        (width, height, 1) into (width, height). A RuntimeError
        is raised when there is more than one channel, unless ignoreMultiChannel=True,
        in which case 'self' is returned.
        '''
        ci = self.channelIndex
        if ci == self.ndim:
            return self
        
        if self.shape[ci] != 1:
            if ignoreMultiChannel:
                return self
            raise RuntimeError("dropChannelAxis(): only allowed when there is a single channel.")
        return self.bindAxis(ci, 0)
    
    def insertChannelAxis(self, order=None):
        '''
        Insert a singleton channel axis.
        This function is for easy transformation of an array shaped 
        (width, height) into (width, height, 1). The 'order' parameter
        determines the position of the new axis: when order is 'F', it
        will become the first axis, otherwise it will become the last
        one. A RuntimeError is raised when the array already contains a
        channel axis.
        '''
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
    
    def view5D(self, order='C'):
        '''
            Create a 5-dimensional view containing the standard tags 'x', 'y', 'z', 't', 'c'
            in the desired 'order' (which can be 'C', 'F', and 'V' with the usual meaning). 
            If 'self' has an axis key that is not among the five admissible keys, an 
            exception is raised. Axes missing in 'self' are added as singleton axes with the appropriate tags. 
        '''
        stdTags = ['x', 'y', 'z', 't', 'c']
        for tag in self.axistags:
            try:
                del stdTags[stdTags.index(tag.key)]
            except:
                raise RuntimeError("VigraArray.view5D(): array contains unsuitable axis key '%s'." % tag.key)
        index = [Ellipsis] + [newaxis(eval('AxisInfo.' + k)) for k in stdTags]
        return self[index].transposeToOrder(order)
    
    def permutationToOrder(self, order):
        '''Create the permutation that would transpose this array into
           an array view with the given order (where order can be 'A',
           'C', 'F', 'V' with the usual meaning). 
        '''
        return list(self.axistags.permutationToOrder(order))
    
    def permutationToNormalOrder(self, types=AxisType.AllAxes):
        '''Create the permutation that would transpose this array to 
           normal order (that is, from the current axis order into 
           ascending order, e.g. 'x y c' into 'c x y'). 
           If 'types' is not 'AxisType.AllAxes', only the axes with the
           desired types are considered.
        '''
        return list(self.axistags.permutationToNormalOrder(types))
    
    def permutationFromNormalOrder(self):
        '''Create the permutation that would transpose an array that is 
           in normal (ascending) order into the axis order of this array.
           (e.g. 'c x y' into 'x y c'). 
        '''
        return list(self.axistags.permutationFromNormalOrder())
    
    def permutationToNumpyOrder(self):
        '''Create the permutation that would transpose this array to 
           numpy order (that is, from the current axis order into 
           descending order, e.g. 'x y c' into 'y x c'). 
        '''
        return list(self.axistags.permutationToNumpyOrder())
    
    def permutationFromNumpyOrder(self):
        '''Create the permutation that would transpose an array that is 
           in numpy (descending) order into the axis order of this array.
           (e.g.  'y x c' into 'x y c'). 
        '''
        return list(self.axistags.permutationFromNumpyOrder())
    
    def permutationToVigraOrder(self):
        '''Create the permutation that would transpose this array to 
           VIGRA order (that is, from the current axis order into 
           ascending spatial order, but with the channel axis at the 
           last position, e.g. 'c x y' into 'x y c'). 
        '''
        return list(self.axistags.permutationToVigraOrder())
    
    def permutationFromVigraOrder(self):
        '''Create the permutation that would transpose an array that is 
           in VIGRA order (ascending spatial order, but with the channel 
           axis at the last position) into the axis order of this array.
           (e.g.  'x y c' into 'c x y'). 
        '''
        return list(self.axistags.permutationFromVigraOrder())
    
    def transposeToOrder(self, order):
        '''
        Get a transposed view onto this array according to the given 'order'.
        Possible orders are:
        
        'A':
            return the array unchanged
        'C':
            transpose to descending axis order (e.g. 'z y x c')
        'F':
            transpose to ascending axis order (e.g. 'c x y z')
        'V':
            transpose to VIGRA order, i.e. ascending spatial axes, but
            the channel axis is last (e.g. 'x y z c')
        '''
        if order == 'A':
            return self
        permutation = self.permutationToOrder(order)
        return self.transpose(permutation)
    
    def transposeToDefaultOrder(self):
        '''Equivalent to self.transposeToOrder(VigraArray.defaultOrder).
        '''
        return self.transposeToOrder(VigraArray.defaultOrder)

    def transposeToNormalOrder(self):
        '''Equivalent to self.transposeToOrder('F').
        '''
        return self.transposeToOrder('F')

    def transposeToVigraOrder(self):
        '''Equivalent to self.transposeToOrder('V').
        '''
        return self.transposeToOrder('V')

    def transposeToNumpyOrder(self):
        '''Equivalent to self.transposeToOrder('C').
        '''
        return self.transposeToOrder('C')

    @property
    def T(self):
        '''
        Equivalent to self.transpose()
        '''
        return self.transpose()

    def __getitem__(self, index):
        '''
        ``array.__getitem__(index)`` implements the indexing operator ``array[index]``.
        In addition to the usual numpy.ndarray indexing functionality, this function
        also updates the axistags of the result array. There are three cases:

            * getitem creates a scalar value => no axistags are required
            * getitem creates an array view => axistags are transferred from the
              corresponding axes of the base array
            * getitem creates a copy of an array (fancy indexing) => all axistags are '?'
         
        If the index contains 'numpy.newaxis', a new singleton axis is inserted at the 
        appropriate position, whose axisinfo is set to '?' (unknown). If the index contains
        'vigra.newaxis(axisinfo)', the singleton axis will get the given axisinfo.
        '''
        try:
            res = numpy.ndarray.__getitem__(self, index)
        except:
            res = numpy.ndarray.__getitem__(self, 
                     map(lambda x: None if isinstance(x, AxisInfo) else x, index))
        if res is not self and hasattr(res, 'axistags'):
            if res.base is self:
                res.axistags = res._transform_axistags(index)
            else:
                res.axistags = res._empty_axistags(res.ndim)
        return res
    
    ###############################################################
    #                                                             #
    #      re-implement ndarray methods to handle axistags        #
    #                                                             #
    ###############################################################    
    
    @_finalize_reduce_result
    @_preserve_doc
    def all(self, axis=None, out=None):
        return numpy.ndarray.all(self, axis, out)

    @_finalize_reduce_result
    @_preserve_doc
    def any(self, axis=None, out=None):
        return numpy.ndarray.any(self, axis, out)

    @_finalize_reduce_result
    @_preserve_doc
    def argmax(self, axis=None, out=None):
        return numpy.ndarray.argmax(self, axis, out)
        
    @_finalize_reduce_result
    @_preserve_doc
    def argmin(self, axis=None, out=None):
        return numpy.ndarray.argmin(self, axis, out)
    
    @_preserve_doc
    def copy(self, order='A'):
        return self.__class__(self, dtype=self.dtype, order=order)
    
    @_finalize_reduce_result
    @_preserve_doc
    def cumsum(self, axis=None, dtype=None, out=None):
        return numpy.ndarray.cumsum(self, axis, dtype, out)

    @_finalize_reduce_result
    @_preserve_doc
    def cumprod(self, axis=None, dtype=None, out=None):
        return numpy.ndarray.cumprod(self, axis, dtype, out)

    @property
    def flat(self):
        return self.transposeToNumpyOrder().view(numpy.ndarray).flat
    
    @_preserve_doc
    def flatten(self, order='C'):
         return self.transposeToNumpyOrder().view(numpy.ndarray).flatten(order)        

    @_finalize_reduce_result
    @_preserve_doc
    def max(self, axis=None, out=None):
        return numpy.ndarray.max(self, axis, out)

    @_finalize_reduce_result
    @_preserve_doc
    def mean(self, axis=None, out=None):
        return numpy.ndarray.mean(self, axis, out)
    
    @_finalize_reduce_result
    @_preserve_doc
    def min(self, axis=None, out=None):
        return numpy.ndarray.min(self, axis, out)
    
    @_preserve_doc
    def nonzero(self):
        res = numpy.ndarray.nonzero(self)
        for k in xrange(len(res)):
            res[k].axistags = copy.copy(self.axistags[k])
        return res

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
    
    @_finalize_reduce_result
    @_preserve_doc
    def prod(self, axis=None, dtype=None, out=None):
        return numpy.ndarray.prod(self, axis, dtype, out)

    @_finalize_reduce_result
    @_preserve_doc
    def ptp(self, axis=None, out=None):
        return numpy.ndarray.ptp(self, axis, out)

    @_preserve_doc
    def ravel(self, order='C'):
        return self.transposeToNumpyOrder().view(numpy.ndarray).ravel(order)        

    @_preserve_doc
    def repeat(self, repeats, axis=None):
        res = numpy.ndarray.repeat(self, repeats, axis)
        if axis is None:
            res.axistags = res._empty_axistags(res.ndim)
        return res        

    @_preserve_doc
    def reshape(self, shape, order='C'):
        res = numpy.ndarray.reshape(self, shape, order)
        res.axistags = res._empty_axistags(res.ndim)
        return res        

    @_preserve_doc
    def resize(self, new_shape, refcheck=True, order=False):
        res = numpy.ndarray.reshape(self, new_shape, refcheck, order)
        res.axistags = res._empty_axistags(res.ndim)
        return res        
            
    @_preserve_doc
    def squeeze(self):
        res = numpy.ndarray.squeeze(self)
        if self.ndim != res.ndim:
            res.axistags = res._copy_axistags()
            for k in xrange(self.ndim-1, -1, -1):
                if self.shape[k] == 1:
                    del res.axistags[k]
        return res        

    @_finalize_reduce_result
    @_preserve_doc
    def std(self, axis=None, dtype=None, out=None, ddof=0):
        return numpy.ndarray.std(self, axis, dtype, out, ddof)

    @_finalize_reduce_result
    @_preserve_doc
    def sum(self, axis=None, dtype=None, out=None):
        return numpy.ndarray.sum(self, axis, dtype, out)
            
    @_preserve_doc
    def swapaxes(self, i, j):
        res = numpy.ndarray.swapaxes(self, i, j)
        res.axistags = res._copy_axistags()
        try:
            res.axistags.swapaxes(i, j)
        except:
            res.axistags[i], res.axistags[j] = res.axistags[j], res.axistags[i]
        return res        
 
    @_preserve_doc
    def take(self, indices, axis=None, out=None, mode='raise'):
        res = numpy.ndarray.take(self, indices, axis, out, mode)
        if axis is None:
            res.axistags = res._empty_axistags(res.ndim)
        return res        
           
    @_preserve_doc
    def transpose(self, *axes):
        res = numpy.ndarray.transpose(self, *axes)
        res.axistags = res._transpose_axistags(*axes)
        return res

    @_finalize_reduce_result
    @_preserve_doc
    def var(self, axis=None, dtype=None, out=None, ddof=0):
        return numpy.ndarray.var(self, axis, dtype, out, ddof)

    ###############################################################
    #                                                             #
    #        reimplement the numerical operators to make          #
    #             sure that array order is preserved              #
    #                                                             #
    ###############################################################    
    
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
    '''
    Factory function for a :class:`~vigra.VigraArray` representing an image (i.e. an array with 
    two spatial axes 'x' and 'y' and optionally a channel axis 'c'). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are not image-like.
    '''
    obj, axistags = _adjustInput(obj, order, 2, 0, axistags, "vigra.Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)
        
def ScalarImage(obj, dtype=numpy.float32, order=None, 
                init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a single-band image (i.e. an 
    array with two spatial axes 'x' and 'y' and no channel axis). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a single-band image.
    '''
    obj, axistags = _adjustInput(obj, order, 2, None, axistags, "vigra.ScalarImage()")
    return VigraArray(obj, dtype, None, init, value, axistags)
        
def Vector2Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 2-band image (i.e. an 
    array with two spatial axes 'x' and 'y' and channel axis 'c' with 2 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 2-band image.
    '''
    obj, axistags = _adjustInput(obj, order, 2, 2, axistags, "vigra.Vector2Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector3Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 3-band image (i.e. an 
    array with two spatial axes 'x' and 'y' and channel axis 'c' with 3 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 3-band image.
    '''
    obj, axistags = _adjustInput(obj, order, 2, 3, axistags, "vigra.Vector3Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector4Image(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 4-band image (i.e. an 
    array with two spatial axes 'x' and 'y' and channel axis 'c' with 4 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 4-band image.
    '''
    obj, axistags = _adjustInput(obj, order, 2, 4, axistags, "vigra.Vector4Image()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def RGBImage(obj, dtype=numpy.float32, order=None, 
             init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a RGB image (i.e. an 
    array with two spatial axes 'x' and 'y' and channel axis 'c' with 3 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for an RGB image.
    '''
    obj, axistags = _adjustInput(obj, order, 2, 3, axistags, "vigra.RGBImage()")
    res = VigraArray(obj, dtype, None, init, value, axistags)
    res.axistags['c'].description = 'RGB'
    return res

#################################################################

def Volume(obj, dtype=numpy.float32, order=None, 
           init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a volume (i.e. an array with 
    three spatial axes 'x', 'y' and 'z' and optionally a channel axis 'c'). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are not volume-like.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 0, axistags, "vigra.Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def ScalarVolume(obj, dtype=numpy.float32, order=None, 
                 init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a single-band volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and no channel axis). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a single-band volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, None, axistags, "vigra.ScalarVolume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector2Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 2-band volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and channel axis 'c' with 2 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 2-band volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 2, axistags, "vigra.Vector2Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector3Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 3-band volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and channel axis 'c' with 3 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 3-band volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 3, axistags, "vigra.Vector3Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector4Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 4-band volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and channel axis 'c' with 4 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 4-band volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 4, axistags, "vigra.Vector4Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def Vector6Volume(obj, dtype=numpy.float32, order=None, 
                  init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing a 6-band volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and channel axis 'c' with 6 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for a 6-band volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 6, axistags, "vigra.Vector6Volume()")
    return VigraArray(obj, dtype, None, init, value, axistags)

def RGBVolume(obj, dtype=numpy.float32, order=None, 
              init=True, value=None, axistags=None):
    '''
    Factory function for a :class:`~vigra.VigraArray` representing an RGB volume (i.e. an 
    array with three spatial axes 'x', 'y' and 'z' and channel axis 'c' with 3 channels). 
    Paramters are interpreted as in the VigraArray constructor, but an exception
    will be raised if the shape or axistags are unsuitable for an RGB volume.
    '''
    obj, axistags = _adjustInput(obj, order, 3, 3, axistags, "vigra.RGBVolume()")
    res = VigraArray(obj, dtype, None, init, value, axistags)
    res.axistags['c'].description = 'RGB'
    return res

#################################################################

class ImagePyramid(list):
    def __init__(self, image, copyImageToLevel = 0, lowestLevel = 0, highestLevel = 0):
        ''' Create a new pyramid.
            The new pyramid levels range from 'lowestLevel' to 'highestLevel' (inclusive),
            and the given 'image' is copied to 'copyImageToLevel'. The images at other
            levels are filled with zeros and sized so that the shape is reduced by half
            when going up (to higher levels), and doubled when going down.
            
            This class can handle multi-channel images, but only when image.channelIndex
            exists and returns image.ndim-1 (i.e. the image must have axistags, and the 
            channel axis must correspond to the last index, as in C- or V-order).
        '''
        if lowestLevel > copyImageToLevel or highestLevel < copyImageToLevel:
            raise ValueError('ImagePyramid(): copyImageToLevel must be between lowestLevel and highestLevel (inclusive)')
        
        list.__init__(self, [image.__class__(image, dtype=image.dtype)])
        self._lowestLevel = copyImageToLevel
        self._highestLevel = copyImageToLevel
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
    
        smooth1 = filters.explictKernel(-1, 1, numpy.array([0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue]))
        smooth2 = filters.explictKernel(-1, 0, numpy.array([0.5, 0.5]));

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

        smooth1 = filters.explictKernel(-1, 1, numpy.array([0.5 - centerValue, 2.0*centerValue, 0.5 - centerValue]))
        smooth2 = filters.explictKernel(-1, 0, numpy.array([0.5, 0.5]));
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
                channelIndex = getattr(image, 'channelIndex', image.ndim)
                if channelIndex < image.ndim:
                    newShape[channelIndex] = image.shape[channelIndex]
                self.append(image.__class__(newShape, dtype=image.dtype))
            self._highestLevel = level
        elif level < self.lowestLevel:
            image = list.__getitem__(self, 0)
            for i in range(self.lowestLevel, level, -1):
                newShape = [2*k-1 for k in image.shape]
                channelIndex = getattr(image, 'channelIndex', image.ndim)
                if channelIndex < image.ndim:
                    newShape[channelIndex] = image.shape[channelIndex]
                self.insert(0, image.__class__(newShape, dtype=image.dtype))
            self._lowestLevel = level
             