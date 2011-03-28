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

import copy
import numpy
    
def preserve_doc(f):
    f.__doc__ = eval('numpy.ndarray.%s.__doc__' % f.__name__)
    return f

class TaggedArray(numpy.ndarray):
    '''
TaggedArray extends numpy.ndarray with an attribute 'axistags'. Any 
axistags object must support the standard sequence interface, and its 
length must match the number of dimensions of the array. Each item in 
the axistags sequence is supposed to provide a description of the 
corresponding array axis. All array functions that change the number or 
ordering of an array's axes (such as transpose() and __getitem__()) are 
overloaded so that they apply the same transformation to the axistags 
object. 

Example:
  >>> axistags = ['x', 'y']
  >>> a = TaggedArray((2,3), axistags=axistags)
  >>> a.axistags
  ['x', 'y']
  >>> a[:,0].axistags
  ['x']
  >>> a[1,...].axistags
  ['y']
  >>> a.transpose().axistags
  ['y', 'x']
  
Except for the new 'axistags' keyword, the 'TaggedArray' constructor is identical to the constructor 
of 'numpy.ndarray'. 
    '''
    def __new__(subtype, shape, dtype=float, buffer=None, offset=0, strides=None, order=None, axistags=None):
        res = numpy.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides, order)
        if subtype is not numpy.ndarray:
            if axistags is None:
                res.axistags = res.empty_axistags(res.ndim)
            else:
                if len(axistags) != res.ndim:
                    raise RuntimeError('TaggedArray(): len(axistags) must match ndim')
                res.axistags = copy.copy(axistags)
        return res
        
    @staticmethod
    def empty_axistags(ndim):
        '''Create an axistags object of the given length with non-informative entries.
        '''
        return [None]*ndim
    
    def copy_axistags(self):
        '''Create a copy of 'self.axistags'. If the array doesn't have axistags, empty_axistags() 
           will be returned.
        '''
        return copy.copy(getattr(self, 'axistags', self.empty_axistags(self.ndim)))
        
    def transpose_axistags(self, axes=None):
        '''Create a copy of 'self.axistags' according to the given axes permutation 
           (internally called in transpose()).
        '''
        axistags = self.empty_axistags(self.ndim)
        if hasattr(self, 'axistags'):
            if axes is None:
                axes = range(self.ndim-1, -1, -1)
            for k in xrange(self.ndim):
                axistags[k] = self.axistags[int(axes[k])]
        return axistags
        
    def transform_axistags(self, index):
        '''Create new axistags from 'self.axistags' according to the given index or 
           slice object (internally called in __getitem__()).
        '''
        # we assume that self.ndim is already set to its new value, whereas
        # self.axistags has just been copied by __array_finalize__
        
        new_axistags = self.empty_axistags(self.ndim)
        
        if hasattr(self, 'axistags'):
            old_axistags = self.axistags
            old_ndim = len(old_axistags)
            new_ndim = len(new_axistags)
            try:
                # make sure that 'index' is a tuple
                len_index = len(index)
            except:
                index = (index,)
                len_index = 1
            len_index -= index.count(numpy.newaxis)
            if len_index < old_ndim and index.count(Ellipsis) == 0:
                index += (Ellipsis,)
                len_index += 1
            
            # how many missing axes are represented by an Ellipsis ?
            len_ellipsis = old_ndim - len_index
            
            knew, kold, kindex = 0, 0, 0
            while knew < new_ndim:
                try:
                    # if index[kindex] is int, the dimension is bound => drop this axis
                    int(index[kindex]) 
                    kold += 1
                    kindex += 1
                except:
                    if index[kindex] is not numpy.newaxis:
                        # copy the tag
                        new_axistags[knew] = old_axistags[kold]
                        kold += 1
                    knew += 1
                    # the first ellipsis represents all missing axes
                    if len_ellipsis > 0 and index[kindex] is Ellipsis:
                        len_ellipsis -= 1
                    else:
                        kindex += 1
        return new_axistags
    
    __array_priority__ = 10.0
    
    def __array_finalize__(self, obj):
        if hasattr(obj, 'axistags'):
            self.axistags = obj.axistags

    @preserve_doc
    def __copy__(self, order = 'C'):
        result = numpy.ndarray.__copy__(self, order)
        result.axistags = result.copy_axistags()
        return result
    
    @preserve_doc
    def __deepcopy__(self, memo):
        result = numpy.ndarray.__deepcopy__(self, memo)
        memo[id(self)] = result
        result.__dict__ = copy.deepcopy(self.__dict__, memo)
        return result
    
    def __repr__(self):
        return "%s(shape=%s, axistags=%s, dtype=%s, data=\n%s)" % \
          (self.__class__.__name__, str(self.shape), str(self.axistags), str(self.dtype), str(self))
          
    @preserve_doc
    def all(self, axis=None, out=None):
        res = numpy.ndarray.all(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res

    @preserve_doc
    def any(self, axis=None, out=None):
        res = numpy.ndarray.any(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res

    @preserve_doc
    def argmax(self, axis=None, out=None):
        res = numpy.ndarray.argmax(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res
        
    @preserve_doc
    def argmin(self, axis=None, out=None):
        res = numpy.ndarray.argmin(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res
    
    @preserve_doc
    def cumsum(self, axis=None, dtype=None, out=None):
        res = numpy.ndarray.cumsum(self, axis, dtype, out)
        if res.ndim != self.ndim:
            res.axistags = res.empty_axistags(res.ndim)
        return res        

    @preserve_doc
    def cumprod(self, axis=None, dtype=None, out=None):
        res = numpy.ndarray.cumprod(self, axis, dtype, out)
        if res.ndim != self.ndim:
            res.axistags = res.empty_axistags(res.ndim)
        return res        

    # FIXME: we should also provide a possibility to determine flattening order by axistags
    #        (the same applies to flat and ravel)
    @preserve_doc
    def flatten(self, order='C'):
        res = numpy.ndarray.flatten(self, order)
        res.axistags = res.empty_axistags(res.ndim)
        return res        

    @preserve_doc
    def max(self, axis=None, out=None):
        res = numpy.ndarray.max(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res

    @preserve_doc
    def mean(self, axis=None, out=None):
        res = numpy.ndarray.mean(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res
    
    @preserve_doc
    def min(self, axis=None, out=None):
        res = numpy.ndarray.min(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res
    
    @preserve_doc
    def nonzero(self):
        res = numpy.ndarray.nonzero(self)
        for k in xrange(len(res)):
            res[k].axistags = copy.copy(self.axistags[k])
        return res

    @preserve_doc
    def prod(self, axis=None, dtype=None, out=None):
        res = numpy.ndarray.prod(self, axis, dtype, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res

    @preserve_doc
    def ptp(self, axis=None, out=None):
        res = numpy.ndarray.ptp(self, axis, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res

    @preserve_doc
    def ravel(self, order='C'):
        res = numpy.ndarray.ravel(self, order)
        res.axistags = res.empty_axistags(res.ndim)
        return res        

    @preserve_doc
    def repeat(self, repeats, axis=None):
        res = numpy.ndarray.repeat(self, repeats, axis)
        if axis is None:
            res.axistags = res.empty_axistags(res.ndim)
        return res        

    @preserve_doc
    def reshape(self, shape, order='C'):
        res = numpy.ndarray.reshape(self, shape, order)
        res.axistags = res.empty_axistags(res.ndim)
        return res        

    @preserve_doc
    def resize(self, new_shape, refcheck=True, order=False):
        res = numpy.ndarray.reshape(self, new_shape, refcheck, order)
        res.axistags = res.empty_axistags(res.ndim)
        return res        
            
    @preserve_doc
    def squeeze(self):
        res = numpy.ndarray.squeeze(self)
        if self.ndim != res.ndim:
            res.axistags = res.copy_axistags()
            for k in xrange(self.ndim-1, -1, -1):
                if self.shape[k] == 1:
                    del res.axistags[k]
        return res        

    @preserve_doc
    def std(self, axis=None, dtype=None, out=None, ddof=0):
        res = numpy.ndarray.std(self, axis, dtype, out, ddof)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        if len(res.shape) == 0:
            res = res.item()
        return res

    @preserve_doc
    def sum(self, axis=None, dtype=None, out=None):
        res = numpy.ndarray.sum(self, axis, dtype, out)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        return res
            
    @preserve_doc
    def swapaxes(self, i, j):
        res = numpy.ndarray.swapaxes(self, i, j)
        res.axistags = res.copy_axistags()
        try:
            res.axistags.swapaxes(i, j)
        except:
            res.axistags[i], res.axistags[j] = res.axistags[j], res.axistags[i]
        return res        
 
    @preserve_doc
    def take(self, indices, axis=None, out=None, mode='raise'):
        res = numpy.ndarray.take(self, indices, axis, out, mode)
        if axis is None:
            res.axistags = res.empty_axistags(res.ndim)
        return res        
           
    @preserve_doc
    def transpose(self, *axes):
        res = numpy.ndarray.transpose(self, *axes)
        res.axistags = res.transpose_axistags(*axes)
        return res

    @preserve_doc
    def var(self, axis=None, dtype=None, out=None, ddof=0):
        res = numpy.ndarray.var(self, axis, dtype, out, ddof)
        if axis is not None:
            res.axistags = res.copy_axistags()
            del res.axistags[axis]
        if len(res.shape) == 0:
            res = res.item()
        return res

    @property
    def T(self):
        return self.transpose()

    def __getitem__(self, index):
        '''x.__getitem__(y) <==> x[y]
         
           In addition to the usual indexing functionality, this function
           also updates the axistags of the result array. There are three cases:
             * getitem creates a scalar value => no axistags are required
             * getitem creates an arrayview => axistags are transferred from the
                                             corresponding axes of the base array,
                                             axes resulting from 'newaxis' get tag 'None'
             * getitem creates a copy of an array (fancy indexing) => all axistags are 'None'
        '''
        res = numpy.ndarray.__getitem__(self, index)
        if res is not self and hasattr(res, 'axistags'):
            if res.base is self:
                res.axistags = res.transform_axistags(index)
            else:
                res.axistags = res.empty_axistags(res.ndim)
        return res

##################################################################

# Python implementation of a more powerful AxisTags object (instead of the dault list)
# We don't need it, because vigranumpycore defines an even better C++ implementation.
class _PyAxisTags(object):
    def __init__(self, *args):
        if len(args) == 0:
            self.tags = []
        elif hasattr(args[0], 'tags'):
            self.tags = copy.copy(args[0].tags)
        elif hasattr(args[0], '__len__'):
            self.tags = list(args[0])
        else:
            try:
                self.tags = [AxisInfo() for k in xrange(args[0])]
            except:
                self.tags = list(args)
    
    def __copy__(self):
        return AxisTags(self)
    
    def __deepcopy__(self, memo):
        result = AxisTags()
        result.tags = copy.deepcopy(self.tags, memo)
        memo[id(self)] = result
        return result
    
    def __repr__(self):
        return ' '.join(map(lambda x: x.key, self.tags))
        
    def __eq__(self, other):
        return self.tags == other.tags
    
    def __ne__(self, other):
        return self.tags != other.tags
    
    def __len__(self):
        return len(self.tags)
    
    def __getitem__(self, index):
        if type(index) is str:
            index = self.index(index)
        return self.tags[index]
    
    def __setitem__(self, index, value):
        if type(index) is str:
            index = self.index(index)
        existing_index = self.index(value.key)
        if existing_index < len(self) and existing_index != index:
            raise RuntimeError("AxisTags.__setitem__(): axis key already exists.")
        self.tags[index] = value
    
    def __delitem__(self, index):
        if type(index) is str:
            index = self.index(index)
        del self.tags[index]
    
    def insert(self, index, value):
        if type(index) is str:
            index = self.index(index)
        if self.index(value.key) < len(self):
            raise RuntimeError("AxisTags.insert(): axis key already exists.")
        self.tags.insert(index, value)
    
    def append(self, value):
        if self.index(value.key) < len(self):
            raise RuntimeError("AxisTags.append(): axis key already exists.")
        self.tags.append(value)
    
    def index(self, key):
        for k in xrange(len(self.tags)):
            if self.tags[k].key == key:
                return k
        return len(self.tags)
   
    @property
    def channelIndex(self):
        return self.index('c')
   
    @property
    def majorNonchannelIndex(self):
        # FIXME: this must be generalized to the case when 'x' is not present.
        return self.index('x')
    
    def axisTypeCount(self, axistype):
        count = 0
        for k in self.tags:
            if k.isType(axistype):
                count += 1
        return count
    
    def permutationToNormalOrder(self):
        return canonicalAxisPermutation(self.tags)
    
    def permutationFromNormalOrder(self):
        return [int(k) for k in numpy.array(map(lambda x: ord(x.key[-1]), self.tags)).argsort().argsort()]
    
    def setChannelDescription(self, description):
        index = self.index('c')
        if index < len(self):
            self.tags[index].description = description
    
    def dropChannelAxis(self):
        index = self.index('c')
        if index < len(self):
            del self.tags[index]
    
    def insertChannelAxis(self):
        index = self.index('c')
        if index < len(self):
            raise RuntimeError("AxisTags.insertChannelDimension(): already have a channel dimension.")
        if defaultOrder == 'F':
            self.tags.insert(0, AxisInfo.c)
        else:
            self.tags.append(AxisInfo.c)
    
    def swapaxes(self, i1, i2):
        self.tags[i1], self.tags[i2] = self.tags[i2], self.tags[i1]
    
    def transpose(self, permutation = None):
        l = len(self)
        if permutation is None:
                permutation = range(l-1, -1, -1)
        result = AxisTags([None]*l)
        for k in xrange(l):
            result.tags[k] = self.tags[permutation[k]]
        return result
    
    def transform(self, index, new_ndim):
        new_axistags = [AxisInfo() for k in xrange(new_ndim)]
        old_axistags = self.tags
        old_ndim = len(old_axistags)

        try:
            # make sure that 'index' is a tuple
            len_index = len(index)
        except:
            index = (index,)
            len_index = 1
        len_index -= index.count(numpy.newaxis)
        if len_index < old_ndim and index.count(Ellipsis) == 0:
            index += (Ellipsis,)
            len_index += 1
        
        # how many missing axes are represented by an Ellipsis ?
        len_ellipsis = old_ndim - len_index
        
        knew, kold, kindex = 0, 0, 0
        while knew < new_ndim:
            try:
                # if index[kindex] is int, the dimension is bound => drop this axis
                int(index[kindex]) 
                kold += 1
                kindex += 1
            except:
                if index[kindex] is not numpy.newaxis:
                    # copy the tag
                    new_axistags[knew] = copy.copy(old_axistags[kold])
                    
                    # adjust the resolution for a possible step in the slice
                    try:
                        new_axistags[knew].resolution *= index[kindex].step
                    except:
                        pass
                        
                    kold += 1
                knew += 1
                # the first ellipsis represents all missing axes
                if len_ellipsis > 0 and index[kindex] is Ellipsis:
                    len_ellipsis -= 1
                else:
                    kindex += 1
        
        return _PyAxisTags(new_axistags)
