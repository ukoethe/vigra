import numpy as np

class InfoArray(np.ndarray):

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
          strides=None, order=None, axistags=None):
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides, order)
        if axistags is None:
            defaulttags = [['x'], ['y','x'], ['z','y','x'], ['z','y','x', 'c']]
            obj.axistags = defaulttags[len(shape)-1]
        else:
            obj.axistags = list(axistags)
        return obj

    def __array_finalize__(self, obj):
        if hasattr(obj, 'axistags'):
            self.axistags = list(obj.axistags)

    def all(self, axis=None, out=None):
        res = np.ndarray.all(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res

    def any(self, axis=None, out=None):
        res = np.ndarray.any(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res

    def argmax(self, axis=None, out=None):
        res = np.ndarray.argmax(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res
        
    def argmin(self, axis=None, out=None):
        res = np.ndarray.argmin(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res
    
    def cumsum(self, axis=None, dtype=None, out=None):
        res = np.ndarray.cumsum(self, axis, dtype, out)
        if res.ndim != self.ndim:
            res.axistags = [None]*res.ndim
        return res        

    def cumprod(self, axis=None, dtype=None, out=None):
        res = np.ndarray.cumprod(self, axis, dtype, out)
        if res.ndim != self.ndim:
            res.axistags = [None]*res.ndim
        return res        

    def flatten(self, order='C'):
        res = np.ndarray.flatten(self, order)
        res.axistags = [None]
        return res        

    def max(self, axis=None, out=None):
        res = np.ndarray.max(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res

    def mean(self, axis=None, out=None):
        res = np.ndarray.mean(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res
    
    def min(self, axis=None, out=None):
        res = np.ndarray.min(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res
    
    def nonzero(self):
        res = np.ndarray.nonzero(self)
        for k in xrange(len(res)):
            res[k].axistags = [self.axistags[k]]
        return res

    def prod(self, axis=None, dtype=None, out=None):
        res = np.ndarray.prod(self, axis, dtype, out)
        if axis is not None:
            del res.axistags[axis]
        return res

    def ptp(self, axis=None, out=None):
        res = np.ndarray.ptp(self, axis, out)
        if axis is not None:
            del res.axistags[axis]
        return res

    def ravel(self, order='C'):
        res = np.ndarray.ravel(self, order)
        res.axistags = [None]
        return res        

    def repeat(self, repeats, axis=None):
        res = np.ndarray.repeat(self, repeats, axis)
        if axis is None:
            res.axistags = [None]*res.ndim
        return res        

    def reshape(self, shape, order='C'):
        res = np.ndarray.reshape(self, shape, order)
        res.axistags = [None]*res.ndim
        return res        

    def resize(self, new_shape, refcheck=True, order=False):
        res = np.ndarray.reshape(self, new_shape, refcheck, order)
        res.axistags = [None]*res.ndim
        return res        
            
    def squeeze(self):
        res = np.ndarray.squeeze(self)
        for k in xrange(self.ndim-1, -1, -1):
            if self.shape[k] == 1:
                del res.axistags[k]
        return res        

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        res = np.ndarray.std(self, axis, dtype, out, ddof)
        if axis is not None:
            del res.axistags[axis]
        if len(res.shape) == 0:
            res = res.item()
        return res

    def sum(self, axis=None, dtype=None, out=None):
        res = np.ndarray.sum(self, axis, dtype, out)
        if axis is not None:
            del res.axistags[axis]
        return res
            
    def swapaxes(self, i, j):
        res = np.ndarray.swapaxes(self, i, j)
        res.axistags[i] = self.axistags[j]
        res.axistags[j] = self.axistags[i]
        return res        
 
    def take(self, indices, axis=None, out=None, mode='raise'):
        res = np.ndarray.take(self, indices, axis, out, mode)
        if axis is None:
            res.axistags = [None]*res.ndim
        return res        
           
    def transpose(self, *axes):
        res = np.ndarray.transpose(self, *axes)
        if len(axes) == 1:
            axes = axes[0]
        if not axes:
            res.axistags.reverse()
        else:
            res.axistags = [self.axistags[k] for k in axes]
        return res

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        res = np.ndarray.var(self, axis, dtype, out, ddof)
        if axis is not None:
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
             * getitem creates a value => no axistags are required
             * getitem creates an arrayview => axistags are transferred from the
                                             corresponding axes of the base array,
                                             axes resulting from 'newaxis' get tag 'None'
             * getitem creates a copy of an array (fancy indexing) => all axistags are 'None'
        '''
        res = np.ndarray.__getitem__(self, index)
        if hasattr(res, 'axistags'):
            # indexing created an array => set default (uninformative) axistags
            lnew = res.ndim
            if res is not self:
                res.axistags = [None]*lnew
            if res.base is self:
                # indexing created an array view => transfer the axistags
                try:
                    # make sure that 'index' is a tuple
                    lindex = len(index)
                except:
                    index = (index,)
                    lindex = 1
                lindex -= index.count(np.newaxis)
                if lindex < self.ndim and index.count(Ellipsis) == 0:
                    index += (Ellipsis,)
                    lindex += 1
                
                # how many missing axes are represented by an Ellipsis ?
                lellipsis = self.ndim - lindex
                
                knew, kold, kindex = 0, 0, 0
                while knew < lnew:
                    try:
                        # if index[kindex] is int, the dimension is bound => drop this axis
                        int(index[kindex]) 
                        kold += 1
                        kindex += 1
                    except:
                        if index[kindex] is np.newaxis:
                            # tag of newaxis is unknown
                            res.axistags[knew] = None
                        else:
                            # tag of existing axis is copied
                            res.axistags[knew] = self.axistags[kold]
                            kold += 1
                        knew += 1
                        # the first ellipsis represents all missing axes
                        if lellipsis > 0 and index[kindex] is Ellipsis:
                            lellipsis -= 1
                        else:
                            kindex += 1
        return res


    for k in ['all', 'any', 'argmax', 'argmin', 'cumsum', 'cumprod', 'flatten', 
               'max', 'mean', 'min', 'nonzero', 'prod', 'ptp', 'ravel', 'repeat', 
               'reshape', 'resize', 'squeeze', 'std', 'sum', 'swapaxes', 'take', 
               'transpose', 'var']:
        exec k + '.__doc__ = np.ndarray.' + k + '.__doc__'

        
def benchmark(expression):
    '''transfer of axistags causes a slowdown by a factor of about 10,
       when getitem returns a value, the slowdown is about 3 (due to Python calls)
    '''
    import timeit, axistags
    reload(axistags)
    repetitions = 100000
    t1 = timeit.Timer(expression, 
         "import numpy, axistags\na = axistags.InfoArray((2,3,4), axistags='zyx', dtype=numpy.uint8)")
    t2 = timeit.Timer(expression, 
         "import numpy, axistags\na = numpy.ndarray((2,3,4), dtype=numpy.uint8)")
    t3 = timeit.Timer(expression, 
         "import numpy, axistags\na = axistags.InfoArray((2,3,4), axistags='zyx', dtype=numpy.uint8).view(numpy.ndarray)")
    print "InfoArray:", t1.timeit(repetitions)/repetitions*1e6,"musec"
    print "ndarray:", t2.timeit(repetitions)/repetitions*1e6,"musec"
    print "InfoArray as ndarray:", t3.timeit(repetitions)/repetitions*1e6,"musec"
