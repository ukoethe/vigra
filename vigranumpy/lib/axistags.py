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

    def swapaxes(self, i, j):
        res = np.ndarray.swapaxes(self, i, j)
        res.axistags[i] = self.axistags[j]
        res.axistags[j] = self.axistags[i]
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
    
    @property
    def T(self):
        return self.transpose()

    def __getitem__(self, index):
        '''There are three cases:
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
