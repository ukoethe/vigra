import numpy as np

class InfoArray(np.ndarray):

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0,
          strides=None, order=None, axistags=None):
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides, order)
        if axistags is None:
            defaulttags = [['x'], ['y','x'], ['z','y','x']]
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
        res = np.ndarray.__getitem__(self, index)
        if hasattr(res, 'axistags') and res.ndim != self.ndim:
            l1 = len(self.shape)
            try:
                l2 = len(index)
            except:
                index = (index, Ellipsis)
                l2 = 2
            k2 = 0
            res.axistags = []
            for k1 in range(l1):
                try:
                    # if index[k2] is int, the dimension is bound => drop this axis
                    int(index[k2]) 
                    k2 += 1
                except:
                    # copy the meaning of the axis
                    res.axistags.append(self.axistags[k1])
                    # the first ellipsis represents all missing axes
                    if l1 != l2 and type(index[k2]) is type(Ellipsis):
                        l2 += 1
                    else:
                        k2 += 1
        return res

def benchmark(expression):
    import timeit, infoarray
    reload(infoarray)
    repetitions = 100000
    t1 = timeit.Timer(expression, 
         "import numpy, infoarray\na = infoarray.InfoArray((2,3,4), axistags='zyx', dtype=numpy.uint8)")
    t2 = timeit.Timer(expression, 
         "import numpy, infoarray\na = numpy.ndarray((2,3,4), dtype=numpy.uint8)")
    t3 = timeit.Timer(expression, 
         "import numpy, infoarray\na = infoarray.InfoArray((2,3,4), axistags='zyx', dtype=numpy.uint8).view(numpy.ndarray)")
    print "InfoArray:", t1.timeit(repetitions)/repetitions*1e6,"musec"
    print "ndarray:", t2.timeit(repetitions)/repetitions*1e6,"musec"
    print "InfoArray as ndarray:", t3.timeit(repetitions)/repetitions*1e6,"musec"
