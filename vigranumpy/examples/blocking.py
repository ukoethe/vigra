import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer
from vigra import blockwise as bw




numpy.random.seed(42)

# input
shape = (100, 500, 500)
blockShape = (50, 50, 50)

data = numpy.random.rand(*shape).astype('float32')


blocking = bw.Blocking3D(shape, blockShape)


resRef = vigra.gaussianSmoothing(data, 2.0)
res = bw.gaussianSmooth(data, 2.0, blocking)

print res-resRef
