import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer
from vigra import blockwise as bw




numpy.random.seed(42)

# input
shape = (1000, 1000, 500)
#blockShape = (50, 50, 50)

data = numpy.random.rand(*shape).astype('float32')

print "make options object"
options = bw.BlockwiseConvolutionOptions3D()
print type(options)

options.stdDev = (1.0,)*3

print "call blockwise filter"

with vigra.Timer("1thread"):
	resRef = vigra.gaussianSmoothing(data, 3.0)
with vigra.Timer("AllThread"):
	res = bw.gaussianSmooth(data, options)

print res-resRef

print numpy.sum(numpy.abs(res-resRef))



resRef = vigra.gaussianSmoothing(data, opt.setSigm(3.0).setNumCores(10).blockShape(64))