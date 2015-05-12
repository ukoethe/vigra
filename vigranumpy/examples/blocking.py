import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer
from vigra import blockwise as bw




numpy.random.seed(42)

# input
shape = (500, 500, 500)

data = numpy.random.rand(*shape).astype('float32')

print "make options object"
options = bw.BlockwiseConvolutionOptions3D()
print type(options)

sigma = 1.0
options.stdDev = (sigma, )*3
options.blockShape = (128, )*3

print "stddev",options.stdDev
print "call blockwise filter"

with vigra.Timer("AllThread"):
	res = bw.gaussianSmooth(data, options)
with vigra.Timer("1thread"):
	resRef = vigra.gaussianSmoothing(data, sigma)


print numpy.sum(numpy.abs(res-resRef))
