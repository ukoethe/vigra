import vigra
from vigra import numpy
from matplotlib import pylab
from time import time
import multiprocessing


path = "12003.jpg"
data = vigra.impex.readImage(path).astype(numpy.float32)
data = vigra.taggedView(100*numpy.random.rand(*data.shape),'xyc').astype('float32') + data
data /= 2.0
vigra.imshow(data)
vigra.show()
cpus = multiprocessing.cpu_count()
policy = vigra.filters.NormPolicy(sigma=10.0, meanDist=300.7, varRatio=0.9)
res = vigra.filters.nonLocalMean2d(data,policy=policy,searchRadius=8,patchRadius=2,nThreads=cpus+1,stepSize=1,verbose=True,sigmaMean=1.0)
res = vigra.taggedView(res,'xyc')
vigra.imshow(res)
vigra.show()
