import vigra
from vigra import numpy
from matplotlib import pylab
from time import time
import multiprocessing


path = "12003.jpg"
data = vigra.impex.readImage(path).astype(numpy.float32)

cpus = multiprocessing.cpu_count()
policy = vigra.filters.RatioPolicy(sigma=5.0, meanRatio=0.7, varRatio=0.5)
res = vigra.filters.nonLocalMean2d(data,policy=policy,searchRadius=4,patchRadius=2,nThreads=cpus+1,stepSize=2,verbose=True,sigmaMean=1.0)
res = vigra.taggedView(res,'xyc')
vigra.imshow(res)
vigra.show()
