import vigra
from vigra import numpy
from matplotlib import pylab
from time import time
import multiprocessing

path = "69015.jpg"
#path = "12074.jpg"
path = "100075.jpg" 
path = "12003.jpg"
data = vigra.impex.readImage(path).astype(numpy.float32)

cpus = multiprocessing.cpu_count()

print "nCpus",cpus

t0 =time()

#for c in range(3):
#    cimg=data[:,:,c]
#    cimg-=cimg.min()
#    cimg/=cimg.max()


iters = 10

#policy = vigra.filters.RatioPolicy(sigma=10.0, meanRatio=0.95, varRatio=0.5)
policy = vigra.filters.NormPolicy(sigma=50.0, meanDist=50, varRatio=0.5)
#data-=100.0
res = vigra.filters.nonLocalMean2d(data,policy=policy,searchRadius=5,patchRadius=1,nThreads=cpus+1,stepSize=2,verbose=True,sigmaMean=10.0)
for i in range(iters-1):
    res = vigra.filters.nonLocalMean2d(res,policy=policy,searchRadius=5,patchRadius=2,nThreads=cpus+1,stepSize=2,verbose=True,sigmaMean=10.0)
t1 = time()

res = vigra.taggedView(res,'xyc')
gma = vigra.filters.gaussianGradientMagnitude(res,4.0)
gmb = vigra.filters.gaussianGradientMagnitude(data,4.0)
#data+=100.0
print t1-t0
imgs  = [data,res,gma,gmb]

for img in imgs:
    for c in range(img.shape[2]):
        cimg=img[:,:,c]
        cimg-=cimg.min()
        cimg/=cimg.max()

f = pylab.figure()
for n, arr in enumerate(imgs):
    arr = arr.squeeze()
    f.add_subplot(1, len(imgs), n)
    pylab.imshow(arr.swapaxes(0,1))

pylab.title('denoised')
pylab.show()
