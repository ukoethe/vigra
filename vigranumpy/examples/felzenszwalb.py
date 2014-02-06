import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab as plt
import math



print "get input"
f       = '100075.jpg'
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 3.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
img     = vigra.filters.gaussianSmoothing(img,0.8)
imgLab  = vigra.colors.transform_RGB2Lab(img)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag)

print gradmag.shape,gradmag.dtype
    
print "get histogram"
hist = vigra.histogram.gaussianHistogram(img,
    minVals=[0.0,0.0,0.0],maxVals=[255.01,255.01,255.01],
    bins=10,sigma=3.0,sigmaBin=1.0).reshape([img.shape[0],img.shape[1],-1])

hsum = numpy.sum(hist,axis=2)
hist/=hsum[:,:,numpy.newaxis]




print "get felzenszwalbSegmentation"
labels = vigra.graphs.felzenszwalbSegmentation(
    nodeFeatures=hist,edgeWeights=0.01*gradmag,k=1.0,
    distanceType="l1")

labels = vigra.analysis.labelImage(labels)

vigra.showSeg(img,labels)
vigra.show()