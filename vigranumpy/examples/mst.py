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
sigma   = 1.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
img     = vigra.filters.gaussianSmoothing(img,sigma)
imgLab  = vigra.colors.transform_RGB2Lab(img)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag)

print gradmag.shape,gradmag.dtype
    

labels = vigra.graphs.minimumSpanningTreeSegmentation(imgLab,nodeNumThreshold=10000)
labels = vigra.analysis.labelImage(labels)

cmap = matplotlib.colors.ListedColormap ( numpy.random.rand ( labels.max()+1,3))
plt.imshow ( numpy.swapaxes(labels,0,1), cmap = cmap)
plt.show()