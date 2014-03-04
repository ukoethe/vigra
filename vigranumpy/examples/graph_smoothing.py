import vigra
import vigra.graphs as vigraph
<<<<<<< HEAD
import numpy
import pylab
=======
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import math

>>>>>>> 95c348d188e7e481cbba03834d574b2acc42cf38

print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 4.0

<<<<<<< HEAD
img          = vigra.impex.readImage(f)
imgLab       = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradmag      = numpy.squeeze(vigra.filters.gaussianGradientMagnitude(imgLabInterpolated, sigma))
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, 10.0, 20)
=======
img          = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab       = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmag      = numpy.squeeze(vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma))
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,20)
>>>>>>> 95c348d188e7e481cbba03834d574b2acc42cf38
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))


gridGraph = vigraph.gridGraph(img.shape[0:2])
<<<<<<< HEAD
gridGraphEdgeIndicator = vigraph.edgeFeaturesFromInterpolatedImage(gridGraph, gradmag)


rag = vigraph.regionAdjacencyGraph(graph=gridGraph, labels=labels, ignoreLabel=0)
=======
gridGraphEdgeIndicator =  vigraph.edgeFeaturesFromInterpolatedImage(gridGraph,gradmag)



rag = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
>>>>>>> 95c348d188e7e481cbba03834d574b2acc42cf38
edgeIndicator  = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)
nodeFeatures   = rag.accumulateNodeFeatures(img)


<<<<<<< HEAD
resultFeatures = vigraph.recursiveGraphSmoothing(rag, nodeFeatures, edgeIndicator,
                                                 gamma=0.0001, edgeThreshold=2.0, scale=1.0, iterations=20)

nodeFeaturesAsImg = rag.projectNodeFeaturesToGridGraph(nodeFeatures)
nodeFeaturesAsImg = vigra.taggedView(nodeFeaturesAsImg, "xyc")

resultImg = rag.projectNodeFeaturesToGridGraph(resultFeatures)
resultImg = vigra.taggedView(resultImg, "xyc")

f = pylab.figure()
ax0 = f.add_subplot(1, 3, 1)
vigra.imshow(img, show=False)
ax0.set_title("raw image")

ax1 = f.add_subplot(1, 3, 2)
vigra.imshow(nodeFeaturesAsImg, show=False)
ax1.set_title("superpixelized image")

ax2 = f.add_subplot(1, 3, 0)
vigra.imshow(resultImg, show=False)
ax2.set_title("smoothed image")

pylab.show()
=======

resultFeatures = vigraph.recursiveGraphSmoothing(rag,nodeFeatures,edgeIndicator,
    gamma=0.0001,edgeThreshold=2.0,scale=1.0,iterations=20)


resultImg = rag.projectNodeFeaturesToGridGraph(resultFeatures)
resultImg = vigra.taggedView(resultImg,"xyc")

rag.show(resultImg,alpha=0.7)

#vigra.imshow(resultImg)
vigra.show()
>>>>>>> 95c348d188e7e481cbba03834d574b2acc42cf38
