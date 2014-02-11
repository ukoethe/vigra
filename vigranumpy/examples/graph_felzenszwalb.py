import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import math


print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 4.0

img          = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab       = vigra.colors.transform_RGB2Lab(img)
gradmag      = numpy.squeeze(vigra.filters.gaussianGradientMagnitude(imgLab,sigma))
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,5)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))


gridGraph      = vigraph.gridGraph(img.shape[0:2])
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)
hyperNodeSizes = vigraph.hyperNodeSizes(rag,gridGraph,labels)
nodeFeatures   = vigraph.graphMap(graph=rag,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(rag,gridGraph,labels,img,nodeFeatures)



edgeWeights = edgeIndicator 
edgeWeights+= 0.1*vigraph.nodeFeatureDistToEdgeWeight(rag,nodeFeatures,metric='l1')


ragLabels   = vigraph.felzenszwalbSegmentation(rag,edgeWeights,k=50)
imageLabels = vigraph.nodeIdsLabels(graph=rag,nodeIds=labels,labels=ragLabels)


vigra.segShow(img,imageLabels)
vigra.show()
