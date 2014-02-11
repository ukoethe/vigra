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
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 3.0

img          = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab       = vigra.colors.transform_RGB2Lab(img)
gradmag      = numpy.squeeze(vigra.filters.gaussianGradientMagnitude(imgLab,sigma))
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,10)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))


gridGraph      = vigraph.gridGraph(img.shape[0:2])
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)
hyperNodeSizes = vigraph.hyperNodeSizes(rag,gridGraph,labels)
nodeFeatures   = vigraph.graphMap(graph=rag,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(rag,gridGraph,labels,img,nodeFeatures)

# visualize node features
imgOut=img.copy()
imageNodeFeatures = vigraph.nodeIdsFeatures(graph=rag,nodeIds=labels,features=nodeFeatures,out=imgOut)
imageNodeFeatures = vigra.taggedView(imageNodeFeatures,"xyc")
vigra.imshow(imageNodeFeatures)
vigra.show()


edgeWeights = edgeIndicator 
edgeWeights+= 0.01*vigraph.nodeFeatureDistToEdgeWeight(rag,nodeFeatures,metric='l1')

## normalize edge weights to [0-1]
#edgeWeights-=edgeWeights[rag.edgeIds()].min()
#edgeWeights/=edgeWeights[rag.edgeIds()].max()
#edgeWeights*=-1.0
#edgeWeights+=1.0
#edgeWeights[numpy.where(edgeWeights<0.8)]=0
#edgeWeights*=0.04
#

#edgeWeights = numpy.exp(-0.1*edgeWeights)
#edgeWeights[numpy.where(edgeWeights<0.75)]=0

print edgeWeights[rag.edgeIds()].min()
print edgeWeights[rag.edgeIds()].max()


resultFeatures = vigraph.recursiveGraphSmoothing(rag,nodeFeatures,edgeWeights,
    gamma=0.0001,edgeThreshold=3.8,scale=1.3,iterations=100)


# visualize node features
imageNodeFeatures = vigraph.nodeIdsFeatures(graph=rag,nodeIds=labels,features=resultFeatures)
imageNodeFeatures = vigra.taggedView(imageNodeFeatures,"xyc")
vigra.imshow(imageNodeFeatures)
vigra.show()