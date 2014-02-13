import vigra
import vigra.graphs as vigraph
import pylab
import numpy
import sys
import matplotlib
import pylab
import matplotlib.pyplot as plt
from   matplotlib.widgets import Button
import math


f       = '100075.jpg'
f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 3.0

print "prepare input"
img                 = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab       	    = vigra.colors.transform_RGB2Lab(img)
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2+1,imgLab.shape[1]*2+1 ])
gradmagInterpolated = vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma)
labels ,nseg 		= vigra.analysis.slicSuperpixels(imgLab,10.0,10)
labels       		= vigra.analysis.labelImage(labels)


print "get rag and grid graph "
gridGraph,rag = vigraph.gridRegionAdjacencyGraph(labels=labels,ignoreLabel=None)


print type(rag)


# get grid graph and edge weights
print "get grid graph edge weights"
gridGraphEdgeWeights =  vigraph.edgeFeaturesFromInterpolatedImage(gridGraph,gradmagInterpolated)

print "get grid graph edge weights and sizes"
ragEdgeSize    = rag.accumulateEdgeSize()
ragEdgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeWeights,acc='mean')

print "get grid graph node features and sizes"
ragNodeSize     = rag.accumulateNodeSize()
ragNodeFeatures = rag.accumulateNodeFeatures(img,acc='mean')


# visualize node features
imgOut=img.copy()
projectedRagNodeFeatures = vigraph.nodeIdsFeatures(graph=rag,nodeIds=labels,features=ragNodeFeatures,out=imgOut)
print projectedRagNodeFeatures.shape
projectedRagNodeFeatures = vigra.taggedView(projectedRagNodeFeatures,"xyc")
vigra.imshow(projectedRagNodeFeatures)
vigra.show()




"""
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
"""