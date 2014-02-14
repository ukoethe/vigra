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
imgLabInterpolated  = vigra.resize(imgLab,[imgLab.shape[0]*2-1,imgLab.shape[1]*2-1 ])
gradmagInterpolated = vigra.filters.gaussianGradientMagnitude(imgLabInterpolated,sigma)
labels ,nseg 		= vigra.analysis.slicSuperpixels(imgLab,10.0,5)
labels       		= vigra.analysis.labelImage(labels)

print "get graph1 and grid graph "
graph0,graph1 = vigraph.gridRegionAdjacencyGraph(labels=labels,ignoreLabel=None)

graph1.show(img)
vigra.show()

print "ragshape ", graph1.shape


# get grid graph and edge weights
print "get grid graph edge weights"
graph0EdgeWeights =  vigraph.edgeFeaturesFromInterpolatedImage(graph0,gradmagInterpolated)

print "get grid graph edge weights and sizes"
graph1EdgeSize    = graph1.accumulateEdgeSize()
graph1EdgeWeights = graph1.accumulateEdgeFeatures(graph0EdgeWeights,acc='mean')

print "get grid graph node features and sizes"
graph1NodeSize     = graph1.accumulateNodeSize()
graph1NodeFeatures = graph1.accumulateNodeFeatures(img,acc='mean')



print "do felzenszwalbSegmentation"
graph1Labels = vigraph.felzenszwalbSegmentation(graph1,graph1EdgeWeights,k=1)
print "get graph 2"
graph2       = vigraph.regionAdjacencyGraph(graph=graph1,labels=graph1Labels,ignoreLabel=None)
graph2EdgeWeights = graph2.accumulateEdgeFeatures(graph1EdgeWeights,acc='mean')


graph1.show(img,graph1Labels)
vigra.show()


graph2.show(img)
vigra.show()

graph2Labels = vigraph.felzenszwalbSegmentation(graph2,graph2EdgeWeights,k=5)
print "get graph 3"
graph3       = vigraph.regionAdjacencyGraph(graph=graph2,labels=graph2Labels,ignoreLabel=None)

graph3.show(img)
vigra.show()

sys.exit(0)







# visualize node features
imgOut=img.copy()
projectedRagNodeFeatures = vigraph.nodeIdsFeatures(graph=graph1,nodeIds=labels,features=ragNodeFeatures,out=imgOut)
print projectedRagNodeFeatures.shape
projectedRagNodeFeatures = vigra.taggedView(projectedRagNodeFeatures,"xyc")
vigra.imshow(projectedRagNodeFeatures)
vigra.show()











"""
edgeIndicator  = vigraph.graphMap(graph=graph1,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(graph1,graph0,hyperEdges,gradmag,edgeIndicator)
hyperNodeSizes = vigraph.hyperNodeSizes(graph1,graph0,labels)
nodeFeatures   = vigraph.graphMap(graph=graph1,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(graph1,graph0,labels,img,nodeFeatures)

# visualize node features
imgOut=img.copy()
imageNodeFeatures = vigraph.nodeIdsFeatures(graph=graph1,nodeIds=labels,features=nodeFeatures,out=imgOut)
imageNodeFeatures = vigra.taggedView(imageNodeFeatures,"xyc")
vigra.imshow(imageNodeFeatures)
vigra.show()


edgeWeights = edgeIndicator 
edgeWeights+= 0.01*vigraph.nodeFeatureDistToEdgeWeight(graph1,nodeFeatures,metric='l1')

## normalize edge weights to [0-1]
#edgeWeights-=edgeWeights[graph1.edgeIds()].min()
#edgeWeights/=edgeWeights[graph1.edgeIds()].max()
#edgeWeights*=-1.0
#edgeWeights+=1.0
#edgeWeights[numpy.where(edgeWeights<0.8)]=0
#edgeWeights*=0.04
#

#edgeWeights = numpy.exp(-0.1*edgeWeights)
#edgeWeights[numpy.where(edgeWeights<0.75)]=0

print edgeWeights[graph1.edgeIds()].min()
print edgeWeights[graph1.edgeIds()].max()


resultFeatures = vigraph.recursiveGraphSmoothing(graph1,nodeFeatures,edgeWeights,
    gamma=0.0001,edgeThreshold=3.8,scale=1.3,iterations=100)


# visualize node features
imageNodeFeatures = vigraph.nodeIdsFeatures(graph=graph1,nodeIds=labels,features=resultFeatures)
imageNodeFeatures = vigra.taggedView(imageNodeFeatures,"xyc")
vigra.imshow(imageNodeFeatures)
vigra.show()
"""