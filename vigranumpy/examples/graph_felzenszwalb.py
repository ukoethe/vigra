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

def labelImageFromRagLabels(ragLabels,spLabels,shape,out):

    for x in range(shape[0]):
        for y in range(shape[1]):
            ragId    = spLabels[x,y]
            out[x,y] = ragLabels[ragId]
    out = vigra.analysis.labelImage(out)
    return out

print "get input"
f       = '100075.jpg'
#f       = '69015.jpg'
f       = '12003.jpg'
sigma   = 4.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
imgLab  = vigra.colors.transform_RGB2Lab(img)
newShape = [img.shape[0]*2-1,img.shape[1]*2-1]
#imgLab     = vigra.resize(imgLab,newShape)
#img    = vigra.resize(img,newShape)
gradmag = vigra.filters.gaussianGradientMagnitude(imgLab,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)
labels ,nseg = vigra.analysis.slicSuperpixels(imgLab,10.0,5)
labels       = numpy.squeeze(vigra.analysis.labelImage(labels))

numLabels = labels.max()
newLabel  = labels.copy()

#vigra.segShow(img,labels)
#vigra.show()




gridGraph      = vigraph.gridGraph(img.shape[0:2])
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
edgeIndicator  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
edgeIndicator  = vigraph.hyperEdgeImageFeatures(rag,gridGraph,hyperEdges,gradmag,edgeIndicator)
hyperNodeSizes = vigraph.hyperNodeSizes(rag,gridGraph,labels)
nodeFeatures   = vigraph.graphMap(graph=rag,item='node',dtype=numpy.float32,channels=3)
nodeFeatures   = vigraph.hyperNodeImageFeatures(rag,gridGraph,labels,imgLab,nodeFeatures)


edgeWeights = edgeIndicator 
edgeWeights+= vigraph.nodeFeatureDistToEdgeWeight(rag,nodeFeatures,metric='l1')


ragLabels      = vigraph.felzenszwalbSegmentation(rag,edgeWeights,k=100)
newLabels=labels.copy()
newLabels      = labelImageFromRagLabels(ragLabels,labels,gradmag.shape,out=newLabels)

vigra.segShow(img,newLabels)
vigra.show()

