import vigra
import vigra.graphs as vigraph
import pylab
import numpy






print "get input"

f       = '100075.jpg'
sigma   = 1.0
img     = vigra.impex.readImage(f)#[0:100,0:100,:]
gradmag = vigra.filters.gaussianGradientMagnitude(img,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)

print "slic"
labels ,nseg = vigra.analysis.slicSuperpixels(vigra.colors.transform_RGB2Lab(img),10.0,15)
labels 		 = numpy.squeeze(vigra.analysis.labelImage(labels))

print "gridGraph"
gridGraph = vigraph.gridGraph(img.shape[0:2])

print "get region adjacency graph"
rag,hyperNodes,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)

print "compute features from hyperNodes and hyperEdges"

print hyperNodes