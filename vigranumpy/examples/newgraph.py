import vigra
import vigra.graphs as vigraph
import pylab
import numpy






print "get input"

f       = '100075.jpg'
sigma   = 1.0
img     = vigra.impex.readImage(f)[0:100,0:100,:]
gradmag = vigra.filters.gaussianGradientMagnitude(img,sigma)
gradmag = numpy.squeeze(gradmag).astype(numpy.float32)


print "slic"
labels ,nseg = vigra.analysis.slicSuperpixels(vigra.colors.transform_RGB2Lab(img),10.0,15)
labels 		 = numpy.squeeze(vigra.analysis.labelImage(labels))






def hierarchicalSuperpixels(labels,edgeIndicator):
	print "gridGraph"
	gridGraph = vigraph.gridGraph(img.shape[0:2])

	print "get region adjacency graph"
	rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)

	print "get the sizes of the hyperedges "
	hyperEdgeSizes=vigraph.hyperEdgeSizes(rag,hyperEdges)

	print "get the mean features (gradmag) for hyperEdges"
	hyperEdgeGradMag  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
	hyperEdgeGradMag  = vigraph.hyperEdgeImageFeatures(rag,hyperEdges,gradmag,hyperEdgeGradMag)

	print "get multichannel out map"
	hyperEdgeImgFeat  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=3)
	hyperEdgeImgFeat  = vigraph.hyperEdgeImageFeatures(rag,hyperEdges,img,hyperEdgeImgFeat)


	print "get hierarchicalClustering W.I.P."

	print "maxEdgeId",rag.maxEdgeId
	print "maxNodeId",rag.maxNodeId

	print "edgeNum",rag.edgeNum
	print "nodeNum",rag.nodeNum


	hc = vigraph.hierarchicalClusteringAdjacencyListGraph(
		rag,
		hyperEdgeGradMag,
		hyperEdgeSizes
	)
	print hc

	hc.cluster()





















print "gridGraph"
gridGraph = vigraph.gridGraph(img.shape[0:2])

print "get region adjacency graph"
rag,hyperEdges = vigraph.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)

print "get the sizes of the hyperedges "
hyperEdgeSizes=vigraph.hyperEdgeSizes(rag,hyperEdges)

print "get the mean features (gradmag) for hyperEdges"
hyperEdgeGradMag  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)
hyperEdgeGradMag  = vigraph.hyperEdgeImageFeatures(rag,hyperEdges,gradmag,hyperEdgeGradMag)

print "get multichannel out map"
hyperEdgeImgFeat  = vigraph.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=3)
hyperEdgeImgFeat  = vigraph.hyperEdgeImageFeatures(rag,hyperEdges,img,hyperEdgeImgFeat)


print "get hierarchicalClustering W.I.P."

print "maxEdgeId",rag.maxEdgeId
print "maxNodeId",rag.maxNodeId

print "edgeNum",rag.edgeNum
print "nodeNum",rag.nodeNum


hc = vigraph.hierarchicalClusteringAdjacencyListGraph(
	rag,
	hyperEdgeGradMag,
	hyperEdgeSizes
)
print hc

hc.cluster()