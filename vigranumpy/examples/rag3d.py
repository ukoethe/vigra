import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer
iEdgeMap = graphs.implicitMeanEdgeMap




numpy.random.seed(42)

# input
shape = [200, 100, 100]
data = numpy.random.rand(*shape).astype(numpy.float32)
data = vigra.taggedView(data,"xyz")

if False:
    labels = numpy.random.randint(5, size=shape[0]*shape[1]*shape[2])
    labels = labels.reshape(shape).astype(numpy.uint32)
    labels = vigra.analysis.labelVolume(labels)
    adjListGraph  = graphs.listGraph()
    gridGraph = graphs.gridGraph(shape)
    rag = graphs.regionAdjacencyGraph(gridGraph, labels)
    rag.writeHDF5("bla.h5", "dset")

else :

    # load the region adjacency graph
    rag = graphs.loadGridRagHDF5("bla.h5","dset")

print rag.labels.shape, rag.labels.dtype ,type(rag.labels)


print "accumulate edge and node features"


edgeCuesMean = rag.accumulateEdgeFeatures( iEdgeMap(rag.baseGraph, data) )
edgeCuesMean = numpy.array([edgeCuesMean, edgeCuesMean]).T

nodeCuesMean = rag.accumulateNodeFeatures(data)
nodeCuesMean = numpy.array([nodeCuesMean, nodeCuesMean]).T



mergeGraph = graphs.mergeGraph(rag)


featureManager = graphs.NeuroDynamicFeatures(rag, mergeGraph)

# assign features
print "edgeCuesShape", edgeCuesMean.shape

featureManager.assignEdgeCues(edgeCuesMean)
featureManager.assignNodeCues(nodeCuesMean)
featureManager.assignEdgeSizes(rag.edgeLengths())
featureManager.assignNodeSizes(rag.nodeSize())

# register the callback mechainsm
featureManager.registerCallbacks()

mgEdge = mergeGraph.edgeFromId(1)

print "edge features", featureManager.getFeatures( mergeGraph.edgeFromId(27885))
mergeGraph.contractEdge(mgEdge)
print "edge features", featureManager.getFeatures( mergeGraph.edgeFromId(27885))


#for edge in mergeGraph.edgeIter():
#    print "edge features", featureManager.getFeatures(edge)


