import vigra
import pylab
import numpy
f       = '69015.jpg'
sigma   = 1.8


class SpecialMap(object):
    def __init__(self):
        pass

    def mergeEdges(self,a,b):
        print "merge Edges",a,b
    def mergeNodes(self,a,b):
        print "merge Nodes",a,b
    def eraseEdge(self,a):
        print "eraseEdge",a




print "get input"
img     = vigra.impex.readImage(f)
gradmag = vigra.filters.gaussianGradientMagnitude(img,sigma)
if False:
    vigra.imshow(gradmag)
    pylab.show()

print "watershed"
labels ,nseg = vigra.analysis.watersheds(gradmag)

print "get cgp"
grid  = vigra.cgp.TopologicalGrid(labels)
graph = vigra.cgp.Cgp(grid)




print "get edge features"
edgeWeights = graph.accumulateCellFeatures(cellType=1,image=vigra.sampling.resize(gradmag[:,:,0],graph.shape),features='Mean')[0]['Mean'].astype(numpy.float32)

ucmWeights  = edgeWeights.copy();



print "get node features"
nodeFeatures = graph.accumulateCellFeatures(cellType=2,image=vigra.sampling.resize(img,graph.shape),features='Mean')[0]['Mean'].astype(numpy.float32)



print edgeWeights

if False :
    vigra.visualize(vigra.sampling.resize(img,graph.shape),graph,edge_data_in=edgeWeights)

print "get merge graph"
mergeGraph = vigra.clustering.MergeGraph(graph.numCells(2),graph.numCells(1))
mergeGraph.setInitalEdges(graph.cell1BoundsArray()-1)

print "get cell sizes"
cell1Size = graph.cellSizes(1).astype(numpy.float32)
cell2Size = graph.cellSizes(2).astype(numpy.float32)

print "get cell size maps"
print cell1Size.shape
size1Map  = vigra.clustering.SumMap0(cell1Size)
size2Map  = vigra.clustering.SumMap0(cell2Size)

print "get maps"
edgeWeightMap     = vigra.clustering.WeightedMeanMap0(edgeWeights,size1Map)
nodeFeatureMaps   = vigra.clustering.WeightedMeanMap1(nodeFeatures,size2Map)
minWeightEdgeMap  = vigra.clustering.MinWeightEdgeMap(mergeGraph,ucmWeights,edgeWeightMap,nodeFeatureMaps)
specialMap        = SpecialMap()
pythonMap         = vigra.clustering.PythonMap(mergeGraph,specialMap)




print "connect edge maps"
mergeGraph.registerMergeEdgeCallBack(size1Map.mergeCallback())
mergeGraph.registerMergeEdgeCallBack(edgeWeightMap.mergeCallback())
mergeGraph.registerMergeEdgeCallBack(minWeightEdgeMap.mergeEdgeCallback())
mergeGraph.registerEraseEdgeCallBack(minWeightEdgeMap.eraseEdgeCallback())

mergeGraph.registerMergeNodeCallBack(pythonMap.mergeNodeCallback())
mergeGraph.registerMergeEdgeCallBack(pythonMap.mergeEdgeCallback())
mergeGraph.registerEraseEdgeCallBack(pythonMap.eraseEdgeCallback())

print "connect node maps"
mergeGraph.registerMergeNodeCallBack(size2Map.mergeCallback())
mergeGraph.registerMergeNodeCallBack(nodeFeatureMaps.mergeCallback())

mergeGraph.mergeParallelEdges()


print "start with ",mergeGraph.numberOfNodes(),"nodes"
while mergeGraph.numberOfNodes()!=1:
    minEdgeLabel = minWeightEdgeMap.minWeightEdgeLabel()
    mergeGraph.mergeRegions(minEdgeLabel)

    if mergeGraph.numberOfNodes()==2000:
        stateOfInitalEdges = mergeGraph.stateOfInitalEdges()

        grid2  = graph.merge2Cells(stateOfInitalEdges.astype(numpy.uint32))
        graph2 = vigra.cgp.Cgp(grid2)

        vigra.visualize(vigra.sampling.resize(img,graph2.shape),graph2)