import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer









# input
shape = [300, 100, 30]
data = numpy.random.rand(*shape).astype(numpy.float32)

with Timer("watershed"):
    seg, nseg = vigra.analysis.watersheds(data)


with Timer("labelVolume"):
    seg = vigra.analysis.labelVolume(seg)


with Timer("findMinMax"):
    minLabel = int(seg.min())
    maxLabel = int(seg.max())

with Timer("ragOptions"):
    ragOptions = graphs.RagOptions(minLabel=minLabel, maxLabel=maxLabel, isDense=True,
                                   nThreads=3,  parallel=True)

with Timer("makeRAG"):

    adjListGraph  = graphs.listGraph()
    gridGraph = graphs.gridGraph(shape)


    affEdges = graphs._makeRag(graph=gridGraph, labels=seg, rag=adjListGraph, options=ragOptions)