import vigra
from vigra import graphs
from vigra import numpy
from vigra import Timer





numpy.random.seed(42)



# input
shape = [100, 100, 100]
data = numpy.random.rand(*shape).astype(numpy.float32)

with Timer("get labels"):
    labels = numpy.random.randint(5, size=shape[0]*shape[1]*shape[2])
    labels = labels.reshape(shape).astype(numpy.uint32)

with Timer("labelVolume"):
    labels = vigra.analysis.labelVolume(labels)


with Timer("findMinMax"):
    minLabel = int(labels.min())
    maxLabel = int(labels.max())

with Timer("ragOptions"):
    ragOptions = graphs.RagOptions(minLabel=minLabel, maxLabel=maxLabel, isDense=True,
                                   nThreads=1,  parallel=False)

with Timer("makeRAG"):

    adjListGraph  = graphs.listGraph()
    gridGraph = graphs.gridGraph(shape)




    rag = graphs.regionAdjacencyGraph(gridGraph, labels)

    rag.writeHdf5("bla.h5", "dset")



with Timer("readRag"):
