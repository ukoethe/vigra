import vigra
from vigra import graphs
from vigra import numpy

# parameter:
filepath = '100075.jpg'   # input image path
sigmaGradMag = 3.0       # sigma Gaussian gradient
superpixelDiameter = 10  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, slicWeight,
                                              superpixelDiameter)
labels = vigra.analysis.labelImage(labels)

# compute gradient on interpolated image
imgLabBig = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradMag = vigra.filters.gaussianGradientMagnitude(imgLabBig, sigmaGradMag)

# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator = graphs.edgeFeaturesFromInterpolatedImage(gridGraph,
                                                                  gradMag)

# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator = graphs.edgeFeaturesFromInterpolatedImage(gridGraph,
                                                                  gradMag)

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

# accumulate edge weights from gradient magnitude
ragEdgeIndicator = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)


affEdges = rag.affiliatedEdges

aVec = affEdges.getEdgeVec(gridGraph, 0)
print aVec


for ei in aVec:
    e = gridGraph.edgeFromId(ei)
    na = gridGraph.u(e)
    nb = gridGraph.v(e)

    print na,nb

    print numpy.array(na.coord())+numpy.array(nb.coord())/2.0