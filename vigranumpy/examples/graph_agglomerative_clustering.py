import vigra
from vigra import graphs
from vigra import numpy

# parameter
filepath = '12003.jpg'
sigmaGRadMag = 2.0
superpixelDiameter = 3
slicWeight = 10.0
beta = 0.5
nodeNumStop = 50

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, slicWeight,
                                              superpixelDiameter)
labels = numpy.squeeze(vigra.analysis.labelImage(labels))

# compute gradient on interpolated image
imgLabBig = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradMag = vigra.filters.gaussianGradientMagnitude(imgLabBig, sigmaGRadMag)

# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator = graphs.edgeFeaturesFromInterpolatedImage(gridGraph,
                                                                  gradMag)

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

# accumulate edge weights from gradient magnitude
edgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)

# accumulate node features from grid graph node map
# which is just a plain image (with channels)
nodeFeatures = rag.accumulateNodeFeatures(imgLab)

# do agglomerativeClustering
labels = graphs.agglomerativeClustering(graph=rag, edgeWeights=edgeWeights,
                                        beta=beta, nodeFeatures=nodeFeatures,
                                        nodeNumStop=nodeNumStop)

# show result
rag.showNested(img, labels)
vigra.show()
