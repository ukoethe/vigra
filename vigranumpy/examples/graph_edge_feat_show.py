import vigra
import vigra.graphs as graphs
import pylab
import numpy
import matplotlib

# parameter:
filepath = '100075.jpg'   # input image path
sigmaGradMag = 3.0       # sigma Gaussian gradient
superpixelDiameter = 20  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight

# load image and convert to LAB
img = vigra.impex.readImage(filepath)

# get super-pixels with slic on LAB image
imgLab = vigra.colors.transform_RGB2Lab(img)
labels, nseg = vigra.analysis.slicSuperpixels(imgLab, slicWeight,
                                              superpixelDiameter)
labels = vigra.analysis.labelImage(labels)

# compute gradient
imgLabBig = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradMag    = vigra.filters.gaussianGradientMagnitude(imgLab, sigmaGradMag)
gradMagBig = vigra.filters.gaussianGradientMagnitude(imgLabBig, sigmaGradMag*2.0)



# get 2D grid graph and  edgeMap for grid graph
# from gradMag of interpolated image
gridGraph = graphs.gridGraph(img.shape[0:2])
gridGraphEdgeIndicator = graphs.edgeFeaturesFromInterpolatedImage(gridGraph,
                                                                  gradMagBig)

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)


# accumulate edge  and ndie weights from gradient magnitude
ragEdgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)
ragNodeWeights = rag.accumulateNodeFeatures(gradMag)


rag.showEdgeFeature(img, ragEdgeWeights)
vigra.show()