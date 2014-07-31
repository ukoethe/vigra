import vigra
import vigra.graphs as graphs

# parameter:
filepath = '12003.jpg'   # input image path
sigmaGradMag = 2.0       # sigma Gaussian gradient
superpixelDiameter = 10  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight
k = 10                   # free parameter in felzenszwalbs method
nodeNumStop = 500        # desired num. nodes in result

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

# get region adjacency graph from super-pixel labels
rag = graphs.regionAdjacencyGraph(gridGraph, labels)

# accumulate edge weights from gradient magnitude
edgeWeights = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)

# do the segmentation with felzenszwalbs method
labels   = graphs.felzenszwalbSegmentation(rag, edgeWeights, 
                                           k=50,nodeNumStop=nodeNumStop)


rag.show(img, labels)
vigra.show()
