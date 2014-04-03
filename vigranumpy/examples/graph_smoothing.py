import vigra
from vigra import graphs

# parameter:
filepath = '12003.jpg'   # input image path
sigmaGradMag = 2.0       # gradient magnitude scale
superpixelDiameter = 10  # super-pixel size
slicWeight = 10.0        # SLIC color - spatial weight
gamma = 0.15             # exp(-gamma * edgeIndicator)
edgeThreshold = 2.5      # values higher are considered as edges
scale = 1.0              # how much smoothing
iterations = 10          # how man smoothing iterations

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
edgeIndicator = rag.accumulateEdgeFeatures(gridGraphEdgeIndicator)

# accumulate node features from grid graph node map
# which is just a plain image (with channels)
nodeFeatures = rag.accumulateNodeFeatures(imgLab)
resultFeatures = graphs.recursiveGraphSmoothing(rag, nodeFeatures,
                                                edgeIndicator,gamma=gamma, 
                                                edgeThreshold=edgeThreshold,
                                                scale=scale,
                                                iterations=iterations)

resultImgLab = rag.projectNodeFeaturesToGridGraph(resultFeatures)
resultImgLab = vigra.taggedView(resultImgLab, "xyc")

vigra.imshow(vigra.colors.transform_Lab2RGB(resultImgLab))
vigra.show()
