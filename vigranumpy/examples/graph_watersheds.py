import vigra
import vigra.graphs as graphs
import pylab


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

# compute gradient
imgLabBig = vigra.resize(imgLab, [imgLab.shape[0]*2-1, imgLab.shape[1]*2-1])
gradMag    = vigra.filters.gaussianGradientMagnitude(imgLab, sigmaGradMag)
gradMagBig = vigra.filters.gaussianGradientMagnitude(imgLabBig, sigmaGradMag*2.0)

vigra.imshow(gradMagBig)
vigra.show()

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

# generate seeds
seeds = graphs.nodeWeightedWatershedsSeeds(rag, ragNodeWeights)

# node weighted watersheds
labelsNodeWeighted  = graphs.nodeWeightedWatersheds(rag, ragNodeWeights, seeds)

# edge weighted watersheds
labelsEdgeWeighted  = graphs.edgeWeightedWatersheds(rag, ragEdgeWeights, seeds)


f = pylab.figure()
ax0 = f.add_subplot(1, 2, 0)
rag.showNested(img, labelsNodeWeighted)
ax0.set_title("node weighted")

ax1 = f.add_subplot(1, 2, 1)
rag.showNested(img, labelsEdgeWeighted)
ax1.set_title("edge weighted")
pylab.show()
