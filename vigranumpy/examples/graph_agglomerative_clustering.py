import vigra
from vigra import graphs
from vigra import numpy
import pylab
# parameter
filepath = '12003.jpg'  # input image path
sigmaGradMag = 2.0      # sigma Gaussian gradient
superpixelDiameter = 10 # super-pixel size
slicWeight = 10.0       # SLIC color - spatial weight
beta = 0.5              # node vs edge weight
nodeNumStop = 50        # desired num. nodes in result


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

# accumulate node features from grid graph node map
# which is just a plain image (with channels)
nodeFeatures = rag.accumulateNodeFeatures(imgLab)

# do agglomerativeClustering
labels = graphs.agglomerativeClustering(graph=rag, edgeWeights=edgeWeights,
                                        beta=beta, nodeFeatures=nodeFeatures,
                                        nodeNumStop=nodeNumStop,wardness=0.0)

# show result
f = pylab.figure()
ax1 = f.add_subplot(2, 2, 1)
vigra.imshow(gradMag,show=False)
ax1.set_title("Input Image")
pylab.axis('off')

ax2 = f.add_subplot(2, 2, 2)
rag.show(img)
ax2.set_title("Over-Segmentation")
pylab.axis('off')

ax3 = f.add_subplot(2, 2, 3)
rag.show(img, labels)
ax3.set_title("Result-Segmentation")
pylab.axis('off')

ax4 = f.add_subplot(2, 2, 4)
rag.showNested(img, labels)
ax4.set_title("Result-Segmentation")
pylab.axis('off')

vigra.show()
