import vigra
from vigra import graphs

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


pathFinder = graphs.ShortestPathPathDijkstra(rag)
pathFinder.run(ragEdgeIndicator,)



# get labels/segmentation for rag
ragLabels = graphs.felzenszwalbSegmentation(rag, ragEdgeIndicator,
                                            k=10, nodeNumStop=1000)

# get more corsair graph from labeled rag
rag2 = graphs.regionAdjacencyGraph(graph=rag, labels=ragLabels)

# accumulate new edge weights
rag2EdgeIndicator = rag2.accumulateEdgeFeatures(ragEdgeIndicator,
                                                acc='mean')

# get labels/segmentation for rag2
rag2Labels = graphs.felzenszwalbSegmentation(rag2, ragEdgeIndicator,
                                             k=20, nodeNumStop=100)

# get more corsair graph from labeled rag2
rag3 = graphs.regionAdjacencyGraph(graph=rag2, labels=rag2Labels)

# visualize results
for g in [rag, rag2, rag3]:
    print g.nodeNum
    g.show(img=img)
    vigra.show()